# Synthetic sub-area inputs: n_muni municipalities of per_muni sub-areas each.
# The first half of the municipalities is drawn (they have sample rows and appear
# in `direct`); the others are out of sample.
.fixture_sub <- function(n_muni = 8L, per_muni = 6L, seed = 11L) {
  withr::local_seed(seed)   # the global random stream is restored on return
  n <- n_muni * per_muni
  pop <- data.table::data.table(
    code = sprintf("%011d", seq_len(n)),
    municipality = rep(seq_len(n_muni), each = per_muni),
    stratum = rep(c("urban", "rural"), length.out = n),
    n_households_census = sample(5:80, n, TRUE),
    n_persons_census    = sample(20:300, n, TRUE),
    z1 = rnorm(n), z2 = rnorm(n), z3 = runif(n))
  drawn <- seq_len(n_muni %/% 2L)
  smp <- pop[municipality %in% drawn,
             .(code, municipality,
               y = pmin(pmax(0.1 + 0.05 * z1 + rnorm(.N, 0, 0.02), 0), 1),
               wt = runif(.N, 50, 500), n_hh = 10L)]
  direct <- data.frame(area_id = drawn, outcome = "migrated_last_12m",
                       theta_direct = runif(length(drawn), .05, .2),
                       var_direct_smoothed = runif(length(drawn), 1e-4, 1e-3))
  list(pop = pop, smp = smp, direct = direct, drawn = drawn, n_muni = n_muni)
}

# Three municipalities with two sub-areas each, small enough to roll up by hand.
# Municipalities 1 and 2 are drawn (sample rows for three of the four sub-areas);
# municipality 3 is not. `share_working_age` is a person share, `z1` a household
# quantity.
.fixture_sub_hand <- function() {
  pop <- data.table::data.table(
    code = c("a1", "a2", "b1", "b2", "c1", "c2"),
    municipality = c(1L, 1L, 2L, 2L, 3L, 3L),
    stratum = rep(c("urban", "rural"), 3),
    n_households_census = c(10L, 30L, 20L, 20L, 10L, 30L),
    n_persons_census    = c(40L, 20L, 20L, 60L, 10L, 50L),
    share_working_age   = c(0.5, 0.3, 0.2, 0.6, 0.1, 0.4),
    z1                  = c(1, 3, 5, 7, 2, 4))
  smp <- data.table::data.table(
    code = c("a1", "a2", "b1"), municipality = c(1L, 1L, 2L),
    y = c(0.2, 0.6, 0.1), wt = c(4, 12, 10), n_hh = 5L)
  direct <- data.frame(area_id = c(1L, 2L), outcome = "migrated_last_12m",
                       theta_direct = c(0.5, 0.1),
                       var_direct_smoothed = c(1e-3, 1e-3))
  list(pop = pop, smp = smp, direct = direct)
}

# Run a fit without povmap's progress output (bootstrap bars and notes on stdout).
.quiet_fit <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(force(expr)))
  out
}

# Department 05, municipality 1. Sub-areas: block U1 (urban, dwellings 1-2),
# block U2 (urban, dwelling 3) and the rural section R (dwellings 4-5, two blocks
# that truncate to the same 11 characters).
.write_cnpv_fixture <- function(root) {
  d <- file.path(root, "Census2018", "05_Antioquia")
  dir.create(d, recursive = TRUE)
  u1 <- "0500110010000010000001"; u2 <- "0500110010000010000002"
  r1 <- "0500130010000010000001"; r2 <- "0500130010000010000099"
  dw <- data.table::data.table(
    U_DPTO = 5L, U_MPIO = 1L,
    UA_CLASE = c(1L, 1L, 1L, 3L, 3L),
    COD_ENCUESTAS = c(1L, 1L, 2L, 3L, 3L),
    U_VIVIENDA = c(1L, 2L, 1L, 1L, 2L))
  mgn <- data.table::copy(dw)[, COD_DANE_ANM := c(u1, u1, u2, r1, r2)]
  hog <- data.table::copy(dw)[, H_AGUA_COCIN := c(1L, 2L, 1L, 1L, 2L)]
  viv <- data.table::copy(dw)[, VA_EE := c(1L, 1L, 2L, 1L, 1L)]
  # persons: (dwelling index, P_EDADR, P_ALFABETA, P_NIVEL_ANOSR, P_TRABAJO, PA1_GRP_ETNIC)
  p <- data.frame(
    dw   = c(1, 1, 2, 2,   3, 3,   4, 4, 5),
    edad = c(4, 14, 3, 6,  5, 20,  8, 8, 2),
    alfa = c(2, 1, 2, 2,   1, 1,   1, 2, 2),
    niv  = c(8, 9, 9, 7,   3, 10,  7, 3, 0),
    trab = c(1, 1, 2, 2,   NA, NA, 1, 2, 2),
    etn  = c(5, 1, 9, 9,   9, 9,   1, 1, 9))
  per <- cbind(dw[p$dw, ],
               data.table::data.table(P_EDADR = p$edad, P_ALFABETA = p$alfa,
                                      P_NIVEL_ANOSR = p$niv, P_TRABAJO = p$trab,
                                      PA1_GRP_ETNIC = p$etn))
  w <- function(x, nm) data.table::fwrite(x, file.path(d, nm))
  w(mgn, "CNPV2018_MGN_05.CSV"); w(per, "CNPV2018_5PER_05.CSV")
  w(hog, "CNPV2018_2HOG_05.CSV"); w(viv, "CNPV2018_1VIV_05.CSV")
  list(u1 = u1, u2 = u2, r = substr(r1, 1, 11))
}

# Minimal ECV design object for direct / direct_validation tests.
.tiny_design_for_direct <- function() {
  structure(list(
    stage1 = list(unit = "municipality", selection = "pps_wor",
                  size_variable = "n_dwellings"),
    stage2 = list(unit = "segment", selection = "systematic_fractional",
                  target_segment_size = 10L, n_per_psu = 2L),
    stage3 = list(unit = "household", selection = "census"),
    strata = data.frame(stratum = c("99_urban", "99_rural"),
                        n_psu   = c(2L, 2L)),
    calibration = list(margins = "stratum")
  ), class = c("ecv_design_spec", "list"))
}

# Fixture for FH-robust / FH-arcsin method tests.
# Generates a 60-area (n_in + n_oos) synthetic direct estimate tibble with
# the columns expected by fit_fh_robust() and fit_fh_arcsin().
.fixture_robust <- function(n_in = 50L, n_oos = 10L, seed = 5042L) {
  withr::local_seed(seed)   # the global random stream is restored on return
  area_id <- as.integer(c(seq_len(n_in), 9000L + seq_len(n_oos)))
  x1 <- rnorm(n_in + n_oos, 0, 1)
  x2 <- rnorm(n_in + n_oos, 0, 1)
  beta <- c(0.05, 0.04, -0.02)
  theta_true <- plogis(beta[1] + beta[2] * x1 + beta[3] * x2)
  vardir <- rgamma(n_in + n_oos, shape = 5, rate = 100)
  theta_direct <- pmin(pmax(theta_true + sqrt(vardir) * rnorm(n_in + n_oos),
                             1e-3), 1 - 1e-3)
  direct <- tibble::tibble(
    area_id              = area_id[seq_len(n_in)],
    theta_direct         = theta_direct[seq_len(n_in)],
    var_direct_smoothed  = vardir[seq_len(n_in)],
    n_eff                = rpois(n_in, 60)
  )
  X_in <- tibble::tibble(municipality = area_id[seq_len(n_in)],
                         x1 = x1[seq_len(n_in)],
                         x2 = x2[seq_len(n_in)])
  X_oos <- tibble::tibble(municipality = area_id[(n_in + 1L):(n_in + n_oos)],
                          x1 = x1[(n_in + 1L):(n_in + n_oos)],
                          x2 = x2[(n_in + 1L):(n_in + n_oos)])
  list(direct = direct, X_in = X_in, X_oos = X_oos, n_in = n_in, n_oos = n_oos)
}

# Build a small synthetic household-level frame for tests.
#
# Mirrors the ECV PSU = municipality + segment structure:
#   - 2 strata (dept 99 x {urban, rural})
#   - 4 municipalities per stratum (8 total)
#   - 6 manzanas per municipality
#   - 5 households per manzana
# Total: 2 * 4 * 6 * 5 = 240 households.
# Each municipality has 30 dwellings, exactly 3 segments of 10 when
# target_size = 10 (no tail segment, so sample sizes are deterministic).
make_synthetic_frame <- function(seed = 42) {
  withr::local_seed(seed)   # the global random stream is restored on return
  strata <- expand.grid(
    dept = 99L,
    area = c("urban", "rural"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  k <- 0L
  for (s in seq_len(nrow(strata))) {
    p_int  <- if (strata$area[s] == "urban") 0.09 else 0.04
    p_intl <- if (strata$area[s] == "urban") 0.03 else 0.01
    for (mun in 1:4) {
      muni_code <- 100L + mun + ifelse(strata$area[s] == "rural", 10L, 0L)
      for (man in 1:6) for (hh in 1:5) {
        k <- k + 1L
        n_persons      <- sample(1:6, 1)
        n_mig_internal <- rbinom(1, n_persons, p_int)
        n_mig_intl     <- rbinom(1, max(n_persons - n_mig_internal, 0L),
                                 p_intl)
        n_mig_total    <- n_mig_internal + n_mig_intl
        rows[[k]] <- data.frame(
          dept                  = strata$dept[s],
          muni                  = muni_code,
          municipality          = strata$dept[s] * 1000L + muni_code,
          area                  = strata$area[s],
          seccion               = sprintf("S%02d_%s", mun, strata$area[s]),
          manzana               = sprintf("M%02d_%02d_%s", mun, man,
                                          strata$area[s]),
          hh_id                 = sprintf("H%05d", k),
          n_persons             = n_persons,
          n_pers_resp           = n_persons,
          n_mig_total           = n_mig_total,
          n_mig_internal        = n_mig_internal,
          n_mig_intl            = n_mig_intl,
          migrated_last_12m     = as.integer(n_mig_total > 0L),
          migrated_internal_12m = as.integer(n_mig_internal > 0L),
          migrated_intl_12m     = as.integer(n_mig_intl > 0L)
        )
      }
    }
  }
  fr <- do.call(rbind, rows)
  fr$stratum <- paste(fr$dept, fr$area, sep = "_")
  tibble::as_tibble(fr)
}

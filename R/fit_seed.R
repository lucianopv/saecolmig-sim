#' Per-fit seeds
#'
#' Every model fit of a replicate runs with its own seed, a function of the
#' fit's identity alone: its configuration, method, outcome and replicate.
#' `sae_fit_seed()` computes it, so a fit can be re-run on the same random
#' stream: pass the value as `seed` to [method_args()] or to the `fit_*()`
#' function directly.
#'
#' seed = 1e7 * config index + 1e6 * method index + 1e4 * outcome index + rep,
#' where each index is read from an explicit map ([sae_seed_config_index],
#' [sae_seed_method_index], [sae_seed_outcome_index]). The indexes are fixed
#' numbers, not positions, so adding a name never moves another; changing a
#' number changes every seed that uses it.
#'
#' A configuration that adds building covariates to a parent configuration
#' carries its parent's config index: `geo_bld` the index of `geo` (2) and
#' `geo_blk_bld` that of `geo_blk` (9). A fit with and a fit without the
#' building covariates then draw the same seed for the same method, outcome and
#' replicate, so their contrast is paired fit by fit (common random numbers).
#'
#' The seed does not depend on how the fits are split into jobs. Re-running a
#' fit with it reproduces the estimates bit for bit only on the same toolchain
#' (package versions, CmdStan, compiler and platform); elsewhere the Bayesian
#' fits agree at Monte Carlo level.
#'
#' @param config,method,outcome character vectors naming the configuration,
#'   method and outcome of each fit
#' @param rep integer vector of replicates, 1-9999
#' @return integer vector of seeds, recycled to the longest argument
#' @examples
#' sae_fit_seed("geo_bld", "bayes_fh_beta", "migrated_intl_12m", 7)
#' @export
sae_fit_seed <- function(config, method, outcome, rep) {
  n <- max(length(config), length(method), length(outcome), length(rep))
  config  <- rep_len(as.character(config), n)
  method  <- rep_len(as.character(method), n)
  outcome <- rep_len(as.character(outcome), n)
  rep     <- rep_len(as.integer(rep), n)
  index_of <- function(x, map, what) {
    bad <- unique(x[!x %in% names(map)])
    if (length(bad)) {
      stop("Unknown ", what, ": ", paste(bad, collapse = ", "), call. = FALSE)
    }
    unname(map[x])
  }
  stopifnot(!anyNA(rep), all(rep >= 1L), all(rep <= 9999L))
  10000000L * index_of(config, sae_seed_config_index, "config") +
    1000000L * index_of(method, sae_seed_method_index, "method") +
    10000L * index_of(outcome, sae_seed_outcome_index, "outcome") + rep
}

#' @rdname sae_fit_seed
#' @format named integer vectors giving the index of each configuration, method
#'   and outcome in the seed formula
#' @export
sae_seed_config_index <- c(
  geo_bld     = 2L,
  geo_blk_bld = 9L)

#' @rdname sae_fit_seed
#' @export
sae_seed_method_index <- c(
  fh             = 1L,
  fh_arcsin      = 2L,
  fh_spatial     = 3L,
  fh_robust      = 4L,
  bayes_fh_beta  = 5L,
  bayes_fh_heavy = 6L,
  xgb_sub        = 7L,
  megb_sub       = 8L)

#' @rdname sae_fit_seed
#' @export
sae_seed_outcome_index <- c(
  migrated_last_12m          = 1L,
  migrated_internal_12m      = 2L,
  migrated_intl_12m          = 3L,
  pers_migrated_12m          = 4L,
  pers_migrated_internal_12m = 5L,
  pers_migrated_intl_12m     = 6L)

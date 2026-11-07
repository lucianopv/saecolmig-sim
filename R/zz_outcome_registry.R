#' Registered simulation outcomes
#'
#' Six estimands derived from PA_VIVIA_1ANO. Household-level outcomes are
#' 0/1 any-member indicators (weighted means; they OVERLAP across types —
#' do not sum them). Person-level outcomes are weighted count ratios
#' (numerator/denominator), additive across types. Person weight equals the
#' household weight (whole-household enumeration, as in the real ECV).
#'
#' Column contract (do not rename):
#' \describe{
#'   \item{outcome}{unique outcome identifier}
#'   \item{level}{"household" or "person"}
#'   \item{indicator}{frame column for household-level 0/1 indicator (NA for
#'     person-level outcomes)}
#'   \item{numerator}{frame column for person-level count (NA for
#'     household-level outcomes)}
#'   \item{denominator}{frame column for person-level denominator (NA for
#'     household-level outcomes)}
#' }
#' @export
sae_outcomes <- tibble::tibble(
  outcome     = c("migrated_last_12m", "migrated_internal_12m",
                  "migrated_intl_12m",
                  "pers_migrated_12m", "pers_migrated_internal_12m",
                  "pers_migrated_intl_12m"),
  level       = c(rep("household", 3L), rep("person", 3L)),
  indicator   = c("migrated_last_12m", "migrated_internal_12m",
                  "migrated_intl_12m", NA, NA, NA),
  numerator   = c(NA, NA, NA,
                  "n_mig_total", "n_mig_internal", "n_mig_intl"),
  denominator = c(NA, NA, NA, rep("n_pers_resp", 3L))
)

#' Look up one outcome's definition as a one-row list
#'
#' @param name character; an entry of \code{sae_outcomes$outcome}
#' @return list with elements outcome, level, indicator, numerator, denominator
#' @export
outcome_spec <- function(name) {
  i <- match(name, sae_outcomes$outcome)
  if (is.na(i)) {
    stop("Unknown outcome: ", name, ". Registered: ",
         paste(sae_outcomes$outcome, collapse = ", "))
  }
  as.list(sae_outcomes[i, ])
}

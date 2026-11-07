#!/usr/bin/env Rscript
# data-raw/build_subarea_pop.R
# The sub-area population frame: one row per census block (urban) or rural
# section, with its census household and person counts. It is the prediction
# universe of the sub-area methods and the `pop_data` side of their povmap
# call, so it covers every sub-area of the country, not only those a replicate
# draws. Sub-area codes are defined in R/subarea_frame.R.
#
#   input   frame/national_household_frame.parquet   (data-raw/build_frame.R)
#   output  aux_tables/subarea_pop.parquet
#
# Prints sub-area counts and quantiles only.

suppressPackageStartupMessages({
  library(saecol2026); library(arrow); library(data.table)
})

t0 <- Sys.time()
frame_path <- file.path(output_path("frame"), "national_household_frame.parquet")
if (!file.exists(frame_path)) stop("missing input: ", frame_path)

fr <- as.data.table(read_parquet(
  frame_path,
  col_select = c("municipality", "area", "manzana", "seccion", "n_persons")))
cat("households:", format(nrow(fr), big.mark = ","), "\n")

pop <- build_subarea_pop(fr)

out <- file.path(output_path("aux_tables", create = TRUE), "subarea_pop.parquet")
write_parquet_safe(tibble::as_tibble(pop), out)

cat("\n--- sub-area frame -----------------------------------------------\n")
cat("sub-areas:            ", format(nrow(pop), big.mark = ","), "\n")
cat("  urban (blocks):     ", format(sum(pop$stratum == "urban"), big.mark = ","), "\n")
cat("  rural (sections):   ", format(sum(pop$stratum == "rural"), big.mark = ","), "\n")
cat("municipalities covered:", uniqueN(pop$municipality), "\n")
cat("households (sum):     ", format(sum(pop$n_households_census), big.mark = ","),
    " (frame:", format(nrow(fr), big.mark = ","), ")\n")
cat("persons (sum):        ", format(sum(pop$n_persons_census), big.mark = ","),
    " (frame:", format(sum(fr$n_persons, na.rm = TRUE), big.mark = ","), ")\n")
cat("\nhouseholds per sub-area (quantiles):\n")
print(round(quantile(pop$n_households_census, c(.05, .25, .5, .75, .95, 1))))
cat("\nsub-areas per municipality (quantiles):\n")
print(round(quantile(pop[, .N, by = municipality]$N, c(.05, .25, .5, .75, .95, 1))))

# Every household lands in exactly one sub-area.
stopifnot(sum(pop$n_households_census) == nrow(fr), uniqueN(pop$code) == nrow(pop))
cat("\nwrote:", out, "\n")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")

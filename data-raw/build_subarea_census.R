#!/usr/bin/env Rscript
# data-raw/build_subarea_census.R
# Census composition covariates per sub-area: the nine shares of the municipal
# matrix (share_urban aside, which the sub-area matrix resolves as `is_urban`),
# grouped by census block or rural section (build_census_subarea()). Records
# reach a sub-area through the MGN key the household frame uses, so a
# household falls in the same block here as in the population frame.
#
#   inputs  Census2018/ (via census_path()), aux_tables/subarea_pop.parquet
#   output  aux_tables/subarea_census.parquet
#           (the population frame's rows and columns plus the nine shares)
#
# The population frame defines the universe: census sub-areas without a
# household in the frame are dropped, and frame sub-areas without census
# records keep NA shares.
#
# Prints sub-area aggregates only.

suppressPackageStartupMessages({
  library(saecol2026); library(arrow); library(data.table)
})

t0 <- Sys.time()
pop_path <- file.path(output_path("aux_tables"), "subarea_pop.parquet")
if (!file.exists(pop_path)) stop("run data-raw/build_subarea_pop.R first: ", pop_path)
pop <- as.data.table(read_parquet(pop_path))

cen <- build_census_subarea()
shares <- grep("^share_", names(cen), value = TRUE)

cat("\n--- coverage against the sub-area frame ---------------------------\n")
cat("sub-areas in frame         :", format(nrow(pop), big.mark = ","), "\n")
cat("sub-areas in census extract:", format(nrow(cen), big.mark = ","), "\n")
cat("frame sub-areas with no census row:", length(setdiff(pop$code, cen$code)), "\n")

m <- merge(pop, cen[, c("code", shares), with = FALSE], by = "code", all.x = TRUE)

chk <- merge(pop[, .(code, n_households_census)], cen[, .(code, n_households_sub)],
             by = "code")
cat("household counts identical, frame vs census extract:",
    sprintf("%.2f%% of %s sub-areas\n",
            100 * mean(chk$n_households_census == chk$n_households_sub),
            format(nrow(chk), big.mark = ",")))

cc <- complete.cases(m[, ..shares])
cat("\ncomplete on all", length(shares), "shares:",
    format(sum(cc), big.mark = ","), sprintf("sub-areas (%.2f%%),", 100 * mean(cc)),
    sprintf("%.2f%% of households,", 100 * sum(m$n_households_census[cc]) /
              sum(m$n_households_census)),
    uniqueN(m$municipality[cc]), "of", uniqueN(m$municipality), "municipalities\n")
cat("\nshare distributions (complete sub-areas):\n")
for (s in shares) {
  v <- m[[s]][cc]
  cat(sprintf("  %-28s min %.3f  q25 %.3f  med %.3f  q75 %.3f  max %.3f\n",
              s, min(v), quantile(v, .25), median(v), quantile(v, .75), max(v)))
}

out <- file.path(output_path("aux_tables", create = TRUE), "subarea_census.parquet")
write_parquet_safe(tibble::as_tibble(m), out)
cat("\nwrote:", out, "\n")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")

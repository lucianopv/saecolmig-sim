# Environment conditions the suite depends on.

# Evaluate `expr` with the warnings whose message matches `pattern` muffled;
# any other warning is passed on. It is for warnings a test cannot expect
# because whether they occur is not fixed: the sampler diagnostics of
# deliberately short test chains (a low effective sample size), which depend on
# the toolchain, and emdi2's warning that sigma2_v sits at the limit of its
# search interval, which fires on some replicates of a loop and not others.
muffle_warnings <- function(expr, pattern) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(pattern, conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

# Skip a test that reads data under SAECOL_DATA_DIR when the variable is unset.
skip_if_no_data_dir <- function() {
  testthat::skip_if(!nzchar(Sys.getenv("SAECOL_DATA_DIR")),
                    "SAECOL_DATA_DIR is not set")
}

# Without CMDSTAN naming an installation, cmdstanr looks for CmdStan in the default
# folder and first probes for WSL by starting `wsl uname` through processx. On Linux
# that start fails, and processx signals its whole process group when it later
# collects the failed handle, which ends the R process with exit status 137 and no
# message. Naming the installation skips the probe.
local({
  if (nzchar(Sys.getenv("CMDSTAN"))) return(invisible(NULL))
  dirs <- Sys.glob(file.path(path.expand("~"), ".cmdstan", "cmdstan-*"))
  if (!length(dirs)) return(invisible(NULL))
  ver <- suppressWarnings(numeric_version(sub("^cmdstan-", "", basename(dirs)),
                                          strict = FALSE))
  Sys.setenv(CMDSTAN = dirs[order(ver, na.last = FALSE)][length(dirs)])
})

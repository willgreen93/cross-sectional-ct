## ============================================================================
## synthetic_truth.R -- the generating parameters for the synthetic scenarios
##
##   Rscript scripts/synthetic_truth.R          # writes data/synthetic_truth.rds
##
## Reads the generating values from the ARRAYS themselves -- each synthetic array is stored
## as list(pars_test, out), so it carries the parameters that produced it.
##
## WHY NOT FROM simulation_study.R
##
## It used to extract the pars_L3/M3/N3/O3 definitions from that script and evaluate them.
## That made the file an independent second record, which sounds like a virtue and was not:
## the definitions are edited by hand between runs, and on 2026-09-05 pars_M3 was edited in
## place while building the P/Q/R scenarios. The script then described P's kinetics under
## M's name, and the file disagreed with every M, N and O chain on six parameters. The
## script drifts; an array cannot drift from itself.
##
## What the file is still FOR, now that it is derived from the arrays: it is a snapshot,
## taken at a known time, of what each array contained. Regenerating an array while chains
## are fitted to the old one leaves the chains fitted to data nobody has a record of -- which
## happened to O on 2026-09-05, when regenerating P/Q/R redrew O's negative-binomial noise
## as a side effect. Comparing a chain's own x$truth against this snapshot catches that.
##
## Re-run this after generating or regenerating any synthetic array, and NOT at any other
## time -- re-running it after a silent regeneration would simply bless the new state.
## ============================================================================

DATA_DIR <- "data"
OUT      <- "data/synthetic_truth.rds"

f <- list.files(DATA_DIR, pattern = "^data_array_vl_.*_N\\.rds$", full.names = TRUE)
if (!length(f)) stop("no data arrays found in ", DATA_DIR)

truth <- list()
skipped <- character(0)
for (ff in f) {
  nm <- sub("^data_array_vl_(.*)_N\\.rds$", "\\1", basename(ff))
  z  <- readRDS(ff)
  ## The four real variants are bare arrays with no generating values -- there is no truth
  ## for real data. Only the synthetic ones are list(pars_test, out).
  if (!is.list(z) || is.null(z$pars_test)) { skipped <- c(skipped, nm); next }
  truth[[nm]] <- z$pars_test
  rm(z); invisible(gc(verbose = FALSE))
}

if (!length(truth)) stop("no synthetic arrays found (none carried $pars_test)")

saveRDS(truth, OUT)

cat(sprintf("%d scenarios written to %s\n", length(truth), OUT))
cat("  ", paste(sort(names(truth)), collapse = " "), "\n")
if (length(skipped))
  cat("   skipped (no generating values, as expected for real data): ",
      paste(sort(skipped), collapse = " "), "\n", sep = "")

## Spot check: each ladder rung shifts every theta by log(10), so rung 3 minus rung 0 should
## be 3 log(10) = 6.91 for any family that has four rungs.
for (fam in c("M", "N", "O", "P", "Q", "R")) {
  a <- truth[[paste0(fam, "0")]]; b <- truth[[paste0(fam, "3")]]
  if (is.null(a) || is.null(b)) next
  cat(sprintf("   %s: theta7 gap rung0 -> rung3 = %.2f (expect %.2f)\n",
              fam, b[["theta7"]] - a[["theta7"]], 3 * log(10)))
}

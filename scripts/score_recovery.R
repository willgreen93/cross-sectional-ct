## ============================================================================
## score_recovery.R -- did the fits recover the parameters that generated the data?
##
##   Rscript scripts/score_recovery.R                     every fit in outputs/
##   Rscript scripts/score_recovery.R outputs C3          one scenario
##   Rscript scripts/score_recovery.R some/other/dir      a different directory
##
## This is the point of the simulation study, so it is worth stating what is and is not being
## claimed. For each parameter the script reports whether the 95% credible interval contains the
## generating value, and the standardised error (posterior mean - truth) / posterior SD. With one
## chain per scenario that is a MARGINAL check, repeated across scenarios and parameters. It is
## not an estimate of coverage: coverage would need many replicates per scenario, which were not
## run, and a single interval either contains the truth or does not.
##
## The generating values come from the fit itself, as x$truth, which input_generator copies from
## the $pars_test vector stored inside the data array. So a fit is always scored against the data
## it was actually fitted to, even if the array has since been regenerated.
##
## The analysis window is taken from x$analysis_from when present and otherwise starts after the
## final covariance adaptation. Draws before it come from a different proposal and are warm-up.
## ============================================================================
args <- commandArgs(trailingOnly = TRUE)
DIR  <- if (length(args) >= 1 && dir.exists(args[1])) args[1] else "outputs"
WANT <- if (length(args) >= 1 && dir.exists(args[1])) args[-1] else args
if (!dir.exists(DIR)) stop("no such directory: ", DIR)

## The parameters the paper reports recovery for: six kinetics, three symptom timing, dispersion.
PARS <- c(a_bar = "eta_a", a_sigma = "omega_a", b_bar = "eta_b", b_sigma = "omega_b",
          vl_max_bar = "vl_peak_bar", vl_max_sigma = "vl_sigma",
          inc1 = "eta_T", inc2 = "omega_T", inc3 = "vl_thresh", rdisp = "xi")

files <- list.files(DIR, pattern = "\\.rds$", full.names = TRUE)
scen  <- sub("^([A-Za-z]+[0-9]*)_.*$", "\\1", basename(files))
keep  <- if (length(WANT)) scen %in% toupper(WANT) else rep(TRUE, length(files))
files <- files[keep]; scen <- scen[keep]
if (!length(files))
  stop("no fits found in ", DIR, if (length(WANT)) paste0(" for ", paste(WANT, collapse = " ")) else "")
## One fit per scenario: if a scenario has several checkpoints, take the furthest advanced.
ord   <- order(scen, -file.info(files)$size)
files <- files[ord]; scen <- scen[ord]
files <- files[!duplicated(scen)]; scen <- scen[!duplicated(scen)]

cat(sprintf("scoring %d fit(s) in %s\n", length(files), DIR))
all <- NULL
for (i in seq_along(files)) {
  o <- readRDS(files[i]); x <- o$x
  if (is.null(x$truth) || all(is.na(x$truth))) {
    cat(sprintf("\n%s: no generating values stored -- not a synthetic fit, skipped\n", scen[i]))
    rm(o, x); invisible(gc(verbose = FALSE)); next
  }
  fl <- max(which(!is.na(o$MCMC_posteriors[, 1])))
  st <- if (!is.null(x$analysis_from)) x$analysis_from else max(x$cov_start) + 1L
  if (st > fl) {
    cat(sprintf("\n%s: only %s draws, analysis window starts at %s -- still warming up\n",
                scen[i], format(fl, big.mark = ","), format(st, big.mark = ",")))
    rm(o, x); invisible(gc(verbose = FALSE)); next
  }
  w <- st:fl
  cat(sprintf("\n== %s == %s draws, scored on %s after warm-up (from %s)\n",
              scen[i], format(fl, big.mark = ","), format(length(w), big.mark = ","),
              format(st, big.mark = ",")))
  cat(sprintf("%-12s %10s %10s %22s %6s %8s\n",
              "parameter", "truth", "post med", "95% CrI", "in?", "z"))
  res <- NULL
  for (p in names(PARS)) {
    if (!p %in% colnames(o$MCMC_output)) next
    d  <- o$MCMC_output[w, p]
    tv <- unname(x$truth[p])
    q  <- unname(quantile(d, c(0.025, 0.5, 0.975)))
    inside <- tv >= q[1] && tv <= q[3]
    z <- (mean(d) - tv) / sd(d)
    cat(sprintf("%-12s %10.4f %10.4f   (%8.4f,%8.4f) %6s %+8.2f\n",
                PARS[[p]], tv, q[2], q[1], q[3], if (inside) "yes" else "NO", z))
    res <- rbind(res, data.frame(scenario = scen[i], parameter = PARS[[p]],
                                 truth = tv, inside = inside, z = z))
  }
  ## The quantities a reader interprets, rather than the sampling-scale parameters.
  nat <- function(v) c(growth = 0.5 + exp(v[["a_bar"]]), decay = 0.25 + exp(v[["b_bar"]]),
                       peak = v[["vl_max_bar"]], onset_med = exp(v[["inc1"]]))
  tn <- nat(x$truth)
  dn <- t(apply(o$MCMC_output[w, c("a_bar", "b_bar", "vl_max_bar", "inc1"), drop = FALSE],
                1, function(r) nat(as.list(setNames(r, c("a_bar","b_bar","vl_max_bar","inc1"))))))
  cat(sprintf("\n%-12s %10s %10s %22s %6s\n", "derived", "truth", "post med", "95% CrI", "in?"))
  for (k in colnames(dn)) {
    q <- unname(quantile(dn[, k], c(0.025, 0.5, 0.975)))
    cat(sprintf("%-12s %10.3f %10.3f   (%8.3f,%8.3f) %6s\n", k, tn[[k]], q[2], q[1], q[3],
                if (tn[[k]] >= q[1] && tn[[k]] <= q[3]) "yes" else "NO"))
  }
  cat(sprintf("\n  %d of %d intervals contain the generating value; largest |z| = %.2f\n",
              sum(res$inside), nrow(res), max(abs(res$z))))
  all <- rbind(all, res)
  rm(o, x, d, dn); invisible(gc(verbose = FALSE))
}

if (!is.null(all) && length(unique(all$scenario)) > 1) {
  cat("\n\n== across all scenarios scored ==\n")
  cat(sprintf("%d of %d intervals contain the generating value (%.1f%%)\n",
              sum(all$inside), nrow(all), 100 * mean(all$inside)))
  miss <- all[!all$inside, ]
  if (nrow(miss)) {
    cat("\nmisses, by parameter:\n")
    print(table(miss$parameter))
    cat("\nmisses, by scenario:\n")
    print(table(miss$scenario))
  }
  cat(sprintf("\nlargest |z|: %.2f (%s, %s)\n", max(abs(all$z)),
              all$scenario[which.max(abs(all$z))], all$parameter[which.max(abs(all$z))]))
}
if (!is.null(all))
  cat("\nCoverage rates were not estimated: one chain per scenario gives a marginal check,\n",
      "not a replicate-based coverage probability.\n", sep = "")

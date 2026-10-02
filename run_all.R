#!/usr/bin/env Rscript
## ============================================================================
## run_all.R -- the entry point for everything in this repository
##
##   Rscript run_all.R check                  verify the environment and the data (minutes)
##   Rscript run_all.R demo                   one short fit, end to end (~20 min)
##   Rscript run_all.R simulations            all 24 reported fits (days; see below)
##   Rscript run_all.R simulations A0 A1       just these scenarios
##   Rscript run_all.R regenerate             rebuild the 24 synthetic datasets
##   Rscript run_all.R longitudinal           refit the priors (needs the cohort data)
##
## Run `check` first. It compiles the C++ kernels and confirms every synthetic dataset
## matches the generating parameters recorded alongside it, which is the whole of what this
## repository can verify without the confidential Pillar 2 data.
##
## WHAT CAN AND CANNOT BE REPRODUCED HERE
##
## The 24 simulation-recovery fits reproduce completely: the datasets are included, they
## carry the parameters that generated them, and the fitting recipe is fixed below.
##
## The four variant fits do NOT. They require the Pillar 2 testing data, which is UKHSA data
## held under agreement and is not ours to publish. The code that fits them is the same code
## used here -- `scripts/fit_model.R wt 1` differs from `scripts/fit_model.R A0 1` only in the
## name of the input array -- but the array itself must be requested from UKHSA.
##
## The longitudinal prior fits need the NBA and ATACCC cohort data, which are third-party
## participant data from the original studies. The two summaries the main model actually
## consumes -- a mean/sd vector and a covariance matrix -- are included, so every simulation
## fit runs without them.
## ============================================================================

STAGES <- c("check", "demo", "simulations", "score", "figures", "regenerate", "longitudinal")

## The reported simulation recipe. These are the settings every one of the 24 fits used;
## fit_model.R reads them from the environment so they are stated once, here.
##   A-C   less variable kinetics profile, onset ladder reading A -> B -> C
##   D-F   more variable kinetics profile, onset ladder reading D -> E -> F
##   0-3   roughly 1,000 / 10,000 / 100,000 / 1,000,000 positive tests
SCENARIOS <- paste0(rep(c("A", "B", "C", "D", "E", "F"), each = 4), 0:3)
RECIPE <- c(FIT_FORM        = "incidence",
            SYMP_DELAY_LIM  = "28",
            FIT_N_ITER      = "250000",
            FIT_COV_START   = "80000,100000,120000,140000",
            FIT_COV_WINDOW  = "20000")
## Everything before the final adaptation at 140,000 is warm-up. Posterior summaries use
## iterations 140,001 onward, which come from a single fixed proposal.
ANALYSIS_FROM <- 140001L

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(trimws(args[1])) else ""
if (!stage %in% STAGES) {
  cat("usage: Rscript run_all.R <stage> [args]\n\nstages:\n")
  cat("  check          verify the environment, the kernels and the synthetic data\n")
  cat("  demo           one short fit end to end, to prove the pipeline runs\n")
  cat("  simulations    the 24 reported simulation-recovery fits\n")
  cat("  score          did the fits recover the known parameters?\n")
  cat("  figures        the manuscript figures that can be built from this repository\n")
  cat("  regenerate     rebuild the synthetic datasets from their generating parameters\n")
  cat("  longitudinal   refit the viral kinetics priors (requires the cohort data)\n")
  quit(status = if (length(args)) 1L else 0L)
}
rest <- args[-1]

ok <- function(...) cat(sprintf("  [ok]   %s\n", sprintf(...)))
no <- function(...) cat(sprintf("  [FAIL] %s\n", sprintf(...)))
nb <- function(...) cat(sprintf("  [note] %s\n", sprintf(...)))
hdr <- function(s) cat(sprintf("\n== %s ==\n", s))
fails <- 0L
fail  <- function(...) { no(...); fails <<- fails + 1L }

## ---------------------------------------------------------------------------
## check
## ---------------------------------------------------------------------------
if (stage == "check") {
  hdr("environment")
  cat(sprintf("  %s on %s\n", R.version.string, R.version$platform))

  need <- c("Rcpp", "dqrng", "BH", "sitmo", "TruncatedNormal", "mvtnorm", "truncnorm")
  opt  <- c(rstan = "longitudinal fits", dplyr = "longitudinal fits",
            ggplot2 = "longitudinal fits", matrixStats = "longitudinal fits",
            posterior = "convergence diagnostics")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) fail("missing required packages: %s", paste(miss, collapse = ", "))
  else ok("all %d required packages present", length(need))
  for (p in names(opt))
    if (!requireNamespace(p, quietly = TRUE))
      nb("%s not installed -- needed only for %s", p, opt[[p]])

  ## MCMC_function compiles both kernels itself, at run time rather than on load, so they are
  ## compiled here explicitly -- otherwise a toolchain problem would not surface until hours
  ## into a fit.
  hdr("compiling the C++ kernels")
  t0 <- proc.time()
  cc <- try({
    Rcpp::sourceCpp("scripts/simulate_onset.cpp")
    Rcpp::sourceCpp("scripts/simulate_incidence.cpp")
  }, silent = TRUE)
  if (inherits(cc, "try-error"))
    fail("a kernel did not compile. OpenMP and a C++11 toolchain are required.\n%s",
         paste("        ", strsplit(as.character(cc), "\n")[[1]], collapse = "\n"))
  else {
    ok("both kernels compiled (%.0f s)", (proc.time() - t0)["elapsed"])
    for (f in c("ct_func2_cpp_thresh_tdist", "ct_func2_cpp_inc"))
      if (exists(f)) ok("%s() exported", f) else fail("%s() not exported by the kernels", f)
  }

  hdr("loading the model code")
  src <- try({
    source("scripts/model.R")
    source("scripts/sampler.R")
  }, silent = TRUE)
  if (inherits(src, "try-error"))
    fail("could not load the model code:\n%s",
         paste("        ", strsplit(as.character(src), "\n")[[1]], collapse = "\n"))
  else
    for (f in c("input_generator", "MCMC_function", "ct_cache_reset"))
      if (exists(f)) ok("%s() available", f) else fail("%s() not defined", f)

  hdr("synthetic datasets")
  miss <- SCENARIOS[!file.exists(sprintf("data/data_array_vl_%s_N.rds", SCENARIOS))]
  if (length(miss)) fail("missing datasets: %s", paste(miss, collapse = " "))
  else ok("all 24 reported datasets present")

  tot <- integer(0)
  for (s in setdiff(SCENARIOS, miss)) {
    z <- readRDS(sprintf("data/data_array_vl_%s_N.rds", s))
    if (is.null(z$out) || is.null(z$pars_test)) { fail("%s: not a list(out, pars_test)", s); next }
    if (!identical(dim(z$out), c(22L, 105L, 7L)))
      fail("%s: dimensions %s, expected 22x105x7", s, paste(dim(z$out), collapse = "x"))
    tot[s] <- round(sum(z$out))
    rm(z)
  }
  if (length(tot) == 24L) {
    ok("every dataset is 22x105x7 and carries its generating parameters")
    for (L in c("A", "B", "C", "D", "E", "F"))
      cat(sprintf("         %s  %s\n", L,
                  paste(sprintf("%s=%s", paste0(L, 0:3),
                                format(tot[paste0(L, 0:3)], big.mark = ",", trim = TRUE)),
                        collapse = "  ")))
  }

  ## The two kinetics profiles must be shared within A-C and within D-F, and must differ
  ## between them; the onset parameters must differ within each triplet. That is the design.
  KIN <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma")
  INC <- c("inc1", "inc2", "inc3")
  pars <- lapply(SCENARIOS, function(s) readRDS(sprintf("data/data_array_vl_%s_N.rds", s))$pars_test)
  names(pars) <- SCENARIOS
  same <- function(a, b, p) isTRUE(all.equal(unname(pars[[a]][p]), unname(pars[[b]][p]), tolerance = 1e-12))
  p1 <- all(vapply(c("A3", "B3", "C3"), same, logical(1), b = "A3", p = KIN))
  p2 <- all(vapply(c("D3", "E3", "F3"), same, logical(1), b = "D3", p = KIN))
  if (p1 && p2 && !same("A3", "D3", KIN)) ok("two distinct kinetics profiles: A-C and D-F")
  else fail("the kinetics profiles are not as the design requires")
  if (!same("A3", "B3", INC) && !same("B3", "C3", INC)) ok("onset varies across each triplet")
  else fail("onset parameters do not vary within a profile")

  hdr("truth snapshot")
  if (!file.exists("data/synthetic_truth.rds")) fail("data/synthetic_truth.rds missing")
  else {
    tr <- readRDS("data/synthetic_truth.rds")
    bad <- character(0)
    for (s in SCENARIOS) {
      t_s <- if (!is.null(tr[[s]])) tr[[s]] else NULL
      if (is.null(t_s)) { bad <- c(bad, sprintf("%s absent", s)); next }
      v <- if (is.list(t_s) && !is.null(t_s$pars_test)) t_s$pars_test else t_s
      sh <- intersect(names(v), names(pars[[s]]))
      if (!length(sh) || !isTRUE(all.equal(unname(v[sh]), unname(pars[[s]][sh]), tolerance = 1e-10)))
        bad <- c(bad, s)
    }
    if (length(bad)) fail("truth snapshot disagrees with the arrays for: %s",
                          paste(bad, collapse = " "))
    else ok("truth snapshot agrees with all 24 arrays")
  }

  hdr("priors consumed by the main model")
  for (f in c("longitudinal_fits/long_post_6.rds", "longitudinal_fits/long_cov_6.rds"))
    if (file.exists(f)) ok("%s (%d bytes)", f, file.info(f)$size) else fail("%s missing", f)

  hdr("model setup on one scenario")
  if (exists("input_generator")) {
    x <- try(input_generator(strain = "A0", data_dir = "data", n_iterations = 1000,
                            test_pop = 1e5, ncores = 1,
                            cov_start = as.integer(strsplit(RECIPE[["FIT_COV_START"]], ",")[[1]]),
                            form = RECIPE[["FIT_FORM"]],
                            symp_delay_lim = as.integer(RECIPE[["SYMP_DELAY_LIM"]])),
              silent = TRUE)
    if (inherits(x, "try-error"))
      fail("input_generator failed on A0:\n%s", as.character(x))
    else {
      ok("A0 set up: data_array %s, %d parameters, form %s, delay window %s",
         paste(dim(x$data_array), collapse = "x"), length(x$parameters), x$form,
         x$symp_delay_lim)
      if (!is.null(x$truth) && !all(is.na(x$truth)))
        ok("generating parameters surfaced as x$truth, so recovery can be scored")
      else fail("x$truth absent -- recovery could not be scored")
    }
  }

  hdr("longitudinal cohort data")
  ## These ship with the repository, so absence is a fault rather than a note. The columns are
  ## checked because the ATACCC file is published reduced to the four the model reads, and a
  ## fuller replacement dropped in by hand should still work.
  COH <- list("data/ct_data_clean.csv" = c("Person.ID", "Date.Index", "log10_GEperML"),
              "data/20210920_ATACCC_L.csv" = c("PersonID", "TestDateIndex", "CtT1", "VacVOC"))
  for (f in names(COH)) {
    if (!file.exists(f)) { fail("%s missing -- it is part of this repository", f); next }
    h <- try(names(read.csv(f, nrows = 1)), silent = TRUE)
    if (inherits(h, "try-error")) { fail("%s could not be read", f); next }
    gap <- setdiff(COH[[f]], h)
    if (length(gap)) fail("%s lacks required column(s): %s", f, paste(gap, collapse = ", "))
    else ok("%s: %d rows, required columns present", f,
            nrow(read.csv(f, colClasses = "character")))
  }

  hdr("optional inputs")
  for (s in c("wt", "alpha", "delta", "omicron")) {
    f <- sprintf("data/data_array_vl_%s_N.rds", s)
    if (file.exists(f)) nb("%s present -- the variant fits can be run", f)
  }

  cat(sprintf("\n%s\n", if (fails == 0L)
    "check passed. `Rscript run_all.R demo` next." else
    sprintf("*** %d check(s) failed ***", fails)))
  quit(status = if (fails == 0L) 0L else 1L)
}

## ---------------------------------------------------------------------------
## helpers for the stages that shell out to fit_model.R
## ---------------------------------------------------------------------------
## The settings are passed by setting them in THIS process and letting the child inherit them.
## system2(env = ) is ignored on Windows, so passing them that way silently ran fit_model.R at
## its 250,000-iteration default instead of the requested schedule.
run_fit <- function(target, chain, ncores, env) {
  old <- Sys.getenv(names(env), unset = NA_character_, names = TRUE)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    set <- old[!is.na(old)]
    if (length(set)) do.call(Sys.setenv, as.list(set))
    if (any(is.na(old))) Sys.unsetenv(names(old)[is.na(old)])
  }, add = TRUE)

  cat(sprintf("\n--- %s chain %d ---\n%s\nRscript scripts/fit_model.R %s %d %d\n",
              target, chain, paste(sprintf("%s=%s", names(env), unlist(env)), collapse = " "),
              target, chain, ncores))
  st <- system2("Rscript", c("scripts/fit_model.R", target, chain, ncores))
  if (!identical(st, 0L)) cat(sprintf("*** %s chain %d exited %s ***\n", target, chain, st))
  st
}

## ---------------------------------------------------------------------------
## demo -- the smallest thing that proves the pipeline works
## ---------------------------------------------------------------------------
if (stage == "demo") {
  target <- if (length(rest)) toupper(rest[1]) else "A0"
  dir.create("outputs", showWarnings = FALSE)
  cat(sprintf("Short fit on %s to prove the pipeline runs end to end.\n", target))
  cat("This is NOT a reproduction: it uses 5,000 iterations against the reported 250,000,\n")
  cat("and test_pop 1e5 against 1e6, so the posterior is not the paper's. Output is written\n")
  cat(sprintf("to outputs/%s_i1_5000.rds (chains checkpoint every 5,000 iterations).\n", target))
  env <- list(FIT_FORM = RECIPE[["FIT_FORM"]], SYMP_DELAY_LIM = RECIPE[["SYMP_DELAY_LIM"]],
              FIT_N_ITER = "5000", FIT_COV_START = "2000,3500", FIT_COV_WINDOW = "1500",
              FIT_TEST_POP = "1e5")
  st <- run_fit(target, 1L, 1L, env)
  f <- sprintf("outputs/%s_i1_5000.rds", target)
  if (file.exists(f)) {
    o <- readRDS(f)
    n <- max(which(!is.na(o$MCMC_posteriors[, 1])))
    cat(sprintf("\nWrote %s: %d draws, final log posterior %.1f.\n", f, n,
                o$MCMC_posteriors[n, 1]))
    cat("Pipeline works. For the reported fits: Rscript run_all.R simulations\n")
  } else {
    cat(sprintf("\n*** no output at %s -- the demo did not complete ***\n", f))
    quit(status = 1L)
  }
  quit(status = if (identical(st, 0L)) 0L else 1L)
}

## ---------------------------------------------------------------------------
## simulations -- the 24 reported fits
## ---------------------------------------------------------------------------
if (stage == "simulations") {
  want <- if (length(rest)) toupper(rest) else SCENARIOS
  bad  <- setdiff(want, SCENARIOS)
  if (length(bad)) stop("unknown scenario: ", paste(bad, collapse = ", "),
                        "\nvalid: ", paste(SCENARIOS, collapse = " "))
  ncores <- as.integer(Sys.getenv("RUN_NCORES", "1"))
  dir.create("outputs", showWarnings = FALSE)
  cat(sprintf("Fitting %d scenario(s) at the reported settings, sequentially.\n", length(want)))
  cat("Each fit is 250,000 iterations and takes several hours on one core, so the full set\n")
  cat("is days of compute in this form. The reported fits were run as independent processes\n")
  cat("in parallel; to do the same, launch one of these per core instead:\n")
  cat(sprintf("  %s Rscript scripts/fit_model.R <scenario> 1 1\n",
              paste(sprintf("%s=%s", names(RECIPE), RECIPE), collapse = " ")))
  cat(sprintf("\nAnalysis window: iterations %s onward.\n\n",
              format(ANALYSIS_FROM, big.mark = ",")))
  env <- as.list(RECIPE)
  bad <- character(0)
  for (s in want) if (!identical(run_fit(s, 1L, ncores, env), 0L)) bad <- c(bad, s)
  cat(sprintf("\n%d of %d completed.%s\n", length(want) - length(bad), length(want),
              if (length(bad)) sprintf(" Failed: %s", paste(bad, collapse = " ")) else ""))
  quit(status = if (length(bad)) 1L else 0L)
}

## ---------------------------------------------------------------------------
## score -- recovery against the known generating parameters
## ---------------------------------------------------------------------------
if (stage == "score") {
  st <- system2("Rscript", c("scripts/score_recovery.R", "outputs", rest))
  quit(status = if (identical(st, 0L)) 0L else 1L)
}

## ---------------------------------------------------------------------------
## figures -- the manuscript figures this repository can build
## ---------------------------------------------------------------------------
if (stage == "figures") {
  cat("Building the manuscript figures that come from the simulation study.\n")
  cat("Three of the paper's figures are hand-drawn schematics, and eleven describe the four\n")
  cat("variants and need the Pillar 2 data. Their code is in scripts/paper_outputs.R and runs\n")
  cat("unchanged for anyone who has that data -- see the header of scripts/make_figures.R.\n")
  st <- system2("Rscript", c("scripts/make_figures.R", rest))
  quit(status = if (identical(st, 0L)) 0L else 1L)
}

## ---------------------------------------------------------------------------
## regenerate -- rebuild the synthetic datasets
## ---------------------------------------------------------------------------
if (stage == "regenerate") {
  fam <- if (length(rest)) toupper(rest) else c("A", "B", "C", "D", "E", "F")
  bad <- setdiff(fam, c("A", "B", "C", "D", "E", "F"))
  if (length(bad)) stop("unknown family: ", paste(bad, collapse = ", "))
  tg <- paste0(rep(fam, each = 4), 0:3)
  fo <- sprintf("data/data_array_vl_%s_N.rds", tg)

  cat("Rebuilding the synthetic datasets from their generating parameter vectors, then\n")
  cat("comparing the result against the shipped arrays byte for byte.\n\n")
  cat("make_scenarios_incidence.R refuses to overwrite an existing array -- regenerating one\n")
  cat("while chains are being fitted to the old copy would leave those chains fitted to data\n")
  cat("nobody has a record of. So the shipped arrays are moved aside first, into\n")
  cat("data/_shipped/, and restored if anything goes wrong.\n\n")

  bk <- "data/_shipped"
  dir.create(bk, showWarnings = FALSE)
  have <- fo[file.exists(fo)]
  md5_before <- tools::md5sum(have)
  ok_move <- file.rename(have, file.path(bk, basename(have)))
  if (!all(ok_move)) {
    cat("*** could not move the shipped arrays aside; nothing has been changed ***\n")
    quit(status = 1L)
  }
  restore <- function() {
    f <- file.path(bk, basename(have))
    file.rename(f[file.exists(f)], have[file.exists(f)])
  }

  st <- system2("Rscript", c("scripts/make_scenarios_incidence.R", fam))
  if (!identical(st, 0L)) {
    cat("\n*** regeneration failed; restoring the shipped arrays ***\n")
    unlink(fo[file.exists(fo)]); restore()
    quit(status = 1L)
  }

  hdr("comparing regenerated against shipped")
  md5_after <- tools::md5sum(have)
  same <- !is.na(md5_after) & md5_after == md5_before[names(md5_after)]
  for (i in seq_along(have))
    cat(sprintf("  %-34s %s\n", basename(have[i]),
                if (isTRUE(same[i])) "identical" else "DIFFERS"))
  nsame <- sum(same, na.rm = TRUE)
  cat(sprintf("\n%d of %d reproduced exactly.\n", nsame, length(have)))
  if (nsame == length(have)) {
    unlink(bk, recursive = TRUE)
    cat("Generation is deterministic; data/_shipped/ removed.\n")
  } else {
    cat(sprintf("The shipped copies are preserved in %s/ for comparison. Regeneration is NOT\n", bk))
    cat("bit-reproducible on this platform -- the arrays differ in their random draws, which\n")
    cat("depend on the RNG stream of the parallel simulation. The generating parameters are\n")
    cat("unaffected, so recovery remains meaningful, but the fits in the paper belong to the\n")
    cat("shipped arrays. Restore them before fitting.\n")
  }

  cat("\nRefreshing the truth snapshot.\n")
  st <- system2("Rscript", c("scripts/synthetic_truth.R"))
  cat(sprintf("\n%s\n", if (identical(st, 0L))
    "Done. Run `Rscript run_all.R check` to confirm the arrays still verify." else
    "*** synthetic_truth.R failed ***"))
  quit(status = if (identical(st, 0L) && nsame == length(have)) 0L else 1L)
}

## ---------------------------------------------------------------------------
## longitudinal -- refit the priors
## ---------------------------------------------------------------------------
if (stage == "longitudinal") {
  miss <- c("data/ct_data_clean.csv", "data/20210920_ATACCC_L.csv")
  miss <- miss[!file.exists(miss)]
  if (length(miss)) {
    cat("The cohort data ships with this repository, so these files should be present.\n")
    cat(sprintf("Missing:\n%s\n\n", paste("  ", miss, collapse = "\n")))
    cat("See the 'Longitudinal priors' section of the README for the columns each file must\n")
    cat("have. The two summaries the main model consumes are committed, so every simulation\n")
    cat("fit runs without this stage.\n")
    quit(status = 1L)
  }
  cat("Refitting the longitudinal model in Stan. Writes long_post_6.rds / long_cov_6.rds\n")
  cat("(pooled NBA + ATACCC, the reported prior) and long_post_ATTAC.rds / long_cov_ATTAC.rds\n")
  cat("(ATACCC only, the cross-platform sensitivity check).\n\n")
  source("scripts/longitudinal_fit.R")
  cat("\nDone.\n")
  quit(status = 0L)
}

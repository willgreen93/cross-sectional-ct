## ============================================================================
## make_scenarios_incidence.R -- the ladder scenarios regenerated under incidence indexing
##
##   Rscript scripts/make_scenarios_incidence.R            # all 24
##   Rscript scripts/make_scenarios_incidence.R V W X      # one or more destination families
##
## Writes data/data_array_vl_{S,T,U,V,W,X}{0,1,2,3}_N.rds: the M/N/O/P/Q/R scenarios
## regenerated with form = "incidence", so incidence is indexed on crossing vl_thresh and the
## delay to symptom onset is convolved against calendar time.
##
##   S <- M   T <- N   U <- O    out-of-prior arm, profile 1 (more variable)
##   V <- P   W <- Q   X <- R    in-prior arm, profile 2 (less variable)
##
## THE GENERATING PARAMETERS ARE NOT RE-CHOSEN. Each source array stores the vector that
## produced it as $pars_test, and that identical vector -- kinetics, onset, testing, rdisp,
## thetas, multiplier -- is pushed through the incidence kernel here. The two studies
## therefore differ in exactly one thing, the indexing, and any difference in recovery is
## attributable to it. It also means V/W/X inherit P/Q/R's in-prior property unchanged, since
## the values have not moved, and that the onset ladder still reads down a column.
##
## What does change is the expected number of tests: the same incidence curve run through a
## different observation model gives somewhat different totals. The thetas are left alone
## rather than recalibrated onto round numbers, because holding the parameters fixed is the
## point; the achieved totals are printed and are what the analysis should quote.
##
## The negative binomial is drawn at WEEKLY resolution, matching equations 10-11. Drawing
## daily and summing gives a variance about sevenfold too small and rdisp then recovers about
## sevenfold under truth -- the same arithmetic as the Figure 3 interval bug (U1).
##
## Generation uses test_pop = 1e7 against the fits' 1e6, so the Monte Carlo error in the data
## is independent of the Monte Carlo error in the likelihood. Sharing it would put the same
## bias on both sides of the comparison and recovery would look good for the wrong reason.
##
## Re-run scripts/synthetic_truth.R afterwards to refresh the truth snapshot.
## ============================================================================
source("scripts/model.R")
Rcpp::sourceCpp("scripts/simulate_onset.cpp")
Rcpp::sourceCpp("scripts/simulate_incidence.cpp")

## TWO GENERATIONS, both derived from the same M-R generating vectors.
##
##   S-X   20-day delay window. The first incidence-indexed set, kept as the archive.
##   A-F   28-day delay window, matching the reported model. The letters are chosen so that
##         the family IS its figure panel: A-C are the less-variable profile (from P, Q, R),
##         D-F the more-variable one (from M, N, O), each reading down the onset ladder.
##
## Only the window differs between them. The priors changed too, but priors play no part in
## generation -- a fixed parameter vector is pushed through the kernel -- so they affect the
## fits alone.
MAP      <- c(S = "M", T = "N", U = "O", V = "P", W = "Q", X = "R",   # 20-day archive
              A = "P", B = "Q", C = "R", D = "M", E = "N", F = "O")   # 28-day, reported
SDL_OF   <- function(dest) if (dest %in% c("A","B","C","D","E","F")) 28L else 20L
RUNGS    <- 0:3
GEN_POP  <- 1e7
GEN_CORES <- 4L
KNOTS    <- c(seq(0, 99, 9)[-1]) / 108      # the ladder's knot rule, as input_generator sets it
POPULATION <- 5.6e7
## Seeded per scenario rather than once at the top, so a family regenerates identically
## whether it is built alone or alongside others. The base differs from the one
## make_scenarios_pqr.R used, so these arrays do not inherit the M-R noise draws, and differs
## again between the two generations so A-F does not inherit S-X's.
SEED_OF  <- function(dest, rung)
  (if (dest %in% c("A","B","C","D","E","F")) 20260925L else 20260917L) +
  100L * match(dest, names(MAP)) + rung

args <- commandArgs(trailingOnly = TRUE)
dest_families <- if (length(args)) toupper(trimws(args)) else names(MAP)
bad <- setdiff(dest_families, names(MAP))
if (length(bad)) stop("unknown destination family: ", paste(bad, collapse = ", "))

ct_cache_reset(GEN_CORES)

for (dest in dest_families) {
  src <- MAP[[dest]]
  SDL <- SDL_OF(dest)
  for (rung in RUNGS) {
    set.seed(SEED_OF(dest, rung))
    f_in  <- sprintf("data/data_array_vl_%s%d_N.rds", src, rung)
    f_out <- sprintf("data/data_array_vl_%s%d_N.rds", dest, rung)
    if (!file.exists(f_in))  stop("no source array ", f_in)
    if (file.exists(f_out))  stop("refusing to overwrite ", f_out)

    z <- readRDS(f_in)
    p <- z$pars_test
    if (is.null(p)) stop(f_in, " carries no generating values")
    ## Never proposed by the sampler; a non-zero l1 makes the kernel draw an extra random
    ## number per individual and shifts the whole stream.
    p[c("l_bar", "l_sigma")] <- 0
    max_day <- dim(z$out)[2]

    ## Expected counts, without the leading undetected class the likelihood drops.
    template <- array(dim = dim(z$out) - c(1L, 0L, 0L),
                      dimnames = lapply(seq_along(dim(z$out)),
                                        function(i) if (i == 1) dimnames(z$out)[[1]][-1] else dimnames(z$out)[[i]]))
    mu <- p_array_func(parameters = p, knots = KNOTS, vk = 1, max_day = max_day,
                       population = POPULATION, data_array = template, test_pop = GEN_POP,
                       ncores = GEN_CORES, form = "incidence", symp_delay_lim = SDL,
                       stoch = 0.5, name = FALSE)[[1]]

    ## Two rows are added back, not one: input_generator strips a leading "vl_0.5" row, and
    ## p_array_func itself drops the lowest remaining bin as the undetected class, so the
    ## expected array is two rows shorter than the stored one. Both carry no observations.
    full <- array(0, dim = dim(z$out), dimnames = dimnames(z$out))
    full[-c(1, 2), , ] <- mu
    full[is.na(full)] <- 0

    wk  <- weekly_aggregator(full)
    obs <- array(rnbinom(length(wk), size = 1 / p[["rdisp"]], mu = wk),
                 dim = dim(wk), dimnames = dimnames(wk))
    out <- array(0L, dim = dim(full), dimnames = dimnames(full))
    out[, seq(1, max_day, 7), ] <- obs

    saveRDS(list(pars_test = p, out = out), f_out)
    cat(sprintf("%s%d (from %s%d, %d-day window): %10s tests, expected %10s | onset %4.1f d | %d days\n",
                dest, rung, src, rung, SDL, format(sum(out), big.mark = ","),
                format(round(sum(full)), big.mark = ","), exp(p[["inc1"]]), max_day))
    rm(z, mu, full, wk, obs, out); invisible(gc(verbose = FALSE))
  }
}

cat("\nrun scripts/synthetic_truth.R next, to refresh data/synthetic_truth.rds\n")

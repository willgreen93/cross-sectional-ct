## ============================================================================
## fit_model.R -- the single fitting script
##
##   Rscript scripts/fit_model.R <target> <chain> [ncores]
##
##   Rscript scripts/fit_model.R wt 1 4          # a variant
##   Rscript scripts/fit_model.R L3 2 2          # a simulation scenario
##
## `target` is one of wt, alpha, delta, omicron, or the name of a synthetic dataset in
## data/synthetic. `chain` is an integer that seeds the starting draw and names the output:
## <target>_<chain>_<iteration>.rds.
##
## SAMPLING SCHEDULE  (one protocol, every fit)
##
##   1       - 100,000   componentwise Metropolis, one parameter per iteration
##   100,000             covariance estimated from iterations 50,000-100,000
##   100,000 - 150,000   block multivariate normal proposal, scaled 2.38^2/d
##   150,000             covariance re-estimated from iterations 100,000-150,000
##   150,000 - 250,000   block proposal, kernel now FIXED -- this is the reported sample
##
## Everything before the final adaptation is warm-up and is discarded. Draws after it come
## from a single fixed kernel, which is what makes them a valid sample and what lets the
## Methods describe the sampler in three sentences.
##
## Both adaptation windows sit entirely after burn-in, which across all four variants has
## been between 5,000 and 33,000 iterations. That matters: a covariance estimated over
## pre-stationary draws is inflated along the directions the chain travelled while
## converging, and the proposal then over-steps in exactly those directions. Alpha ran
## 100,000 iterations at 6.5% acceptance and ESS 5 on a_bar from that mistake.
##
## STARTING POINTS
##
## Self-contained -- nothing is inherited from a previous fit, so a chain can be reproduced
## from this repository alone. Each chain is dispersed from the others, which is what makes
## agreement between them evidence rather than an artefact of a shared start.
## ============================================================================

## ---- settings ---------------------------------------------------------------
## Overridable via environment variables so one-off schedules (the simulation study has run
## with adaptation at 80,000/20,000 rather than 100,000/50,000) don't require editing this
## file -- and so every other target keeps the documented default without having to remember
## to change it back. Set FIT_N_ITER / FIT_COV_START / FIT_COV_WINDOW before launching:
##   $env:FIT_COV_START = "80000,120000"; $env:FIT_COV_WINDOW = "20000"
N_ITER     <- as.integer(Sys.getenv("FIT_N_ITER", "250000"))
COV_START  <- as.integer(strsplit(Sys.getenv("FIT_COV_START", "100000,150000"), ",")[[1]])
COV_WINDOW <- as.integer(Sys.getenv("FIT_COV_WINDOW", "50000"))
## 1e6 is the value every reported fit used; FIT_TEST_POP exists so run_all.R's smoke test can
## drop it for speed. Lowering it raises the Monte Carlo error in the likelihood, so a run with
## a different value is a pipeline check, not a reproduction.
TEST_POP   <- as.numeric(Sys.getenv("FIT_TEST_POP", "1e6"))
VARIANTS   <- c("wt", "alpha", "delta", "omicron")
## The synthetic arrays live alongside the real ones as data/data_array_vl_<scenario>_N.rds.
SYNTH_DIR  <- "data"

## The generating parameters come from the array file itself -- each synthetic array is saved
## as a list with the parameters that produced it in $pars_test, and input_generator surfaces
## them as x$truth. They are used for one thing only: pinning the gauge coefficient at its
## true value, so recovery can be checked without re-gauging.
##
## NOT from data/synthetic_truth.rds, which re-derives them by evaluating the definitions in
## simulation_study.R. That is a second record of the same thing and it has already drifted
## once: editing pars_M3 and pars_O3 after the arrays were generated left the file three
## onset parameters out of date, with nothing to announce it. The array cannot drift from
## itself.

THETA_SD  <- 1     # spline coefficients: N(prior mean, THETA_SD)
MULT_SD   <- 0.1   # multiplier: N(1, MULT_SD)
## Starts are drawn to SPAN the known modes rather than to sit near either: the slow-growth mode
## has an onset near 15 days and the fast one near 4, and a chain started under both (as the
## longitudinal prior's growth of ~0.83 is) slides into whichever is nearer. The mode a chain
## reaches is then evidence rather than an artefact of where it began.
## The symptom-timing parameters start from their priors (see the draw below); only growth and
## its spread are drawn to straddle the modes, since that is where the two solutions separate.

## ---- arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: Rscript scripts/fit_model.R <target> <chain> [ncores]")

target <- trimws(args[1])
chain  <- as.integer(args[2])
ncores <- if (length(args) >= 3) as.integer(args[3]) else 1
if (is.na(chain)  || chain  < 1) stop("chain must be a positive integer")
if (is.na(ncores) || ncores < 1) stop("ncores must be a positive integer")

## form selects which model is fitted. "thresh_tdist" is the one the paper reports and stays
## the default, so every existing invocation is unaffected. "incidence" indexes incidence on
## crossing detectability and convolves a skew-normal delay to symptom onset, which makes the
## viral load composition at a fixed test delay depend on the epidemic -- the mechanism the
## reported model cannot represent.
FORM <- if (length(args) >= 4) trimws(args[4]) else Sys.getenv("FIT_FORM", "thresh_tdist")
if (!FORM %in% c("thresh_tdist", "incidence"))
  stop("form must be thresh_tdist or incidence, got: ", FORM)
SYMP_DELAY_LIM <- as.integer(Sys.getenv("SYMP_DELAY_LIM", "20"))

is_variant <- tolower(target) %in% VARIANTS
if (is_variant) target <- tolower(target)

## Distinct per chain and reproducible from the arguments alone.
seed <- sum(utf8ToInt(target)) * 1000L + chain

## ---- setup ------------------------------------------------------------------
source("scripts/model.R")
source("scripts/sampler.R")
set.seed(seed)

## One call for both: input_generator already defaults data_dir to "data", where the real
## and synthetic arrays sit side by side, and branches internally on whether the strain is
## one of the four variants. (An `x <- if (...) ... else ...` split across lines parses as
## complete at the first branch and errors on the `else`.)
x <- input_generator(strain = target, data_dir = SYNTH_DIR, n_iterations = N_ITER,
                     test_pop = TEST_POP, ncores = ncores, cov_start = COV_START,
                     form = FORM, symp_delay_lim = SYMP_DELAY_LIM)

x$cov_window    <- COV_WINDOW
x$proposal_type <- "MH"        # block proposal takes over at the first adaptation
## Incidence-form chains are tagged i1, i2, ... so they cannot overwrite the reported
## thresh_tdist chains, which are tagged 1, 2, 3 and written to the same directory.
x$gene          <- if (FORM == "incidence") paste0("i", chain) else as.character(chain)

## Where the analysis window opens: the draw after the final covariance adaptation, so every
## retained draw comes from one fixed proposal. Recorded on the chain itself because everything
## downstream -- score_recovery.R, and the recovery and precision figures -- needs it, and the
## alternative is for each of them to re-derive it. Deriving it from cov_start is only correct
## for a chain that ran the whole schedule from iteration 1, which is exactly what this script
## produces; a continuation has a different history and carries the value written by whichever
## script joined it.
x$analysis_from <- as.integer(max(COV_START) + 1L)

## ---- starting point -----------------------------------------------------------
nmv <- names(x$parameters)
vl  <- 1:6
inc <- match(c("inc1", "inc2", "inc3"), nmv)
th  <- grep("^theta", nmv)
pin <- x$ignored_spline_par
lb  <- x$lower_bound
lnm <- intersect(names(lb), nmv)

## Reject rather than clip. A symptom threshold above the peak leaves nobody detectable; a
## collapsed inc2 starts the chain inside the degenerate solution where the onset
## distribution has no width; and an onset outside 1-20 days is not a state worth starting
## a chain from whatever the priors permit.
## inc1 and inc2 mean the same thing in both forms, so the delay checks are shared. The
## difference is inc3: a symptom threshold under thresh_tdist, which must sit below the peak or
## nobody is detectable; unused and pinned under the incidence form, where the delay must
## instead fit inside the convolution window.
## Braced deliberately: a top-level `x <- if (...) ... else ...` split across lines parses as
## complete at the first branch and errors on the else.
viable <- {
  if (FORM == "incidence")
    ## Same checks, plus one: the delay is convolved over a finite window, so a median onset
    ## near symp_delay_lim would have most of its mass truncated away.
    function(v) v[["inc3"]] < v[["vl_max_bar"]] - 1 && v[["inc2"]] > 0.05 &&
                exp(v[["inc1"]]) > 1 && exp(v[["inc1"]]) < SYMP_DELAY_LIM - 4
  else
    function(v) v[["inc3"]] < v[["vl_max_bar"]] - 1 && v[["inc2"]] > 0.05 &&
                exp(v[["inc1"]]) > 1 && exp(v[["inc1"]]) < 20
}

## ---- the incidence scale, computed once ------------------------------------------
## None of this depends on the draw, so it sits outside the retry loop -- an apply() over
## the data array on all 500 attempts would be pure waste.
daily    <- apply(x$data_array[-1, , , drop = FALSE], 2, sum, na.rm = TRUE)
knot_day <- pmax(round(c(1, x$knots * x$max_day)), 1)

## Each knot takes the count from the nearest day that actually has data, rather than from
## its own day floored at 0.5. Knots falling in the spline lead-in -- the model only starts
## filling at start_j = 35 -- would otherwise all return the floor, and log(0.5 / 5.6e7) =
## -18.5 for every scenario regardless of its incidence. That made the pinned coefficient
## identical across rungs differing by three orders of magnitude, the opposite of what this
## initialisation is for.
pos <- which(daily > 0)
if (!length(pos)) stop("no observed counts in ", target, "'s data array")
obs <- daily[vapply(pmin(knot_day, length(daily)),
                    function(d) pos[which.min(abs(pos - d))], integer(1))]

## ---- where to pin the gauge ------------------------------------------------------
## One spline coefficient is held at its starting value for the whole run, because that is
## what breaks the degeneracy between the thetas and the multiplier: incidence is
## population * exp(spline(theta)) * multiplier, so adding a constant to every theta and
## dividing the multiplier by exp(it) leaves it unchanged.
##
## For a SYNTHETIC scenario it is pinned at the generating value. The gauge is not
## estimated from the data -- it is a convention -- so fixing it at truth adds no
## information about the epidemic or the kinetics; it only chooses which point on the
## orbit gets reported. What it buys is a direct recovery check: the fitted thetas and the
## multiplier can be compared to truth without first re-gauging them, and the multiplier
## should land near 1 rather than absorbing an arbitrary constant.
##
## For the REAL variants there is no truth, so the value comes from the data. That keeps
## this script usable on both without ever reaching for an answer it should not have.
pin_name <- names(x$parameters)[pin]
pin_val  <- log(obs[match(pin, th)] / x$population)
pin_src  <- "data"
if (!is_variant) {
  tr <- x$truth
  if (!is.null(tr) && !all(is.na(tr)) && pin_name %in% names(tr)) {
    pin_val <- tr[[pin_name]]
    pin_src <- "truth"
  } else {
    warning("no generating value for ", target, "'s ", pin_name,
            " in the array's own record; pinning from the data instead", call. = FALSE)
  }
}

## FIT_INIT_KIN_DIR / FIT_INIT_KIN_STRAIN: take the kinetics and the onset delay from ANOTHER
## variant's fitted posterior, leaving the spline, the multiplier and rdisp to the usual draw.
##
## For a target whose own chains are suspected of sharing a mode, the choice of starting region
## is the whole experiment, and the two obvious candidates are both compromised: the
## longitudinal prior sits at growth ~0.83, below both Alpha modes, and Alpha's own posterior
## begs the question. Another variant's posterior is neither -- it is a plausible region for
## SARS-CoV-2 kinetics chosen without reference to this variant's likelihood.
##
## Drawn jointly from one posterior sample rather than marginally, so the correlations between
## growth, decay, peak and their between-individual spreads are those of a real solution. Only
## the kinetics and eta_T/omega_T transfer: the epidemic curve belongs to this variant.
KIN_DIR    <- Sys.getenv("FIT_INIT_KIN_DIR", "")
KIN_STRAIN <- Sys.getenv("FIT_INIT_KIN_STRAIN", "")
kin_draws  <- NULL
if (nzchar(KIN_DIR)) {
  if (!nzchar(KIN_STRAIN)) stop("FIT_INIT_KIN_DIR needs FIT_INIT_KIN_STRAIN")
  kf <- list.files(KIN_DIR, sprintf("^%s_[^_]+(_[0-9]+)?\\.rds$", KIN_STRAIN), full.names = TRUE)
  if (!length(kf)) stop("no chains for ", KIN_STRAIN, " in ", KIN_DIR)
  ks  <- readRDS(sort(kf)[1])
  kfl <- max(which(!is.na(ks$MCMC_posteriors)))
  kb  <- if (!is.null(ks$x$analysis_from)) as.integer(ks$x$analysis_from) else floor(kfl / 2)
  kin_draws <- ks$MCMC_output[kb:kfl, intersect(c(names(x$parameters)[vl], "inc1", "inc2"),
                                                colnames(ks$MCMC_output)), drop = FALSE]
  rm(ks); invisible(gc(verbose = FALSE))
  cat(sprintf("kinetics start drawn from %s in %s (%s draws)\n", KIN_STRAIN, KIN_DIR,
              format(nrow(kin_draws), big.mark = ",")))
}

for (att in 1:500) {
  p <- x$parameters                                    # input_generator's own prior draw

  ## Viral kinetics from the longitudinal prior -- the only place in the start that uses
  ## external information, and it is the same prior the model itself uses.
  p[vl] <- rnorm(6, x$prior_mean[vl], x$prior_sds[vl])

  ## The symptom-timing parameters likewise start from their own priors, which the tightened
  ## priors now make sensible -- the old N(1, 1) on inc1 admitted delays of months, so the
  ## start had to be hand-set. Drawing from the prior means the start tracks the prior if it
  ## changes and cannot quietly favour a region: the previous draw pinned vl_thresh at 0.10,
  ## within a hair of where it is always estimated, so a threshold found at the bottom of its
  ## range was not evidence of anything. inc3 is drawn on the natural scale and converted
  ## after the viability check below.
  ## vl_thresh is drawn from its prior RESTRICTED TO BELOW THE LIMIT OF DETECTION, unlike the
  ## other two. This is a statement about the sampler, not about the posterior. The conditional
  ## posterior for vl_thresh has very steep walls once the other parameters have adapted around
  ## it -- from a settled state, moving it alone by 0.5 costs hundreds of log units, in either
  ## direction -- and escaping requires vl_thresh, growth and the onset delay to move together,
  ## which componentwise sampling cannot do. A chain started above the detection limit therefore
  ## stays there for the whole run rather than exploring: on 2026-09-23, eight of twelve chains
  ## started from the unrestricted prior were still locked at vl_thresh 3.5-4.8 at 60,000
  ## iterations, roughly 240 log units below the solution wild type had already found.
  ##
  ## Starting below the detection limit is where the biology points -- innate sensing is
  ## cell-intrinsic and begins when virus is unmeasurable in a swab -- and where every fit to
  ## date has settled. The restriction applies ONLY to the starting draw; the prior is
  ## untouched, and the chain is free to move above the limit if the data support it.
  VT_START_MAX <- as.numeric(sub("^vl_", "", rownames(x$data_array)[2]))   # the limit of detection
  p[inc] <- c(rnorm(1, x$prior_mean[["inc1"]], x$prior_sds[["inc1"]]),
              truncnorm::rtruncnorm(1, a = 0, mean = x$prior_mean[["inc2"]], sd = x$prior_sds[["inc2"]]),
              truncnorm::rtruncnorm(1, a = 0, b = VT_START_MAX,
                                    mean = x$prior_mean[["inc3"]], sd = x$prior_sds[["inc3"]]))

  ## Optional: confine the random start to the region that avoids the secondary mode.
  ## U5 found that every fresh Alpha chain starting with growth above 1.0 reached the
  ## better-supported mode (9 of 9), and 15 of 17 starting with onset below 5.5 days did, while
  ## draws from the longitudinal prior sit at growth ~0.83 -- below that line. Starts stay
  ## random (uniform within the bounds, so chains are still dispersed) but are drawn only
  ## where U5 says the good mode is reached.
  ##   FIT_INIT_GROWTH_MIN  lower bound on growth 0.5 + exp(a_bar); a_bar drawn over +1 on
  ##                        the log scale above the implied minimum
  ##   FIT_INIT_ONSET_MAX   upper bound on median onset exp(inc1) in days; inc1 drawn between
  ##                        3 days and this bound
  ## One posterior draw, replacing both blocks above. Applied before the confinement options
  ## below, so those still win if they are also set.
  if (!is.null(kin_draws)) {
    kd <- kin_draws[sample.int(nrow(kin_draws), 1L), ]
    p[names(kd)] <- kd
  }

  if (nzchar(Sys.getenv("FIT_INIT_GROWTH_MIN"))) {
    a_min <- log(as.numeric(Sys.getenv("FIT_INIT_GROWTH_MIN")) - 0.5)
    p[["a_bar"]] <- runif(1, a_min, a_min + 1)
  }
  if (nzchar(Sys.getenv("FIT_INIT_ONSET_MAX"))) {
    p[["inc1"]] <- runif(1, log(3), log(as.numeric(Sys.getenv("FIT_INIT_ONSET_MAX"))))
  }

  p[th] <- log(obs / x$population) + rnorm(length(th), 0, THETA_SD)
  p[["multiplier"]] <- rnorm(1, 1, MULT_SD)

  ## The pinned coefficient, taken from the data WITHOUT jitter.
  ##
  ## Incidence is population * exp(spline(theta)) * multiplier, so adding a constant to
  ## every theta and dividing the multiplier by exp(it) leaves it unchanged. One
  ## coefficient is held fixed to break that degeneracy, and it must be held at the same
  ## value in every chain of a target -- otherwise chains sit on different points of one
  ## orbit and disagree about parameters while describing an identical epidemic. A
  ## data-derived value satisfies that: it depends only on the array, which every chain
  ## shares, so it is identical across chains while also being in the right place.
  ##
  ## This was previously the prior mean of -9, which suits the four variants but is wrong
  ## by seven log units for the low-incidence synthetic rungs -- and being pinned it never
  ## recovers, so the multiplier compensates for the whole run.
  p[pin]            <- pin_val
  ## Not proposed, and a non-zero l1 makes the kernel draw an extra random number per
  ## simulated individual, shifting the whole stream.
  p[x$ignored_pars] <- 0

  bad <- is.finite(lb[lnm]) & p[lnm] <= lb[lnm]
  if (any(bad)) p[lnm][bad] <- lb[lnm][bad] + 1e-3
  if (viable(p)) break
}
if (att == 500) stop("no viable start for ", target, " chain ", chain)

## The draw above works in vl_thresh; the sampler works in its log when x$log_vthresh is set.
## Converted here, after viable(), so every check above reads the quantity it names.
if (isTRUE(x$log_vthresh)) p[["inc3"]] <- log(max(p[["inc3"]], 1e-3))
x$parameters <- p

## ---- optional: start from a converged fit of the other form -----------------------
## FIT_INIT_DIR names a directory of fitted chains; the final state of this target's first
## chain there replaces the random draw above.
##
## Why this is safe across forms: the hybrid incidence model shares its entire parameter
## vector with thresh_tdist -- same kinetics, same eta_T and omega_T, same vl_thresh -- so
## every element transfers with its meaning intact. It would NOT have been safe against the
## detectability-anchored version, where inc1-inc3 meant something else.
##
## Why it is worth doing: the likelihood is multimodal (a slow-growth mode with inflated
## a_sigma), there is one chain per variant here, and no second chain to cross-check against.
## Starting from a settled fit avoids spending the run escaping a bad mode.
##
## What it costs: if the estimates do not move, that is weak evidence of agreement, because a
## chain that never left its starting point looks the same as one that had nowhere to go. Only
## movement is strong evidence under this design. Say so when writing it up.
##
## One expected mismatch: the spline is ONSET incidence under thresh_tdist and THRESHOLD
## CROSSING incidence here, so the transferred coefficients start the epidemic curve about a
## delay -- roughly six days -- too late. The componentwise phase has to shift it, and the
## early likelihood will look poor until it does.
INIT_DIR <- Sys.getenv("FIT_INIT_DIR", "")
if (nzchar(INIT_DIR)) {
  cand <- list.files(INIT_DIR, sprintf("^%s_[^_]+(_[0-9]+)?\\.rds$", target), full.names = TRUE)
  if (!length(cand)) stop("FIT_INIT_DIR set but no chain for ", target, " in ", INIT_DIR)
  src <- sort(cand)[1]
  s   <- readRDS(src)
  sfl <- max(which(!is.na(s$MCMC_posteriors)))
  init <- s$MCMC_output[sfl, ]
  ## inc3 may be on a different scale in the source than in this fit -- chains written before
  ## the reparameterisation hold vl_thresh, later ones hold its log. Transferring it untouched
  ## would start the chain at exp(0.01) or log(2.1) without any error being raised.
  src_log <- isTRUE(s$x$log_vthresh)
  if (src_log && !isTRUE(x$log_vthresh)) init[["inc3"]] <- exp(init[["inc3"]])
  if (!src_log && isTRUE(x$log_vthresh)) init[["inc3"]] <- log(max(init[["inc3"]], 1e-3))
  rm(s); invisible(gc(verbose = FALSE))

  common <- intersect(names(init), names(x$parameters))
  if (length(common) != length(x$parameters))
    message(sprintf("  init: %d of %d parameters matched by name; the rest keep their draw",
                    length(common), length(x$parameters)))
  x$parameters[common] <- init[common]
  x$parameters[x$ignored_pars] <- 0
  cat(sprintf("initialised from %s (iteration %s)\n", basename(src), format(sfl, big.mark = ",")))
}

## ---- report and run -------------------------------------------------------------
geom <- function(v) {
  a <- 0.5 + exp(v[["a_bar"]]); b <- 0.25 + exp(v[["b_bar"]])
  vl0 <- as.numeric(sub("^vl_", "", rownames(x$data_array)[2]))
  sprintf("growth %.2f | decay %.2f | peak %.2f | a_sigma %.2f | onset %.2f d | detectable %.1f d",
          a, b, v[["vl_max_bar"]], v[["a_sigma"]], exp(v[["inc1"]]),
          (v[["vl_max_bar"]] - vl0) * (a + b) / (a * b))
}

## The form is named in the header because it is set by an environment variable as often as
## by an argument, so the log is the only record of which model a chain actually fitted.
cat(sprintf("%s | chain %d (seed %d) | form %s | fitted from %s upward | %d core%s | %s\n",
            target, chain, seed, FORM, rownames(x$data_array)[2], ncores,
            if (ncores > 1) "s" else "", R.version.string))
cat(sprintf("  schedule    MH to %s, adapt at %s over the preceding %s, sample from %s\n",
            format(COV_START[1], big.mark = ","),
            paste(format(COV_START, big.mark = ","), collapse = " and "),
            format(COV_WINDOW, big.mark = ","),
            format(COV_START[length(COV_START)], big.mark = ",")))
cat(sprintf("  start       %s\n", geom(x$parameters)))
## Printed because it is held for the entire run and never proposed: if it starts in the
## wrong place the multiplier spends the whole chain compensating, and nothing downstream
## complains. Worth a glance before committing a day of compute.
cat(sprintf("  pinned      %s = %.3f (from %s) | multiplier %.3f | %d thetas from the data\n",
            pin_name, x$parameters[[pin]], pin_src,
            x$parameters[["multiplier"]], length(th)))
## x$gene, not the chain number: an incidence-form chain is written as <target>_i<chain>_*.rds,
## so printing <target>_<chain>_*.rds sent anyone looking for the output to a filename that does
## not exist.
cat(sprintf("  %s iterations | test_pop %.0e | analysis window from %s | writing outputs/%s_%s_*.rds | started %s\n",
            format(N_ITER, big.mark = ","), TEST_POP,
            format(x$analysis_from, big.mark = ","), target, x$gene, format(Sys.time())))

t0 <- proc.time()
MCMC_function(x)
cat(sprintf("\n%s chain %d finished %s -- %.2f h\n", target, chain, format(Sys.time()),
            (proc.time() - t0)["elapsed"] / 3600))

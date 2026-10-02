
## Scenarios named <variant>sim are synthetic data generated at that variant's own
## posterior and fitted with that variant's real priors -- see input_generator. Defined once
## here because input_generator tests it in three places and a fifth variant name added in
## only two of them would fail silently.
VARIANT_SIM_RE <- "^(wt|alpha|delta|omicron)sim[0-9]*$"

.ct_cache <- new.env(parent = emptyenv())

ct_cache_reset <- function(ncores, seed = 20260813L){
  ## The simulation seeds must be the same on every call -- that is what makes the
  ## likelihood deterministic. But drawing them must NOT disturb the caller's random
  ## stream. This used to call set.seed() directly, which reset R's global RNG every time:
  ## MCMC_function calls this once at startup, so every chain's proposal sequence began
  ## from seed 20260813 no matter what the caller had set, and any code calling it in a
  ## loop got identical "random" draws after the first.
  if(!is.null(seed)){
    old <- if(exists(".Random.seed", envir = .GlobalEnv))
             get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    .ct_cache$pseeds <- unlist(dqrng::generateSeedVectors(ncores, 1))
    if(is.null(old)) rm(".Random.seed", envir = .GlobalEnv)
    else             assign(".Random.seed", old, envir = .GlobalEnv)
  } else {
    .ct_cache$pseeds <- unlist(dqrng::generateSeedVectors(ncores, 1))
  }
  .ct_cache$store  <- list()
  .ct_cache$hits   <- 0
  .ct_cache$misses <- 0
  invisible(.ct_cache$pseeds)
}

ct_cache_stats <- function(){
  c(hits = .ct_cache$hits, misses = .ct_cache$misses)
}

## A handful of entries is enough: a rejected viral kinetics proposal must not evict the
## current state's matrix, which the following theta proposals still need.
## form selects the kernel, and with it what the trajectory clock measures.
##
##   "thresh_tdist" (default, the model reported in the paper) -- ct_func2_cpp_thresh_tdist.
##     Rows are days since SYMPTOM ONSET. inc1-inc3 enter the kernel as eta_T, omega_T and
##     vl_thresh, so the onset delay is integrated inside the C++.
##
##   "incidence" -- ct_func2_cpp_inc. Rows are days since crossing vl_thresh, the SAME origin,
##     but the kernel stops there: the delay to symptom onset is convolved in R by
##     p_array_func against incidence indexed on threshold crossing. inc1 and inc2 keep their
##     meaning as eta_T and omega_T, and inc3 is vl_thresh in both.
##
## The two forms therefore share a parameter vector entirely and differ in exactly one
## structural respect -- whether incidence is indexed on symptom onset or on crossing the
## threshold -- which is what makes their estimates comparable parameter by parameter.
ct_sim_cached <- function(parameters, max_inf, vl_range, test_pop, ncores, vk, stoch,
                          form = "thresh_tdist", max_entries = 4L){
  if(is.null(.ct_cache$pseeds)) ct_cache_reset(ncores)

  key <- paste(c(form, sprintf("%.17g", c(parameters[1:11], max_inf, vl_range, test_pop,
                                          ncores, vk, stoch))), collapse="|")

  hit <- .ct_cache$store[[key]]
  if(!is.null(hit)){
    .ct_cache$hits <- .ct_cache$hits + 1
    return(hit)
  }

  res <- if(form == "incidence")
    ## vl_thresh (inc3, parameter 11) is appended as p[5] in R, read as p[4] in C++.
    ct_func2_cpp_inc(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)],
                     l=parameters[c(7,8)],
                     p=c(max_inf-1, vl_range[1], vl_range[2], test_pop, parameters[11]),
                     nthreads=ncores, pseeds=.ct_cache$pseeds, vg=6, vk, vt=1, stoch=stoch)
  else
    ct_func2_cpp_thresh_tdist(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)],
                                   l=parameters[c(7,8)], t=parameters[c(9,10,11)],
                                   p=c(max_inf-1, vl_range[1], vl_range[2], test_pop),
                                   nthreads=ncores, pseeds=.ct_cache$pseeds, vk, vt=1, stoch=stoch)

  .ct_cache$misses <- .ct_cache$misses + 1
  .ct_cache$store[[key]] <- res
  if(length(.ct_cache$store) > max_entries) .ct_cache$store <- .ct_cache$store[-1]

  return(res)
}
## ---------------------------------------------------------------------------------

p_array_func <- function(parameters, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5, name=F){
  if(length(dim(data_array))==3){
    max_test_day <- dim(data_array)[3]-1
    max_inf <- 2*(max_test_day+14)+1
  }
  if(length(dim(data_array))==2) max_inf <- 30
    
  inf_lims <- floor(max_inf/2)
  
  ### CHANGE THIS
  vl_vec <- as.numeric(substr(rownames(data_array), 4, nchar(rownames(data_array))))
  vl_range <- vl_vec[c(1,length(vl_vec))]
  vl_vec[1]="negative"
  
  theta_start <- which(names(parameters)=="theta1")
  thetas <- thetas_generator(parameters)[theta_start:(length(parameters)-1)]
  
  ## EXTRACT PARAMS
  inc <- parameters[c('inc1','inc2','inc3')]
  
  test_param <- c(parameters[grepl("test",names(parameters))],test6=1)
  
  #test_delay <- c(test_param, test6=max(0,1-sum(test_param)))
  
  test_delay <- setNames(vector(length=length(test_param)), names(test_param))
  test_delay[1] = test_param[1]
  for (i in 2:length(test_delay)) test_delay[i] <- prod(1-test_param[1:(i-1)])*test_param[i]
  
  days <- 0:max_inf
  infecteds <- infecteds_generator(parameters, knots, population, max_day, form) # generate infecteds
  
  #plot(infecteds)
  
  p_ct_tau_mat <- ct_sim_cached(parameters, max_inf, vl_range, test_pop, ncores, vk, stoch, form)

  dimnames(p_ct_tau_mat) <- list(0:(max_inf-1), rev(vl_vec))
  start_j = 36

  if(length(dim(data_array))==2) return(p_ct_tau_mat)

  else{
    ## CALCULATE p_ct_tau
    p_array <- array(NA, dim=c(length(vl_vec)-1, max_day, max_test_day+1))

    if(form == "incidence"){
      ## Incidence is indexed on crossing DETECTABILITY, and symptom onset follows after a
      ## delay tau with distribution inc_period. An individual tested k days after onset whose
      ## delay was tau therefore became detectable on day j-k-tau and is read off the kernel at
      ## trajectory time k+tau. Summing over tau is the convolution:
      ##
      ##   N[v,j,k] = test_delay[k] * sum_tau  infecteds[j-k-tau] * f(tau) * q[v, k+tau]
      ##
      ## This is what makes the form phase-sensitive. The weight on tau carries the calendar
      ## through infecteds[j-k-tau], so when incidence is rising, short delays are
      ## over-represented among a given day's onsets and the viral load composition at a FIXED
      ## k shifts with the epidemic. The thresh_tdist branch below has no such sum: its
      ## incidence is already indexed on onset, so the calendar only scales the counts.
      ##
      ## LOGNORMAL delay from crossing vl_thresh to symptom onset -- the same distribution
      ## thresh_tdist draws per individual inside the kernel, with the same parameters
      ## (inc1 = eta_T, inc2 = omega_T), discretised by rounding to the nearest whole day.
      ##
      ## The archived implementation used a skew normal in days with a third skew parameter.
      ## A skew normal has mass below zero, which the discretisation silently dropped rather
      ## than renormalised, by an amount that grew with the scale; the third parameter earned
      ## little; and a days-scale location could not be compared with thresh_tdist's log-scale
      ## eta_T without translation. The lognormal removes all three problems at once.
      inc_period <- stats::plnorm(0:symp_delay_lim + 0.5, meanlog=inc[1], sdlog=inc[2]) -
                    stats::plnorm(pmax(0:symp_delay_lim - 0.5, 0), meanlog=inc[1], sdlog=inc[2])
      ## Renormalise: mass beyond symp_delay_lim is truncated, and the sampler is free to
      ## propose a delay wide enough for that to matter.
      inc_period <- inc_period / sum(inc_period)

      ## The kernel must span k + tau for the largest of each.
      need <- max_test_day + symp_delay_lim + 1
      if(nrow(p_ct_tau_mat) < need)
        stop("max_inf gives ", nrow(p_ct_tau_mat), " trajectory rows but the incidence form ",
             "needs ", need, " (max_test_day ", max_test_day, " + symp_delay_lim ",
             symp_delay_lim, ")")

      for(j in start_j:max_day) for(k in 0:max_test_day){
        p_array[,j,k+1] <-
          t(c(infecteds[(j-k):(j-k-symp_delay_lim)] * inc_period * test_delay[k+1])) %*%
          p_ct_tau_mat[(k+1):(k+symp_delay_lim+1), (ncol(p_ct_tau_mat)-1):1]
      }
      incubation <- inc_period

    } else {
      # generate an array where the row pertains to ct value, the column to calendar day and the slice to days since onset
      for(j in start_j:max_day){
        p_array[,j,] <- t(c(infecteds[j:(j-max_test_day)]*test_delay) * p_ct_tau_mat[1:(max_test_day+1),(ncol(p_ct_tau_mat)-1):1])
      }
      incubation <- NA
    }

    dimnames(p_array) = list(as.numeric(vl_vec[-1]),1:max_day,0:max_test_day)

    return(list(p_array, p_ct_tau_mat, incubation))
  }
}

infecteds_generator <- function(parameters, knots, population, max_day, form){
  #thetas <- thetas_generator(parameters)
  thetas <- parameters[grepl("theta", names(parameters))]
  infecteds <- population*exp(stats::splinefun(x=c(0,knots*max_day), y=thetas)(1:(max_day)))*parameters['multiplier']
  if(form=="peak") infecteds <- c(rep(0,14),infecteds)
  
  names(infecteds) = 1:length(infecteds)
  
  return(infecteds)
}

thetas_generator <- function(parameters){
  start <- which(names(parameters) == "theta1")
  end <- length(parameters)-1
  
  #parameters[start:end] <- parameters[start:end]*parameters[length(parameters)]
  #parameters[length(parameters)] <- 1
  
  return(parameters)
}

input_generator <- function(n_iterations=1e6, strain="all", test_pop=1e7, ncores=8, proposal_type="MH", cov_matrix=NULL, population=5.6e7,
                            parameters=NA, cov_start=150000, max_test_day=6, gene="N", tag="", data_dir="data",
                            form="thresh_tdist", symp_delay_lim=20, log_vthresh=FALSE){
  
  
  if(strain %in% c("alpha", "delta", "omicron", "wt")){
    data_array <- readRDS(file.path(data_dir, paste0("data_array_vl_",strain,"_N",".rds")))[,,1:(max_test_day+1)]
    data_array <- data_array[,-1,]
    truth <- NA
  }
  ## <variant>sim -- synthetic data generated at that variant's own posterior.
  ##
  ## Deliberately NOT matched by the [LMNOPQR] pattern below, so none of the synthetic
  ## widening applies: these fits use the variants' real longitudinal priors, their real knot
  ## spacing and their real array geometry. That is the entire point of the design. The
  ## ladder scenarios ask whether the method works when the priors are weak; these ask
  ## whether it recovers the values the paper actually reports, under the priors the paper
  ## actually used -- which is the regime a reviewer cares about and the one where the
  ## posterior sits 3-16 prior SDs from the prior mean.
  ##
  ## The array is stored in the synthetic list form (out + pars_test) rather than the
  ## variants' bare-array form, because the generating values have to travel with it.
  if(grepl(VARIANT_SIM_RE, strain)){
    data_array_raw <- readRDS(file.path(data_dir, paste0("data_array_vl_",strain,"_N",".rds")))
    data_array <- data_array_raw$out[,,1:(max_test_day+1)]
    truth <- data_array_raw$pars_test
  }
  if(grepl("^[A-FLMNOPQRSTUVWX][0-3]a?$", strain)){
    data_array_raw <- readRDS(file.path(data_dir, paste0("data_array_vl_",strain,"_N",".rds")))
    data_array <- data_array_raw$out[,,1:(max_test_day+1)]
    truth <- data_array_raw$pars_test
  }  
  if("vl_0.5" %in% rownames(data_array)) data_array <- data_array[-which(rownames(data_array)=="vl_0.5"),,]
  #if("vl_1" %in% rownames(data_array)) data_array <- data_array[-which(rownames(data_array)=="vl_1"),,]
  
  array_dim_names <- dimnames(data_array)
  
  n_vls <- dim(data_array)[1]
  vls <- as.numeric(substr(array_dim_names[[1]], 4, 4+nchar(array_dim_names[[1]])-3))
  
  n_days <- dim(data_array)[2]
  dates <- array_dim_names[[2]]
  days <- 1:n_days
  
  n_test_days <- dim(data_array)[3]
  test_days <- as.numeric(substr(array_dim_names[[3]], 5, 5))
  
  if(grepl("^[A-FLMNOPQRSTUVWX][0-3]a?$", strain))     knots <- c(seq(0,99,9)[-1])/108
  if(strain=="all")                         knots <- c(seq(0,490,7)[-1], seq(497,n_days,7))/n_days
  if(strain=="omicron")                     knots <- c(seq(0,n_days-7,7))[-1]/n_days
  if(strain %in% c("alpha", "delta", "wt")) knots <- seq(0,n_days,14)[-1]/n_days
  ## Each <variant>sim inherits its parent's knot rule. The sim array has the parent's day
  ## count by construction, so this reproduces the parent's knots exactly -- which is what
  ## makes the spline dimension, and therefore the parameter vector, identical between the
  ## real fit and its simulated twin.
  if(grepl(VARIANT_SIM_RE, strain)){
    parent <- sub("sim$", "", strain)
    knots  <- if(parent == "omicron") c(seq(0,n_days-7,7))[-1]/n_days
              else                    seq(0,n_days,14)[-1]/n_days
  }
  
  theta_names <- paste0("theta", seq(1:(length(knots)+1)))
  
  #priors_file <- "longitudinal_fits/longitudinal_posteriors.rds"
  priors_file <- "longitudinal_fits/long_post_6.rds"
  priors <- readRDS(priors_file)
  vl_prior_means <- c(priors$mean, "l_sigma"=0)
  vl_prior_sds <- c(priors$sd , "l_sigma"=0)
  
  #priors_cov_file <- "longitudinal_fits/longitudinal_cov.rds"
  priors_cov_file <- "longitudinal_fits/long_cov_6.rds"
  prior_cov <- readRDS(priors_cov_file)
  
  ## vl_max_bar's prior is the tightest in the model -- SD 0.100 from the longitudinal fit,
  ## against variant posteriors sitting 6 to 16 prior SDs above its mean. It is singled out
  ## for widening because the peak enters every duration as (vl_peak - vl_0)/rate, so holding
  ## it too low forces compensating error into the growth and decay rates rather than
  ## absorbing it locally.
  ##
  ## Scaling the row AND column by k, with the diagonal thereby scaled by k^2, raises this
  ## one variable's SD by k while leaving every correlation unchanged.
  ##
  ## VL_PEAK_WIDEN extends the same treatment to the <variant>sim scenarios, where it is a
  ## deliberate experiment rather than a default: does loosening this one prior improve
  ## recovery of the OTHER kinetics parameters? Unset, those fits use the variants' real
  ## prior, which is the point of that study.
  vl_widen <- if(grepl("^[A-FLMNOPQRSTUVWX][0-3]a?$", strain)) 2 else
              if(grepl(VARIANT_SIM_RE, strain))
                as.numeric(Sys.getenv("VL_PEAK_WIDEN", "1")) else 1
  if(vl_widen != 1){
    vl_prior_sds["vl_max_bar"] <- vl_prior_sds["vl_max_bar"] * vl_widen
    prior_cov["vl_max_bar", ]  <- prior_cov["vl_max_bar", ] * vl_widen
    prior_cov[, "vl_max_bar"]  <- prior_cov[, "vl_max_bar"] * vl_widen
    message(sprintf("  vl_max_bar prior SD widened %gx to %.3f", vl_widen,
                    vl_prior_sds[["vl_max_bar"]]))
  }

  ## KIN_WIDEN -- the real-data prior sensitivity analysis. Widens the WHOLE kinetics block
  ## for the four real variants, which VL_PEAK_WIDEN deliberately does not do.
  ##
  ## Why the block and not vl_max_bar alone. Measured on the reported chains, the ratio of
  ## prior SD to posterior SD is 1.72 for b_bar and 1.98 for b_sigma against 4.06 for
  ## vl_max_bar -- so the DECAY priors constrain the posterior more than the peak prior does,
  ## and vl_max_bar is only the fourth tightest in absolute SD (0.100, behind b_sigma 0.074,
  ## vl_max_sigma 0.079 and b_bar 0.087). Singling out the peak would be an eyeballed choice
  ## that a reviewer repeating this check would immediately question.
  ##
  ## inc1-inc3 are deliberately EXCLUDED: their prior/posterior SD ratios are 18-424, so the
  ## likelihood already dominates them by one to two orders of magnitude and widening them
  ## changes nothing while adding a choice to defend.
  ##
  ## Scaling rows AND columns by k scales the [kin,kin] block by k^2 and the cross terms by
  ## k, so every SD in the block rises by k and every correlation is unchanged -- equivalent
  ## to Sigma' = D Sigma D with D = diag(k). The cross-block terms are zero here anyway.
  ##
  ## 2x matches the width used for the L/M/N/O/P/Q/R simulation ladder, so the two analyses
  ## speak to each other: at the prior width where simulation shows the method recovers
  ## truth, this measures how far the real-data estimates move.
  kin_widen <- if(strain %in% c("wt","alpha","delta","omicron"))
                 as.numeric(Sys.getenv("KIN_WIDEN", "1")) else 1
  if(kin_widen != 1){
    kin <- intersect(c("a_bar","a_sigma","b_bar","b_sigma","vl_max_bar","vl_max_sigma"),
                     names(vl_prior_sds))
    vl_prior_sds[kin] <- vl_prior_sds[kin] * kin_widen
    prior_cov[kin, ]  <- prior_cov[kin, ] * kin_widen
    prior_cov[, kin]  <- prior_cov[, kin] * kin_widen
    message(sprintf("  kinetics prior SDs widened %gx: %s", kin_widen,
                    paste(sprintf("%s %.4f", kin, vl_prior_sds[kin]), collapse = ", ")))
  }
  
  test_prior_means_raw <- setNames((apply(data_array, 3, sum)/sum(data_array)), paste0("test",seq(0,max_test_day)))
  test_prior_means <- setNames(vector(length=length(test_prior_means_raw)-1), names(test_prior_means_raw[-length(test_prior_means_raw)]))
  
  test_prior_means[1] <- test_prior_means_raw[1]
  for(i in 2:length(test_prior_means)) test_prior_means[i] = test_prior_means_raw[i]/prod(1-test_prior_means[1:(i-1)])
  
  test_prior_sds <-   setNames(rep(0.5, max_test_day), paste0("test",seq(0,max_test_day-1)))
  
  spline_prior_means <- setNames(rep(-9,length(theta_names)),theta_names)
  spline_prior_sds <- setNames(c(2, 2, rep(9,length(theta_names)-4), 1, 1),theta_names)

  ## Synthetic scenarios get a flat spline prior instead of the tight ends.
  ##
  ## The mean of -9 is a statement about incidence in a population of 5.6e7 during a real UK
  ## wave, and for the variants it is close to right -- their pinned coefficients come out at
  ## -8.4 to -8.6. The simulation ladder deliberately moves incidence down three orders of
  ## magnitude while holding the population fixed, so the same statement becomes false by up
  ## to 7 log units, and an SD of 1 at the end knots turns a false statement into a hard
  ## constraint. Evaluated at the generating values the log-prior penalty is 16.6 at the
  ## richest rung and 126.1 at the poorest -- twice the 60.7-unit gap between alpha's two
  ## likelihood modes, which was enough to trap chains for an entire run.
  ##
  ## Worse, the penalty grows quadratically as incidence falls, so it is confounded with the
  ## thing the ladder exists to measure: some of the apparent loss of precision at small N is
  ## the prior becoming more wrong, not the data becoming less informative.
  ##
  ## The end regularisation is dropped rather than re-centred. Re-centring on the pinned
  ## coefficient would import the generating value into every theta, and the spline ends
  ## being weakly identified is a result worth showing, not one to shrink away.
  if(grepl("^[A-FLMNOPQRSTUVWX][0-3]a?$", strain))
    spline_prior_sds <- setNames(rep(9, length(theta_names)), theta_names)
  
  spline_multiplier_mean <- c(multiplier=1)
  spline_multiplier_sd <- c(multiplier=10)

  ## Both forms give inc1 and inc2 the same meaning -- the log-location and log-scale of the
  ## delay to symptom onset -- so they share a prior, and the two fits are comparable without
  ## translation. inc3 is vl_thresh under thresh_tdist and unused under the incidence form,
  ## where it is pinned by ignored_pars below; its prior is then a constant either way.
  ## Symptom-timing priors, all on the NATURAL scale of the quantity they describe -- inc3's
  ## prior is a statement about vl_thresh itself even when the sampler works with its log.
  ##
  ## inc1 (eta_T, the log median delay from crossing vl_thresh to symptom onset): median 3 days,
  ##   95% range 0.9-9.7. Anchored on published incubation-period estimates, which centre near
  ##   5 days for the ancestral variant and are shorter for Delta and Omicron, less a latent
  ##   period of roughly 2 days between infection and detectable virus. Wide, because the delay
  ##   starts at the fitted threshold rather than at infection and would be much shorter were
  ##   that threshold estimated close to peak.
  ## inc2 (omega_T): centred on the log-scale spread of published lognormal incubation fits,
  ##   about 0.42-0.5. At 0.45 an individual's delay spans a factor of ~2.4 either side of the
  ##   median; the old N(0.5, 0.5) admitted a factor of 20, and lost 16% of its mass to the
  ##   truncation at zero against 1.2% here.
  ## inc3 (vl_thresh): half-normal, i.e. monotonically decreasing from zero, so the prior
  ##   encodes an ordering -- lower is more plausible -- rather than a location, since this
  ##   quantity has never been measured in swab units. Innate sensing is cell-intrinsic, so the
  ##   response that leads to symptoms is initiated when very few cells are infected, at viral
  ##   loads below anything an assay can measure; 55% of the mass lies below the limit of
  ##   detection. The scale puts 95% below 4.5 log10, well under a peak of 8-9: a threshold
  ##   near peak would mean symptoms beginning only as viral load peaks, in everyone.
  ##   It also keeps every ladder scenario's generating value (1.3-4.0) inside the central 95%.
  incu_prior_means <- c(inc1=log(3), inc2=0.45, inc3=0)
  incu_prior_sds   <- c(inc1=0.6,    inc2=0.20, inc3=2)

  ## ONSET_WIDEN scales the symptom-timing priors -- the counterpart to KIN_WIDEN above, and
  ## the other half of the real-data prior sensitivity check. The four real variants only, for
  ## the same reason: the synthetic scenarios already run under widened priors by construction.
  ##
  ## Only the SD moves. The means stay, so the widened fit tests whether the data can hold the
  ## posterior in place against a prior that has stopped pinning it, rather than testing what
  ## happens when the prior is moved somewhere else.
  ##
  ## Two things to watch in the result, both consequences of the reasoning in the synthetic
  ## block below. Widening inc1 opens the long-onset region where the spurious mode lives --
  ## slow growth, onset near 15 days, a_sigma above 2.3 -- so a continuation that drifts there
  ## has found that mode, not moved within the good one, and should be discarded rather than
  ## reported. And doubling inc3's half-normal scale from 2 to 4 moves its 95% point from 4.5
  ## to 7.8 log10, which is no longer comfortably below a peak of 8-9: of the three, this is
  ## the prior that genuinely stops being informative, so it is the one whose posterior moving
  ## would mean the most.
  onset_widen <- if(strain %in% c("wt","alpha","delta","omicron"))
                   as.numeric(Sys.getenv("ONSET_WIDEN", "1")) else 1
  if(onset_widen != 1){
    incu_prior_sds <- incu_prior_sds * onset_widen
    message(sprintf("  onset prior SDs widened %gx: %s", onset_widen,
                    paste(sprintf("%s %.4f", names(incu_prior_sds), incu_prior_sds),
                          collapse = ", ")))
  }

  ## Synthetic scenarios: the six VIRAL KINETICS priors are widened.
  ##
  ## This is what makes the simulation study a test of the METHOD rather than of the prior.
  ## R1.4b asks whether the approach works without informative priors, and E1 asks for
  ## recovery that is not an artefact of them; under the variants' priors, a parameter
  ## landing on its generating value could mean the likelihood found it or simply that the
  ## prior was already sitting there. The SDs come from the longitudinal fit and run 0.07 to
  ## 0.22, tight enough to dominate.
  ##
  ## Note this compounds with the vl_max_bar doubling above, so that one parameter ends at
  ## 2 * VL_PRIOR_WIDEN times the variants' SD while the other five are at VL_PRIOR_WIDEN.
  ##
  ## Only these six. The onset priors are deliberately left alone: widening inc1 and inc2
  ## opens up the long-onset region, which is exactly where the spurious mode lives (slow
  ## growth, onset around 15 days, a_sigma above 2.3), and only 2 of 10 dispersed alpha
  ## chains escaped it when it was reachable. The remaining priors -- rdisp N(1, 10) against
  ## a truth near 0.05, inc3 N(4, 4), the test hazards at SD 0.5 on probabilities -- are
  ## already weak relative to their own scale, so scaling them would change little.
  ##
  ## The covariance is scaled by the SQUARE so the CORRELATION structure of the longitudinal
  ## posterior survives. Scaling only vl_prior_sds would leave prior_cov -- which is what
  ## prior() actually evaluates for parameters 1:7 -- untouched, and the two records of the
  ## same prior would silently disagree.
  ##
  ## The thetas are excluded here: they are already flat at SD 9 above, which against a mean
  ## of -9 is a coefficient of variation of 1.
  VL_PRIOR_WIDEN <- 2
  if(grepl("^[A-FLMNOPQRSTUVWX][0-3]a?$", strain)){
    vl_prior_sds <- vl_prior_sds * VL_PRIOR_WIDEN
    prior_cov    <- prior_cov * VL_PRIOR_WIDEN^2
  }

  prior_mean <- c(vl_prior_means, incu_prior_means, test_prior_means, rdisp=1,  spline_prior_means, spline_multiplier_mean)
  prior_sds <-  c(vl_prior_sds,   incu_prior_sds,   test_prior_sds,   rdisp=10, spline_prior_sds,   spline_multiplier_sd)
  
  ignored_spline_par <- length(prior_mean) - ceiling(length(theta_names)/2)
  
  ## 7 and 8 are l_bar and l_sigma, the degenerate plateau parameters, held fixed in every fit.
  ## inc3 is free in both forms: it is vl_thresh in both, and under the incidence form it sets
  ## the origin of the trajectory clock that incidence is indexed on.
  ignored_pars <- c(7, 8)
  
  if(is.na(parameters)[1]){
    parameters <- c(prior_mean)
    parameters[-ignored_spline_par] <- rnorm(n=length(parameters)-1, mean=parameters[-ignored_spline_par], sd=abs(parameters[-ignored_spline_par])*0.1)
  }

  shift <- unname(prior_mean[ignored_spline_par] - parameters[ignored_spline_par])
  if(shift != 0){
    parameters[theta_names]  <- parameters[theta_names] + shift
    parameters['multiplier'] <- parameters['multiplier'] * exp(-shift)
  }
  parameters[ignored_spline_par] <- prior_mean[ignored_spline_par]

  parameters[ignored_pars] <- 0

  MCMC_sds         <- abs(parameters)/200
  MCMC_sds['inc3'] <- 0.05

  ## LOG PARAMETERISATION OF vl_thresh -- available, but OFF by default. See the warning below.
  ##
  ## With log_vthresh the sampler works with theta = log(vl_thresh) and the kernels are passed
  ## exp(theta), so positivity is automatic and the hard bound at zero disappears -- about 11%
  ## of block proposals were rejected there, and inc3's tail ESS sat below its bulk ESS, the
  ## signature of a boundary.
  ##
  ## WHY IT IS OFF. The target is provably unchanged (the Jacobian below makes it an exact
  ## reparameterisation, verified to 1e-12), but the SAMPLER behaves differently and worse.
  ## The Jacobian tilts the density in theta so that its mode sits at vl_thresh = sd rather
  ## than at zero, and in the first few thousand iterations -- while the spline and kinetics
  ## are still poor and the likelihood barely constrains the threshold -- chains drift up to
  ## 3.5-5. The other parameters then adapt around that value and the chain is trapped: at
  ## such a state, moving vl_thresh alone to 0.02 costs 560 log units for alpha and 7,642 for
  ## delta, because escape needs vl_thresh, growth and the onset delay to move together and
  ## componentwise sampling moves one at a time. Eight of twelve chains were stuck this way
  ## on 2026-09-23, against none under the natural scale, where the bound at zero keeps early
  ## excursions short. The boundary rejections are the cheaper problem.
  ##
  ## The PRIOR IS UNCHANGED: prior() evaluates the truncated normal above at exp(theta) and adds
  ## theta as the Jacobian, so this is a reparameterisation and not a different model. That also
  ## keeps the left tail proper -- the likelihood is flat below the limit of detection, so a flat
  ## prior on theta would let the chain drift to minus infinity, while the e^theta factor does not.
  ##
  ## Everything outside the sampler stays on the natural scale: MCMC_function exponentiates
  ## before the likelihood, and the draw helpers in paper_outputs.R convert on read, so figures,
  ## tables and the recovery comparisons against generating values are unaffected.
  if(log_vthresh){
    parameters['inc3'] <- log(max(parameters[['inc3']], 1e-3))
    ## A step of 0.05 is right for a parameter of order 1, not for one ranging over several
    ## log units; the componentwise phase would otherwise crawl along theta.
    MCMC_sds['inc3']   <- 0.5
  }

  ## The incidence form reads the kernel at trajectory time k + tau, so it needs
  ## max_test_day + symp_delay_lim + 1 rows; 43 covers 6 + 20 + 1 = 27 with room to spare.
  ## symp_delay_lim is NA under thresh_tdist, where no convolution happens.
  max_inf <- 43
  if(form != "incidence") symp_delay_lim <- NA
  
  spline_pars <- names(parameters)[which(grepl("theta",names(parameters)))]
  
  lower_bound_vl <- c(a_bar=-Inf, a_sigma=0.05, b_bar=-Inf, b_sigma=0.05, vl_max_bar=-Inf, vl_max_sigma=0, l_bar=-Inf, l_sigma=0)
  ## inc3's bound is the natural-scale one and is applied inside prior() after exponentiating;
  ## on the log scale the parameter itself is unbounded.
  lower_bound_inc <- c(inc1= -Inf, inc2=0, inc3=if(log_vthresh) -Inf else 0)
  lower_bound_test <- c(test0=0, test1=0, test2=0, test3=0, test4=0, test5=0, rdisp=-Inf)
  lower_bound_splines <- setNames(rep(-Inf, length(spline_pars)), spline_pars)
  
  lower_bound <- c(lower_bound_vl, lower_bound_inc, lower_bound_test, lower_bound_splines)
  
  prior_cov_final <- matrix(0, nrow=length(parameters), ncol=length(parameters))
  colnames(prior_cov_final) <- names(parameters)
  rownames(prior_cov_final) <- names(parameters)
  
  prior_cov_final[1:7,1:7] <- prior_cov
  for(i in 9:length(parameters)) prior_cov_final[i,i] <- prior_sds[i]^2
  
  input <- list(n_iterations=n_iterations,
                parameters=parameters,
                knots=knots,
                data_array=data_array,
                test_pop=test_pop,
                MCMC_sds=MCMC_sds,
                prior_mean=prior_mean,
                prior_sds=prior_sds,
                lower_bound=lower_bound,
                prior_cov_final=prior_cov_final,
                ncores=ncores,
                proposal_type=proposal_type,
                cov_matrix = cov_matrix,
                population = population,
                tag=tag,
                ignored_spline_par = ignored_spline_par,
                strain=strain,
                form=form,
                ## Objects written before this existed sampled vl_thresh on the natural scale;
                ## every reader treats a missing flag as FALSE.
                log_vthresh=log_vthresh,
                max_inf=max_inf,
                symp_delay_lim=symp_delay_lim,
                max_day=n_days,
                cov_start=cov_start,
                gene=gene,
                n_vls = n_vls,
                vls=vls,
                n_days=n_days, 
                days=days,
                dates=dates,
                n_test_days=n_test_days,
                test_days=test_days,
                n_parameters=length(parameters),
                ignored_pars=ignored_pars,
                truth=truth)
  
  return(input)
}


## prior function
## log_vthresh: parameters['inc3'] holds log(vl_thresh). The prior is still the truncated normal
## on vl_thresh itself, so it is evaluated at exp(inc3) with the Jacobian |d vl_thresh/d inc3| =
## exp(inc3) added on the log scale. Sampling is reparameterised; the model is not.
prior <- function(parameters, prior_mean, prior_cov_final, lower_bound, form, ignored_pars,
                  log_vthresh=FALSE){
  spline_pars <- names(parameters)[which(grepl("theta",names(parameters)))]
  #parameters[spline_pars] <- parameters[spline_pars]*parameters['multiplier']
  
  vl_pars <- parameters[1:8][-ignored_pars]
  vl_pars_names <- names(vl_pars)
  
  last_par <- length(parameters)-1
  
  non_vl <- parameters[9:last_par]
  lb_non <- lower_bound[9:last_par]
  jacobian <- 0
  if(log_vthresh){
    i3 <- which(names(non_vl) == "inc3")
    jacobian   <- non_vl[[i3]]            # log|d vl_thresh / d log(vl_thresh)| = log(vl_thresh)
    non_vl[i3] <- exp(non_vl[[i3]])
    lb_non[i3] <- 0                       # the natural-scale prior keeps its truncation at zero
  }

  prior_val_vl <- TruncatedNormal::dtmvnorm(x=vl_pars, mu=prior_mean[vl_pars_names], sigma=prior_cov_final[vl_pars_names,vl_pars_names], lb=lower_bound[vl_pars_names], log=T)
  prior_val_non_vl <- log(truncnorm::dtruncnorm(x=non_vl, a=lb_non, mean=prior_mean[9:last_par], sd=sqrt(diag(prior_cov_final)[9:last_par])))
  prior_val <- sum(prior_val_vl, prior_val_non_vl) + jacobian
  
  return(prior_val)
}

## likelihood function
likelihood_function3 <- function(parameters, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim){
  rdisp <- parameters['rdisp']
  p_array <- p_array_func(parameters, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5)[[1]]
  
  p_array <- weekly_aggregator(p_array)
  
  index_start <- which(is.na(p_array[1,,1])==F)[1]
  
  likelihood <- sum(dnbinom(x=data_array[-1,index_start:ncol(p_array),], size=1/rdisp, mu=p_array[,index_start:ncol(p_array),], log=T))
  
  return(likelihood)
}

likelihood_REACT <- function(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stop=F){
  p_array <- p_array_func(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5, name=T)
  
  infecteds <- round(infecteds_generator(parameters, knots, population, max_day, form))
  
  N_matrix <- matrix(NA, nrow=nrow(data_array), ncol=ncol(data_array))
  rownames(N_matrix) <- c(rev(rownames(data_array)[-1]),"negative") 
  colnames(N_matrix) <- colnames(data_array)
  
  for(j in 31:ncol(N_matrix)) N_matrix[,j] <- t(p_array) %*% infecteds[j:(j-29)]
  
  N_matrix[nrow(N_matrix),] <- population-colSums(N_matrix[1:(nrow(N_matrix)-1),])
  p_matrix <- N_matrix/population
  
  if(stop==T) return(p_matrix)
  
  likelihood <- sum(data_array[nrow(data_array):1,(31:max_day)] * log(p_matrix[,31:max_day]))
  
  return(likelihood)
}

weekly_aggregator <- function(daily_array){
  dimnames2 <- seq(3.5, dim(daily_array)[2]-3.5, 7)
  
  if(dim(daily_array)[2] %% 7 != 0) print("non integer number of weeks", quote=F)
  
  if(length(dim(daily_array))==2){
    weekly_array <- array(dim=c(dim(daily_array)[1], ceiling(dim(daily_array)[2]/7)))
    for(j in 1:ncol(weekly_array)) weekly_array[,j] <- apply(daily_array[,(7*(j-1)+1):(7*j)],c(1),sum)
    dimnames(weekly_array) <- list(dimnames(daily_array)[[1]], dimnames2)
  } 
  
  if(length(dim(daily_array))==3){
    weekly_array <- array(dim=c(dim(daily_array)[1], ceiling(dim(daily_array)[2]/7), dim(daily_array)[3]))
    for(j in 1:ncol(weekly_array)) weekly_array[,j,] <- apply(daily_array[,(7*(j-1)+1):(7*j),],c(1,3),sum)
    dimnames(weekly_array) <- list(dimnames(daily_array)[[1]], dimnames2, dimnames(daily_array)[[3]])
  }
  
  return(weekly_array)
}

output_list <- function(strain, form, corr, flat, recentness){
  files <- list.files(paste0(strain))
  form_filter = grepl(paste0("N_",form,"_Cov"), files) | grepl(paste0("N_",form,"_MH"), files)
  #corr_filter = grepl(paste0("corr",corr), files)
  flat_filter = grepl("flat", files)
  vg_filter = grepl(paste0("vg=",corr), files)
  mult_filter = grepl("mult", files)
  check_filter <- !grepl("checking", files)
  #vg_filter = !grepl("vg", files)
  
  if(flat==0) flat_filter <- !flat_filter
  
  filtered_files <- files[which(form_filter*flat_filter*vg_filter*mult_filter*check_filter==1)]
  
  return(filtered_files)
}

output_finder <- function(strain, form, corr, flat, recentness=1){
  filtered_files <- output_list(strain, form, corr, flat, recentness)
  
  #recentness <- (length(filtered_files))
  
  if(length(filtered_files)==0){
    print("No file") 
    return("No file")
  } 
  info <- file.info(paste0(strain,"/",filtered_files))
  
  times <- info$ctime
  file_name <- paste0(strain, "/", filtered_files[which(times==sort(times, decreasing=T)[recentness])])
  print(file_name, quote=F)
  return(readRDS(file_name))
}

output_finder2 <- function(strain){
  
  files <- list.files("outputs/")
  strain_file <- files[grepl(strain, files) * grepl("Nc", files)] 
  
  
  if(length(filtered_files)==0){
    print("No file") 
    return("No file")
  } 
  
  return(readRDS(strain_file))
}

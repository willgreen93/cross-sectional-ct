## MCMC function
MCMC_function <- function(x){
  n_iterations <- x$n_iterations
  parameters <- x$parameters
  knots <- x$knots
  population <- x$population
  data_array_raw <- x$data_array
  test_pop <- x$test_pop
  MCMC_sds <- x$MCMC_sds
  prior_mean <- x$prior_mean
  prior_sds <- x$prior_sds
  ncores <- x$ncores
  proposal_type <- x$proposal_type
  cov_matrix <- x$cov_matrix
  tag <- x$tag
  strain <- x$strain
  max_inf <- x$max_inf
  symp_delay_lim <- x$symp_delay_lim
  max_day <- x$max_day
  ignored_spline_par <- x$ignored_spline_par
  cov_start <- x$cov_start
  gene <- x$gene
  days <- x$days
  dates <- x$dates
  n_days <- x$n_days
  vls <- x$vls
  n_vls <- x$n_vls
  n_test_days <- x$n_test_days
  test_days <- x$test_days
  n_parameters <- x$n_parameters
  ignored_pars <- x$ignored_pars 
  vg <- x$vg
  prior_cov_final <- x$prior_cov_final
  lower_bound <- x$lower_bound
  ## Taken from the fitted object, defaulting to the form the paper reports so that every
  ## existing chain and every object written before x$form existed behaves exactly as before.
  form <- if(is.null(x$form)) "thresh_tdist" else x$form
  ## Whether inc3 holds log(vl_thresh). Absent in objects written before the reparameterisation,
  ## which sampled it on the natural scale.
  log_vth <- isTRUE(x$log_vthresh)
  truth <- x$truth

  ## Both kernels are compiled: the incidence form calls ct_func2_cpp_inc, the threshold form
  ## ct_func2_cpp_thresh_tdist, and ct_sim_cached dispatches on form. Compiling both costs a
  ## few seconds once and removes a failure mode where the wrong kernel is loaded for the form.
  Rcpp::sourceCpp("scripts/simulate_onset.cpp")
  Rcpp::sourceCpp("scripts/simulate_incidence.cpp")

  source("scripts/model.R")

  ## Fix the simulation seeds for the lifetime of this run and clear any cached matrix
  ## left over from a previous one.
  ct_cache_reset(ncores)

  if(length(dim(data_array_raw))==3) data_array <- weekly_aggregator(daily_array=data_array_raw)
  if(length(dim(data_array_raw))==2) data_array <- data_array_raw
  
  ## posterior function
  posterior <- function(parameters, prior_mean, prior_cov_final, lower_bound, knots, vk, max_day, population, data_array, test_pop, ncores, form="thresh_tdist", symp_delay_lim, ignored_pars){
    a <- prior(parameters, prior_mean, prior_cov_final, lower_bound, form, ignored_pars, log_vth)

    if(a==-Inf) return(-Inf)

    ## The kernels take vl_thresh on its natural scale, so the transform happens here and
    ## nowhere else: p_array_func, ct_sim_cached and the C++ are untouched by the
    ## reparameterisation, and the cache is keyed on the natural-scale vector as before.
    pars_nat <- parameters
    if(log_vth) pars_nat[["inc3"]] <- exp(pars_nat[["inc3"]])

    if(length(dim(data_array))==3) b <- likelihood_function3(pars_nat, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim)
    if(length(dim(data_array))==2) b <- likelihood_REACT(pars_nat, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim)
    
    if(is.nan(b)) b <- -Inf
    
    return(c(a+b,b))
  }
  
  ## proposal function
  proposal <- function(parameters, iteration, MCMC_sds, proposal_type, cov_matrix, ignored_pars, ignored_spline_par){
    if(proposal_type == "MH"){
      param_interest <- (iteration-1) %% length(parameters) + 1
      
      new_param <- rnorm(n=1, mean=parameters[param_interest], sd=MCMC_sds[param_interest])
      parameters[param_interest] <- new_param
    }
    
    if(proposal_type == "Cov"){
      output <- c(mvtnorm::rmvnorm(n=1, mean=parameters, sigma=cov_matrix))
      names(output) <- names(parameters)
      output[ignored_spline_par] <- parameters[ignored_spline_par]
      output[ignored_pars] <- parameters[ignored_pars]
      #if(test_form=='empirical') output[c('test1', 'test2')] <- parameters[c('test1', 'test2')]
      
      parameters <- output
    }
    
    return(parameters)
  }
  
  MCMC_output <- matrix(nrow = n_iterations + 1, ncol = n_parameters)
  MCMC_output[1,] <- parameters
  colnames(MCMC_output) <- names(parameters)
  
  vk <- 1
  
  p <- proc.time()
  #test_vec <- matrix(nrow=10, ncol=2)
  #for(i in 1:nrow(test_vec)){
  #  test_vec[i,1] <- posterior(parameters+0.0001*i, prior_mean, prior_cov_final, lower_bound, knots, vk=vk, max_day, population, data_array, test_pop=1e6, ncores, form, symp_delay_lim, ignored_pars)[1]
  #  test_vec[i,2] <- posterior(parameters+0.0001*i, prior_mean, prior_cov_final, lower_bound, knots, vk=vk, max_day, population, data_array, test_pop=1e7, ncores, form, symp_delay_lim, ignored_pars)[1]
  #}
  current_posterior <- posterior(parameters, prior_mean, prior_cov_final, lower_bound, knots, vk=vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars)
  current_posterior
  print((proc.time()-p)["elapsed"])
  
  MCMC_posteriors <- matrix(nrow=n_iterations + 1, ncol=1)
  MCMC_likelihoods <- matrix(nrow=n_iterations + 1, ncol=1)
  
  MCMC_posteriors[1] <- current_posterior[1]
  MCMC_likelihoods[1] <- current_posterior[2]
  
  Acceptances <- vector(length=n_iterations)
  time_vec <- vector(length=n_iterations)
  
  pr <- 100

  ## Draws used for each covariance adaptation. Defaulted rather than required so that
  ## fits saved before this argument existed still run.
  cov_window <- if (is.null(x$cov_window)) 50000L else x$cov_window

  ## Block timer. Reported at the TOP of the iteration rather than the bottom so that the
  ## pinned-parameter `next` below cannot skip the reset -- if it did, the following block
  ## would silently cover 200 iterations and read twice as slow.
  p_blk <- proc.time()
  i_blk <- 0L

  for(i in 1:n_iterations){
    #if(i > 150000) test_pop <- 1e7
    p <- proc.time()
    if(i %% pr == 0){
      cat("\n")
      ## Cost of the block just finished, not of one sampled iteration. A single iteration
      ## is a poor guide: cache hits and misses differ by an order of magnitude, and one
      ## print in fifty lands on the save iteration and carries the saveRDS cost with it.
      blk <- (proc.time() - p_blk)[["elapsed"]]
      nb  <- i - 1L - i_blk
      ## Name the chain on every block. Each chain has its own log, so this is redundant
      ## there -- but it makes any view that mixes them (a status listing, a concatenated
      ## tail, a screen you scrolled past the header on) self-identifying.
      cat(sprintf("%s_%s | iterations %s-%s | %.1f s | %.3f s/iter | %.1f h remaining of %s\n",
                  strain, gene,
                  format(i_blk + 1L, big.mark = ","), format(i - 1L, big.mark = ","),
                  blk, blk / nb, blk / nb * (n_iterations - i) / 3600,
                  format(n_iterations, big.mark = ",")))
      p_blk <- proc.time()
      i_blk <- i - 1L
    }
    param_interest <- (i-1) %% n_parameters + 1
    if(i %% pr == 0) print(i)
    if(i %% pr == 0) print(param_interest)
    
    ## Switching the proposal and adapting the covariance are decisions about the CHAIN, so
    ## they must happen before the per-parameter skip below -- which calls next() and would
    ## otherwise carry them off with it. When iteration cov_start landed on a pinned
    ## parameter the adaptation was never run: the following iteration switched to "Cov"
    ## with cov_matrix still NULL, and rmvnorm died on isSymmetric(NULL). Every simulation
    ## scenario had n_parameters = 31 and 100000 %% 31 = 25, the pinned spline coefficient,
    ## so all twelve were destined to fail at exactly the same iteration.
    if(i >= min(cov_start)) proposal_type <- "Cov"

    ## Adapt at every point in cov_start -- a scalar still works -- over the preceding
    ## cov_window draws. The window must lie entirely after burn-in: estimating the
    ## covariance over pre-stationary draws inflates it along whatever directions the chain
    ## travelled while converging, and the 2.38^2/d scaling then over-steps in exactly those
    ## directions. That is what mistuned alpha, which ran 100,000 iterations at 6.5%
    ## acceptance and reached ESS 5 on a_bar before it was found.
    if(i %in% cov_start)
      cov_matrix <- cov(MCMC_output[max(1L, i - cov_window):i,])*(2.38^2)/n_parameters

    ## Belt and braces: any other route into the block proposal without a covariance would
    ## fail the same way, and failing 100,000 iterations in is expensive.
    if(proposal_type=="Cov" && is.null(cov_matrix))
      cov_matrix <- cov(MCMC_output[max(1L, i - cov_window):i,])*(2.38^2)/n_parameters

    ## One parameter at a time during the componentwise phase; the pinned coefficients are
    ## simply not proposed. In "Cov" mode every parameter moves together and proposal()
    ## holds the pinned ones fixed itself, so this branch correctly stops applying.
    if(proposal_type=="MH" && (i %% n_parameters %in% c(ignored_spline_par,ignored_pars))){
      if(i %% pr == 0) print("Rejected", quote=F)
      MCMC_output[i+1,] <- MCMC_output[i,]
      Acceptances[i] <- NA
      MCMC_posteriors[i+1] <- MCMC_posteriors[i]
      MCMC_likelihoods[i+1] <- MCMC_likelihoods[i]
      next
    }

    #if((proposal_type == "Cov") & (i == 100000)) cov_matrix <- cov(MCMC_output[(i-50000):i,])*(2.38^2)/n_parameters
    
    current_parameters <- MCMC_output[i,]
    current_posterior <- MCMC_posteriors[i]
    
    proposed_parameters <- proposal(parameters=current_parameters,iteration=i,MCMC_sds=MCMC_sds,proposal_type=proposal_type, cov_matrix=cov_matrix, ignored_pars=ignored_pars, ignored_spline_par=ignored_spline_par)
    if(i %% pr == 0) print(c("proposed parameters = ", proposed_parameters), quote=F)
    
    posterior_calc <- posterior(parameters=proposed_parameters, prior_mean, prior_cov_final, lower_bound, knots, vk, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars)
    
    proposed_posterior <- posterior_calc[1]
    proposed_likelihood <- posterior_calc[2]
    
    if(i %% pr == 0)print(c("proposed posterior = ", proposed_posterior), quote=FALSE)
    if(i %% pr == 0)print(c("proposed likelihood = ", proposed_likelihood), quote=FALSE)
    
    likelihood_ratio <- ifelse(proposed_posterior == -Inf || is.nan(proposed_posterior), 0, exp(proposed_posterior-current_posterior))
    
    if(runif(1) < likelihood_ratio){
      if(i %% pr == 0) print("Accepted", quote=F)
      MCMC_output[i+1,] <- proposed_parameters
      MCMC_posteriors[i+1] <- proposed_posterior
      MCMC_likelihoods[i+1] <- proposed_likelihood
      Acceptances[i] <- 1
    }
    
    else{
      if(i %% pr == 0) print("Rejected", quote=F)
      MCMC_output[i+1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      ## The current state's posterior used to be recomputed with probability 0.01 to
      ## refresh the Monte Carlo noise in the likelihood estimate. With seeds now fixed
      ## the likelihood is deterministic, so re-evaluating returns the same value and the
      ## correction is redundant.
      MCMC_posteriors[i+1] <- MCMC_posteriors[i]
      MCMC_likelihoods[i+1] <- MCMC_likelihoods[i]
    }
    
    if(i %% pr == 0)print(c("Last accepted = ", max(which(Acceptances==1))), quote=F)
    if(i %% pr == 0)print(c("Acceptance rate =", round(sum(Acceptances[1:i], na.rm=T)/i,2)), quote=F)
    
    save_freq <- 5000
    
    if(i %% save_freq == 0){
      if (interactive()) plot(MCMC_posteriors[1:i,], type="l")
      Acceptances_vec <- vector(length=n_parameters)
      for(j in 1:length(Acceptances_vec)){
        Acceptances_param <- Acceptances[seq(j,i,n_parameters)]
        n_accept_param <- sum(Acceptances_param, na.rm=T)
        n_it_param <- length(Acceptances_param[is.na(Acceptances_param)==FALSE])
        Acceptances_vec[j] <- n_accept_param/n_it_param
      }
      names(Acceptances_vec) <- names(parameters)
      print(Acceptances_vec, quote=F)
      
      old_file_name <- paste0("outputs/",strain,"_",gene,"_",format(i-save_freq,scientific=F),".rds")
      file_name <- paste0("outputs/",strain,"_",gene,"_",format(i,scientific=F),".rds")
      
      output <- list(MCMC_output=MCMC_output, MCMC_posteriors=MCMC_posteriors, MCMC_likelihoods=MCMC_likelihoods, Acceptances=Acceptances, cov_matrix=cov_matrix, time_vec=time_vec, x=x)
      saveRDS(output, file=file_name)

      unlink(old_file_name)

      ## Optional live trace, off unless x$trace_file is set, so the real analysis path
      ## is unaffected. Wrapped in try() -- a plotting failure must not kill a long run.
      if(!is.null(x$trace_file)) try({
        tp <- if(is.null(x$trace_pars)) c("a_bar","b_bar","vl_max_bar","inc1","inc3","rdisp","multiplier") else x$trace_pars
        tp <- intersect(tp, colnames(MCMC_output))
        png(x$trace_file, width=1200, height=850, res=110)
        op <- par(mfrow=c(2, ceiling((length(tp)+1)/2)), mar=c(3.4,3.6,2.2,0.8), mgp=c(2.1,0.6,0))
        plot(MCMC_posteriors[1:i], type="l", xlab="iteration", ylab="log posterior",
             main=sprintf("log posterior (iter %s)", format(i, big.mark=",")), col="grey20")
        if(!is.null(cov_start) && i > cov_start) abline(v=cov_start, col="steelblue", lty=3)
        for(pn in tp){
          plot(MCMC_output[1:i, pn], type="l", xlab="iteration", ylab="", main=pn, col="grey20")
          if(!is.null(x$truth) && pn %in% names(x$truth))
            abline(h=x$truth[[pn]], col="firebrick", lwd=2, lty=2)
          if(!is.null(cov_start) && i > cov_start) abline(v=cov_start, col="steelblue", lty=3)
        }
        par(op); dev.off()
      }, silent=TRUE)
    }
    ## Per-iteration time is still recorded for later analysis; the block line above is what
    ## gets printed, so the single-iteration print here would only add noise.
    time_vec[i] <- ((proc.time()-p)["elapsed"])
  }
}


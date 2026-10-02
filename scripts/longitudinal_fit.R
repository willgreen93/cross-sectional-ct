library(ggplot2)
library(dplyr)
library(rstan)
library(matrixStats)

## Run from the repository root; every path below is relative to it. This previously began
## setwd("Q:/ClusterPaper"), which was one author's machine and made the script unrunnable
## anywhere else.
if (!file.exists("data/ct_data_clean.csv"))
  stop("run this from the repository root: data/ct_data_clean.csv not found from ", getwd())

## KISSLER DATA FROM HERE https://github.com/gradlab/CtTrajectories/blob/main/data/ct_dat_clean.csv
data1 <- read.csv("data/ct_data_clean.csv")
N1 <- nrow(data1) # number of data points
t1 <- data1$Date.Index # time points in data
individual1 <- as.numeric(as.factor(data1$Person.ID)) # individual in study #1
K1 <- length(unique(individual1)) # number of individuals
s1 <- rep(1, N1) # assign study number
voc1 <- rep(0,N1) # variant status
PersonID1 <- individual1 # individual number in aggregated data
vl1 <- data1$log10_GEperML # viral load in data

ggplot(data1 %>% filter(log10_GEperML>2.66), aes(x=Date.Index, y=log10_GEperML)) +
  geom_point(aes(color=Novel.Persistent.Infection)) +
  facet_wrap(~Person.ID, scales="free") +
  theme_bw()

## FERGUSON LALVANI DATA 
data2 <- read.csv("data/20210920_ATACCC_L.csv") %>% mutate(vl=log10(133.33*exp((37.933-CtT1)/1.418)))
N2 <- nrow(data2) # number of data points
t2 <- data2$TestDateIndex # time points in data
individual2 <- as.numeric(as.factor(data2$PersonID)) # individual in study #2
K2 <- length(unique(individual2)) # number of individuals in study 2
s2 <- rep(2, N2) # assign study number 2
voc2 <- data2$VacVOC # variant status
PersonID2 <- individual2+K1 # individual number in aggregated data
vl2 <- data2$vl # viral load in data

floor_vl <- min(data2$vl, na.rm = TRUE)
d <- data2 %>% mutate(detected = vl > floor_vl + 1e-9)

p1 <- ggplot(d, aes(TestDateIndex, vl)) +
  geom_line(aes(group = PersonID), colour = "grey65", alpha = 0.3, linewidth = 0.3) +
  geom_point(aes(shape = detected), size = 1, alpha = 0.5) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  stat_summary(fun = median, geom = "line", colour = "firebrick", linewidth = 1.1) +
  geom_hline(yintercept = floor_vl, linetype = 3) +
  labs(x = NULL, y = expression(log[10]~"viral load"),
       subtitle = "open circles = at the detection floor; median includes them") +
  theme_bw()

p2 <- d %>% group_by(TestDateIndex) %>%
  summarise(prop = mean(detected), n = n()) %>%
  ggplot(aes(TestDateIndex, prop)) +
  geom_col(fill = "grey60") +
  labs(x = "days", y = "proportion detectable") +
  theme_bw()

ggpubr::ggarrange(p1, p2, ncol = 1, heights = c(3, 1), align = "v")

## Combined data sets
N_c <- sum(N1,N2)
K_c <- sum(K1, K2)
## ct_c <- c(ct1, ct2) removed: ct1 and ct2 were never defined anywhere, and ct_c was never
## read, so this line only ever errored on a clean top-to-bottom run. The model works from
## viral load (vl_c), not Ct, so nothing downstream wanted it.
t_c <- c(t1,t2)
s_c <- c(s1,s2)
voc_c <- c(voc1, voc2)
individual_c <- c(individual1, individual2)
PersonID_c <- c(PersonID1, PersonID2)
vl_c <- c(vl1, vl2)

df_c_raw <- data.frame(individual=individual_c, vl=vl_c, t=t_c, s=s_c, voc=voc_c, PersonID=PersonID_c, study=s_c) %>% mutate(vl_min=c(min(vl1), min(vl2))[study]) 
duplicated_rows <- which(duplicated(df_c_raw[, c("individual", "t", "PersonID", "study")])==TRUE)-1
df_c = df_c_raw[-duplicated_rows,]

df_ATTAC_raw <- data.frame(individual=individual2, vl=vl2, t=t2, s=s2, voc=voc2, PersonID=PersonID2, study=s2) %>% mutate(vl_min=c(min(vl2))) 
duplicated_ATTAC_rows <- which(duplicated(df_ATTAC_raw[, c("individual", "t", "PersonID", "study")])==TRUE)-1
df_ATTAC_c = df_ATTAC_raw[-duplicated_ATTAC_rows,]


parameter_list_generator <- function(vl_df, formula){
  three_dps <- vl_df %>% dplyr::filter(vl != vl_min) %>% group_by(PersonID) %>% summarise(non_LOD_count = n()) %>% dplyr::filter(non_LOD_count >= 3) %>% dplyr::select(PersonID) %>% pull() 
  #growth_data <- vl_df %>% group_by(PersonID) %>% filter(vl != vl_min, t<0) %>% select(PersonID) %>% pull() %>% unique()
  #initial_lod <- vl_df %>% group_by(PersonID) %>% filter(t == min(t), vl==vl_min) %>% select(PersonID) %>% pull()
  
  #IDs_include <- Reduce(intersect, list(three_dps,growth_data,initial_lod))
  #IDs_include <- Reduce(intersect, list(three_dps,growth_data))
  IDs_include <- Reduce(intersect, list(three_dps))
  
  df_filtered <- vl_df #%>% #filter(PersonID %in% IDs_include) %>% 
  #mutate(PersonID = as.numeric(factor(PersonID)))
  #df_filtered <- vl_df %>% mutate(PersonID = as.numeric(factor(PersonID)))
  
  N <- nrow(df_filtered)
  K <- max(df_filtered$PersonID)
  t <- df_filtered$t
  individual <- df_filtered$PersonID
  vl <- df_filtered$vl
  vl_min <- df_filtered$vl_min
  study <- df_filtered$study
  
  k = formula
  #if(formula=="flat-corr0") k=4
  #if(formula=="flat-corr1") k=5
  #if(formula=="flat-corr2") k=6
  
  return(list(vl_df_filtered=df_filtered, N=N, K=K, formula=formula, k=k, t=t, individual=individual, vl=vl, vl_min=vl_min, study=study))
}

model_multiple_vl <- stan_model("longitudinal_fits/longitudinal_model.stan")

fitting_func <- function(stan_model, vl_df, formula, iter, chains, cores, seed=1){
  p <- proc.time()
  input_list <- parameter_list_generator(vl_df=vl_df, formula=formula)
  
  fit <- sampling(stan_model, input_list, iter=iter, chains=chains, cores=cores, seed=5, control=list(adapt_delta=0.99), refresh = 20)
  
  combined <- list(fit=fit, vl_df=input_list$vl_df_filtered, k=input_list$k, K=input_list$K, formula=formula, N=input_list$N, iter=iter, chains=chains, cores=cores, seed=1)
  
  output <- list()
  
  output$mean <- summary(fit)$summary[c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",  "l1"),c("mean")]
  output$sd   <- summary(fit)$summary[c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",  "l1"),c("sd")]
  
  combined_samples <- data.frame(extract(fit, c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma", "l1")))
  
  covariance_matrix <- cov(as.matrix(combined_samples))
  
  #out <- list(combined=combined, output=output, cov_matrix=cov_matrix)
  
  #saveRDS(out, paste0("longitudinal_fits/long_fit_", formula, ".rds"))
  
  saveRDS(combined, paste0("longitudinal_fits/long_fit_",formula,".rds"))
  saveRDS(output, file=paste0("longitudinal_fits/long_post_",formula,".rds"))
  saveRDS(covariance_matrix, paste0("longitudinal_fits/long_cov_",formula,".rds"))
  
  print(proc.time()-p)["elapsed"]
  return(combined)        
}

reject_IDS <- c(11, 18, 19, 20, 21, 22, 23, 24, 30, 
                34, 35, 38, 41, 43, 50, 51, 54, 59, 
                60, 61, 66, 67, 79, 138, 149, 167,
                197, 201, 228)

fitting_func_ATTAC <- function(stan_model, vl_df, formula, iter, chains, cores, seed=1){
  p <- proc.time()
  input_list <- parameter_list_generator(vl_df=vl_df, formula=formula)
  
  fit <- sampling(stan_model, input_list, iter=iter, chains=chains, cores=cores, seed=5, control=list(adapt_delta=0.99), refresh = 20)
  
  combined <- list(fit=fit, vl_df=input_list$vl_df_filtered, k=input_list$k, K=input_list$K, formula=formula, N=input_list$N, iter=iter, chains=chains, cores=cores, seed=1)
  
  output <- list()
  
  output$mean <- summary(fit)$summary[c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",  "l1"),c("mean")]
  output$sd   <- summary(fit)$summary[c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",  "l1"),c("sd")]
  
  combined_samples <- data.frame(extract(fit, c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma", "l1")))
  
  covariance_matrix <- cov(as.matrix(combined_samples))
  
  #out <- list(combined=combined, output=output, cov_matrix=cov_matrix)
  
  #saveRDS(out, paste0("longitudinal_fits/long_fit_", formula, ".rds"))
  
  saveRDS(combined, paste0("longitudinal_fits/long_fit_ATTAC.rds"))
  saveRDS(output, file=paste0("longitudinal_fits/long_post_ATTAC.rds"))
  saveRDS(covariance_matrix, paste0("longitudinal_fits/long_cov_ATTAC.rds"))
  
  print(proc.time()-p)["elapsed"]
  return(combined)        
}

## The two fitting frames are built BEFORE they are used. They were previously defined below
## the fit_ATTAC call, which only ever worked because the script was run out of order in an
## interactive session; sourced top to bottom it failed on a missing vl_df_ATTAC.
vl_df <- df_c %>% dplyr::filter(!PersonID %in% reject_IDS) %>% mutate(facet=floor(PersonID/10), num=PersonID%%10)
vl_df_ATTAC <- df_ATTAC_c %>% dplyr::filter(!PersonID %in% reject_IDS) %>% mutate(facet=floor(PersonID/10), num=PersonID%%10)

## ATACCC only. Writes long_post_ATTAC.rds / long_cov_ATTAC.rds, the cross-platform
## sensitivity check (the prior refitted without the NBA cohort).
fit_ATTAC <- fitting_func_ATTAC(stan_model=model_multiple_vl, vl_df=vl_df_ATTAC, formula=6, iter=2000, chains=4, cores=4, seed=1)

ggplot(vl_df, aes(x=t, y=vl, color=factor(num))) +
  geom_line() +
  facet_wrap(~factor(facet), scales="free") +
  theme_bw() +
  theme(panel.grid=element_blank())

fit8.6 <- fitting_func(stan_model=model_multiple_vl, vl_df=vl_df, formula=6, iter=2000, chains=4, cores=4, seed=1)
fit8.7 <- fitting_func(stan_model=model_multiple_vl, vl_df=vl_df, formula=7, iter=2000, chains=4, cores=4, seed=1)

vl_func <- function(a, b, tmax, t, l=0, vl_max, k, vl_min){
  if(k==1) out <- ifelse(t < tmax, vl_max+exp(a)*(t-tmax), vl_max-exp(b)*(t-tmax))
  if(k==2) out <- ifelse(t < tmax, vl_max+exp(a)*(t-tmax), vl_max-b/exp(a)*(t-tmax))
  if(k==3) out <- ifelse(t < tmax, vl_max*exp(a)+exp(a)*(t-tmax), vl_max*exp(a)-b/exp(a)*(t-tmax))
  if(k==4) out <- ifelse(t < tmax, vl_max+vl_max/exp(a)*(t-tmax), vl_max-exp(b)/vl_max*(t-tmax))
  
  if(k==6) out <- ifelse(t < tmax, vl_max+(0.5+exp(a))*(t-tmax), vl_max-(0.25+exp(b))*(t-tmax))
  if(k==7) out <- ifelse(t < tmax-0.5*exp(l), vl_max+(0.5+exp(a))*(t-tmax+0.5*exp(l)),
                         ifelse(t < tmax+0.5*exp(l), vl_max,
                                vl_max-(0.25+exp(b))*(t-tmax-0.5*exp(l))))
  #if(k==4) out <- ifelse(t < tmax, vl_max+exp(a)*(t-tmax), ifelse(t < tmax+l, vl_max, vl_max*exp(a)-exp(b)*(t-tmax-l)))
  #if(k==5) out <- ifelse(t < tmax, vl_max+exp(a)*(t-tmax), ifelse(t < tmax+l, vl_max, vl_max-b/exp(a)*(t-tmax-l)))
  #if(k==6) out <- ifelse(t < tmax, vl_max*exp(a)+exp(a)*(t-tmax), ifelse(t < tmax+l, vl_max*exp(a), vl_max*exp(a)-b/exp(a)*(t-tmax-l)))
  
  return(pmax(vl_min,out))
} 

output_plot <- function(fitting_output, nsamples){
  stan_fit <- fitting_output$fit
  fit_input <- fitting_output$vl_df
  iter <- fitting_output$iter/2
  indiv <- fitting_output$K
  unique_IDS = fit_input$PersonID %>% unique()
  
  samples <- extract(fitting_output$fit)
  
  a_samples <- (as.vector(samples$a_bar) + as.matrix(samples$a_raw)*as.vector(samples$a_sigma))[seq(iter/nsamples,iter, length.out = nsamples),]
  b_samples <- (as.vector(samples$b_bar) + as.matrix(samples$b_raw)*as.vector(samples$b_sigma))[seq(iter/nsamples,iter, length.out = nsamples),]
  l_samples <- as.vector(samples$l1)
                
  vl_max_samples <- (as.vector(samples$vl_max_bar) + as.matrix(samples$vl_max_raw)*as.vector(samples$vl_max_sigma))[seq(iter/nsamples,iter, length.out = nsamples),]
  tmax_samples <- (as.matrix(samples$tmax))[seq(iter/nsamples,iter, length.out = nsamples),]
  
  for(k in unique_IDS){
    vl_min <- fitting_output$vl_df %>% filter(PersonID==k) %>% dplyr::select(vl_min) %>% unique() %>% pull()
    individual_hold <- matrix(nrow=301, ncol=nsamples)
    
    for(i in 1:nsamples) individual_hold[,i] <- vl_func(a=a_samples[i,k], b=b_samples[i,k], tmax=tmax_samples[i,k], t=seq(-10,20,0.1), l=l_samples[k], vl_max=vl_max_samples[i,k], k=fitting_output$k, vl_min=vl_min)
    
    individual_range_raw <- data.frame(t=seq(-10,20,0.1),
                                       rowQuantiles(individual_hold, probs=c(0.025, 0.5, 0.975)),
                                       PersonID=k) %>% magrittr::set_colnames(c("t", "lower", "median", "upper", "PersonID")) 
    
    if(k==1) individual_range <- individual_range_raw
    else individual_range <- rbind(individual_range, individual_range_raw)
  }
  
  plotting_output <- left_join(individual_range, fit_input %>% dplyr::select(PersonID, t, vl), by=c("PersonID", "t")) 
  
  p <- ggplot(plotting_output, aes(x=t)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
    geom_line(aes(y=median)) +
    geom_point(aes(y=vl)) +
    facet_wrap(~PersonID, nrow=9) +
    theme_bw() +
    labs(y=bquote(log[10](vl)), x="Time from peak (days)") + 
    theme(axis.line = element_line(colour = "black"),
          #panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position="none",
          strip.background = element_blank())#,
          #strip.text.x = element_blank())
  
  ggsave(filename=paste0("longitudinal_fits/viral_kinetics_fits_",fitting_output$k,".png"), plot = p, device = png, width=7.5, height=5)
  
  print(p)
  
  #ggsave(filename=paste0("viral_kinetics_long_params",fitting_output$k,".png"), plot=q, device = png, width=7.5, height=5)

}

output_plot(fitting_output=fit8.6, nsamples=100)
output_plot(fitting_output=fit8.7, nsamples=100)
#output_plot(fitting_output=fit8.5, nsamples=100)
#output_plot(fitting_output=fit8.6, nsamples=100)

cov_prior_generator <- function(flat, corr){
  files <- list.files("longitudinal_fits/")
  fit_files <- grepl("fit", files)
  flatness <- grepl("flat", files)
  correlation <- grepl(paste0("corr",corr), files)
  if(flat==0) flatness = !flatness
  
  correct_file <- paste0("longitudinal_fits/",files[which(fit_files*flatness*correlation==1)])
  print(correct_file)
  
  fitting_output <- readRDS(correct_file)$fit
  
  combined_samples <- data.frame(extract(fitting_output, c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma", "l_bar", "l_sigma")))
  
  covariance_matrix <- cov(as.matrix(combined_samples))
  
  flat_text <- ifelse(flat==0, "", ifelse(flat==1, "flat-"))
  
  saveRDS(covariance_matrix, paste0("longitudinal_fits/longitudinal_cov_",flat_text,"corr",corr,"_final.rds"))
  
  return(covariance_matrix)
}
#
#cov_prior_generator(flat=0, corr=0)
#cov_prior_generator(flat=1, corr=0)
#cov_prior_generator(flat=0, corr=1)
#cov_prior_generator(flat=1, corr=1)
#cov_prior_generator(flat=0, corr=2)
#cov_prior_generator(flat=1, corr=2)

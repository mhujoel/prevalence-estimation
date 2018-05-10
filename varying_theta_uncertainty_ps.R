library(rjags)
library(rmeta)
library(truncnorm)
library(lme4)
require(ggplot2)
require(gridExtra)
#  --------------------------------- DATA GENERATION --------------------------------- #
simulationPopulation_fixedTheta_ratio <- function(ns,p_ratio,theta){
  # simulate one population with a given fixed theta (population prevalence)
  # ratio of probabilities of being in the study given carrier status (pplus)
  # over the probability given non carrier status (pminus). 
  # ns gives the study size
  # output: a vector with 4 numbers:
  # number of carriers, number in the study, 
  # the prevalence used to generate the study, the ratio: p+/p-
  alpha <- p_ratio*(theta/(p_ratio*theta+1-theta))
  sample_gen <- rbinom(ns,size=1,prob=alpha)
  num <- length(sample_gen)
  carriers <- sum(sample_gen)
  return(c(carriers,num,theta,p_ratio))
}

studyData_fixed <- function(N,study_size,seed,theta, pRatioRange,pRatio_fixed = TRUE, pRatio_value){
  # overview: a function that will generate data from N studies
  # inputs:
  # N: number of studies
  # study_size: size of each study we want
  # seed: the seed we want (allows for replication)
  # theta: one theta value <- the population prevalence 
  # pRatioRange: the range p+/p- can be in our data generations
  # pRatio_fixed: TRUE: given ratio is for each study otherwise we generate random p_ratios
  # pRatio_value: the ratio if pRatio_fixed == TRUE
  # output: a data frame with the following
  # carriers: a vector containing the carriers in each study
  # study_sizes: a vector with the number in each study
  # p_ratios: a vector with the p+/p- ratio for each study
  # theta_s: a vector with the true prevalence value for each study 
  # N: number of studies
  set.seed(seed)
  carriers_sample <- c()
  p_ratios <- c()
  theta_s <- c()
  study_sizes <- c()
  for (i in 1:N) {
    pRatio_val <- ifelse(pRatio_fixed, pRatio_value[i], 
                         runif(n = 1, min = pRatioRange[1], max = pRatioRange[2]))
    y = simulationPopulation_fixedTheta_ratio(ns = study_size[i], p_ratio = pRatio_val,theta = theta)
    carriers_sample <- append(carriers_sample, y[1])
    study_sizes <- append(study_sizes, y[2])
    p_ratios <- append(p_ratios, y[4])
    theta_s <- append(theta_s, y[3])
  }
  return(list("carriers" = carriers_sample, "study_sizes" = study_sizes, "p_ratios"=p_ratios, "theta_s" = theta_s,"N"=N))
}

studyData_varyingTheta <- function(N,study_size,seed,theta, pRatio_values){
  # overview: a function that will generate data from N studies
  # inputs:
  # N: number of studies
  # study_size: size of each study we want
  # seed: the seed we want (allows for replication)
  # theta: a vector of N thetas <- the population prevalence underlying each study
  # pRatio_values: the p_ratios for the studies
  # output: a data frame with the following
  # carriers: a vector containing the carriers in each study
  # study_sizes: a vector with the number in each study
  # p_ratios: a vector with the p+/p- ratio for each study
  # theta_s: a vector with the true prevalence value for each study 
  # N: number of studies
  set.seed(seed)
  carriers_sample <- c()
  p_ratios <- c()
  theta_s <- c()
  study_sizes <- c()
  for (i in 1:N) {
    y = simulationPopulation_fixedTheta_ratio(ns = study_size[i], p_ratio =  pRatio_values[i], theta = theta[i])
    carriers_sample <- append(carriers_sample, y[1])
    study_sizes <- append(study_sizes, y[2])
    p_ratios <- append(p_ratios, y[4])
    theta_s <- append(theta_s, y[3])
  }
  return(list("carriers" = carriers_sample, "study_sizes" = study_sizes, "p_ratios"=p_ratios, "theta_s" = theta_s,"N"=N))
}
#  --------------------------------- RANDOM THETA - MISSPECIFICATION --------------------------------- #

simulation_varying_theta_varyingPratio_missppecified <- function(seed.val=as.numeric(Sys.time()),theta_simulation_vals,
                                                                 num_sims=100,num_studies=10,
                                                                 study_sizes_vec=rep(5000,num_studies),
                                                                 p_Ratios_vec = rep(30,num_studies), 
                                                                 p_ratios_se_implementation = 1,
                                                                 mean_misspecification_factor=rep(1,num_studies),
                                                                 n.adapt.bayes=5000,n.iter.bayes=10000){
  set.seed(seed.val)
  
  bayesian.theta.beta.fixed <- rep(NA,num_sims)
  ci.bayesian.beta.fixed <- matrix(NA,nrow = num_sims,ncol=2)

  bayesian.theta.uniform.fixed <- rep(NA,num_sims)
  ci.bayesian.uniform.fixed <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian.theta.logitNormal.fixed <- rep(NA,num_sims)
  ci.bayesian.logitNormal.fixed <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian.theta.beta.uncertain <- rep(NA,num_sims)
  ci.bayesian.beta.uncertain <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian.theta.uniform.uncertain <- rep(NA,num_sims)
  ci.bayesian.uniform.uncertain <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian.theta.logitNormal.uncertain <- rep(NA,num_sims)
  ci.bayesian.logitNormal.uncertain <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian_model_beta ="model {
  for( i in 1 : Num ) 
  {
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i]))); 
  theta[i] ~ dbeta(a,b);
  }
  a = (((1-mu)/sigmaSq) - (1/mu))*(mu^2);
  b = a*((1/mu)-1);
  mu ~ dunif(0,0.01);
  sigmaSq ~ dunif(0,mu*(1-mu));
}
"

bayesian_model_uniform ="model {
for( i in 1 : Num ) 
{
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i]))); 
  theta[i] ~ dunif(mean-0.5*width,mean+0.5*width);
}
mean ~ dunif(0,0.01);
width ~ dunif(0,2*mean);
}
"

bayesian_model_transformed ="model {
          for( i in 1 : Num )
{
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])));
  theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
  logittheta[i] ~ dnorm(mu,tau);
}
mu ~ dunif(-50,0);
tau ~ dunif(1,50);
}
"

  bayesian_model_beta_pratio_uncertainty ="model {
  for( i in 1 : Num ) 
  {
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i]))); 
  p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0,);
  theta[i] ~ dbeta(a,b);
  }
  a = (((1-mu)/sigmaSq) - (1/mu))*(mu^2);
  b = a*((1/mu)-1);
  mu ~ dunif(0,0.01);
  sigmaSq ~ dunif(0,mu*(1-mu));
}
"
bayesian_model_uniform_pratio_uncertainty ="model {
for( i in 1 : Num ) 
{
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])));
  p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0,);
  theta[i] ~ dunif(mean-0.5*width,mean+0.5*width);
}
mean ~ dunif(0,0.01);
width ~ dunif(0,2*mean);
}
"

bayesian_model_transformed_pratio_uncertainty ="model {
          for( i in 1 : Num )
{
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])));
  p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0,);
  theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
  logittheta[i] ~ dnorm(mu,tau);
}
mu ~ dunif(-50,0);
tau ~ dunif(1,50);
}
"


for (j in 1:num_sims){
  s = seed.val + j # set our seed
  
  x = studyData_varyingTheta(N=num_studies,study_size=study_sizes_vec,seed=s,
                                           theta=theta_simulation_vals, pRatio_values = p_Ratios_vec)
  
  seed.values = abs(.Random.seed[1:3])+s
  
  data_fixed=list(x=x$carriers,n=x$study_sizes, p_ratio=p_Ratios_vec*mean_misspecification_factor, Num=x$N)
  initial_values <- list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper"))
  
  # BAYESIAN BETA PRIOR - FIXED
  model_bayesian_beta_fixed =jags.model(textConnection(bayesian_model_beta), data=data_fixed, n.chains = 3,
                                  inits=initial_values, quiet = TRUE, n.adapt= 5000)
  mcmc_output_beta_fixed = coda.samples(model_bayesian_beta_fixed,c("theta"), n.iter = 10000,thin = 10,progress.bar="none")
  
  median_chains_bayesian_beta_fixed <- c(apply(mcmc_output_beta_fixed[[1]][,1:num_studies],1,median),
                                   apply(mcmc_output_beta_fixed[[2]][,1:num_studies],1,median),
                                   apply(mcmc_output_beta_fixed[[3]][,1:num_studies],1,median))
  bayesian.theta.beta.fixed[j] = mean(median_chains_bayesian_beta_fixed)
  ci.bayesian.beta.fixed[j,] = quantile(median_chains_bayesian_beta_fixed, probs = c(0.025,0.975))
  
  # BAYESIAN UNIFORM APPROACH - FIXED
  model_bayesian_uni_fixed =jags.model(textConnection(bayesian_model_uniform), data=data_fixed, n.chains = 3,
                                 inits=initial_values,quiet = TRUE,n.adapt= 5000)
  mcmc_output_uni_fixed = coda.samples(model_bayesian_uni_fixed,c("theta"), n.iter = 10000, thin = 10,progress.bar="none")
  
  median_chains_bayesian_uni_fixed <- c(apply(mcmc_output_uni_fixed[[1]][,1:num_studies],1,median),
                                  apply(mcmc_output_uni_fixed[[2]][,1:num_studies],1,median),
                                  apply(mcmc_output_uni_fixed[[3]][,1:num_studies],1,median))

  bayesian.theta.uniform.fixed[j] = mean(median_chains_bayesian_uni_fixed)
  ci.bayesian.uniform.fixed[j,] = quantile(median_chains_bayesian_uni_fixed, probs = c(0.025,0.975))
  # BAYESIAN LOGIT NORMAl APPROACH - FIXED
  model_bayesian_logitN_fixed =jags.model(textConnection(bayesian_model_transformed), data=data_fixed, n.chains = 3,
                                       inits=initial_values,quiet = TRUE,n.adapt= 5000)
  mcmc_output_logitN_fixed = coda.samples(model_bayesian_logitN_fixed,c("theta"), n.iter = 10000, thin = 10,progress.bar="none")
  
  median_chains_bayesian_logitN_fixed <- c(apply(mcmc_output_logitN_fixed[[1]][,1:num_studies],1,median),
                                        apply(mcmc_output_logitN_fixed[[2]][,1:num_studies],1,median),
                                        apply(mcmc_output_logitN_fixed[[3]][,1:num_studies],1,median))
  
  bayesian.theta.logitNormal.fixed[j] = mean(median_chains_bayesian_logitN_fixed)
  ci.bayesian.logitNormal.fixed[j,] = quantile(median_chains_bayesian_logitN_fixed, probs = c(0.025,0.975))
  # DATA WITH UNCERTAINTY
  se_vec <- rep(p_ratios_se_implementation,num_studies)
  data=list(x=x$carriers,n=x$study_sizes, p_ratio_mean=p_Ratios_vec*mean_misspecification_factor,
            p_ratio_sd=se_vec, Num=x$N)
  # BAYESIAN BETA PRIOR - VARIANCE
  model_bayesian_beta =jags.model(textConnection(bayesian_model_beta_pratio_uncertainty), data=data, n.chains = 3,
                                  inits=initial_values, quiet = TRUE,n.adapt= n.adapt.bayes)
  mcmc_output_beta = coda.samples(model_bayesian_beta,c("theta"), n.iter = n.iter.bayes , thin = 10,progress.bar="none")
  
  median_chains_bayesian_beta <- c(apply(mcmc_output_beta[[1]][,1:num_studies],1,median),
                              apply(mcmc_output_beta[[2]][,1:num_studies],1,median),
                              apply(mcmc_output_beta[[3]][,1:num_studies],1,median))
  # store what we want to keep
  bayesian.theta.beta.uncertain[j] = mean(median_chains_bayesian_beta)
  ci.bayesian.beta.uncertain[j,] = quantile(median_chains_bayesian_beta, probs = c(0.025,0.975))
  # BAYESIAN UNIFORM APPROACH - VARIANCE
  model_bayesian_uni =jags.model(textConnection(bayesian_model_uniform_pratio_uncertainty), data=data, n.chains = 3,
                                 inits=initial_values, quiet = TRUE,n.adapt= n.adapt.bayes)
  mcmc_output_uni = coda.samples(model_bayesian_uni,c("theta"), n.iter = n.iter.bayes,
                                 thin = 10,progress.bar="none")
  
  median_chains_bayesian_uni <- c(apply(mcmc_output_uni[[1]][,1:num_studies],1,median),
                                   apply(mcmc_output_uni[[2]][,1:num_studies],1,median),
                                   apply(mcmc_output_uni[[3]][,1:num_studies],1,median))
  # store what we want to keep
  bayesian.theta.uniform.uncertain[j] = mean(median_chains_bayesian_uni)
  ci.bayesian.uniform.uncertain[j,] = quantile(median_chains_bayesian_uni, probs = c(0.025,0.975))
  # BAYESIAN LOGIT NORMAL APPROACH - VARIANCE
  model_bayesian_logitN =jags.model(textConnection(bayesian_model_transformed_pratio_uncertainty), data=data, n.chains = 3,
                                 inits=initial_values, quiet = TRUE,n.adapt= n.adapt.bayes)
  mcmc_output_logitN = coda.samples(model_bayesian_logitN,c("theta"), n.iter = n.iter.bayes,
                                 thin = 10,progress.bar="none")
  
  median_chains_bayesian_logitN <- c(apply(mcmc_output_logitN[[1]][,1:num_studies],1,median),
                                  apply(mcmc_output_logitN[[2]][,1:num_studies],1,median),
                                  apply(mcmc_output_logitN[[3]][,1:num_studies],1,median))
  # store what we want to keep
  bayesian.theta.logitNormal.uncertain[j] = mean(median_chains_bayesian_logitN)
  ci.bayesian.logitNormal.uncertain[j,] = quantile(median_chains_bayesian_logitN, probs = c(0.025,0.975))
  
}
# true median of the prevalences:
true_underlying_median <- median(theta_simulation_vals)
# FIXED
outside.bayesian.beta.fixed <- length(which(true_underlying_median > ci.bayesian.beta.fixed[,2]))+ length(which(true_underlying_median < ci.bayesian.beta.fixed[,1]))
coverage.bayesian.beta.fixed  <- (num_sims-outside.bayesian.beta.fixed)/num_sims
outside.bayesian.uniform.fixed  <- length(which(true_underlying_median > ci.bayesian.uniform.fixed[,2]))+ length(which(true_underlying_median < ci.bayesian.uniform.fixed[,1]))
coverage.bayesian.uniform.fixed  <- (num_sims-outside.bayesian.uniform.fixed)/num_sims
outside.bayesian.logitN.fixed  <- length(which(true_underlying_median > ci.bayesian.logitNormal.fixed[,2]))+ length(which(true_underlying_median < ci.bayesian.logitNormal.fixed[,1]))
coverage.bayesian.logitN.fixed  <- (num_sims-outside.bayesian.logitN.fixed)/num_sims
# UNCERTAINTY
outside.bayesian.beta.uncertain <- length(which(true_underlying_median > ci.bayesian.beta.uncertain[,2]))+ length(which(true_underlying_median < ci.bayesian.beta.uncertain[,1]))
coverage.bayesian.beta.uncertain <- (num_sims-outside.bayesian.beta.uncertain)/num_sims
outside.bayesian.uniform.uncertain <- length(which(true_underlying_median > ci.bayesian.uniform.uncertain[,2]))+ length(which(true_underlying_median < ci.bayesian.uniform.uncertain[,1]))
coverage.bayesian.uniform.uncertain <- (num_sims-outside.bayesian.uniform.uncertain)/num_sims
outside.bayesian.logitN.uncertain <- length(which(true_underlying_median > ci.bayesian.logitNormal.uncertain[,2]))+ length(which(true_underlying_median < ci.bayesian.logitNormal.uncertain[,1]))
coverage.bayesian.logitN.uncertain <- (num_sims-outside.bayesian.logitN.uncertain)/num_sims

return(list(trueMedian = true_underlying_median,
            beta_coverage_fixed = coverage.bayesian.beta.fixed, beta_coverage_uncertain = coverage.bayesian.beta.uncertain,
            uniform_coverage_fixed = coverage.bayesian.uniform.fixed, uniform_coverage_uncertain = coverage.bayesian.uniform.uncertain, 
            logitN_coverage_fixed = coverage.bayesian.logitN.fixed, logitN_coverage_uncertain = coverage.bayesian.logitN.uncertain, 
            ci_beta_fixed = ci.bayesian.beta.fixed, ci_beta_uncertain= ci.bayesian.beta.uncertain,  
            theta.beta_fixed = bayesian.theta.beta.fixed, theta.beta_uncertain = bayesian.theta.beta.uncertain,
            ci_uniform_fixed = ci.bayesian.uniform.fixed, ci_uniform_uncertain = ci.bayesian.uniform.uncertain, 
            theta.uniform_fixed = bayesian.theta.uniform.fixed,theta.uniform_uncertain = bayesian.theta.uniform.uncertain,
            ci_logitN_fixed = ci.bayesian.logitNormal.fixed, ci_logitN_uncertain = ci.bayesian.logitNormal.uncertain, 
            theta.logitN_fixed = bayesian.theta.logitNormal.fixed,theta.logitN_uncertain = bayesian.theta.logitNormal.uncertain))
}

num_studies = 10
num_simulations = 100

set.seed(20030314)
theta_simuls_larger <- runif(n = num_studies,min = 0.0009,max = 0.0011)
mean_vec = seq(10,30,length.out = num_studies)

set.seed(19610527)
theta_simuls_smaller <- runif(n = num_studies,min = 0.00009,max = 0.00011)
mean_vec_smaller_vec = seq(30,50,length.out = num_studies)

# Uncertainty in the p_s estimates
uncertainty_simulation <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 74,theta_simulation_vals=theta_simuls_larger,
                                                                               p_Ratios_vec = mean_vec)
saveRDS(uncertainty_simulation,file="varying_theta_scenario4.rds")
# redo with the smaller theta values to get their respective plots:
uncertainty_simulation <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 511,theta_simulation_vals=theta_simuls_smaller,
                                                                               p_Ratios_vec = mean_vec_smaller_vec)
saveRDS(uncertainty_simulation,file="varying_theta_scenario6.rds")
#  -------- COVERAGE: DIFFERENT DEGREES OF MISSPECIFICATION ----------
little_vec = rep(0.1,num_studies)
one_above <- (mean_vec+2.5*little_vec)/(mean_vec)
one_below <- (mean_vec-2.5*little_vec)/(mean_vec)
two_above <- (mean_vec+5*little_vec)/(mean_vec)
two_below <- (mean_vec-5*little_vec)/(mean_vec)
three_above <- (mean_vec+7.5*little_vec)/(mean_vec)
three_below <- (mean_vec-7.5*little_vec)/(mean_vec)
four_above <- (mean_vec+10*little_vec)/(mean_vec)
four_below <- (mean_vec-10*little_vec)/(mean_vec)
# no misspecification is scenario 4: above
# misspecify 0.25 above the mean
one_unit_above <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 690,theta_simulation_vals=theta_simuls_larger,
                                                                     p_Ratios_vec = mean_vec,
                                                     mean_misspecification_factor=one_above)
saveRDS(one_unit_above,file="varying_theta_scenario4_misspecified1above.rds")
# misspecify 0.5 above the mean
two_unit_above <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 1230,theta_simulation_vals=theta_simuls_larger,
                                                                     p_Ratios_vec = mean_vec, 
                                                     mean_misspecification_factor=two_above)
saveRDS(two_unit_above,file="varying_theta_scenario4_misspecified2above.rds")
# misspecify 0.75 above the mean : 
three_unit_above <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 573,theta_simulation_vals=theta_simuls_larger,
                                                                       p_Ratios_vec = mean_vec,
                                                                       mean_misspecification_factor=three_above)
saveRDS(three_unit_above,file="varying_theta_scenario4_misspecified3above.rds")
# misspecify 1 above the mean
four_unit_above <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 314,theta_simulation_vals=theta_simuls_larger,
                                                                       p_Ratios_vec = mean_vec, 
                                                                       mean_misspecification_factor=four_above)
saveRDS(four_unit_above,file="varying_theta_scenario4_misspecified4above.rds")
# misspecify 0.25 below the mean
one_unit_below <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 313,theta_simulation_vals=theta_simuls_larger,
                                                                     p_Ratios_vec = mean_vec, 
                                                     mean_misspecification_factor=one_below)
saveRDS(one_unit_below,file="varying_theta_scenario4_misspecified1below.rds")
# misspecify 0.5 below the mean
two_unit_below <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 69,theta_simulation_vals=theta_simuls_larger,
                                                                     p_Ratios_vec = mean_vec, 
                                                     mean_misspecification_factor=two_below)
saveRDS(two_unit_below,file="varying_theta_scenario4_misspecified2below.rds")
# misspecify 0.75 below the mean
three_unit_below <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 19,theta_simulation_vals=theta_simuls_larger,
                                                                       p_Ratios_vec = mean_vec,
                                                                       mean_misspecification_factor=three_below)
saveRDS(three_unit_below,file="varying_theta_scenario4_misspecified3below.rds")
# misspecify 1 below the mean
four_unit_below <- simulation_varying_theta_varyingPratio_missppecified(seed.val = 822,theta_simulation_vals=theta_simuls_larger,
                                                                       p_Ratios_vec = mean_vec, 
                                                                       mean_misspecification_factor=four_below)
saveRDS(four_unit_below,file="varying_theta_scenario4_misspecified4below.rds")





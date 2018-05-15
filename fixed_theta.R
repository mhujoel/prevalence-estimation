library(rjags)
library(rmeta)
library(truncnorm)
library(lme4)
require(ggplot2)
require(gridExtra)
library(boot)
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



#  --------------------------------- FREQUENTIST THEORY --------------------------------- #
likelihood.fcn <- function(theta,p.ratios,study.sizes,carriers){
  # function used to generate the value of the likelihood ratio for a given theta,
  # a vector of the ascertainment probability ratios, the study sizes and the carriers
  # it returns the derivative of the log likelihood 
  alphas = p.ratios * theta / ((p.ratios * theta) + 1 - theta)
  deriv.alpha = p.ratios * (p.ratios*theta + 1-theta)^(-2)
  deriv.log.likeli = sum(carriers*(1/alphas)*deriv.alpha-(study.sizes-carriers)*(1/(1-alphas))*deriv.alpha)
  return(deriv.log.likeli)
} 

mle.maximizer <- function(carriers,study.sizes,p.ratios){
  # function that returns the root of the derivative of the log likelihood: or equivalently,
  # the theta which maximizes the likelihood
  mle = uniroot(likelihood.fcn, p.ratios = p.ratios, study.sizes=study.sizes,
                carriers = carriers, lower=0.0000000001, upper=0.99,tol=.Machine$double.eps)$root
  return(mle)
}

fisher.info <- function(theta,p.ratios,study.sizes,carriers){
  # function that calculates the fisher information
  alphas = p.ratios * theta / ((p.ratios * theta) + 1 - theta)
  deriv.alpha = p.ratios * (p.ratios*theta + 1-theta)^(-2)
  fisher.info <- sum(study.sizes*deriv.alpha^2*(alphas*(1-alphas))^(-1))
  return(fisher.info)
}

ci.mle.fisher <- function(theta,fisher.information,alpha.val=0.05){
  # function that given the fisher information constructs a log-transformed
  # CI
  log.variance <- 1/(fisher.information*theta^2)
  log.ci.hi <- log(theta) + qnorm(1-(alpha.val/2))*sqrt(log.variance)
  log.ci.lo <- log(theta) + qnorm(alpha.val/2)*sqrt(log.variance)
  
  hi.ci <- exp(log.ci.hi)
  lo.ci <- exp(log.ci.lo)
  
  return(c(lo.ci,hi.ci))
}

ci.mle.transform <- function(theta,var,alpha.val=0.05){
  # function that given the variance constructs a log-transformed
  # CI
  log.variance <- var*(1/theta^2)
  log.ci.hi <- log(theta) + qnorm(1-(alpha.val/2))*sqrt(log.variance)
  log.ci.lo <- log(theta) + qnorm(alpha.val/2)*sqrt(log.variance)
  
  hi.ci <- exp(log.ci.hi)
  lo.ci <- exp(log.ci.lo)
  
  return(c(lo.ci,hi.ci))
}
#  --------------------------------- FIXED THETA --------------------------------- #
#  --------------------------------- FIGURE 1 --------------------------------- #
simulation_fixed_theta <- function(seed.val=as.numeric(Sys.time()),theta_true,num_sims=100,num_studies=10,
                                   study_sizes_vec=rep(5000,num_studies),p_Ratios_vec = rep(30,num_studies),
                                   n.adapt.bayes = 5000,n.iter.bayes = 10000){
  # seed.val: the seed, allows replicatability
  # theta_true: the true underlying theta from which our study data is drawn
  # num_sims: number of simulations we wish to conduct 
  #           default is 100
  # num_studies: number of studies per simulations
  #             default is 10
  # study_sizes_vec: a vector the length of the number of studies containing the study sizes wanted per study
  #             default is same sample size and 5,000
  # p_Ratios_vec: the vector of underlying ascertainment probability ratios for each study
  #           default is 30 for each study
  set.seed(seed.val)
  
  mle.theta <- rep(NA,num_sims)
  ci.mle <- matrix(NA,nrow = num_sims,ncol=2)
  bayesian.theta <- rep(NA,num_sims)
  ci.bayesian <- matrix(NA,nrow = num_sims,ncol=2)
  
  bayesian_model ="model {
  for( i in 1 : Num ) 
  {
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta/(p_ratio[i]*theta+(1-theta))); 
  }
  theta ~ dunif(0,1);
}
"

for (j in 1:num_sims){
  s = seed.val + j
  # generate study data
  x = studyData_fixed(N=num_studies, study_size = study_sizes_vec, seed=s,
                      theta=theta_true, pRatio_fixed = TRUE, pRatio_value = p_Ratios_vec)
  # Frequentist Approach 
  mle.theta[j] = mle.maximizer(x$carriers,x$study_sizes,x$p_ratios)
  fish.info <- fisher.info(mle.theta[j],x$p_ratios,x$study_sizes,x$carriers)
  ci.mle[j,] = ci.mle.fisher(mle.theta[j],fish.info)
  # Bayesian Approach 
  seed.values = abs(.Random.seed[1:3])+s
  data=list(x=x$carriers,n=x$study_sizes, p_ratio=x$p_ratios, Num=x$N)
  initial_values <- list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper"))
  
  model_bayesian =jags.model(textConnection(bayesian_model), data=data, n.chains = 3,
                             inits= initial_values, quiet = TRUE, n.adapt = n.adapt.bayes)
  
  mcmc_output = coda.samples(model_bayesian,c("theta"), n.iter = n.iter.bayes,
                              thin = 10,progress.bar="none")
  
  summary <- summary(mcmc_output)
  bayesian.theta[j] = summary$statistics[1]
  ci.bayesian[j,] = summary$quantiles[c(1,5)]
}

df.mle <- data.frame(x =1:num_sims, MLE = mle.theta, L.mle = ci.mle[,1], U.mle = ci.mle[,2])
df.bayesian <- data.frame(x =1:num_sims, BayesianEst = bayesian.theta, L.bayes = ci.bayesian[,1], U.bayes = ci.bayesian[,2])

return(c(df.mle,df.bayesian))
}
# In order to generate Figure 1: one uses the above function with theta_true = 0.001 or 0.0001 and the seed value 11273014.
theta_low = 0.001
theta_extra_low = 0.0001
# use function to generate simulation data
fig_1_low <- simulation_fixed_theta(seed.val=11,theta_true=theta_low)
fig_1_extra_low <- simulation_fixed_theta(seed.val=27,theta_true=theta_extra_low)
saveRDS(fig_1_low,file="Data/fixed_theta_scenario1.rds")
saveRDS(fig_1_extra_low,file="Data/fixed_theta_scenario2.rds")



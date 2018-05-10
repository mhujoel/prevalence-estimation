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

#  --------------------------------- FREQUENTIST THEORY --------------------------------- #
mle.one.study <- function(carriers,study.sizes,p.ratios){
  # return the mle for the one study
  mle = carriers / (study.sizes*p.ratios - p.ratios + 1)
  return(mle)
}

fisher.info <- function(theta,p.ratios,study.sizes,carriers){
  # function that calculates the fisher information
  alphas = p.ratios * theta / ((p.ratios * theta) + 1 - theta)
  deriv.alpha = p.ratios * (p.ratios*theta + 1-theta)^(-2)
  fisher.info <- sum(study.sizes*deriv.alpha^2*(alphas*(1-alphas))^(-1))
  return(fisher.info)
}

likelihood_added_beta <- function(theta_and_hyper, p_studies, num_carriers, study_sizes) {
  # theta will be a vector, as will y 
  hyper = theta_and_hyper[1:2]
  theta = theta_and_hyper[-c(1,2)]
  num_Studies <- length(p_studies)
  alpha_studies <- (p_studies*theta)/((p_studies*theta)+1-theta) 
  if (min(theta_and_hyper) > 0 & max(theta) < 1 & min(alpha_studies) > 0){
    logL <- (-num_Studies)*lbeta(hyper[1], hyper[2]) +
      sum(num_carriers*log(alpha_studies)+
            (study_sizes-num_carriers)*log(1-alpha_studies)+
            (hyper[1]-1)*log(theta)+
            (hyper[2]-1)*log(1-theta))
    negLogL <- -1 * logL
  } else{
    negLogL <- Inf
  }
  return(negLogL)
}


# function to obtain median from RE-beta and CI from the data 
bootstrap_re_beta <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample 
  num_s <- length(d$carriers)
  theta.start = rep(0,num_s)
  for (i in 1:length(theta.start)){
    if (d$carriers[i] > 0){
      theta.start[i] = mle.one.study(carriers = d$carriers[i],study.sizes = d$study_sizes[i],
                                     p.ratios = d$p_ratios[i]) 
    } else {
      theta.start[i] = 0.000001
    }
  }
  if(max(theta.start)==min(theta.start)){
    theta.start[which.max(theta.start)] = max(theta.start) + 0.00001
  }
  theta_start_hyper = c(1,1,theta.start)
  mle.optim <- optim(par = theta_start_hyper, 
                     fn = likelihood_added_beta, 
                     p_studies=d$p_ratios, 
                     num_carriers=d$carriers,
                     study_sizes = d$study_sizes)
  estimate_median <- median(mle.optim$par[-c(1,2)])
  return(estimate_median)
} 

# function to obtain weighted median and CI from the data 
bootstrap_median <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample 
  estimate <- meta.summaries(d=d$yj,se = sqrt(d$vj),method = "random")
  estimate_median <- matrixStats::weightedMedian(estimate$effects,estimate$weights)
  return(estimate_median)
} 

#  --------------------------------- RANDOM THETA --------------------------------- #
#  ----------------------FIGURE 2/3/4 and TABLE 1 --------------------------------- #

simulation_varying_theta <- function(seed.val=as.numeric(Sys.time()),theta_simulation_vals,num_sims=100,num_studies=10,
                                     study_sizes_vec=rep(5000,num_studies),p_Ratios_vec = rep(30,num_studies)){
  set.seed(seed.val)
  
  mle.theta.re.beta <- rep(NA,num_sims)
  ci.mle.re.beta <- matrix(NA,nrow = num_sims,ncol=2)
  mle.theta.all.study.estimates.re.beta <- matrix(NA,nrow=num_sims,ncol=num_studies)
  
  bayesian.theta.beta <- rep(NA,num_sims)
  ci.bayesian.beta <- matrix(NA,nrow = num_sims,ncol=2)
  bayesian.theta.beta.all.study.estimates <- matrix(NA,nrow=num_sims,ncol=num_studies)
  
  bayesian.theta.uniform <- rep(NA,num_sims)
  ci.bayesian.uniform <- matrix(NA,nrow = num_sims,ncol=2)
  bayesian.theta.uniform.all.study.estimates <- matrix(NA,nrow=num_sims,ncol=num_studies)
  
  bayesian.theta.transformed <- rep(NA,num_sims)
  ci.bayesian.transformed <- matrix(NA,nrow = num_sims,ncol=2)
  bayesian.theta.transformed.all.study.estimates <- matrix(NA,nrow=num_sims,ncol=num_studies)
  
  zero.study <- rep(NA,num_sims)
  zero.study.num <- rep(NA,num_sims)
  
  bayesian_model_beta ="model {
  for( i in 1 : Num ) 
  {
  x[i] ~ dbin(alpha[i],n[i]);
  alpha[i] = p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i]))); 
  theta[i] ~ dbeta(a+0.01,b+0.01);
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

for (j in 1:num_sims){
  s = seed.val + j # set our seed
  x = studyData_varyingTheta(N=num_studies, study_size = study_sizes_vec, seed=s, 
                             theta=theta_simulation_vals, pRatio_values = p_Ratios_vec) 
  # MLE/ Frequentist approach: RE - Beta
  theta.start = rep(0,num_studies)
  for (i in 1:length(theta.start)){
    if (x$carriers[i] > 0){
    theta.start[i] = mle.one.study(carriers = x$carriers[i],study.sizes = x$study_sizes[i],p.ratios = x$p_ratios[i])
    } else {
      theta.start[i] = 0.000001
    }
  }
  theta.start.hyper = c(1,1,theta.start)
  mle.optim <- optim(par = theta.start.hyper, 
                     fn = likelihood_added_beta, 
                     p_studies=x$p_ratios, 
                     num_carriers=x$carriers,
                     study_sizes = x$study_sizes)
  mle.theta.all.study.estimates.re.beta[j,] = mle.optim$par[-c(1,2)] # store individual theta study estimates
  # bootstrapping with 1000 replications
  data_mle_re <- as.data.frame(cbind(x$carriers,x$study_sizes,x$p_ratios))
  colnames(data_mle_re) <- c("carriers","study_sizes","p_ratios")
  boot_strap_median_results_re <- boot(data=data_mle_re, statistic=bootstrap_re_beta, R=1000)
  mle.theta.re.beta[j] = boot_strap_median_results_re$t0
  ci.mle.re.beta[j,] <- quantile(boot_strap_median_results_re$t, probs = c(0.025,0.975))
  
  # make seed values for initialization and initialize data
  seed.values = abs(.Random.seed[1:3])+s
  data=list(x=x$carriers,n=x$study_sizes, p_ratio=x$p_ratios, Num=x$N)
  initial_values <- list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                         list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper"))
  
  # BAYESIAN BETA PRIOR
  model_bayesian_beta =jags.model(textConnection(bayesian_model_beta), data=data, n.chains = 3,
                                  inits=initial_values,
                                  quiet = TRUE, n.adapt= 5000)
  mcmc_output_beta = coda.samples(model_bayesian_beta,c("theta"), n.iter = 10000,
                                  thin = 10,progress.bar="none")
  
  median_chains_bayesian_beta <- c(apply(mcmc_output_beta[[1]][,1:num_studies],1,median),
                                   apply(mcmc_output_beta[[2]][,1:num_studies],1,median),
                                   apply(mcmc_output_beta[[3]][,1:num_studies],1,median))
  
  median_est_bayesian_beta <- mean(median_chains_bayesian_beta)
  median_est_CI_bayesian_beta <- quantile(median_chains_bayesian_beta, probs = c(0.025,0.975))
  
  summary_beta <- summary(mcmc_output_beta)
  # store what we want to keep
  bayesian.theta.beta.all.study.estimates[j,] = summary_beta$statistics[1:(num_studies),1] # each study estimate
  bayesian.theta.beta[j] = median_est_bayesian_beta
  ci.bayesian.beta[j,] = median_est_CI_bayesian_beta
  
  # BAYESIAN UNIFORM APPROACH
  model_bayesian_uni =jags.model(textConnection(bayesian_model_uniform), data=data, n.chains = 3,
                                 inits=initial_values,
                                 quiet = TRUE,n.adapt= 5000)
  mcmc_output_uni = coda.samples(model_bayesian_uni,c("theta"), n.iter = 10000,
                                 thin = 10,progress.bar="none")
  
  median_chains_bayesian_uni <- c(apply(mcmc_output_uni[[1]][,1:num_studies],1,median),
                                  apply(mcmc_output_uni[[2]][,1:num_studies],1,median),
                                  apply(mcmc_output_uni[[3]][,1:num_studies],1,median))
  
  median_est_bayesian_uni <- mean(median_chains_bayesian_uni)
  median_est_CI_bayesian_uni <- quantile(median_chains_bayesian_uni, probs = c(0.025,0.975))
  
  summary_uni <- summary(mcmc_output_uni)
  # store what we want to keep
  bayesian.theta.uniform.all.study.estimates[j,] = summary_uni$statistics[1:(num_studies),1]
  bayesian.theta.uniform[j] = median_est_bayesian_uni
  ci.bayesian.uniform[j,] = median_est_CI_bayesian_uni
  
  # BAYESIAN TRANSFORMED APPROACH
  model_bayesian_transformed =jags.model(textConnection(bayesian_model_transformed), data=data, n.chains = 3,
                                 inits=initial_values,
                                 quiet = TRUE,n.adapt= 5000)
  mcmc_output_transformed = coda.samples(model_bayesian_transformed,c("theta"), n.iter = 10000,
                                 thin = 10,progress.bar="none")
  
  median_chains_bayesian_transformed <- c(apply(mcmc_output_transformed[[1]][,1:num_studies],1,median),
                                  apply(mcmc_output_transformed[[2]][,1:num_studies],1,median),
                                  apply(mcmc_output_transformed[[3]][,1:num_studies],1,median))
  
  median_est_bayesian_transformed <- mean(median_chains_bayesian_transformed)
  median_est_CI_bayesian_transformed <- quantile(median_chains_bayesian_transformed, probs = c(0.025,0.975))
  
  summary_transformed <- summary(mcmc_output_transformed)
  # store what we want to keep
  bayesian.theta.transformed.all.study.estimates[j,] = summary_transformed$statistics[1:(num_studies),1]
  bayesian.theta.transformed[j] = median_est_bayesian_transformed
  ci.bayesian.transformed[j,] = median_est_CI_bayesian_transformed
  
  # study with 0 carriers?
  zero.study[j]=ifelse(0 %in% x$carriers, 1 , 0)
  zero.study.num[j]=length(which(x$carriers == 0))
}
return(list(mle.theta.all.study.estimates.re.beta,mle.theta.re.beta,ci.mle.re.beta,
            bayesian.theta.beta.all.study.estimates,bayesian.theta.beta,ci.bayesian.beta,
            bayesian.theta.uniform.all.study.estimates,bayesian.theta.uniform,ci.bayesian.uniform,
            bayesian.theta.transformed.all.study.estimates,bayesian.theta.transformed,ci.bayesian.transformed,
            zero.study,zero.study.num))
}
#  ------------------------------------------------------------------------#
# realistic simulation
number_studies = 6
study_sizes_vector <- c(1000,125,1000,200,500,50)
p_ratios_vector <- c(10,5,5,55,55,20)
set.seed(74822)
theta_simulation_values <- runif(n=number_studies, min = 0.0009,max = 0.0011)
theta_varying_simul_data <- simulation_varying_theta(seed.val = 61,
                                                     theta_simulation_vals=theta_simulation_values,
                                                     study_sizes_vec = study_sizes_vector,
                                                     p_Ratios_vec = p_ratios_vector,
                                                     num_studies = number_studies)
saveRDS(theta_varying_simul_data,file="varying_theta_scenario12.rds")

set.seed(3141231)
theta_simulation_values <- runif(n=number_studies, min = 0.00009,max = 0.00011)
theta_varying_simul_data <- simulation_varying_theta(seed.val = 2103,
                                                     theta_simulation_vals=theta_simulation_values,
                                                     study_sizes_vec = study_sizes_vector,
                                                     p_Ratios_vec = p_ratios_vector,
                                                     num_studies = number_studies)
saveRDS(theta_varying_simul_data,file="varying_theta_scenario13.rds")
#  ------------------------------------------------------------------------#



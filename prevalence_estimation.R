# needed librarys
library(rmeta)
library(rjags)
library(boot)
#  --------------------------------- MLE THEORY --------------------------------- #
mle.one.study <- function(carriers,study.sizes,p.ratios,vec_pratio_type){
  # return the mle for the one study
  mle = ifelse(vec_pratio_type == 1,
               carriers / (p.ratios*(study.sizes-carriers) + carriers),
               carriers/(study.sizes*p.ratios) )
  return(mle)
}

likelihood.fcn <- function(theta,p.ratios,study.sizes,carriers,vec_pratio_type){
  alphas = ifelse(vec_pratio_type == 1,
                  p.ratios * theta / ((p.ratios * theta) + 1 - theta),
                  p.ratios* theta)
  deriv.alpha = ifelse(vec_pratio_type == 1,
                       p.ratios * (p.ratios*theta + 1-theta)^(-2),
                       p.ratios)
  deriv.log.likeli = sum(carriers*(1/alphas)*deriv.alpha-(study.sizes-carriers)*(1/(1-alphas))*deriv.alpha)
  return(deriv.log.likeli)
}

mle.maximizer <- function(carriers,study.sizes,p.ratios,vec_pratio_type){
  mle = uniroot(likelihood.fcn, p.ratios = p.ratios, study.sizes=study.sizes,
                vec_pratio_type=vec_pratio_type,
                carriers = carriers, lower=0.0000000001, upper=0.1,tol=.Machine$double.eps)$root
  return(mle)
}
fisher.info <- function(theta,p.ratios,study.sizes,carriers,vec_pratio_type){
  alphas = ifelse(vec_pratio_type == 1,
                  p.ratios * theta / ((p.ratios * theta) + 1 - theta),
                  p.ratios* theta)
  deriv.alpha = ifelse(vec_pratio_type == 1,
                       p.ratios * (p.ratios*theta + 1-theta)^(-2),
                       p.ratios)
  fisher.info <- sum(study.sizes*deriv.alpha^2*(alphas*(1-alphas))^(-1))
  return(fisher.info)
}

double.deriv <- function(p,theta,sum_x,n,vec_pratio_type){
  alpha = ifelse(vec_pratio_type == 1,
                  p * theta / ((p * theta) + 1 - theta),
                  p* theta)
  alpha.dot = ifelse(vec_pratio_type == 1,
                       p * (p*theta + 1-theta)^(-2),
                       p)
  alpha.ddot = ifelse(vec_pratio_type == 1,
                      -2*p*(p-1)*((1+(theta*p)-theta)^(-3)),
                      0)
  double.deriv.val = -sum_x*alpha.dot^2*alpha^(-2) +
    sum_x*alpha.ddot*alpha^(-1) -
    (n-sum_x)*alpha.dot^2*(1-alpha)^2-
    (n-sum_x)*(1-alpha)^(-1)*alpha.ddot
  return(double.deriv.val)
}

ci.mle.fisher <- function(theta,fisher.information,alpha.val=0.05){
  log.variance <- 1/(fisher.information*theta^2)
  log.ci.hi <- log(theta) + qnorm(1-(alpha.val/2))*sqrt(log.variance)
  log.ci.lo <- log(theta) + qnorm(alpha.val/2)*sqrt(log.variance)

  hi.ci <- exp(log.ci.hi)
  lo.ci <- exp(log.ci.lo)

  return(c(lo.ci,hi.ci))
}

ci.mle.transform <- function(theta,var,alpha.val=0.05){
  log.variance <- var*(1/theta^2)
  log.ci.hi <- log(theta) + qnorm(1-(alpha.val/2))*sqrt(log.variance)
  log.ci.lo <- log(theta) + qnorm(alpha.val/2)*sqrt(log.variance)

  hi.ci <- exp(log.ci.hi)
  lo.ci <- exp(log.ci.lo)

  return(c(lo.ci,hi.ci))
}

likelihood_added_beta <- function(theta_and_hyper, p_studies, num_carriers, study_sizes,vec_pratio_type) {
  hyper = theta_and_hyper[1:2]
  theta = theta_and_hyper[-c(1,2)]
  num_Studies <- length(p_studies)
  alpha_studies <- ifelse(vec_pratio_type == 1,
                          (p_studies*theta)/((p_studies*theta)+1-theta) ,
                          p_studies* theta)
  if (min(theta_and_hyper) > 0 & max(theta) < 1 
      & min(alpha_studies) > 0 &  max(alpha_studies) < 1){
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
  	if(d$carriers[i] > 0){
    theta.start[i] = mle.maximizer(carriers = d$carriers[i],
                                   study.sizes = d$study_sizes[i],
                                   p.ratios = d$p_ratios[i],
                                   vec_pratio_type = d$p_ratio_type[i])
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
                     vec_pratio_type = d$p_ratio_type,
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

#' Prevalence Estimation Function
#'
#' This function combines estimates from multiple studies through a meta-analysis,
#' incorporating study-specific ascertainment mechanisms into a joint likelihood function.
#' There are various methods for estimation of the consensus estimate of prevalence using
#' both a frequentist and a Bayesian based approach. The inputs required are the number of
#' carriers in the studies, the study sizes, as well as the ascertainment probability ratio
#' between carriers and non carriers. The function outputs study specific estimates, an
#' overall mean prevalence estimate and confidence interval.
#'
#' @param carriers A vector containing the number of carriers in each study.
#' @param study_sizes A vector containing the number of individuals in each study.
#' @param p_ratio_vals A vector containing the ascertainment probability ratios for each study
#' (either the true assumed value or the mean of a distribution -- see p_ratios_unknown).
#' @param p_ratio_type A vector containing whether the values in p_ratio_vals
#' is p_+/p_- or p_+/c - a 1 indicates you have entered the p_+/p_- ratio whereas
#' a 0 indicates you have entered p_+/c
#' @param frequentist A boolean whether the estimation method should be Frequentist or Bayesian. Defaults to FALSE.
#'@param frequentist_method If Frequentist method is desired, which frequentist method: "RE-Beta" or "WM". Defaults to "RE-Beta".
#' @param bayesian_method If Bayesian method is desired, which Bayesian method: "uniform" or "beta" or "logitNormal" or "logitT3". Defaults to "logitNormal".
#' @param random_prevalence Boolean whether the underlying prevalence estimates are considered as varying
#' from study to study. Defaults to TRUE.
#' @param seed.value The seed value to be used by the Bayesian methods. Defaults to as.numeric(Sys.time()).
#' @param bootstrap_R Number of Bootstrap resamples for the estimate of the variance of the weighted median for
#' the frequentist approach, default is 5000
#' @param p_ratios_unknown Boolean for whether we are treating the ascertainment probability ratios as unknown variables
#' with a given variance. If TRUE, then the values in p_ratio_vals represent the mean of this distribution.
#' Defaults to FALSE
#' @param p_ratios_sd If p_ratios_unknown is TRUE, then the standard deviations of the ascertainment probability ratios.
#' @param chain_convergence_check If TRUE, the JAGs Model is plotted so chain convergence for the mean can be verified.
#' This check is only valid for Bayesian approaches. Default is FALSE.
#' @param double_derivative_check If TRUE, the double derivative of the likelihood evaluated at the MLE is printed so it can be
#' verified that it is negative (and thus a maximum was acheived). This check is only valid for the frequentist approach.
#' If only one underlying value of theta is assumed, we check that the sum of the double derivative of the log likelihoods for the
#' studies is negative at the MLE, otherwise for the WM approaxh we check that each individual study theta is negative. No current
#' check for the frequentist RE-Beta approach. Default is FALSE.
#' @param n.adapt.bayes Parameter that allows the user to increase the number of iterations used for adaptation. The default is 5,000.
#' @param n.iter.bayes Parameter that allows the user to increase the number of iterations monitored. The default is 10,000.
#' @param thin.value Parameter that allows the user to specify thinning value for bayesian chains
#' @keywords
#' @export
#' @examples
#' prevalence_estimate(c(2,3),study_sizes = c(100,200),p_ratio_vals = c(30,50),frequentist =FALSE,random_prevalence = TRUE,p_ratios_unknown = TRUE,p_ratios_sd = c(4,4),bayesian_method = "uniform")
prevalence_estimate <- function(carriers,study_sizes,
	p_ratio_type,p_ratio_vals,
	frequentist=FALSE, frequentist_method="RE-Beta",
	bayesian_method="logitNormal",
	random_prevalence=TRUE,seed.value= as.numeric(Sys.time()),
	p_ratios_unknown = FALSE,p_ratios_sd,
	chain_convergence_check=FALSE,double_derivative_check = FALSE,
	n.adapt.bayes = 5000,n.iter.bayes = 10000, bootstrap_R= 1000, thin.value=10){
  if (frequentist == TRUE & p_ratios_unknown == TRUE){
    stop("Frequentist Method cannot handle uncertainty in P_Ratio at this time")
  }
  if (length(carriers) == 1 & random_prevalence == TRUE & frequentist==TRUE){
    stop("With only one study in the Frequentist paradigm it must be assumed there is one underlying prevalence")
  }
  set.seed(seed.value)
  num_studies = length(carriers)
  study_estimates <- rep(NA,num_studies)
  if(frequentist == TRUE){
    if(random_prevalence == TRUE){
      if (frequentist_method=="RE-Beta"){ # Frequentist approach: RE - Beta
        theta.start = rep(0,num_studies)
        for (i in 1:length(theta.start)){
        	if(carriers[i] > 0){
          theta.start[i] = mle.one.study(carriers = carriers[i],study.sizes = study_sizes[i],
                                         p.ratios = p_ratio_vals[i], vec_pratio_type = p_ratio_type[i])
          } else {
          	theta.start[i] = 0.000001
          }
        }
        theta.start.hyper = c(1,1,theta.start)
        mle.optim <- optim(par = theta.start.hyper,  fn = likelihood_added_beta, 
                           p_studies=p_ratio_vals, vec_pratio_type = p_ratio_type,
                           num_carriers=carriers,study_sizes = study_sizes)
        study_estimates = mle.optim$par[-c(1,2)] # store individual theta study estimates
        # bootstrapping with 5000 replications (or what user specifies)
        data_mle_re <- as.data.frame(cbind(carriers,study_sizes,p_ratio_vals,p_ratio_type))
        colnames(data_mle_re) <- c("carriers","study_sizes","p_ratios","p_ratio_type")
        boot_strap_median_results_re <- boot(data=data_mle_re, statistic=bootstrap_re_beta, R=bootstrap_R)
        prevalence_estimate = boot_strap_median_results_re$t0
        prevalence_CI <- quantile(boot_strap_median_results_re$t, probs = c(0.025,0.975))
      } else{ # Frequentist approach : Weighted Median
      vj <- rep(NA,num_studies)
      for (i in 1:num_studies){
        study_estimates[i] <- mle.maximizer(carriers = carriers[i],study.sizes = study_sizes[i],
                                            p.ratios = p_ratio_vals[i], vec_pratio_type = p_ratio_type[i])
        vj[i] <- 1/(fisher.info(theta = study_estimates[i],p.ratios =p_ratio_vals[i],
                                study.sizes = study_sizes[i],carriers = carriers[i], 
                                vec_pratio_type = p_ratio_type[i]))
      }
      # bootstrapping with 5000 replications (or what user specifies)
      yj = study_estimates
      data_mle <- as.data.frame(cbind(yj,vj))
      boot_strap_median_results <- boot(data=data_mle, statistic=bootstrap_median, R=bootstrap_R)
      prevalence_estimate <- boot_strap_median_results$t0
      prevalence_CI <- quantile(boot_strap_median_results$t,probs=c(0.025,0.975))
      } } else {
      prevalence_estimate = mle.maximizer(carriers=carriers,study.sizes=study_sizes,
                                          p.ratios=p_ratio_vals,vec_pratio_type=p_ratio_type)
      study_estimates = rep(prevalence_estimate,num_studies)
      fisher_information <- fisher.info(theta=prevalence_estimate,p.ratios=p_ratio_vals,
                                        study.sizes=study_sizes,carriers=carriers,
                                        vec_pratio_type=p_ratio_type)
      prevalence_CI = ci.mle.fisher(prevalence_estimate, fisher_information)
    }
    if (double_derivative_check == TRUE){
      double.derivative.per.study <- rep(NA,num_studies)
      for (i in 1:num_studies){
        double.derivative.per.study[i] = double.deriv(p=p_ratio_vals[i],theta=study_estimates[i],sum_x=carriers[i],n=study_sizes[i])
      }
      ifelse(random_prevalence == TRUE, 
             ifelse(frequentist_method == "RE-Beta",
                    print("No Current Check for RE-Beta Approach") ,
                    print(paste("Individual study log double derivative value: ",double.derivative.per.study))),
             print(paste("Sum of study log double derivative values: ",sum(double.derivative.per.study))))
    }
    prevalence_list <- list("prevalence_estimate"= prevalence_estimate,
                            "prevalence_CI" = prevalence_CI,
                            "study_estimates" = study_estimates)
  } else {
    set.seed(seed.value)
    data=list(x=carriers,n=study_sizes, p_ratio=p_ratio_vals, Num=num_studies, Type = p_ratio_type)
    seed.values = abs(.Random.seed[1:3])
    
    if(random_prevalence==FALSE){
      bayesian_model ="model {
      for( i in 1 : Num )
      {
      x[i] ~ dbin(alpha[i],n[i]);
      alpha[i] = Type[i]*p_ratio[i]*(theta/(p_ratio[i]*theta+(1-theta))) + (1-Type[i])*p_ratio[i]*theta;
      }
      theta ~ dunif(0,1);
      }
      "
      initial_values <- list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                             list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                             list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper"))
      
      model_bayesian =jags.model(textConnection(bayesian_model), data=data, n.chains = 3,
                                 inits= initial_values,
                                 quiet = TRUE,n.adapt = n.adapt.bayes)
      mcmc_output = coda.samples(model_bayesian,c("theta"), n.iter = n.iter.bayes, 
                                 thin = thin.value,progress.bar="none")
      summary <- summary(mcmc_output)
      
      if(chain_convergence_check){
        plot(mcmc_output[,c("theta")])
      }
      
      # store what we want to keep
      prevalence_estimate = summary$statistics[1]
      prevalence_CI = summary$quantiles[c(1,5)]
      study_estimates = rep(prevalence_estimate,num_studies)
    } else if (random_prevalence==TRUE & p_ratios_unknown == FALSE){
      if(bayesian_method=="logitNormal"){
        variable_names = c("theta","mu","tau")
        bayesian_model ="model {
          for( i in 1 : Num )
          {
          x[i] ~ dbin(alpha[i],n[i]);
          alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
          theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
          logittheta[i] ~ dnorm(mu,tau);
          }
          mu ~ dunif(-50,0);
          tau ~ dunif(1,50);
          }
          "
      }else if(bayesian_method == "logitT3"){
        variable_names = c("theta","mu","tau")
        bayesian_model ="model {
        for( i in 1 : Num )
        {
        x[i] ~ dbin(alpha[i],n[i]);
        alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
        theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
        logittheta[i] ~ dt(mu,tau,3);
        }
        mu ~ dunif(-50,0);
        tau ~ dunif(1,50);
      }
        "
      }
      else if(bayesian_method=="beta"){
        variable_names = c("theta","mu","sigmaSq")
        bayesian_model ="model {
          for( i in 1 : Num )
          {
          x[i] ~ dbin(alpha[i],n[i]);
          alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
          theta[i] ~ dbeta(a,b);
          }
          a = (((1-mu)/sigmaSq) - (1/mu))*(mu^2);
          b = a*((1/mu)-1);
          mu ~ dunif(0,0.01);
          sigmaSq ~ dunif(0,mu*(1-mu));
          }
          "
      } else {
        variable_names = c("theta","mu","length")
        bayesian_model ="model {
        for( i in 1 : Num )
        {
        x[i] ~ dbin(alpha[i],n[i]);
        alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
        theta[i] ~ dunif(mu-0.5*length,mu+0.5*length);
        }
        mu ~ dunif(0,0.01);
        length ~ dunif(0,2*mu);
        }
        "
      }
      model_bayesian =jags.model(textConnection(bayesian_model), data=data, n.chains = 3,
                                 inits= list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                                             list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                                             list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper")),
                                 quiet = TRUE, n.adapt= n.adapt.bayes)
      mcmc_output = coda.samples(model_bayesian,variable.names = variable_names, n.iter = n.iter.bayes,n.chains=3, thin = thin.value,
                                 progress.bar="none")
      summary <- summary(mcmc_output)
      median_chains_bayesian <- c(apply(mcmc_output[[1]][,(1:num_studies)+2],1,median),
                                       apply(mcmc_output[[2]][,(1:num_studies)+2],1,median),
                                       apply(mcmc_output[[3]][,(1:num_studies)+2],1,median))
      median_est_bayesian <- mean(median_chains_bayesian)
      median_est_CI_bayesian <- quantile(median_chains_bayesian, probs = c(0.025,0.975))
      
      if(chain_convergence_check){
        plot(mcmc_output,density=FALSE)
        
        a = apply(mcmc_output[[1]][,(1:num_studies)+2],1,median)
        b = apply(mcmc_output[[2]][,(1:num_studies)+2],1,median)
        c = apply(mcmc_output[[3]][,(1:num_studies)+2],1,median)
        plot(a, ylim=range(a, b, c), col='black',type = "l",
             ylab="Median Value",lwd=1)
        lines(b, col='red',type = "l",lwd=1)
        lines(c, col='green',type = "l",lwd=1)
      }

      study_estimates = summary$statistics[(1:num_studies)+2,1]
      prevalence_estimate = median_est_bayesian
      prevalence_CI = median_est_CI_bayesian

    } else {
      data=list(x=carriers,n=study_sizes,p_ratio_mean=p_ratio_vals,p_ratio_sd=p_ratios_sd ,Num=num_studies, Type = p_ratio_type)
      if(bayesian_method=="logitNormal"){
        variable_names = c("theta","mu","tau")
        bayesian_model ="model {
        for( i in 1 : Num )
        {
        x[i] ~ dbin(alpha[i],n[i]);
        alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
        p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0, );
        theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
        logittheta[i] ~ dnorm(mu,tau);
        }
        mu ~ dunif(-50,0);
        tau ~ dunif(1,50);
      }
        "
    } else if(bayesian_method == "logitT3"){
      variable_names = c("theta","mu","tau")
      bayesian_model ="model {
      for( i in 1 : Num )
      {
      x[i] ~ dbin(alpha[i],n[i]);
      alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
      p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0, );
      theta[i] <- exp(logittheta[i])*(1+exp(logittheta[i]))^(-1);
      logittheta[i] ~ dt(mu,tau,3);
      }
      mu ~ dunif(-50,0);
      tau ~ dunif(1,50);
    }
      "
      } else if(bayesian_method=="beta"){
          variable_names = c("theta","mu","sigmaSq")
          bayesian_model = "model {
          for( i in 1 : Num )
          {
          x[i] ~ dbin(alpha[i],n[i]);
          alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
          p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0, );
          theta[i] ~ dbeta(a,b);
          }
          a = (((1-mu)/sigmaSq) - (1/mu))*(mu^2);
          b = a*((1/mu)-1);
          mu ~ dunif(0,0.01);
          sigmaSq ~ dunif(0,mu*(1-mu));
          }
          "
      } else{
        variable_names = c("theta","mu","length")
        bayesian_model ="model {
        for( i in 1 : Num )
        {
        x[i] ~ dbin(alpha[i],n[i]);
        alpha[i] = Type[i]*p_ratio[i]*(theta[i]/(p_ratio[i]*theta[i]+(1-theta[i])))+ (1-Type[i])*p_ratio[i]*theta[i];
        p_ratio[i] ~ dnorm(p_ratio_mean[i],1/p_ratio_sd[i]^2)T(0, );
        theta[i] ~ dunif(mu-0.5*length,mu+0.5*length);
        }
        mu ~ dunif(0,0.01);
        length ~ dunif(0,2*mu);
        }
        "
      }
      model_bayesian =jags.model(textConnection(bayesian_model), data=data, n.chains = 3,
                                 inits= list(list(.RNG.seed=seed.values[1],.RNG.name="base::Super-Duper"),
                                            list(.RNG.seed=seed.values[2],.RNG.name="base::Super-Duper"),
                                            list(.RNG.seed=seed.values[3],.RNG.name="base::Super-Duper")),
                                 quiet = TRUE, n.adapt= n.adapt.bayes)
      mcmc_output = coda.samples(model_bayesian,variable.names = variable_names, n.iter = n.iter.bayes, thin = thin.value,
                                 progress.bar="none")
      summary <- summary(mcmc_output)
      median_chains_bayesian <- c(apply(mcmc_output[[1]][,(1:num_studies)+2],1,median),
                                  apply(mcmc_output[[2]][,(1:num_studies)+2],1,median),
                                  apply(mcmc_output[[3]][,(1:num_studies)+2],1,median))
      median_est_bayesian <- mean(median_chains_bayesian)
      median_est_CI_bayesian <- quantile(median_chains_bayesian, probs = c(0.025,0.975))
      
      study_estimates = summary$statistics[(1:num_studies)+2,1]
      prevalence_estimate = median_est_bayesian
      prevalence_CI = median_est_CI_bayesian
      
      if(chain_convergence_check){
        plot(mcmc_output,density=FALSE)
        
        a = apply(mcmc_output[[1]][,(1:num_studies)+2],1,median)
        b = apply(mcmc_output[[2]][,(1:num_studies)+2],1,median)
        c = apply(mcmc_output[[3]][,(1:num_studies)+2],1,median)
        plot(a, ylim=range(a, b, c), col='black',type = "l",
             ylab="Median Value",lwd=1)
        lines(b, col='red',type = "l",lwd=1)
        lines(c, col='green',type = "l",lwd=1)
      }
    }
    prevalence_list <- list("prevalence_estimate"= prevalence_estimate,
                            "prevalence_CI" = prevalence_CI,
                            "study_estimates" = study_estimates)
  }
  return(prevalence_list)
}

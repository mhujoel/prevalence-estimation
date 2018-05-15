require(ggplot2)
require(gridExtra)
library(truncnorm)
# open generated data
fixed_theta_scenario1 = readRDS(file = "Data/fixed_theta_scenario1.rds")
fixed_theta_scenario2 = readRDS(file = "Data/fixed_theta_scenario2.rds")
varied_theta_scenario3 = readRDS(file = "Data/varying_theta_scenario3.rds")
varied_theta_scenario4 = readRDS(file = "Data/varying_theta_scenario4.rds")
varied_theta_scenario5 = readRDS(file = "Data/varying_theta_scenario5.rds")
varied_theta_scenario6 = readRDS(file = "Data/varying_theta_scenario6.rds")
varied_theta_scenario7 = readRDS(file = "Data/varying_theta_scenario7.rds")
varied_theta_scenario8 = readRDS(file = "Data/varying_theta_scenario8.rds")
varied_theta_scenario9 = readRDS(file = "Data/varying_theta_scenario9.rds")
varied_theta_scenario10 = readRDS(file = "Data/varying_theta_scenario10.rds")
varied_theta_scenario11 = readRDS(file = "Data/varying_theta_scenario11.rds")
varied_theta_scenario12 = readRDS(file = "Data/varying_theta_scenario12.rds")
varied_theta_scenario13 = readRDS(file = "Data/varying_theta_scenario13.rds")

# true theta:
theta_low = 0.001
theta_extra_low = 0.0001
# varying theta:
num_studies = 10
set.seed(20030314) # low
theta_simulation_values_scenario3 <- runif(n = num_studies,min = 0.0009,max = 0.0011)
set.seed(19610527) # extra low
theta_simulation_values_scenario5 <- runif(n = num_studies,min = 0.00009,max = 0.00011)
set.seed(19881231)
theta_simulation_values_scenario7 <- rtruncnorm(n=10, a=0, b=1, mean = 0.001, sd = 0.0001)
set.seed(05111961) 
theta_simulation_values_scenario8 <- rtruncnorm(n=10, a=0, b=1, mean = 0.001, sd = 0.001)
set.seed(07041994) 
theta_simulation_values_scenario9 <- runif(n=10, min = 0.0001,max = 0.0019)
set.seed(08221994) # low prevalence
theta_simulation_values_scenario10 <- runif(n = 5,min = 0.0009,max = 0.0011)
set.seed(4325567) #extremely low prevalence 
theta_simulation_values_scenario11 <- runif(n = 5,min = 0.00009,max = 0.00011)
# realistic 
set.seed(74822)
theta_simulation_values_scenario12 <- runif(n=6, min = 0.0009,max = 0.0011)
set.seed(3141231)
theta_simulation_values_scenario13 <- runif(n=6, min = 0.00009,max = 0.00011)


# Bias
bias_fixed_theta_low_mle <- mean(fixed_theta_scenario1$MLE-theta_low)
bias_fixed_theta_low_bayesian <- mean(fixed_theta_scenario1$BayesianEst-theta_low)
bias_fixed_theta_extra_low_mle <- mean(fixed_theta_scenario2$MLE-theta_extra_low)
bias_fixed_theta_extra_low_bayesian <- mean(fixed_theta_scenario2$BayesianEst-theta_extra_low)
# MSE
mse_fixed_theta_low_mle <- mean( (fixed_theta_scenario1$MLE-theta_low)^2)
mse_fixed_theta_low_bayesian <- mean( (fixed_theta_scenario1$MLE-theta_low)^2)
mse_fixed_theta_extra_low_mle <- mean( (fixed_theta_scenario2$MLE-theta_extra_low)^2)
mse_fixed_theta_extra_low_bayesian <- mean( (fixed_theta_scenario2$MLE-theta_extra_low)^2)

output <- function(scenario,num_simulations=100,true_underlying_median){
  theta_vals <- paste("theta_simulation_values_scenario",scenario,sep="")
  file_name <- paste("varied_theta_scenario",scenario,sep="")
  median_scenario <- median(get(theta_vals))
  
  biasWM <- mean(get(file_name)[[2]] - median_scenario)
  biasREBETA <- mean(get(file_name)[[5]] - median_scenario)
  biasBayesianBeta <- mean(get(file_name)[[8]] - median_scenario)
  biasBayesianUni <- mean(get(file_name)[[11]] - median_scenario)
  biasBayesianLogitN <- mean(get(file_name)[[14]] - median_scenario)
  
  mseWM <- mean((get(file_name)[[2]] - median_scenario)^2)
  mseREBETA <- mean((get(file_name)[[5]] - median_scenario)^2)
  mseBayesianBeta <- mean((get(file_name)[[8]] - median_scenario)^2)
  mseBayesianUni <- mean((get(file_name)[[11]] - median_scenario)^2)
  mseBayesianLogitN <- mean((get(file_name)[[14]] - median_scenario)^2)
  
  outside.mle.wm <- length(which(median_scenario > get(file_name)[[3]][,2]))+ length(which(median_scenario < get(file_name)[[3]][,1]))
  coverageWM <- (num_simulations-outside.mle.wm)/num_simulations
  
  outside.mle.re.beta <- length(which(median_scenario > get(file_name)[[6]][,2]))+ length(which(median_scenario < get(file_name)[[6]][,1]))
  coverageREBETA <- (num_simulations-outside.mle.re.beta)/num_simulations
  
  outside.bayesian.beta <- length(which(median_scenario > get(file_name)[[9]][,2]))+ length(which(median_scenario < get(file_name)[[9]][,1]))
  coverageBayesianBeta <- (num_simulations-outside.bayesian.beta)/num_simulations
  
  outside.bayesian.uniform <- length(which(median_scenario > get(file_name)[[12]][,2]))+ length(which(median_scenario < get(file_name)[[12]][,1]))
  coverageBayesianUni <- (num_simulations-outside.bayesian.uniform)/num_simulations
  
  outside.bayesian.logitN <- length(which(median_scenario > get(file_name)[[15]][,2]))+ length(which(median_scenario < get(file_name)[[15]][,1]))
  coverageBayesianLogitN <- (num_simulations-outside.bayesian.logitN)/num_simulations
  
  bias_vector = c(biasWM,biasREBETA,biasBayesianBeta,biasBayesianUni,biasBayesianLogitN)
  return(data.frame(bias_scenario = bias_vector,
                    relative_bias = bias_vector/true_underlying_median,
                    mse_scenario = c(mseWM,mseREBETA,mseBayesianBeta,mseBayesianUni,mseBayesianLogitN),
                    coverage_scenario = c(coverageWM,coverageREBETA,coverageBayesianBeta,coverageBayesianUni,coverageBayesianLogitN)))
}

realistic_output <- function(file_name,num_simulations=100,true_underlying_median){
  theta_varying_simul_data = get(file_name)
  # how the data is stored is a list:
  mle.theta.re.beta = theta_varying_simul_data[[2]]
  ci.mle.re.beta = theta_varying_simul_data[[3]]
  bayesian.theta.beta = theta_varying_simul_data[[5]]
  ci.bayesian.beta = theta_varying_simul_data[[6]]
  bayesian.theta.uniform = theta_varying_simul_data[[8]]
  ci.bayesian.uniform = theta_varying_simul_data[[9]]
  bayesian.theta.transformed = theta_varying_simul_data[[11]]
  ci.bayesian.transformed = theta_varying_simul_data[[12]]
  # proportion of simulation replicates that have at least 1 study with 0 carriers
  zero_studies = theta_varying_simul_data[[13]]
  zero_studies_number = theta_varying_simul_data[[14]]
  prop <- sum(zero_studies)/num_simulations
  mean <- mean(zero_studies_number)
  # BIAS #
  biasREBeta = mean(mle.theta.re.beta-true_underlying_median)
  biasBayesianBeta = mean(bayesian.theta.beta-true_underlying_median)
  biasBayesianUni = mean(bayesian.theta.uniform-true_underlying_median)
  biasBayesianTransformed = mean(bayesian.theta.transformed-true_underlying_median)
  # MSE #
  mseREBeta = mean((mle.theta.re.beta-true_underlying_median)^2)
  mseBayesianBeta = mean((bayesian.theta.beta-true_underlying_median)^2)
  mseBayesianUni = mean((bayesian.theta.uniform-true_underlying_median)^2)
  mseBayesianTransformed = mean((bayesian.theta.transformed-true_underlying_median)^2)
  # COVERAGE #
  outside.mle.re.beta <- length(which(true_underlying_median > ci.mle.re.beta[,2]))+ length(which(true_underlying_median < ci.mle.re.beta[,1]))
  coverage.mle.re.beta <- (num_simulations-outside.mle.re.beta)/num_simulations
  outside.bayesian.beta <- length(which(true_underlying_median > ci.bayesian.beta[,2]))+ length(which(true_underlying_median < ci.bayesian.beta[,1]))
  coverage.bayesian.beta <- (num_simulations-outside.bayesian.beta)/num_simulations
  outside.bayesian.uniform <- length(which(true_underlying_median > ci.bayesian.uniform[,2]))+ length(which(true_underlying_median < ci.bayesian.uniform[,1]))
  coverage.bayesian.uniform <- (num_simulations-outside.bayesian.uniform)/num_simulations
  outside.bayesian.transformed <- length(which(true_underlying_median > ci.bayesian.transformed[,2]))+ length(which(true_underlying_median < ci.bayesian.transformed[,1]))
  coverage.bayesian.transformed <- (num_simulations-outside.bayesian.transformed)/num_simulations
  #
  bias_vector = c(biasREBeta,biasBayesianBeta,biasBayesianUni,biasBayesianTransformed)
  # 
  print(summary(zero_studies_number))
  return(data.frame(bias_scenario = bias_vector,
                    relative_bias = bias_vector/true_underlying_median,
                    prop_zero_carriers = prop,
                    mean_zero_carriers_study = mean,
                    mse_scenario = c(mseREBeta,mseBayesianBeta,mseBayesianUni,mseBayesianTransformed),
                    coverage_scenario = c(coverage.mle.re.beta,coverage.bayesian.beta,coverage.bayesian.uniform,coverage.bayesian.transformed)))
}
output_misspecification <- function(file_name,num_simulations=100){
  median_scenario = get(file_name)[[1]]
  # Bayesian Beta, Bayesian Beta + Uncertainty, Bayesian Uni, Bayesian Uni + Uncertainty,
  # Coverage
  coverage = c(get(file_name)[[2]],get(file_name)[[3]],get(file_name)[[4]],get(file_name)[[5]],
               get(file_name)[[6]],get(file_name)[[7]])
  # Bias
  biasBayesianBeta <- mean(get(file_name)[[10]] - median_scenario)
  biasBayesianBetaUncertain <- mean(get(file_name)[[11]] - median_scenario)
  biasBayesianUni <- mean(get(file_name)[[14]] - median_scenario)
  biasBayesianUniUncertain <- mean(get(file_name)[[15]] - median_scenario)
  biasBayesianLogitN <- mean(get(file_name)[[18]] - median_scenario)
  biasBayesianLogitNUncertain <- mean(get(file_name)[[19]] - median_scenario)
  bias_vector = c(biasBayesianBeta,biasBayesianBetaUncertain,
                  biasBayesianUni,biasBayesianUniUncertain,
                  biasBayesianLogitN,biasBayesianLogitNUncertain)
  # MSE
  mseBayesianBeta <- mean((get(file_name)[[10]] - median_scenario)^2)
  mseBayesianBetaUncertain <- mean((get(file_name)[[11]] - median_scenario)^2)
  mseBayesianUni <- mean((get(file_name)[[14]] - median_scenario)^2)
  mseBayesianUniUncertain <- mean((get(file_name)[[15]] - median_scenario)^2)
  mseBayesianLogitN <- mean((get(file_name)[[18]] - median_scenario)^2)
  mseBayesianLogitNUncertain <- mean((get(file_name)[[19]] - median_scenario)^2)
  # return Coverage, Bias, MSE
  return(data.frame(bias_scenario = bias_vector,
                    relative_bias = bias_vector/median_scenario,
                    mse_scenario = c(mseBayesianBeta,mseBayesianBetaUncertain,
                                     mseBayesianUni,mseBayesianUniUncertain,
                                     mseBayesianLogitN,mseBayesianLogitNUncertain),
                    coverage_scenario = coverage))
}

scenario_3 = output(3,true_underlying_median=median(theta_simulation_values_scenario3))
scenario_4 = output_misspecification("varied_theta_scenario4")
scenario_5 = output(5,true_underlying_median=median(theta_simulation_values_scenario5))
scenario_6 = output_misspecification("varied_theta_scenario6")
scenario_7 = output(7,true_underlying_median=median(theta_simulation_values_scenario7))
scenario_8 = output(8,true_underlying_median=median(theta_simulation_values_scenario8))
scenario_9 = output(9,true_underlying_median=median(theta_simulation_values_scenario9))
scenario_10 = output(10,true_underlying_median=median(theta_simulation_values_scenario10))
scenario_11 = output(11,true_underlying_median=median(theta_simulation_values_scenario11))

scenario_12 = realistic_output(file_name = "varied_theta_scenario12",
                               true_underlying_median=median(theta_simulation_values_scenario12))
scenario_12$prop_zero_carriers 
scenario_12$mean_zero_carriers_study 

scenario_13 = realistic_output(file_name = "varied_theta_scenario13",
                               true_underlying_median=median(theta_simulation_values_scenario13))
scenario_13$prop_zero_carriers 
scenario_13$mean_zero_carriers_study 


# logitNormal: Circle 
# beta:  triangle
# uniform: square

# "mle" : "darkblue" // 8
# "bayesian": "black" // 16
# "mle-WM": "#F8766D" // 4
# "mle-RE-Beta": ="#FF61CC" // 13
# "bayesian-Beta": "gold" // 2
# "Bayesian-Uniform": "#00A9FF" // 0
# "Bayesian-LogitNormal": "chocolate1" // 1
# "bayesian-Beta+Uncertainty": "gold" // 17
# "Bayesian-Uniform+Uncertainty": "#00A9FF" // 15
# "Bayesian-LogitNormal+Uncertainty": "chocolate1" // 16

# Coverage
coverage_mle_low <- (length(fixed_theta_scenario1$x) - length(which(fixed_theta_scenario1$L.mle >= theta_low)) - length(which(fixed_theta_scenario1$U.mle <= theta_low)))/length(fixed_theta_scenario1$x)
coverage_bayesian_low <- (length(fixed_theta_scenario1$x) - length(which(fixed_theta_scenario1$L.bayes >= theta_low)) - length(which(fixed_theta_scenario1$U.bayes <= theta_low)))/length(fixed_theta_scenario1$x)
coverage_mle_extra_low <- (length(fixed_theta_scenario2$x) - length(which(fixed_theta_scenario2$L.mle >= theta_extra_low)) - length(which(fixed_theta_scenario2$U.mle <= theta_extra_low)))/length(fixed_theta_scenario2$x)
coverage_bayesian_extra_low <- (length(fixed_theta_scenario2$x) - length(which(fixed_theta_scenario2$L.bayes >= theta_extra_low)) - length(which(fixed_theta_scenario2$U.bayes <= theta_extra_low)))/length(fixed_theta_scenario2$x)


# output order: WM, REBeta, BayesianBeta, BayesianUni, BayesianTransformed
# realistic_output order: REBeta, BayesianBeta, BayesianUni, BayesianTransformed
# output_misspecification: BayesianBeta, BayesianBetaUncertain,BayesianUni,BayesianUniUncertain,
#                     BayesianLogitN,BayesianLogitNUncertain
# 4,6: misspecifcation, 12,13: realistic

relative_bias <- c(bias_fixed_theta_low_mle/0.001,
                  bias_fixed_theta_low_bayesian/0.001,
                  bias_fixed_theta_extra_low_mle/0.0001,
                  bias_fixed_theta_extra_low_bayesian/0.0001,
                  scenario_3$relative_bias[c(1,2,5)],scenario_4$relative_bias[c(5,6)],
                  scenario_5$relative_bias[c(1,2,5)],scenario_6$relative_bias[c(5,6)],
                  scenario_7$relative_bias[c(1,2,5)],scenario_8$relative_bias[c(1,2,5)],
                  scenario_9$relative_bias[c(1,2,5)],scenario_10$relative_bias[c(1,2,5)],
                  scenario_11$relative_bias[c(1,2,5)],
                  0,scenario_12$relative_bias[c(1,4)],
                  0,scenario_13$relative_bias[c(1,4)])

mse <- c(mse_fixed_theta_low_mle,mse_fixed_theta_low_bayesian,
         mse_fixed_theta_extra_low_mle,mse_fixed_theta_extra_low_bayesian,
         scenario_3$mse_scenario[c(1,2,5)],scenario_4$mse_scenario[c(5,6)],
         scenario_5$mse_scenario[c(1,2,5)],scenario_6$mse_scenario[c(5,6)],
         scenario_7$mse_scenario[c(1,2,5)],scenario_8$mse_scenario[c(1,2,5)],
         scenario_9$mse_scenario[c(1,2,5)],scenario_10$mse_scenario[c(1,2,5)],
         scenario_11$mse_scenario[c(1,2,5)],
         0,scenario_12$mse_scenario[c(1,4)],
         0,scenario_13$mse_scenario[c(1,4)])
coverage <- c(coverage_mle_low,coverage_bayesian_low,
              coverage_mle_extra_low,coverage_bayesian_extra_low,
              scenario_3$coverage_scenario[c(1,2,5)],scenario_4$coverage_scenario[c(5,6)],
              scenario_5$coverage_scenario[c(1,2,5)],scenario_6$coverage_scenario[c(5,6)],
              scenario_7$coverage_scenario[c(1,2,5)],scenario_8$coverage_scenario[c(1,2,5)],
              scenario_9$coverage_scenario[c(1,2,5)],scenario_10$coverage_scenario[c(1,2,5)],
              scenario_11$coverage_scenario[c(1,2,5)],
              0.95,scenario_12$coverage_scenario[c(1,4)],
              0.95,scenario_13$coverage_scenario[c(1,4)])

# output order: WM, REBeta, BayesianBeta, BayesianUni, BayesianTransformed
# realistic_output order: REBeta, BayesianBeta, BayesianUni, BayesianTransformed
# output_misspecification: BayesianBeta, BayesianBetaUncertain,BayesianUni,BayesianUniUncertain,
#                     BayesianLogitN,BayesianLogitNUncertain
# 4,6: misspecifcation, 12,13: realistic

estimation <- c("mle","bayesian","mle","bayesian",
                rep(c("mle-WM","mle-RE-Beta","BayesianTransformed"),1),
                rep(c("BayesianTransformed","BayesianTransformed+Uncertainty"),1),
                rep(c("mle-WM","mle-RE-Beta","BayesianTransformed"),1),
                rep(c("BayesianTransformed","BayesianTransformed+Uncertainty"),1),
                rep(c("mle-WM","mle-RE-Beta","BayesianTransformed"),5),
                rep(c("mle-WM","mle-RE-Beta","BayesianTransformed"),2))
scenario <- c(1,1,2,2,rep(3,3),rep(4,2),rep(5,3),rep(6,2),rep(7,3),rep(8,3),rep(9,3),
              rep(10,3),rep(11,3),rep(12,3),rep(13,3))
truth <- c(rep(1,29),0,1,1,0,1,1)
simulation_scenarios <- data.frame(cbind(estimation,scenario,
                                         relative_bias,mse,
                                         coverage,truth)) 

simulation_scenarios$scenario <- as.numeric(as.character(simulation_scenarios$scenario))
simulation_scenarios$relative_bias <- as.numeric(as.character(simulation_scenarios$relative_bias))
simulation_scenarios$mse <- as.numeric(as.character(simulation_scenarios$mse))
simulation_scenarios$coverage <- as.numeric(as.character(simulation_scenarios$coverage))

percent_bias_plot <- ggplot(simulation_scenarios,aes(factor(scenario))) +
  geom_bar(stat = "identity",aes(x=factor(scenario),y=relative_bias,
               fill=factor(estimation),color=factor(truth)),position = "dodge") +
  theme_bw()  + ylab("Relative Bias") + xlab("Scenario Number") + 
  theme(legend.position="none") +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("mle","bayesian","mle-WM","mle-RE-Beta",
                             "BayesianTransformed","BayesianTransformed+Uncertainty"),
                    labels = c("Fixed MLE","Fixed Bayesian","Frequentist-WM","Frequentist-RE-Beta", 
                               "Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty"),
                    values= c("mle"="#b04d88","bayesian"="#46bf97",
                              "mle-WM"="#b44c48","mle-RE-Beta"="#729d4a",
                              "BayesianTransformed"="#746bba",
                              "BayesianTransformed+Uncertainty"="#ba8c3c"))  +
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

mse_plot <- ggplot(simulation_scenarios,aes(factor(scenario))) +
  geom_bar(stat = "identity",aes(x=factor(scenario),y=sqrt(mse),
                                   fill=factor(estimation),color=factor(truth)),
             position = "dodge") +
  theme_bw() + ylab(expression(sqrt(MSE))) + xlab("Scenario Number") + theme(legend.position="none") +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("mle","bayesian","mle-WM","mle-RE-Beta",
                             "BayesianTransformed","BayesianTransformed+Uncertainty"),
                    labels = c("Fixed MLE","Fixed Bayesian","Frequentist-WM","Frequentist-RE-Beta", 
                               "Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty"),
                    values= c("mle"="#b04d88","bayesian"="#46bf97",
                              "mle-WM"="#b44c48","mle-RE-Beta"="#729d4a",
                              "BayesianTransformed"="#746bba",
                              "BayesianTransformed+Uncertainty"="#ba8c3c")) +
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

coverage_plot <- ggplot(simulation_scenarios,aes(factor(scenario))) +
  geom_bar(stat = "identity",aes(x=factor(scenario),y=100*(coverage-0.95),
                                 fill=factor(estimation),color=factor(truth)),
             position = "dodge") +
  theme_bw() + ylab("Coverage - 95%") + xlab("Scenario Number") + theme(legend.position="right")  +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("mle","bayesian","mle-WM","mle-RE-Beta",
                             "BayesianTransformed","BayesianTransformed+Uncertainty"),
                    labels = c("Fixed MLE","Fixed Bayesian","Frequentist-WM","Frequentist-RE-Beta", 
                               "Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty"),
                    values= c("mle"="#b04d88","bayesian"="#46bf97",
                              "mle-WM"="#b44c48","mle-RE-Beta"="#729d4a",
                              "BayesianTransformed"="#746bba",
                              "BayesianTransformed+Uncertainty"="#ba8c3c")) +
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)


simulation_plot_percent_bias <- grid.arrange(percent_bias_plot,mse_plot,coverage_plot,layout_matrix = rbind(c(1,2), c(3,3)))
ggsave(simulation_plot_percent_bias,filename = "Figures/Bar_Chart_Scenarios1-13.pdf",width=11,height = 8,units = "in")


# Realistic Simulations:
scenario_12$prop_zero_carriers # 74%
scenario_13$prop_zero_carriers # 100%

# Realistic Scenarios: 4/6 + misspecification in 4
varied_theta_scenario4 = readRDS(file = "varying_theta_scenario4.rds")
varied_theta_scenario4_1above = readRDS(file = "varying_theta_scenario4_misspecified1above.rds")
varied_theta_scenario4_2above = readRDS(file = "varying_theta_scenario4_misspecified2above.rds")
varied_theta_scenario4_3above = readRDS(file = "varying_theta_scenario4_misspecified3above.rds")
varied_theta_scenario4_4above = readRDS(file = "varying_theta_scenario4_misspecified4above.rds")
varied_theta_scenario4_1below = readRDS(file = "varying_theta_scenario4_misspecified1below.rds")
varied_theta_scenario4_2below = readRDS(file = "varying_theta_scenario4_misspecified2below.rds")
varied_theta_scenario4_3below = readRDS(file = "varying_theta_scenario4_misspecified3below.rds")
varied_theta_scenario4_4below = readRDS(file = "varying_theta_scenario4_misspecified4below.rds")
# No misspecification
misspecification_scenario4_none = output_misspecification("varied_theta_scenario4")
# misspecification
misspecification_scenario4_1above = output_misspecification("varied_theta_scenario4_1above")
misspecification_scenario4_2above = output_misspecification("varied_theta_scenario4_2above")
misspecification_scenario4_3above = output_misspecification("varied_theta_scenario4_3above")
misspecification_scenario4_4above = output_misspecification("varied_theta_scenario4_4above")
misspecification_scenario4_1below = output_misspecification("varied_theta_scenario4_1below")
misspecification_scenario4_2below = output_misspecification("varied_theta_scenario4_2below")
misspecification_scenario4_3below = output_misspecification("varied_theta_scenario4_3below")
misspecification_scenario4_4below = output_misspecification("varied_theta_scenario4_4below")

# output_misspecification: BayesianBeta, BayesianBetaUncertain,BayesianUni,BayesianUniUncertain,
#                     BayesianLogitN,BayesianLogitNUncertain

relative_bias_misspecification <- c(0,misspecification_scenario4_4below$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_3below$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_2below$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_1below$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_none$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_1above$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_2above$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_3above$relative_bias[c(5,6)],0,
                                    0,misspecification_scenario4_4above$relative_bias[c(5,6)],0)
mse_misspecification <- c(0,misspecification_scenario4_4below$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_3below$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_2below$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_1below$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_none$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_1above$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_2above$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_3above$mse_scenario[c(5,6)],0,
                          0,misspecification_scenario4_4above$mse_scenario[c(5,6)],0)
coverage_misspecification <- c(.95,misspecification_scenario4_4below$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_3below$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_2below$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_1below$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_none$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_1above$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_2above$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_3above$coverage_scenario[c(5,6)],0.95,
                               0.95,misspecification_scenario4_4above$coverage_scenario[c(5,6)],0.95)

# estimation_misspecification <- rep(c("BayesianTransformed","BayesianTransformed+Uncertainty"),9)
# scenario_misspecification <- c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2))
# simulation_scenarios_misspecification <- data.frame(cbind(estimation = estimation_misspecification,
#                                                           scenario = scenario_misspecification,
#                                          relative_bias = relative_bias_misspecification,
#                                          mse = mse_misspecification,
#                                          coverage = coverage_misspecification)) 
# simulation_scenarios_misspecification$scenario <- as.numeric(as.character(simulation_scenarios_misspecification$scenario))
# simulation_scenarios_misspecification$relative_bias <- as.numeric(as.character(simulation_scenarios_misspecification$relative_bias))
# simulation_scenarios_misspecification$mse <- as.numeric(as.character(simulation_scenarios_misspecification$mse))
# simulation_scenarios_misspecification$coverage <- as.numeric(as.character(simulation_scenarios_misspecification$coverage))

estimation_misspecification <- rep(c("a","BayesianTransformed","BayesianTransformed+Uncertainty","c"),9)
scenario_misspecification <- c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),
                               rep(8,4),rep(9,4))
simulation_scenarios_misspecification <- data.frame(cbind(estimation = estimation_misspecification,
                                                          scenario = scenario_misspecification,
                                                          relative_bias = relative_bias_misspecification,
                                                          mse = mse_misspecification,
                                                          coverage = coverage_misspecification)) 
simulation_scenarios_misspecification$scenario <- as.numeric(as.character(simulation_scenarios_misspecification$scenario))
simulation_scenarios_misspecification$relative_bias <- as.numeric(as.character(simulation_scenarios_misspecification$relative_bias))
simulation_scenarios_misspecification$mse <- as.numeric(as.character(simulation_scenarios_misspecification$mse))
simulation_scenarios_misspecification$coverage <- as.numeric(as.character(simulation_scenarios_misspecification$coverage))

simulation_scenarios_misspecification$true_value = rep(c(0,1,1,0),9)

 

relative_bias_plot_misspecification <- ggplot(simulation_scenarios_misspecification,aes(factor(scenario))) +
  geom_bar(stat = "identity",aes(x=factor(scenario),y=relative_bias,
                                               fill=factor(estimation),color=factor(true_value)),position="dodge") +
  theme_bw() +   geom_hline(yintercept = 0,color="grey",alpha=0.5) + ylab("Relative Bias") + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="none") +
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  scale_fill_manual(name = "Estimation Technique",
                      breaks=c("a","BayesianTransformed","BayesianTransformed+Uncertainty","c"),
                      labels = c("a","Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty","c"),
                      values= c("a"="white","BayesianTransformed"="#746bba",
                                "BayesianTransformed+Uncertainty"="#ba8c3c",
                                "c"="white")) +
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) 

mse_plot_misspecification <- ggplot(simulation_scenarios_misspecification,aes(factor(scenario))) +
  geom_bar(stat = "identity",aes(x=factor(scenario),y=sqrt(mse),
                                 fill=factor(estimation),
                                 color=factor(true_value)),position="dodge") +
  theme_bw() +   geom_hline(yintercept = 0,color="grey",alpha=0.5) + ylab(expression(sqrt(MSE))) + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="none") +
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("a","BayesianTransformed","BayesianTransformed+Uncertainty","c"),
                    labels = c("a","Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty","c"),
                    values= c("a"="white","BayesianTransformed"="#746bba",
                              "BayesianTransformed+Uncertainty"="#ba8c3c",
                              "c"="white"))+
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) 

simulation_scenarios_misspecification$estimation2<-factor(simulation_scenarios_misspecification$estimation,
                                    labels=c("a","Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty","c"))

coverage_plot_misspecification <- ggplot(simulation_scenarios_misspecification,aes(factor(scenario))) +
  geom_bar(stat = "identity",position="dodge",
           aes(x=factor(scenario),y=(coverage-0.95)*100,fill=factor(estimation2),color=factor(true_value))) +
  theme_bw() + ylab("Coverage - 95%") + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="right") +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("Bayesian-LogitNormal","Bayesian-LogitNormal + Uncertainty"),
                    values= c("a"="white","Bayesian-LogitNormal"="#746bba",
                              "Bayesian-LogitNormal + Uncertainty"="#ba8c3c",
                              "c"="white"))+
  scale_color_manual(breaks=c(1,0),values=c("white","black"),guide = FALSE) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) + 
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

misspecification_plot <- grid.arrange(relative_bias_plot_misspecification,mse_plot_misspecification,coverage_plot_misspecification,layout_matrix = rbind(c(1,2), c(3,3)))
ggsave(misspecification_plot,filename = "Figures/Bars_Misspecification_Scenario_Graphics.pdf",width=11,height = 8,units = "in")


##### ----------------- Appendix Plots ------------------ ####
relative_bias_beta_uni <- c(scenario_3$relative_bias[c(3,4)],
                            scenario_4$relative_bias[c(1:4)],
                   scenario_5$relative_bias[c(3,4)],scenario_6$relative_bias[c(1:4)],
                   scenario_7$relative_bias[c(3,4)],scenario_8$relative_bias[c(3,4)],
                   scenario_9$relative_bias[c(3,4)],scenario_10$relative_bias[c(3,4)],
                   scenario_11$relative_bias[c(3,4)],scenario_12$relative_bias[c(2,3)], 
                   scenario_13$relative_bias[c(2,3)])
mse_beta_uni <- c(scenario_3$mse_scenario[c(3,4)],scenario_4$mse_scenario[c(1:4)],
         scenario_5$mse_scenario[c(3,4)],scenario_6$mse_scenario[c(1:4)],
         scenario_7$mse_scenario[c(3,4)],scenario_8$mse_scenario[c(3,4)],
         scenario_9$mse_scenario[c(3,4)],scenario_10$mse_scenario[c(3,4)],
         scenario_11$mse_scenario[c(3,4)],scenario_12$mse_scenario[c(2,3)],
         scenario_13$mse_scenario[c(2,3)])
coverage_beta_uni <- c(scenario_3$coverage_scenario[c(3,4)],
                       scenario_4$coverage_scenario[c(1:4)],
              scenario_5$coverage_scenario[c(3,4)],
              scenario_6$coverage_scenario[c(1:4)],
              scenario_7$coverage_scenario[c(3,4)],
              scenario_8$coverage_scenario[c(3,4)],
              scenario_9$coverage_scenario[c(3,4)],
              scenario_10$coverage_scenario[c(3,4)],
              scenario_11$coverage_scenario[c(3,4)],
              scenario_12$coverage_scenario[c(2,3)],
              scenario_13$coverage_scenario[c(2,3)])

estimation_beta_uni <- c(rep(c("BayesianBeta","BayesianUni"),1),
                rep(c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),1),
                rep(c("BayesianBeta","BayesianUni"),1),
                rep(c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),1),
                rep(c("BayesianBeta","BayesianUni"),7))
scenario_beta_uni <- c(rep(3,2),rep(4,4),rep(5,2),rep(6,4),rep(7,2),rep(8,2),rep(9,2),rep(10,2),rep(11,2),rep(12,2),rep(13,2))
simulation_scenarios_beta_uni <- data.frame(cbind(estimation=estimation_beta_uni,scenario=scenario_beta_uni,
                                                  relative_bias=relative_bias_beta_uni,
                                                  mse=mse_beta_uni, coverage= coverage_beta_uni)) 

simulation_scenarios_beta_uni$scenario <- as.numeric(as.character(simulation_scenarios_beta_uni$scenario))
simulation_scenarios_beta_uni$relative_bias <- as.numeric(as.character(simulation_scenarios_beta_uni$relative_bias))
simulation_scenarios_beta_uni$mse <- as.numeric(as.character(simulation_scenarios_beta_uni$mse))
simulation_scenarios_beta_uni$coverage <- as.numeric(as.character(simulation_scenarios_beta_uni$coverage))

percent_bias_plot_beta_uni <- ggplot(simulation_scenarios_beta_uni,
                                     aes(factor(scenario))) +
  geom_bar(stat="identity",aes(x=factor(scenario),
                               y=relative_bias,fill=factor(estimation)),position="dodge",color="black") +
  theme_bw() + ylab("Relative Bias") + xlab("Scenario Number") + theme(legend.position="none") +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

mse_plot_beta_uni <- ggplot(simulation_scenarios_beta_uni,aes(factor(scenario))) +
  geom_bar(stat="identity",aes(x=factor(scenario),y=sqrt(mse),
                               fill=factor(estimation)),position = "dodge",color="black") +
  theme_bw() + ylab(expression(sqrt(MSE))) + xlab("Scenario Number") + theme(legend.position="none") +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

coverage_plot_beta_uni <- ggplot(simulation_scenarios_beta_uni,aes(factor(scenario))) +
  geom_bar(stat="identity",aes(x=factor(scenario),y=100*(coverage-0.95),
                               fill=factor(estimation)),position = "dodge",color="black") +
  theme_bw() + ylab("Coverage - 95%") + xlab("Scenario Number") + theme(legend.position="right")   +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)


simulation_plot_percent_bias_beta_uni <- grid.arrange(percent_bias_plot_beta_uni,mse_plot_beta_uni,
                                             coverage_plot_beta_uni,layout_matrix = rbind(c(1,2), c(3,3)))
ggsave(simulation_plot_percent_bias_beta_uni,filename = "Figures/Bars_Scenarios1-13_Supp.pdf",width=11,height = 8,units = "in")







relative_bias_misspecification_beta_uni <- c(misspecification_scenario4_4below$relative_bias[c(1:4)],misspecification_scenario4_3below$relative_bias[c(1:4)],
                                    misspecification_scenario4_2below$relative_bias[c(1:4)],misspecification_scenario4_1below$relative_bias[c(1:4)],
                                    misspecification_scenario4_none$relative_bias[c(c(1:4))],
                                    misspecification_scenario4_1above$relative_bias[c(1:4)],misspecification_scenario4_2above$relative_bias[c(1:4)],
                                    misspecification_scenario4_3above$relative_bias[c(1:4)],misspecification_scenario4_4above$relative_bias[c(1:4)])
mse_misspecification_beta_uni <- c(misspecification_scenario4_4below$mse_scenario[c(1:4)],misspecification_scenario4_3below$mse_scenario[c(1:4)],
                          misspecification_scenario4_2below$mse_scenario[c(1:4)],misspecification_scenario4_1below$mse_scenario[c(1:4)],
                          misspecification_scenario4_none$mse_scenario[c(1:4)],
                          misspecification_scenario4_1above$mse_scenario[c(1:4)],misspecification_scenario4_2above$mse_scenario[c(1:4)],
                          misspecification_scenario4_3above$mse_scenario[c(1:4)],misspecification_scenario4_4above$mse_scenario[c(1:4)])
coverage_misspecification_beta_uni <- c(misspecification_scenario4_4below$coverage_scenario[c(1:4)],misspecification_scenario4_3below$coverage_scenario[c(1:4)],
                               misspecification_scenario4_2below$coverage_scenario[c(1:4)],misspecification_scenario4_1below$coverage_scenario[c(1:4)],
                               misspecification_scenario4_none$coverage_scenario[c(1:4)],
                               misspecification_scenario4_1above$coverage_scenario[c(1:4)],misspecification_scenario4_2above$coverage_scenario[c(1:4)],
                               misspecification_scenario4_3above$coverage_scenario[c(1:4)],misspecification_scenario4_4above$coverage_scenario[c(1:4)])

estimation_misspecification_beta_uni <- rep(c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),9)
scenario_misspecification_beta_uni <- c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,4),rep(9,4))
simulation_scenarios_misspecification_beta_uni <- data.frame(cbind(estimation = estimation_misspecification_beta_uni,
                                                          scenario = scenario_misspecification_beta_uni,
                                                          relative_bias = relative_bias_misspecification_beta_uni,
                                                          mse = mse_misspecification_beta_uni,
                                                          coverage = coverage_misspecification_beta_uni)) 
simulation_scenarios_misspecification_beta_uni$scenario <- as.numeric(as.character(simulation_scenarios_misspecification_beta_uni$scenario))
simulation_scenarios_misspecification_beta_uni$relative_bias <- as.numeric(as.character(simulation_scenarios_misspecification_beta_uni$relative_bias))
simulation_scenarios_misspecification_beta_uni$mse <- as.numeric(as.character(simulation_scenarios_misspecification_beta_uni$mse))
simulation_scenarios_misspecification_beta_uni$coverage <- as.numeric(as.character(simulation_scenarios_misspecification_beta_uni$coverage))

relative_bias_plot_misspecification_beta_uni <- ggplot(simulation_scenarios_misspecification_beta_uni,aes(factor(scenario))) +
  geom_bar(stat="identity",color="black",aes(x=factor(scenario),y=relative_bias,
                               fill=factor(estimation)),position = "dodge") +
  theme_bw()  + ylab("Relative Bias") + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="none") +
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

mse_plot_misspecification_beta_uni <- ggplot(simulation_scenarios_misspecification_beta_uni,aes(factor(scenario))) +
  geom_bar(stat="identity",color="black",aes(x=factor(scenario),y=sqrt(mse),fill=factor(estimation)),position = "dodge") +
  theme_bw() + ylab(expression(sqrt(MSE))) + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="none") +
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

coverage_plot_misspecification_beta_uni <- ggplot(simulation_scenarios_misspecification_beta_uni,aes(factor(scenario))) +
  geom_bar(stat="identity",color="black",aes(x=factor(scenario),y=100*(coverage-.95),fill=factor(estimation)),position = "dodge") +
  theme_bw() + ylab("Coverage - 95%") + xlab("Misspecification of Ascertainment Ratio") + theme(legend.position="right")   +
  scale_x_discrete(labels = c('-1','-0.75','-0.5',"-0.25","0","+0.25","+0.5","+0.75","+1")) +
  scale_fill_manual(name = "Estimation Technique",
                    breaks=c("BayesianBeta","BayesianBetaUncertain","BayesianUni","BayesianUniUncertain"),
                    labels = c("Bayesian-Beta","Bayesian-Beta + Uncertainty","Bayesian-Uniform","Bayesian-Uniform + Uncertainty"),
                    values= c("BayesianBeta"="gold","BayesianBetaUncertain"="goldenrod1",
                              "BayesianUni"="#00A9FF","BayesianUniUncertain"="navy")) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)) + 
  geom_hline(yintercept = 0,color="grey",alpha=0.5)

misspecification_plot_beta_uni <- grid.arrange(relative_bias_plot_misspecification_beta_uni,
                                      mse_plot_misspecification_beta_uni,
                                      coverage_plot_misspecification_beta_uni,layout_matrix = rbind(c(1,2), c(3,3)))
ggsave(misspecification_plot_beta_uni,filename = "Figures/Bar_Misspecification_Scenario_Graphics_Supp.pdf",width=11,height = 8,units = "in")

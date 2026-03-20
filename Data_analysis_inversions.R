rm(list = ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

Tmax = 30000
log_interval = 5e3
q = 5   # Scaling factor
h = log_interval/q
n_rep = 5



replace_NAN <- function(vec) {
  for (i in 1:length(vec)) {
    if (vec[i] == "NAN") {
      vec[i] <- 0
    }
  }
  vec
}

##   ========= Sensitivity analyses ========== 

## 1)  Mutation rate

Mu = c("1.0e-06","5.0e-06","1.0e-05","5.0e-05")

growth_rate_I_matrix_Mu = matrix(nrow = 5, ncol = 4)
growth_rate_N_matrix_Mu = matrix(nrow = 5, ncol = 4)
initial_effective_dominance_Mu = 1:4
initial_effective_dominance_matrix_Mu = matrix(nrow = 5, ncol = 4)

overall_growth_I_Mu = 1:5
overall_growth_N_Mu = 1:5
mean_growth_rates_I_Mu = 1:4
mean_growth_rates_N_Mu = 1:4
mean_growth_I_Mu = 1:4
mean_growth_N_Mu = 1:4
mean_growth_rates_per_mutation_I_Mu = 1:4
mean_growth_rates_per_mutation_N_Mu = 1:4

inversions_Mu <- data.frame(tick = integer(),
                 freq_I=numeric(), 
                 freq_II=numeric(),
                 freq_IN=numeric(),
                 freq_NN=numeric(),
                 number_haplosomes_I = integer(),
                 number_haplosomes_N = integer(),
                 number_mutations_I = character(),
                 number_mutations_N = numeric(),
                 marginal_fitness_inversion_I = character(),
                 marginal_fitness_inversion_N = character(),
                 mean_fitness_inversion_II = character(),
                 mean_fitness_inversion_IN = character(),
                 mean_fitness_inversion_NN = numeric(),
                 fitness_load_I = character(),
                 fitness_load_N = character(),
                 freq_homozygous_mut_II = character(),
                 freq_homozygous_mut_IN = character(),
                 freq_homozygous_mut_NN = character(),
                 mean_fitness_global_II = character(),
                 mean_fitness_global_IN = character(),
                 mean_fitness_global_NN = numeric(),
                 covariance_mutation_inversion_I = character(),
                 covariance_mutation_inversion_N = character(),
                 covariance_out_of_marginal_fitness_I = character(),
                 covariance_out_of_marginal_fitness_N = character(),
                 effective_selection = character(),
                 effective_selection_formula = character(),
                 effective_dominance = character(),
                 gamma = character(),
                 dominance_variance = character(),
                 Rep = integer(),
                 Mu = character(),
                 number_mutations_no_NANs_I = numeric(),
                 growth_rate_I = numeric(),
                 growth_rate_N = numeric(),
                 growth_rate_per_mutation_I = numeric(),
                 growth_rate_per_mutation_N = numeric(),
                 stringsAsFactors=FALSE) 



for (i in (1:length(Mu))) {
  for (Rep in (1:n_rep)) {
    File_Mu = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ", "0.05" ," Mu ",Mu[i]," Rep ",Rep," NoAnalysis"," .csv", sep ="")
    
    aux_inversions_Mu <- read.csv(File_Mu, stringsAsFactors = FALSE)
    
    n_I_Mu = length(aux_inversions_Mu$number_mutations_I)
    n_N_Mu = length(aux_inversions_Mu$number_mutations_N)
    
      aux_inversions_Mu <- aux_inversions_Mu %>%
        mutate(Rep = Rep,
               Mu = as.numeric(Mu[i]),
               number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
               growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_Mu] - number_mutations_no_NANs_I[1:(n_I_Mu-1)])/h,NA),
               growth_rate_N = c((number_mutations_N[2:n_N_Mu] - number_mutations_N[1:(n_N_Mu-1)])/h,NA),
               growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
               growth_rate_per_mutation_N = growth_rate_N/number_mutations_N)
      
      inversions_Mu <- inversions_Mu %>%
        rbind(aux_inversions_Mu)
  
  }}


n_Mu = 11  ## first n after inversion

summary_by_Rep <- inversions_Mu %>% 
  group_by(Mu,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_Mu:length(growth_rate_I)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_Mu:length(growth_rate_N)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[n_Mu:length(growth_rate_I)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[n_Mu:length(growth_rate_N)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)] - number_mutations_no_NANs_I[n_Mu])/(h*(length(number_mutations_no_NANs_I)-n_Mu)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)] - number_mutations_N[n_Mu])/(h*(length(number_mutations_N)-n_Mu))) %>% 
  arrange(Mu)

long <- summary_by_Rep %>% 
  pivot_longer(
    cols = c(mean_growth_I_by_Rep,mean_growth_N_by_Rep),
    names_to = "arrangement",
    values_to = "mean_growth_by_Rep"
  )

print(long)
  
summary <- summary_by_Rep %>% 
  group_by(Mu) %>% 
  summarize(mean_growth_I = c(mean(mean_growth_I_by_Rep,na.rm = T),mean(mean_growth_N_by_Rep , na.rm = T)),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T)) %>% 
  arrange(Mu)

  #   # initial_effective_dominance_matrix_Mu[Rep,i] = mean(as.numeric(inversions_Mu$effective_dominance[11:12]),na.rm = TRUE)
  # # initial_effective_dominance_Mu[i] = mean(initial_effective_dominance_matrix_Mu[,i], na.rm = TRUE)

par(mfrow = c(2,2))

plot(x = as.numeric(Mu), y = summary$mean_growth_I, col = "black", xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(-0.05,0.25))
points(x = Mu, y = summary$mean_growth_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = Mu, y = summary$overall_growth_I , col = "black", xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean mutation growth rate after inversion",ylim = c(-0.05,0.25))
points(x = Mu, y = summary$overall_growth_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = Mu, y = summary$mean_growth_rates_per_mutation_I , col = "black", xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation after inversion")
points(x = Mu, y = summary$mean_growth_rates_per_mutation_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


# plot(x = Mu, y = initial_effective_dominance,col = "red", xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation rates",ylim = c(0,2))


summary_by_Rep

par(mfrow = c(2,3))

summary_by_Rep %>%  ggplot(aes(x = as.factor(Mu),y=mean_growth_I_by_Rep)) +
  geom_boxplot()

summary_by_Rep %>% ggplot(aes(x = as.factor(Mu),y = mean_growth_N_by_Rep)) +
  geom_boxplot()

summary_by_Rep %>% ggplot(aes(x = as.factor(Mu), y = overall_growth_I_by_Rep)) +
  geom_boxplot()

summary_by_Rep %>% ggplot(aes(x = as.factor(Mu), y = overall_growth_N_by_Rep)) +
  geom_boxplot()

summary_by_Rep %>% ggplot(aes(x = as.factor(Mu), y = mean_growth_rates_per_mutation_I_by_Rep)) +
  geom_boxplot()

summary_by_Rep %>% ggplot(aes(x = as.factor(Mu), y = mean_growth_rates_per_mutation_N_by_Rep)) +
  geom_boxplot()

# boxplot(x = summary_by_Rep$mean_growth_I_by_Rep, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# boxplot(x = summary_by_Rep$mean_growth_N_by_Rep, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# 
# boxplot(x = growth_rate_per_mutation_I_matrix_Mu, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = growth_rate_per_mutation_N_matrix_Mu, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = initial_effective_dominance_matrix,border = "red", names = Mu, xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance after inversion")

## 2) Inversion benefit

S_INV = c("0.01","0.02","0.05","0.1")

growth_rate_I_matrix_SINV = matrix(nrow = 5, ncol = 4)
growth_rate_N_matrix_SINV = matrix(nrow = 5, ncol = 4)
initial_effective_dominance_SINV = 1:4
initial_effective_dominance_matrix_SINV = matrix(nrow = 5, ncol = 4)

overall_growth_I_SINV = 1:5
overall_growth_N_SINV = 1:5
mean_growth_rates_I_SINV = 1:4
mean_growth_rates_N_SINV = 1:4
mean_growth_I_SINV = 1:4
mean_growth_N_SINV = 1:4

for (i in (1:4)) {
  for (Rep in (1:5)) {
    File_SINV = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ",S_INV[i]," Mu ","1.0e-05"," Rep ",Rep," InversionBenefitSensitivityAnalysis"," .csv", sep ="")
    
    inversions_SINV <- read.csv(File_SINV, stringsAsFactors = FALSE)
    
    
    n_mut_per_inv_no_NANs_SINV = inversions_SINV$number_mutations_I
    for (j in 1:length(n_mut_per_inv_no_NANs_SINV)) {
      if (inversions_SINV$number_mutations_I[j] == "NAN") {
        n_mut_per_inv_no_NANs_SINV[j] <- as.numeric(0)
      }
    }
    initial_effective_dominance_matrix_SINV[Rep,i] = mean(as.numeric(inversions_SINV$effective_dominance[11:12]),na.rm = TRUE)
    n_I_SINV = length(inversions_SINV$number_mutations_I)
    n_N_SINV = length(inversions_SINV$number_mutations_N)
    growth_rate_I_SINV = (as.numeric(n_mut_per_inv_no_NANs_SINV[2:n_I_SINV]) - as.numeric(n_mut_per_inv_no_NANs_SINV[1:(n_I_SINV-1)]))/h
    growth_rate_N_SINV = (as.numeric(inversions_SINV$number_mutations_N[2:n_N_SINV]) - as.numeric(inversions_SINV$number_mutations_N[1:(n_N_SINV-1)]))/h
    overall_growth_I_SINV[Rep] = (as.numeric(n_mut_per_inv_no_NANs_SINV[n_I_SINV]) - as.numeric(n_mut_per_inv_no_NANs_SINV[11]))/(h*(n_I_SINV-11))
    overall_growth_N_SINV[Rep] = (as.numeric(inversions_SINV$number_mutations_N[n_N_SINV]) - as.numeric(n_mut_per_inv_no_NANs_SINV[11]))/(h*(n_N_SINV-11))
    
    n_SINV = 11  ## first n after inversion
    
    growth_rate_I_matrix_SINV[Rep,i] <- mean(growth_rate_I_SINV[n_SINV:length(growth_rate_I_SINV)])
    growth_rate_N_matrix_SINV[Rep,i] <- mean(growth_rate_N_SINV[n_SINV:length(growth_rate_N_SINV)])
  }
  mean_growth_rates_I_SINV[i] = mean(growth_rate_I_matrix_SINV[,i],na.rm = TRUE)
  mean_growth_rates_N_SINV[i] = mean(growth_rate_N_matrix_SINV[,i],na.rm = TRUE)
  mean_growth_I_SINV[i] = mean(overall_growth_I_SINV,na.rm = TRUE)
  mean_growth_N_SINV[i] = mean(overall_growth_N_SINV,na.rm = TRUE)
  initial_effective_dominance_SINV[i] = mean(initial_effective_dominance_matrix_SINV[,i], na.rm = TRUE)
}

par(mfrow = c(2,2))

plot(x = S_INV, y = mean_growth_rates_I_SINV , col = "black", xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion", ylim = c(0,0.25))
points(x = S_INV, y = mean_growth_rates_N_SINV, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = S_INV, y = mean_growth_I_SINV , col = "black", xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean mutation growth rate after inversion",ylim = c(0,0.25))
points(x = S_INV, y = mean_growth_N_SINV, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

# plot(x = S_INV, y = initial_effective_dominance, col = "red", xlab = "Inversion benefit", ylab = "Initial effective dominance", main = "Initial effective dominance for different inversion benefits")

boxplot(x = growth_rate_I_matrix_SINV, border = "black", names = S_INV, xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(0,0.08))
boxplot(x = growth_rate_N_matrix_SINV, border = "blue", names = S_INV, xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(0,0.08))

# boxplot(x = initial_effective_dominance_matrix, border = "red", names = S_INV, xlab = "Inversion benefit", ylab = "Initial effective dominance", main = "Initial effective dominance for different inversion benefits")



## 3) Deleterious mutation dominance

H = c("0.0","0.1","0.25","0.5")

growth_rate_I_matrix_H = matrix(nrow = 5, ncol = 4)
growth_rate_N_matrix_H = matrix(nrow = 5, ncol = 4)
initial_effective_dominance_H = 1:4
initial_effective_dominance_matrix_H = matrix(nrow = 5, ncol = 4)

overall_growth_I_H = 1:5
overall_growth_N_H = 1:5
mean_growth_rates_I_H = 1:4
mean_growth_rates_N_H = 1:4
mean_growth_I_H = 1:4
mean_growth_N_H = 1:4
crashed_simulations = rep(0,times = 4)

for (i in (1:4)) {
  for (Rep in (1:5)) {
    File_H = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ",H[i], " S_INV ","0.05"," Mu ","1.0e-05"," Rep ",Rep," MutationDominanceSensitivityAnalysis"," .csv", sep ="")
    
    inversions_H <- read.csv(File_H, stringsAsFactors = FALSE)
    
    n_mut_per_inv_no_NANs_H = inversions_H$number_mutations_I
    for (j in 1:length(n_mut_per_inv_no_NANs_H)) {
      if (inversions_H$number_mutations_I[j] == "NAN") {
        n_mut_per_inv_no_NANs_H[j] <- as.numeric(0)
      }
    }
    initial_effective_dominance_matrix_H[Rep,i] = mean(as.numeric(inversions_H$effective_dominance[11:12]),na.rm = TRUE)
    n_I_H = length(inversions_H$number_mutations_I)
    n_N_H = length(inversions_H$number_mutations_N)
    growth_rate_I_H = (as.numeric(n_mut_per_inv_no_NANs_H[2:n_I_H]) - as.numeric(n_mut_per_inv_no_NANs_H[1:(n_I_H-1)]))/h
    growth_rate_N_H = (as.numeric(inversions_H$number_mutations_N[2:n_N_H]) - as.numeric(inversions_H$number_mutations_N[1:(n_N_H-1)]))/h
    overall_growth_I_H[Rep] = (as.numeric(n_mut_per_inv_no_NANs_H[n_I_H]) - as.numeric(n_mut_per_inv_no_NANs_H[11]))/(h*(n_I_H-11))
    overall_growth_N_H[Rep] = (as.numeric(inversions_H$number_mutations_N[n_N_H]) - as.numeric(n_mut_per_inv_no_NANs_H[11]))/(h*(n_N_H-11))
    
    n_H = 11  ## first n after inversion
    
    growth_rate_I_matrix_H[Rep,i] <- mean(growth_rate_I_H[n_H:length(growth_rate_I_H)])
    growth_rate_N_matrix_H[Rep,i] <- mean(growth_rate_N_H[n_H:length(growth_rate_N_H)])
    if () {
      crashed_simulations[i] = crashed_simulations[i] + 1
    }
  }
  mean_growth_rates_I_H[i] = mean(growth_rate_I_matrix_H[,i],na.rm = TRUE)
  mean_growth_rates_N_H[i] = mean(growth_rate_N_matrix_H[,i],na.rm = TRUE)
  mean_growth_I_H[i] = mean(overall_growth_I_H,na.rm = TRUE)
  mean_growth_N_H[i] = mean(overall_growth_N_H,na.rm = TRUE)
  initial_effective_dominance_H[i] = mean(initial_effective_dominance_matrix_H[,i], na.rm = TRUE)
}

par(mfrow = c(2,2))

plot(x = H, y = mean_growth_rates_I_H , col = "black", xlab = "Mutation dominance", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion", ylim = c(0,0.25))
points(x = H, y = mean_growth_rates_N_H, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = H, y = mean_growth_I_H , col = "black", xlab = "Mutation dominance", ylab = "Mutation growth rates", main = "Mean mutation growth rate after inversion",ylim = c(0,0.25))
points(x = H, y = mean_growth_N_H, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

# plot(x = H, y = initial_effective_dominance, col = "red", xlab = "Mutation dominance", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation dominances")


data <- data.frame(
  name = c(rep("0.0",5),rep("0.1",5),rep("0.25",5),rep("0.5",5)),
  value = c(growth_rate_I_matrix,growth_rate_N_matrix)
)

data %>% 
  ggplot( aes(x = name, y = value, fill = name)) +
   geom_boxplot()

df = data.frame()


boxplot(x = growth_rate_I_matrix_H, border = "black", names = H, xlab = "Mutation dominance", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion (I)",ylim = c(0,0.08))
boxplot(x = growth_rate_N_matrix_H, border = "blue", names = H, xlab = "Mutation dominance", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion (N)",ylim = c(0,0.08))

# boxplot(x = initial_effective_dominance_matrix, border = "red", names = H, xlab = "Mutation dominance", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation dominances")



for (i in (1:4)) {
  for (Rep in (1:5)) {
    plot(x = inversions$tick, y = inversions$dominance_variance, col = "green", type = "l", xlab = "Generations", ylab = "Dominance variance", main = "Evolution of the dominance variance",xlim = c(0,Tmax))
  }
}

plot(x = inversions$tick, y = inversions$effective_dominance, col = "red", type = "l", xlab = "Generations", ylab = "Effective dominance", main = "Evolution of the effective dominance",xlim = c(0,Tmax))








## ===========      Plotting all results  ===============


# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ 90%Inv .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithoutCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ HigherMutRate .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Before_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ After_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Recoded_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ N_1000_WithCrossingOver .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Without_migration .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ HigherSelectionMutation .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ S_INV = 0.1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ S_INV = 0.05 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ NoInversion Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ MoreTime Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ No Inversion_HigherSelectionMutation (SINV=0.2,Mu=1e-5) Rep 5 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ No Inversion_HigherSelectionMutation (SINV=0.1,Mu=1e-5) Rep 5 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Test_fitness Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Degenerescence_behaviour Rep 3 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Degenerescence_behaviour MoreTime (S=0.05) Rep 2 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Degenerescence_behaviour (N=4000) Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Mutation_Selection equilibrium (N=5000) (Mu=1e-5,S=-0.01) Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Inversion_frequency (N=2000)_Degenerescence_behaviour MoreTime (S=0.05, Mu=1e-5) Rep 1 (Seed = 9151250013170199950) Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Mutation_Selection equilibrium (N=5000) (Mu=1e-5,S=-0.01) Rep 1 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Freq Equilibrium 0.2-0.4 Rep 1 .csv"
File = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0"," S_INV ","0.05", " Mu ","1.0e-05 ","Rep ","1"," NoAnalysis .csv", sep ="")

inversions <- read.csv(File, stringsAsFactors = FALSE)

n_mut_per_inv_no_NANs = inversions$number_mutations_I
for (i in 1:length(n_mut_per_inv_no_NANs)) {
  if (inversions$number_mutations_I[i] == "NAN") {
    n_mut_per_inv_no_NANs[i] <- as.numeric(0)
  }
}

marginal_fitness_I = inversions$marginal_fitness_inversion_I
for (i in 1:length(marginal_fitness_I)) {
  if (marginal_fitness_I[i] == "NAN") {
    marginal_fitness_I[i] <- as.numeric(0)
  }
}

marginal_fitness_N = inversions$marginal_fitness_inversion_N
for (i in 1:length(marginal_fitness_N)) {
  if (marginal_fitness_N[i] == "NAN") {
    marginal_fitness_N[i] <- as.numeric(0)
  }
}

fitness_load_I_no_NANs = inversions$fitness_load_I
for (i in 1:length(fitness_load_I_no_NANs)) {
  if (fitness_load_I_no_NANs[i] == "NAN") {
    fitness_load_I_no_NANs[i] <- as.numeric(0)
  }
}

fitness_load_N_no_NANs = inversions$fitness_load_N
for (i in 1:length(fitness_load_N_no_NANs)) {
  if (fitness_load_N_no_NANs[i] == "NAN") {
    fitness_load_N_no_NANs[i] <- as.numeric(0)
  }
}




inverted_max_mut <- max(as.numeric(n_mut_per_inv_no_NANs))
inverted_max_fitness <- max(as.numeric(marginal_fitness_I))
non_inverted_max_mut <- max(as.numeric(inversions$number_mutations_N))
non_inverted_max_fitness <- max(as.numeric(marginal_fitness_N))
max = max(inverted_max_mut,non_inverted_max_mut)
max_fitness = max(inverted_max_fitness,non_inverted_max_fitness)


# Calculating the growth rate of mutations


n_I = length(inversions$number_mutations_I)
n_N = length(inversions$number_mutations_N)
growth_rate_I = (as.numeric(n_mut_per_inv_no_NANs [2:n_I]) - as.numeric(n_mut_per_inv_no_NANs [1:(n_I-1)]))/h
growth_rate_N = (as.numeric(inversions$number_mutations_N[2:n_N]) - as.numeric(inversions$number_mutations_N[1:(n_N-1)]))/h
growth_rate_per_mutation_I = growth_rate_I[11:(n_I-1)]/as.numeric(n_mut_per_inv_no_NANs[11:(n_I-1)])
growth_rate_per_mutation_N = growth_rate_N[11:(n_N-1)]/as.numeric(inversions$number_mutations_N[11:(n_N-1)])
mean(growth_rate_I)
mean(growth_rate_N)
mean(growth_rate_per_mutation_I,na.rm = T)
mean(growth_rate_per_mutation_N, na.rm = T)

dev.off()



par(mfrow = c(4,4))

plot(x = inversions$tick, y = inversions$freq_I, col = "red", type = "l", xlab = "Generations", ylab = "Frequency of the inversion", main = "Inversion frequency evolution",xlim = c(0,Tmax))

plot(x = inversions$tick, y = inversions$freq_II, col = "black",  type = "l", xlab = "Generations", ylab = "Arrangement haplotype proportions", main = "Evolution of the arrangement haplotypes",xlim = c(0,Tmax),ylim = c(0,1.05))
lines(x = inversions$tick, y = inversions$freq_IN, col = "violet", type = "l")
lines(x = inversions$tick, y = inversions$freq_NN, col = "blue", type = "l")
legend("topleft", legend=c("NN","IN","II"),
       col=c("blue", "violet","black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$mean_fitness_mutation_II, col = "black",  type = "l", xlab = "Generations", ylab = "Mean fitness caused by mutation", main = "Evolution of the mean fitness caused by mutations in arrangements",xlim = c(0,Tmax),ylim = c(0,1.05))
lines(x = inversions$tick, y = inversions$mean_fitness_mutation_IN, col = "violet", type = "l")
lines(x = inversions$tick, y = inversions$mean_fitness_mutation_NN, col = "blue", type = "l")
legend("topleft", legend=c("NN","IN","II"),
       col=c("blue", "violet","black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$fitness_load_I, col = "black", type = "l", xlab = "Generations", ylab = "Fitness load", main = "Evolution of the fitness load of haplosomes",xlim = c(0,Tmax))
lines(x = inversions$tick, y = inversions$fitness_load_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = inversions$number_mutations_I, col = "black", type = "l", xlab = "Generations", ylab = "Number of mutations in the inversion area", main = "Evolution of the number of mutations in the inversion area",xlim = c(0,Tmax),ylim = c(0,max*1.2))
lines(x = inversions$tick, y = inversions$number_mutations_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

T = length(inversions$tick)
plot(x = inversions$tick[1:(T-1)], y = growth_rate_I, col = "black", type = "l", xlab = "Generations", ylab = "Mutation growth rates", main = "Evolution of the mutation growth rates",xlim = c(0,Tmax),ylim = c(-0.05,0.2))
lines(x = inversions$tick[1:(T-1)], y = growth_rate_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$marginal_fitness_inversion_I, col = "black", type = "l", xlab = "Generations", ylab = "Marginal fitness of haplosomes", main = "Evolution of the marginal fitness of haplosomes",xlim = c(0,Tmax),ylim = c(1,max_fitness+0.005))
lines(x = inversions$tick, y = inversions$marginal_fitness_inversion_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = inversions$freq_homozygous_mut_II, col = "black",  type = "l", xlab = "Generations", ylab = "Homozygote frequency of mutations", main = "Evolution of the homozygote frequency of mutations",xlim = c(0,Tmax),ylim = c(0,1))
lines(x = inversions$tick, y = inversions$freq_homozygous_mut_IN, col = "violet",  type = "l")
lines(x = inversions$tick, y = inversions$freq_homozygous_mut_NN, col = "blue",  type = "l")
legend("topleft", legend=c("NN","IN","II"),
       col=c("blue", "violet","black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$mean_fitness_global_II, col = "black",  type = "l", xlab = "Generations", ylab = "Global fitness", main = "Evolution of the global fitness of arrangements",xlim = c(0,Tmax),ylim = c(0.1,1.5))
lines(x = inversions$tick, y = inversions$mean_fitness_global_IN, col = "violet",  type = "l")
lines(x = inversions$tick, y = inversions$mean_fitness_global_NN, col = "blue",  type = "l")
legend("topleft", legend=c("NN","IN","II"),
       col=c("blue", "violet","black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = as.numeric(inversions$marginal_fitness_inversion_I)*(1-as.numeric(fitness_load_I_no_NANs)), col = "black",  type = "l", xlab = "Generations", ylab = "Product of marginal fitnesses", main = "Evolution of the product of marginal fitnesses",xlim = c(0,Tmax),ylim = c(0.1,1.5))
lines(x = inversions$tick, y = as.numeric(inversions$marginal_fitness_inversion_N)*(1-as.numeric(fitness_load_N_no_NANs)), col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = inversions$covariance_mutation_inversion_I, col = "black", type = "l", xlab = "Generations", ylab = "Covariance mutation-inversion", main = "Evolution of the covariance mutation-inversion of haplosomes",xlim = c(0,Tmax),ylim = c(-0.02,0.02))
lines(x = inversions$tick, y = inversions$covariance_mutation_inversion_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = as.numeric(inversions$mean_fitness_global_II)*inversions$freq_I + as.numeric(inversions$mean_fitness_global_IN)*(1-inversions$freq_I), col = "black", type = "l", xlab = "Generations", ylab = "Marginal fitness", main = "Evolution of the marginal fitness",xlim = c(0,Tmax),ylim = c(0.1,1.5))
lines(x = inversions$tick, y = as.numeric(inversions$mean_fitness_global_IN)*inversions$freq_I + as.numeric(inversions$mean_fitness_global_NN)*(1-inversions$freq_I), col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

# plot(x = inversions$tick, y = inversions$covariance_out_of_marginal_fitness_I, col = "black", type = "l", xlab = "Generations", ylab = "Proportion of fitness explained by covariance", main = "Evolution of the proportion of fitness explained by covariance",xlim = c(0,Tmax),ylim = c(-0.02,0.02))
# lines(x = inversions$tick, y = inversions$covariance_out_of_marginal_fitness_N, col = "blue", type = "l")
# legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
#        col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = inversions$effective_selection, col = "orange", type = "l", xlab = "Generations", ylab = "Effective selection", main = "Evolution of the effective selection",xlim = c(0,Tmax))
plot(x = inversions$tick, y = inversions$effective_selection_formula, col = "red", type = "l", xlab = "Generations", ylab = "Effective selection with formula", main = "Evolution of the effective selection (obtained with formula)",xlim = c(0,Tmax))
plot(x = inversions$tick, y = inversions$gamma, col = "brown", type = "l", xlab = "Generations", ylab = "Gamma (overdominance)", main = "Evolution of gamma (overdominance)",xlim = c(0,Tmax))


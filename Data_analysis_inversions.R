rm(list = ls())
library(dplyr)

Tmax = 30000
log_interval = 5e3
q = 5   # Scaling factor
h = log_interval/q

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
File = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0"," S_INV ","0.05", " Mu ","1.0e-05 ","Rep ","1"," .csv", sep ="")

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
mean(growth_rate_I)
mean(growth_rate_N)

dev.off()

## Plotting summary statistics for mutation rate and S_INV
Mu = c("1.0e-06","5.0e-06","1.0e-05","5.0e-05")

growth_rate_I_matrix = matrix(nrow = 5, ncol = 4)
growth_rate_N_matrix = matrix(nrow = 5, ncol = 4)
initial_effective_dominance = 1:4
initial_effective_dominance_matrix = matrix(nrow = 5, ncol = 4)

overall_growth_I = 1:5
overall_growth_N = 1:5
mean_growth_rates_I = 1:4
mean_growth_rates_N = 1:4
mean_growth_I = 1:4
mean_growth_N = 1:4

for (i in (1:4)) {
  for (Rep in (1:5)) {
    File = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_s_0001_h_0_sinv_005_mig_0005_freqinv_005_scaling_5.csv Mu ",Mu[i]," Rep ",Rep," .csv", sep ="")
    
    inversions <- read.csv(File, stringsAsFactors = FALSE)
    
    n_mut_per_inv_no_NANs = inversions$number_mutations_I
    for (j in 1:length(n_mut_per_inv_no_NANs)) {
      if (inversions$number_mutations_I[j] == "NAN") {
        n_mut_per_inv_no_NANs[j] <- as.numeric(0)
      }
    }
    initial_effective_dominance_matrix[Rep,i] = mean(as.numeric(inversions$effective_dominance[11:12]),na.rm = TRUE)
    n_I = length(inversions$number_mutations_I)
    n_N = length(inversions$number_mutations_N)
    growth_rate_I = (as.numeric(n_mut_per_inv_no_NANs[2:n_I]) - as.numeric(n_mut_per_inv_no_NANs[1:(n_I-1)]))/h
    growth_rate_N = (as.numeric(inversions$number_mutations_N[2:n_N]) - as.numeric(inversions$number_mutations_N[1:(n_N-1)]))/h
    overall_growth_I[Rep] = (as.numeric(n_mut_per_inv_no_NANs[n_I]) - as.numeric(n_mut_per_inv_no_NANs[11]))/(h*(n_I-11))
    overall_growth_N[Rep] = (as.numeric(inversions$number_mutations_N[n_N]) - as.numeric(n_mut_per_inv_no_NANs[11]))/(h*(n_N-11))
    
    n = 11  ## first n after inversion
    
    growth_rate_I_matrix[Rep,i] <- mean(growth_rate_I[n:length(growth_rate_I)])
    growth_rate_N_matrix[Rep,i] <- mean(growth_rate_N[n:length(growth_rate_N)])
  }
  mean_growth_rates_I[i] = mean(growth_rate_I_matrix[,i])
  mean_growth_rates_N[i] = mean(growth_rate_N_matrix[,i])
  mean_growth_I[i] = mean(overall_growth_I)
  mean_growth_N[i] = mean(overall_growth_N)
  initial_effective_dominance[i] = mean(initial_effective_dominance_matrix[,i], na.rm = TRUE)
}

par(mfrow = c(2,3))

plot(x = Mu, y = mean_growth_rates_I , col = "black", xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(-0.05,0.25))
points(x = Mu, y = mean_growth_rates_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = Mu, y = mean_growth_I , col = "black", xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean mutation growth rate after inversion",ylim = c(-0.05,0.25))
points(x = Mu, y = mean_growth_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = Mu, y = initial_effective_dominance,col = "red", xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation rates",ylim = c(0,2))

boxplot(x = growth_rate_I_matrix, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
boxplot(x = growth_rate_N_matrix, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))

boxplot(x = initial_effective_dominance_matrix,border = "red", names = Mu, xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance after inversion")

## S_INV

S_INV = c("0.01","0.02","0.05","0.1")

initial_effective_dominance = 1:4
initial_effective_dominance_matrix = matrix(nrow = 5, ncol = 4)

for (i in (1:4)) {
  for (Rep in (1:5)) {
    File = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ",S_INV[i]," Mu ","1.0e-05"," Rep ",Rep," .csv", sep ="")
    
    inversions <- read.csv(File, stringsAsFactors = FALSE)
    
    n_mut_per_inv_no_NANs = inversions$number_mutations_I
    for (j in 1:length(n_mut_per_inv_no_NANs)) {
      if (inversions$number_mutations_I[j] == "NAN") {
        n_mut_per_inv_no_NANs[j] <- as.numeric(0)
      }
    }
    initial_effective_dominance_matrix[Rep,i] = mean(as.numeric(inversions$effective_dominance[11:12]),na.rm = TRUE)
    n_I = length(inversions$number_mutations_I)
    n_N = length(inversions$number_mutations_N)
    growth_rate_I = (as.numeric(n_mut_per_inv_no_NANs[2:n_I]) - as.numeric(n_mut_per_inv_no_NANs[1:(n_I-1)]))/h
    growth_rate_N = (as.numeric(inversions$number_mutations_N[2:n_N]) - as.numeric(inversions$number_mutations_N[1:(n_N-1)]))/h
    overall_growth_I[Rep] = (as.numeric(n_mut_per_inv_no_NANs[n_I]) - as.numeric(n_mut_per_inv_no_NANs[11]))/(h*(n_I-11))
    overall_growth_N[Rep] = (as.numeric(inversions$number_mutations_N[n_N]) - as.numeric(n_mut_per_inv_no_NANs[11]))/(h*(n_N-11))
    
    n = 11  ## first n after inversion
    
    growth_rate_I_matrix[Rep,i] <- mean(growth_rate_I[n:length(growth_rate_I)])
    growth_rate_N_matrix[Rep,i] <- mean(growth_rate_N[n:length(growth_rate_N)])
  }
  mean_growth_rates_I[i] = mean(growth_rate_I_matrix[,i],na.rm = TRUE)
  mean_growth_rates_N[i] = mean(growth_rate_N_matrix[,i],na.rm = TRUE)
  mean_growth_I[i] = mean(overall_growth_I,na.rm = TRUE)
  mean_growth_N[i] = mean(overall_growth_N,na.rm = TRUE)
  initial_effective_dominance[i] = mean(initial_effective_dominance_matrix[,i], na.rm = TRUE)
}

par(mfrow = c(2,3))

plot(x = S_INV, y = mean_growth_rates_I , col = "black", xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion", ylim = c(0,0.25))
points(x = S_INV, y = mean_growth_rates_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = S_INV, y = mean_growth_I , col = "black", xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean mutation growth rate after inversion",ylim = c(0,0.25))
points(x = S_INV, y = mean_growth_N, col = "blue")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = S_INV, y = initial_effective_dominance, col = "red", xlab = "Inversion benefit", ylab = "Initial effective dominance", main = "Initial effective dominance for different inversion benefits")

boxplot(x = growth_rate_I_matrix, border = "black", names = S_INV, xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(0,0.08))
boxplot(x = growth_rate_N_matrix, border = "blue", names = S_INV, xlab = "Inversion benefit", ylab = "Mutation growth rates", main = "Mean of mutation growth rates after inversion",ylim = c(0,0.08))

boxplot(x = initial_effective_dominance_matrix, border = "red", names = S_INV, xlab = "Inversion benefit", ylab = "Initial effective dominance", main = "Initial effective dominance for different inversion benefits")




for (i in (1:4)) {
  for (Rep in (1:5)) {
    plot(x = inversions$tick, y = inversions$dominance_variance, col = "green", type = "l", xlab = "Generations", ylab = "Dominance variance", main = "Evolution of the dominance variance",xlim = c(0,Tmax))
  }
}

plot(x = inversions$tick, y = inversions$effective_dominance, col = "red", type = "l", xlab = "Generations", ylab = "Effective dominance", main = "Evolution of the effective dominance",xlim = c(0,Tmax))


## Plotting all results

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
plot(x = inversions$tick, y = inversions$effective_dominance, col = "red", type = "l", xlab = "Generations", ylab = "Effective dominance", main = "Evolution of the effective dominance",xlim = c(0,Tmax))

plot(x = inversions$tick, y = inversions$dominance_variance, col = "green", type = "l", xlab = "Generations", ylab = "Dominance variance", main = "Evolution of the dominance variance",xlim = c(0,Tmax))

y = matrix(data = c(1,1.2,1.3,0,0.2,0.3),nrow = 3)
y = c(1,1.2,1.3)

boxplot(y ~ )
 
rm(list = ls())
library(dplyr)

Tmax = 30000

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
File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_K1_20e3_K2_20e3_chrom_length_3e6_inv_length_2e6_recrate_1e8_mutrate_1e8_s_0001_h_0_sinv_005_mig_0005_freqinv_005_scaling_5.csv Rep 1 .csv"

inversions <- read.csv(File, stringsAsFactors = FALSE)

n_mut_per_inv_no_NANs = inversions$n_mut_per_inv_haplosome
for (i in 1:length(n_mut_per_inv_no_NANs)) {
  if (inversions$n_mut_per_inv_haplosome[i] == "NAN") {
    n_mut_per_inv_no_NANs[i] <- as.numeric(0)
  }
}

marginal_fitness_I = inversions$inv_marginal_fitness_inv_inversion
for (i in 1:length(marginal_fitness_I)) {
  if (marginal_fitness_I[i] == "NAN") {
    marginal_fitness_I[i] <- as.numeric(0)
  }
}

marginal_fitness_N = inversions$non_inv_marginal_fitness_inversion
for (i in 1:length(marginal_fitness_N)) {
  if (marginal_fitness_N[i] == "NAN") {
    marginal_fitness_N[i] <- as.numeric(0)
  }
}




inverted_max_mut <- max(as.numeric(n_mut_per_inv_no_NANs))
inverted_max_fitness <- max(as.numeric(marginal_fitness_I))
non_inverted_max_mut <- max(as.numeric(inversions$n_mut_per_non_inv_haplosome))
non_inverted_max_fitness <- max(as.numeric(marginal_fitness_N))
max = max(inverted_max_mut,non_inverted_max_mut)
max_fitness = max(inverted_max_fitness,non_inverted_max_fitness)

## Plotting all results

par(mfrow = c(3,3))

plot(x = inversions$tick, y = inversions$inv_freq, col = "red", type = "l", xlab = "Generations", ylab = "Frequency of the inversion", main = "Inversion frequency evolution",xlim = c(0,Tmax))

plot(x = inversions$tick, y = inversions$mean_fitness_mutation_II, col = "black",  type = "l", xlab = "Generations", ylab = "Mean fitness caused by mutation", main = "Evolution of the mean fitness caused by mutations in arrangements",xlim = c(0,Tmax),ylim = c(0,1.05))
lines(x = inversions$tick, y = inversions$mean_fitness_mutation_IN, col = "violet", type = "l")
lines(x = inversions$tick, y = inversions$mean_fitness_mutation_NN, col = "blue", type = "l")
legend("topleft", legend=c("NN","IN","II"),
       col=c("blue", "violet","black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$fitness_load_I, col = "black", type = "l", xlab = "Generations", ylab = "Fitness load", main = "Evolution of the fitness load of haplosomes",xlim = c(0,Tmax),ylim = c(0,0.4))
lines(x = inversions$tick, y = inversions$fitness_load_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)

plot(x = inversions$tick, y = inversions$n_mut_per_inv_haplosome, col = "black", type = "l", xlab = "Generations", ylab = "Number of mutations in the inversion area", main = "Evolution of the number of mutations in the inversion area",xlim = c(0,Tmax),ylim = c(0,max*1.2))
lines(x = inversions$tick, y = inversions$n_mut_per_non_inv_haplosome, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


plot(x = inversions$tick, y = inversions$inv_marginal_fitness_inv_inversion, col = "black", type = "l", xlab = "Generations", ylab = "Marginal fitness of haplosomes", main = "Evolution of the marginal fitness of haplosomes",xlim = c(0,Tmax),ylim = c(1,max_fitness+0.005))
lines(x = inversions$tick, y = inversions$non_inv_marginal_fitness_inversion, col = "blue", type = "l")
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

plot(x = inversions$tick, y = inversions$covariance_mutation_inversion_I, col = "black", type = "l", xlab = "Generations", ylab = "Covariance mutation-inversion", main = "Evolution of the covariance mutation-inversion of haplosomes",xlim = c(0,Tmax),ylim = c(-0.02,0.02))
lines(x = inversions$tick, y = inversions$covariance_mutation_inversion_N, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


rm(list = ls())
library(dplyr)
library()

## File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ 90%Inv .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithoutCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ HigherMutRate .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ Test_m2 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Before_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ After_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Recoded_recombination_function .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ N_1000_WithCrossingOver .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ Without_migration .csv"
File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ HigherSelectionMutation .csv"



inversions <- read.csv(File, stringsAsFactors = FALSE)

for (i in 1:length(inversions$n_mut_per_inv_haplosome)) {
  if (inversions$n_mut_per_inv_haplosome[i] == "NAN") {
    inversions$n_mut_per_inv_haplosome[i] <- as.numeric(0)
  } else {
    inversions$n_mut_per_inv_haplosome[i] <- as.numeric(inversions$n_mut_per_inv_haplosome[i] )
  }
}

inverted_max_mut <- max(as.numeric(inversions$n_mut_per_inv_haplosome))
non_inverted_max_mut <- max(as.numeric(inversions$n_mut_per_non_inv_haplosome))
max = max(inverted_max_mut,non_inverted_max_mut)

par(mfrow = c(1,2))
plot(x = inversions$tick, y = inversions$inv_freq, col = "red", type = "l", xlab = "Generations", ylab = "Frequency of the inversion", main = "Inversion frequency evolution",xlim = c(0,200000))
plot(x = inversions$tick, y = inversions$n_mut_per_inv_haplosome, col = "black", type = "l", xlab = "Generations", ylab = "Number of mutations in the inversion area", main = "Evolution of the number of mutations in the inversion area",xlim = c(0,200000),ylim = c(0,max*1.2))
lines(x = inversions$tick, y = inversions$n_mut_per_non_inv_haplosome, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


rm(list = ls())
library(dplyr)

## File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ 90%Inv .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithoutCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ HigherMutRate .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ WithCrossover .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ Test_m2 .csv"
File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ N=5000 .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/freq_ N_1000_WithCrossingOver .csv"

inversions <- read.csv(File, stringsAsFactors = FALSE)
inversions$n_mut_per_inv_haplosome[is.nan(inversions$n_mut_per_inv_haplosome)] <- 0

par(mfrow = c(1,2))
plot(x = inversions$tick, y = inversions$inv_freq, col = "red", type = "l", xlab = "Generations", ylab = "Frequency of the inversion", main = "Inversion frequency evolution")
plot(x = inversions$tick, y = inversions$n_mut_per_inv_haplosome, col = "black", type = "l", xlab = "Generations", ylab = "Number of mutations in the inversion area", main = "Evolution of the number of mutations in the inversion area")
lines(x = inversions$tick, y = inversions$n_mut_per_non_inv_haplosome, col = "blue", type = "l")
legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
       col=c("blue", "black"), lty=1:2, cex=0.8)


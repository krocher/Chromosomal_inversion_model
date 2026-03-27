rm(list = ls())
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(cowplot)

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

Mu = c("1.0e-06","5.0e-06","1.0e-05","5.0e-05","0.0001","0.0005")

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
                 extinction_time = character(),
                 Rep = integer(),
                 Mu = character(),
                 number_mutations_no_NANs_I = numeric(),
                 growth_rate_I = numeric(),
                 growth_rate_N = numeric(),
                 growth_rate_per_mutation_I = numeric(),
                 growth_rate_per_mutation_N = numeric(),
                 decrease_rate = numeric(),
                 decrease_rate_per_mutation = numeric(),
                 stringsAsFactors=FALSE) 



for (i in (1:length(Mu))) {
  for (Rep in (1:n_rep)) {
    File_Mu = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ", "0.05" ," Mu ",Mu[i]," Rep ",Rep," MutationRateSensitivityAnalysis"," .csv", sep ="")
    
    aux_inversions_Mu <- read.csv(File_Mu, stringsAsFactors = FALSE)
    
    n_I_Mu = length(aux_inversions_Mu$number_mutations_I)
    n_N_Mu = length(aux_inversions_Mu$number_mutations_N)
    n_inv_Mu = length(aux_inversions_Mu$freq_I)
    
      aux_inversions_Mu <- aux_inversions_Mu %>%
        mutate(Rep = Rep,
               Mu = as.numeric(Mu[i]),
               number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
               growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_Mu] - number_mutations_no_NANs_I[1:(n_I_Mu-1)])/h,NA),
               growth_rate_N = c((number_mutations_N[2:n_N_Mu] - number_mutations_N[1:(n_N_Mu-1)])/h,NA),
               growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
               growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
               decrease_rate = c(abs(freq_I[2:n_inv_Mu] - freq_I[1:(n_inv_Mu-1)])/h,NA),
               decrease_rate_per_mutation = decrease_rate/freq_I)
      
      inversions_Mu <- inversions_Mu %>%
        rbind(aux_inversions_Mu)
  
  }}


n_Mu = 11  ## first n after inversion

summary_by_Rep_Mu <- inversions_Mu %>% 
  group_by(Mu,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_Mu:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_Mu:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_Mu+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_Mu+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_Mu])/(h*((length(number_mutations_no_NANs_I)-1)-n_Mu)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_Mu])/(h*(length(number_mutations_N)-1-n_Mu)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[n_Mu:(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[n_Mu:(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time)) %>% 
  arrange(Mu) %>% 
  mutate(Mu = as.factor(Mu))


summary_by_Rep_I_Mu <- summary_by_Rep_Mu[c("Mu","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted")]
summary_by_Rep_N_Mu <- summary_by_Rep_Mu[c("Mu","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted")]

summary_by_Rep_I_Mu  <- summary_by_Rep_I_Mu %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_Mu <- summary_by_Rep_N_Mu %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_Mu = rbind(summary_by_Rep_I_Mu,summary_by_Rep_N_Mu) 
pivot_summary_by_Rep_Mu <- pivot_summary_by_Rep_Mu %>% 
  mutate(Mu = as.factor(Mu))

summary_Mu <- summary_by_Rep_Mu %>% 
  group_by(Mu) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(Mu)

  #   # initial_effective_dominance_matrix_Mu[Rep,i] = mean(as.numeric(inversions_Mu$effective_dominance[11:12]),na.rm = TRUE)
  # # initial_effective_dominance_Mu[i] = mean(initial_effective_dominance_matrix_Mu[,i], na.rm = TRUE)

# plot(x = Mu, y = initial_effective_dominance,col = "red", xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation rates",ylim = c(0,2))

pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=overall_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Overall growth rate for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Overall mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates per mutation for different mutation rates") +
 labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_Mu %>% ggplot(aes(x = Mu,y=mean_decrease_rate_per_mutation_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean decrease rates per mutation for different mutation rates") +
  labs(x = "Mutation rate", y = "Mean mutation decrease rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_Mu %>% ggplot(aes(x = Mu,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different mutation rates") +
  labs(x = "Mutation rate", y = "Extinction time") +
  theme_cowplot(12)



# boxplot(x = summary_by_Rep$mean_growth_I_by_Rep, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# boxplot(x = summary_by_Rep$mean_growth_N_by_Rep, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# 
# boxplot(x = growth_rate_per_mutation_I_matrix_Mu, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = growth_rate_per_mutation_N_matrix_Mu, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = initial_effective_dominance_matrix,border = "red", names = Mu, xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance after inversion")

## 2) Inversion benefit

S_INV = c("0.01","0.02","0.05","0.1")

inversions_S_INV <- data.frame(tick = integer(),
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
                            extinction_time = character(),
                            Rep = integer(),
                            S_INV = character(),
                            number_mutations_no_NANs_I = numeric(),
                            growth_rate_I = numeric(),
                            growth_rate_N = numeric(),
                            growth_rate_per_mutation_I = numeric(),
                            growth_rate_per_mutation_N = numeric(),
                            decrease_rate = numeric(),
                            decrease_rate_per_mutation = numeric(),
                            stringsAsFactors=FALSE) 



for (i in (1:length(S_INV))) {
  for (Rep in (1:n_rep)) {
    File_S_INV = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ",S_INV[i] ," Mu ","1.0e-05"," Rep ",Rep," InversionBenefitSensitivityAnalysis"," .csv", sep ="")
    
    aux_inversions_S_INV <- read.csv(File_S_INV, stringsAsFactors = FALSE)
    
    n_I_S_INV = length(aux_inversions_S_INV$number_mutations_I)
    n_N_S_INV = length(aux_inversions_S_INV$number_mutations_N)
    n_inv_S_INV = length(aux_inversions_S_INV$freq_I)
    
    aux_inversions_S_INV <- aux_inversions_S_INV %>%
      mutate(Rep = Rep,
             S_INV = as.numeric(S_INV[i]),
             number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
             growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_S_INV] - number_mutations_no_NANs_I[1:(n_I_S_INV-1)])/h,NA),
             growth_rate_N = c((number_mutations_N[2:n_N_S_INV] - number_mutations_N[1:(n_N_S_INV-1)])/h,NA),
             growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
             growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
             decrease_rate = c(abs(freq_I[2:n_inv_S_INV] - freq_I[1:(n_inv_S_INV-1)])/h,NA),
             decrease_rate_per_mutation = decrease_rate/freq_I)
    
    inversions_S_INV <- inversions_S_INV %>%
      rbind(aux_inversions_S_INV)
    
  }}


n_S_INV = 11  ## first n after inversion

summary_by_Rep_S_INV <- inversions_S_INV %>% 
  group_by(S_INV,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_S_INV:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_S_INV:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[n_S_INV:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[n_S_INV:(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_S_INV])/(h*((length(number_mutations_no_NANs_I)-1)-n_S_INV)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_S_INV])/(h*(length(number_mutations_N)-1-n_S_INV)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[n_S_INV:(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[n_S_INV:(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time)) %>% 
  arrange(S_INV) %>% 
  mutate(S_INV = as.factor(S_INV))


summary_by_Rep_I_S_INV <- summary_by_Rep_S_INV[c("S_INV","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted")]
summary_by_Rep_N_S_INV <- summary_by_Rep_S_INV[c("S_INV","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted")]

summary_by_Rep_I_S_INV  <- summary_by_Rep_I_S_INV %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_S_INV <- summary_by_Rep_N_S_INV %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_S_INV = rbind(summary_by_Rep_I_S_INV,summary_by_Rep_N_S_INV) 
pivot_summary_by_Rep_S_INV <- pivot_summary_by_Rep_S_INV %>% 
  mutate(S_INV = as.factor(S_INV))

summary_S_INV <- summary_by_Rep_S_INV %>% 
  group_by(S_INV) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(S_INV)

#   # initial_effective_dominance_matrix_Mu[Rep,i] = mean(as.numeric(inversions_Mu$effective_dominance[11:12]),na.rm = TRUE)
# # initial_effective_dominance_Mu[i] = mean(initial_effective_dominance_matrix_Mu[,i], na.rm = TRUE)

# plot(x = Mu, y = initial_effective_dominance,col = "red", xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation rates",ylim = c(0,2))

pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=overall_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Overall growth rate for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Overall mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates per mutation for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_S_INV %>% ggplot(aes(x = S_INV,y=mean_decrease_rate_per_mutation_by_Rep)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean decrease rates per mutation for different mutation rates") +
  labs(x = "Mutation rate", y = "Mean mutation decrease rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_S_INV %>% ggplot(aes(x = S_INV,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean decrease rates per mutation for different mutation rates") +
  labs(x = "Mutation rate", y = "Extinction time")+
  theme_cowplot(12)

# boxplot(x = initial_effective_dominance_matrix, border = "red", names = S_INV, xlab = "Inversion benefit", ylab = "Initial effective dominance", main = "Initial effective dominance for different inversion benefits")



## 3) Deleterious mutation dominance

H = c("0.0","0.1","0.25","0.5")

inversions_H <- data.frame(tick = integer(),
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
                            extinction_time = character(),
                            Rep = integer(),
                            H = character(),
                            number_mutations_no_NANs_I = numeric(),
                            growth_rate_I = numeric(),
                            growth_rate_N = numeric(),
                            growth_rate_per_mutation_I = numeric(),
                            growth_rate_per_mutation_N = numeric(),
                            decrease_rate = numeric(),
                            decrease_rate_per_mutation = numeric(),
                            stringsAsFactors=FALSE) 



for (i in (1:length(H))) {
  for (Rep in (1:n_rep)) {
    File_H = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ",H[i], " S_INV ", "0.05" ," Mu ","1.0e-05"," Rep ",Rep," MutationDominanceSensitivityAnalysis"," .csv", sep ="")
    
    aux_inversions_H <- read.csv(File_H, stringsAsFactors = FALSE)
    
    n_I_H = length(aux_inversions_H$number_mutations_I)
    n_N_H = length(aux_inversions_H$number_mutations_N)
    n_inv_H = length(aux_inversions_H$freq_I)
    
    aux_inversions_H <- aux_inversions_H %>%
      mutate(Rep = Rep,
             H = as.numeric(H[i]),
             number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
             growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_H] - number_mutations_no_NANs_I[1:(n_I_H-1)])/h,NA),
             growth_rate_N = c((number_mutations_N[2:n_N_H] - number_mutations_N[1:(n_N_H-1)])/h,NA),
             growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
             growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
             decrease_rate = c(abs(freq_I[2:n_inv_H] - freq_I[1:(n_inv_H-1)])/h,NA),
             decrease_rate_per_mutation = decrease_rate/freq_I)
    
    inversions_H <- inversions_H %>%
      rbind(aux_inversions_H)
    
  }}


n_H = 11  ## first n after inversion

summary_by_Rep_H <- inversions_H %>% 
  group_by(H,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_H:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_H:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[n_H:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[n_H:(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_H])/(h*((length(number_mutations_no_NANs_I)-1)-n_H)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_Mu])/(h*(length(number_mutations_N)-1-n_Mu)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[n_H:(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[n_H:(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time)) %>% 
  arrange(H) %>% 
  mutate(H = as.factor(H))


summary_by_Rep_I_H <- summary_by_Rep_H[c("H","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted")]
summary_by_Rep_N_H <- summary_by_Rep_H[c("H","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted")]

summary_by_Rep_I_H  <- summary_by_Rep_I_H %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_H <- summary_by_Rep_N_H %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_H = rbind(summary_by_Rep_I_H,summary_by_Rep_N_H) 
pivot_summary_by_Rep_H <- pivot_summary_by_Rep_H %>% 
  mutate(H = as.factor(H))

summary_H <- summary_by_Rep_H %>% 
  group_by(H) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(H)

#   # initial_effective_dominance_matrix_Mu[Rep,i] = mean(as.numeric(inversions_Mu$effective_dominance[11:12]),na.rm = TRUE)
# # initial_effective_dominance_Mu[i] = mean(initial_effective_dominance_matrix_Mu[,i], na.rm = TRUE)

# plot(x = Mu, y = initial_effective_dominance,col = "red", xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance for different mutation rates",ylim = c(0,2))

pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=overall_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Overall growth rate for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Overall mutation growth rates") +
  theme_cowplot(12)

pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates per mutation for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_H %>% ggplot(aes(x = H,y=mean_decrease_rate_per_mutation_by_Rep)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean decrease rates per mutation for different mutation rates") +
  labs(x = "Mutation rate", y = "Mean mutation decrease rates per mutation") +
  theme_cowplot(12)

summary_by_Rep_H %>% ggplot(aes(x = H,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean decrease rates per mutation for different mutation rates") +
  labs(x = "Mutation rate", y = "Extinction time") +
  theme_cowplot(12)



# boxplot(x = summary_by_Rep$mean_growth_I_by_Rep, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# boxplot(x = summary_by_Rep$mean_growth_N_by_Rep, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates", main = "Mean of mutation growth rates for different mutation rates",ylim = c(0,0.1))
# 
# boxplot(x = growth_rate_per_mutation_I_matrix_Mu, border = "black", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = growth_rate_per_mutation_N_matrix_Mu, border = "blue", names = Mu, xlab = "Mutation rate", ylab = "Mutation growth rates per mutation", main = "Mean of mutation growth rates per mutation for different mutation rates",ylim = c(0,0.1))
# boxplot(x = initial_effective_dominance_matrix,border = "red", names = Mu, xlab = "Mutation rate", ylab = "Initial effective dominance", main = "Initial effective dominance after inversion")




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
n_inv = length(inversions$inv_freq)
growth_rate_I = (as.numeric(n_mut_per_inv_no_NANs [2:n_I]) - as.numeric(n_mut_per_inv_no_NANs [1:(n_I-1)]))/h
growth_rate_N = (as.numeric(inversions$number_mutations_N[2:n_N]) - as.numeric(inversions$number_mutations_N[1:(n_N-1)]))/h
growth_rate_per_mutation_I = growth_rate_I[11:(n_I-1)]/as.numeric(n_mut_per_inv_no_NANs[11:(n_I-1)])
growth_rate_per_mutation_N = growth_rate_N[11:(n_N-1)]/as.numeric(inversions$number_mutations_N[11:(n_N-1)])
decrease_rate = abs((as.numeric(inversions$inv_freq[11:n_inv]) - as.numeric(inversions$inv_freq[11:(n_inv-1)]))/h)
decrease_rate_per_mutation = decrease_rate/as.numeric(inversions$inv_freq)
mean(growth_rate_I)
mean(growth_rate_N)
mean(growth_rate_per_mutation_I,na.rm = T)
mean(growth_rate_per_mutation_N, na.rm = T)

dev.off()



par(mfrow = c(1,2))

plot(x = inversions$tick, y = inversions$freq_I, col = "red", type = "l", xlab = "Generations", ylab = "Frequency of the inversion", main = "Inversion frequency evolution",xlim = c(0,Tmax), ylim = c(0,1.3))

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

s <- as.numeric(inversions$effective_selection)
m <- 0.005
gamma <- as.numeric(inversions$gamma)

## Equilibria ( when s > 0 )

eq_freq_0 = rep(0,23)
eq_freq_1 = (-2*s-3*gamma+m*gamma+sqrt(gamma^2+m^2*gamma^2+2*m*(2*s^2+gamma*(4+gamma)+s*(2+4*gamma))))/(2*(m-1)*(s + 2*gamma))
eq_freq_2 = -(2*s+3*gamma-m*gamma+sqrt(gamma^2+m^2*gamma^2+2*m*(2*s^2+gamma*(4+gamma)+s*(2+4*gamma))))/(2*(m-1)*(s + 2*gamma))

lines(x = inversions$tick, y = eq_freq_0, col = "grey", type = "l", xlab = "Generations", ylab = "QLE frequency", main = "Evolution of QLE inversion frequency",xlim = c(0,Tmax), ylim = c(0,1.5))
lines(x = inversions$tick, y = eq_freq_1, col = "pink")
lines(x = inversions$tick, y = eq_freq_2, col = "green")

legend("topleft", legend=c("Equilibrium frequency 0", "Equilibrium frequency 1", "Equilibrium frequency 2"),
       col=c("grey", "pink","green"), lty=1:2, cex=0.8)

rm(list = ls())
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)

Tmax = 30000
log_interval = 5e3
q = 5   # Scaling factor
h = log_interval/q
n_rep = 2
inv_length = 2e3

three_cols = c("black","violet","blue")
two_cols = c("black","blue")
four_cols = c("grey","pink","green","red")



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
                 stringsAsFactors = FALSE) 



for (i in (1:length(Mu))) {
  for (Rep in (1:n_rep)) {
    File_Mu = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ", "0.05" ," Mu ",Mu[i]," Rep ",Rep," MutationRateSensitivityAnalysis"," 1"," .csv", sep ="")
    
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
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_Mu+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_Mu+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_Mu:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_Mu:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(Mu) %>% 
  mutate(Mu = as.factor(Mu))


summary_by_Rep_I_Mu <- summary_by_Rep_Mu[c("Mu","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_Mu <- summary_by_Rep_Mu[c("Mu","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_Mu  <- summary_by_Rep_I_Mu %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_Mu <- summary_by_Rep_N_Mu %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
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


## Plotting summary statistics

mean_growth_Mu <- pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_Mu <- pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different mutation rates") +
  labs(fill = "Arrangement", x = "Mutation rate", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_Mu <- pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different mutation rates") +
#   labs(fill = "Arrangement", x = "Mutation rate", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_Mu <- pivot_summary_by_Rep_Mu %>%  ggplot(aes(x = Mu,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different mutation rates") +
#  labs(fill = "Arrangement", x = "Mutation rate", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_Mu <- summary_by_Rep_Mu %>% ggplot(aes(x = Mu,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different mutation rates") +
  labs(x = "Mutation rate", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)



# mean_decrease_per_mutation_Mu <- summary_by_Rep_Mu %>% ggplot(aes(x = Mu,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different mutation rates") +
#   labs(x = "Mutation rate", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_Mu <- summary_by_Rep_Mu %>% ggplot(aes(x = Mu,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different mutation rates") +
  labs(x = "Mutation rate", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_Mu,proportion_mutation_Mu,mean_decrease_Mu,extinction_time_Mu,nrow = 4)


## 2) Inversion benefit

S_INV = c("0.001","0.005","0.01","0.05","0.1")

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
                            stringsAsFactors = FALSE) 



for (i in (1:length(S_INV))) {
  for (Rep in (1:n_rep)) {
    File_S_INV = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0", " S_INV ",S_INV[i] ," Mu ","1.0e-05"," Rep ",Rep," InversionBenefitSensitivityAnalysis"," 2 ","ReverseMutations F"," .csv", sep ="")
    
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
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_S_INV+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_S_INV+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_S_INV])/(h*((length(number_mutations_no_NANs_I)-1)-n_S_INV)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_S_INV])/(h*(length(number_mutations_N)-1-n_S_INV)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_S_INV+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_S_INV+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_S_INV:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_S_INV:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(S_INV) %>% 
  mutate(S_INV = as.factor(S_INV))


summary_by_Rep_I_S_INV <- summary_by_Rep_S_INV[c("S_INV","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_S_INV <- summary_by_Rep_S_INV[c("S_INV","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_S_INV  <- summary_by_Rep_I_S_INV %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_S_INV <- summary_by_Rep_N_S_INV %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
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


## Plotting summary statistics

mean_growth_S_INV <- pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different inversion benefits") +
  labs(fill = "Arrangement", x = "Inversion benefit", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_S_INV <- pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different inversion benefits") +
  labs(fill = "Arrangement", x = "Inversion benefit", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_S_INV <- pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different inversion benefits") +
#   labs(fill = "Arrangement", x = "Inversion benefit", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_S_INV <- pivot_summary_by_Rep_S_INV %>%  ggplot(aes(x = S_INV,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different inversion benefits") +
#  labs(fill = "Arrangement", x = "Inversion benefit", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_S_INV <- summary_by_Rep_S_INV %>% ggplot(aes(x = S_INV,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different inversion benefits") +
  labs(x = "Inversion benefit", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)

# mean_decrease_per_mutation_S_INV <- summary_by_Rep_S_INV %>% ggplot(aes(x = S_INV,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different inversion benefits") +
#   labs(x = "Inversion benefit", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_S_INV <- summary_by_Rep_S_INV %>% ggplot(aes(x = S_INV,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different inversion benefits") +
  labs(x = "Inversion benefit", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_S_INV,proportion_mutation_S_INV,mean_decrease_S_INV,extinction_time_S_INV,nrow = 4)


## 3) Deleterious mutation dominance

H = c("0.0","0.01","0.05","0.1","0.25","0.5")

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
                            stringsAsFactors = FALSE) 



for (i in (1:length(H))) {
  for (Rep in (1:n_rep)) {
    File_H = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ",H[i], " S_INV ", "0.05" ," Mu ","1.0e-05"," Rep ",Rep," MutationDominanceSensitivityAnalysis"," 2"," ReverseMutations F"," .csv", sep ="")
    
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
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_H+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_H+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_H])/(h*((length(number_mutations_no_NANs_I)-1)-n_H)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_H])/(h*(length(number_mutations_N)-1-n_H)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_H+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_H+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_H:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_H:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(H) %>% 
  mutate(H = as.factor(H))


summary_by_Rep_I_H <- summary_by_Rep_H[c("H","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_H <- summary_by_Rep_H[c("H","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_H  <- summary_by_Rep_I_H %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_H <- summary_by_Rep_N_H %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
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


## Plotting summary statistics

mean_growth_H <- pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation dominances") +
  labs(fill = "Arrangement", x = "Mutation dominance", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_H <- pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different mutation dominances") +
  labs(fill = "Arrangement", x = "Mutation dominance", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_H <- pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different mutation dominances") +
#   labs(fill = "Arrangement", x = "Mutation dominance", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_H <- pivot_summary_by_Rep_H %>%  ggplot(aes(x = H,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different mutation dominances") +
#  labs(fill = "Arrangement", x = "Mutation dominance", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_H <- summary_by_Rep_H %>% ggplot(aes(x = H,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different mutation dominances") +
  labs(x = "Mutation dominance", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)



# mean_decrease_per_mutation_H <- summary_by_Rep_H %>% ggplot(aes(x = H,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different mutation dominances") +
#   labs(x = "Mutation dominance", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_H <- summary_by_Rep_H %>% ggplot(aes(x = H,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different mutation dominances") +
  labs(x = "Mutation dominance", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_H,proportion_mutation_H,mean_decrease_H,extinction_time_H,nrow = 4)


## 4) Initial inversion proportion (at T_INV)

I = c("0.02","0.05","0.1","0.25","0.5","0.75")

inversions_I <- data.frame(tick = integer(),
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
                           I = character(),
                           number_mutations_no_NANs_I = numeric(),
                           growth_rate_I = numeric(),
                           growth_rate_N = numeric(),
                           growth_rate_per_mutation_I = numeric(),
                           growth_rate_per_mutation_N = numeric(),
                           decrease_rate = numeric(),
                           decrease_rate_per_mutation = numeric(),
                           stringsAsFactors = FALSE) 



for (i in (1:length(I))) {
  for (Rep in (1:n_rep)) {
    File_I = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_R_1e-5_S_0001_m_0005_scaling_5 I ",I[i]," H ","0.0", " S_INV ", "0.05" ," Mu ","1.0e-05"," Rep ",Rep," InversionProportionSensitivityAnalysis"," 2"," ReverseMutations F"," .csv", sep ="")
    
    aux_inversions_I <- read.csv(File_I, stringsAsFactors = FALSE)
    
    n_I_I = length(aux_inversions_I$number_mutations_I)
    n_N_I = length(aux_inversions_I$number_mutations_N)
    n_inv_I = length(aux_inversions_I$freq_I)
    
    aux_inversions_I <- aux_inversions_I %>%
      mutate(Rep = Rep,
             I = as.numeric(I[i]),
             number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
             growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_I] - number_mutations_no_NANs_I[1:(n_I_I-1)])/h,NA),
             growth_rate_N = c((number_mutations_N[2:n_N_I] - number_mutations_N[1:(n_N_I-1)])/h,NA),
             growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
             growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
             decrease_rate = c(abs(freq_I[2:n_inv_I] - freq_I[1:(n_inv_I-1)])/h,NA),
             decrease_rate_per_mutation = decrease_rate/freq_I)
    
    inversions_I <- inversions_I %>%
      rbind(aux_inversions_I)
    
  }}


n_I = 11  ## first n after inversion

summary_by_Rep_I <- inversions_I %>% 
  group_by(I,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_I:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_I:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_I+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_I+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_I])/(h*((length(number_mutations_no_NANs_I)-1)-n_I)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_I])/(h*(length(number_mutations_N)-1-n_I)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_I+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_I+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_I:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_I:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(I) %>% 
  mutate(I = as.factor(I))


summary_by_Rep_I_I <- summary_by_Rep_I[c("I","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_I <- summary_by_Rep_I[c("I","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_I  <- summary_by_Rep_I_I %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_I <- summary_by_Rep_N_I %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_I = rbind(summary_by_Rep_I_I,summary_by_Rep_N_I) 
pivot_summary_by_Rep_I <- pivot_summary_by_Rep_I %>% 
  mutate(I = as.factor(I))

summary_I <- summary_by_Rep_I %>% 
  group_by(I) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(I)


## Plotting summary statistics

mean_growth_I <- pivot_summary_by_Rep_I %>%  ggplot(aes(x = I,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different initial inversion proportions") +
  labs(fill = "Arrangement", x = "Initial inversion proportion", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_I <- pivot_summary_by_Rep_I %>%  ggplot(aes(x = I,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different initial inversion proportions") +
  labs(fill = "Arrangement", x = "Initial inversion proportion", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_I <- pivot_summary_by_Rep_I %>%  ggplot(aes(x = I,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different initial inversion proportions") +
#   labs(fill = "Arrangement", x = "Initial inversion proportion", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_I <- pivot_summary_by_Rep_I %>%  ggplot(aes(x = I,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different initial inversion proportions") +
#  labs(fill = "Arrangement", x = "Initial inversion proportion", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_I <- summary_by_Rep_I %>% ggplot(aes(x = I,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different initial inversion proportions") +
  labs(x = "Initial inversion proportion", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)



# mean_decrease_per_mutation_I <- summary_by_Rep_I %>% ggplot(aes(x = I,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different initial inversion proportions") +
#   labs(x = "Initial inversion proportion", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_I <- summary_by_Rep_I %>% ggplot(aes(x = I,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different initial inversion proportions") +
  labs(x = "Initial inversion proportion", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_I,proportion_mutation_I,mean_decrease_I,extinction_time_I,nrow = 4)


## 5) Mutation selection coefficient

S = c("-0.0005","-0.001","-0.005","-0.01","-0.05","-0.1")

inversions_S <- data.frame(tick = integer(),
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
                           S = character(),
                           number_mutations_no_NANs_I = numeric(),
                           growth_rate_I = numeric(),
                           growth_rate_N = numeric(),
                           growth_rate_per_mutation_I = numeric(),
                           growth_rate_per_mutation_N = numeric(),
                           decrease_rate = numeric(),
                           decrease_rate_per_mutation = numeric(),
                           stringsAsFactors = FALSE) 



for (i in (1:length(S))) {
  for (Rep in (1:n_rep)) {
    File_S = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_R_1e-5_m_0005_scaling_5 S ",S[i]," I ","0.05"," H ","0.0", " S_INV ", "0.05" ," Mu ","1.0e-05"," Rep ",Rep," MutationSelectionCoefficientSensitivityAnalysis"," 2"," ReverseMutations F"," .csv", sep ="")
    
    aux_inversions_S <- read.csv(File_S, stringsAsFactors = FALSE)
    
    n_I_S = length(aux_inversions_S$number_mutations_I)
    n_N_S = length(aux_inversions_S$number_mutations_N)
    n_inv_S = length(aux_inversions_S$freq_I)
    
    aux_inversions_S <- aux_inversions_S %>%
      mutate(Rep = Rep,
             S = as.numeric(S[i]),
             number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
             growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_S] - number_mutations_no_NANs_I[1:(n_I_S-1)])/h,NA),
             growth_rate_N = c((number_mutations_N[2:n_N_S] - number_mutations_N[1:(n_N_S-1)])/h,NA),
             growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
             growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
             decrease_rate = c(abs(freq_I[2:n_inv_S] - freq_I[1:(n_inv_S-1)])/h,NA),
             decrease_rate_per_mutation = decrease_rate/freq_I)
    
    inversions_S <- inversions_S %>%
      rbind(aux_inversions_S)
    
  }}


n_S = 11  ## first n after inversion

summary_by_Rep_S <- inversions_S %>% 
  group_by(S,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_S:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_S:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_S+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_S+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_S])/(h*((length(number_mutations_no_NANs_I)-1)-n_S)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_S])/(h*(length(number_mutations_N)-1-n_S)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_S+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_S+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_S:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_S:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(S) %>% 
  mutate(S = as.factor(S))


summary_by_Rep_I_S <- summary_by_Rep_S[c("S","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_S <- summary_by_Rep_S[c("S","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_S  <- summary_by_Rep_I_S %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_S <- summary_by_Rep_N_S %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_S = rbind(summary_by_Rep_I_S,summary_by_Rep_N_S) 
pivot_summary_by_Rep_S <- pivot_summary_by_Rep_S %>% 
  mutate(S = as.factor(S))

summary_S <- summary_by_Rep_S %>% 
  group_by(S) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(S)


## Plotting summary statistics

mean_growth_S <- pivot_summary_by_Rep_S %>%  ggplot(aes(x = S,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different mutation selection coefficients") +
  labs(fill = "Arrangement", x = "Mutation selection coefficient", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_S <- pivot_summary_by_Rep_S %>%  ggplot(aes(x = S,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different mutation selection coefficients") +
  labs(fill = "Arrangement", x = "Mutation selection coefficient", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_S <- pivot_summary_by_Rep_S %>%  ggplot(aes(x = S,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different mutation selection coefficients") +
#   labs(fill = "Arrangement", x = "Mutation selection coefficient", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_S <- pivot_summary_by_Rep_S %>%  ggplot(aes(x = S,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different mutation selection coefficients") +
#  labs(fill = "Arrangement", x = "Mutation selection coefficient", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_S <- summary_by_Rep_S %>% ggplot(aes(x = S,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different mutation selection coefficients") +
  labs(x = "Mutation selection coefficient", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)



# mean_decrease_per_mutation_S <- summary_by_Rep_S %>% ggplot(aes(x = S,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different mutation selection coefficients") +
#   labs(x = "Mutation selection coefficient", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_S <- summary_by_Rep_S %>% ggplot(aes(x = S,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different mutation selection coefficients") +
  labs(x = "Mutation selection coefficient", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_S,proportion_mutation_S,mean_decrease_S,extinction_time_S,nrow = 4)



## 6) Recombination rate

R = c("1.0e-08","1.0e-07","1.0e-06","1.0e-05")

inversions_R <- data.frame(tick = integer(),
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
                           R = character(),
                           number_mutations_no_NANs_I = numeric(),
                           growth_rate_I = numeric(),
                           growth_rate_N = numeric(),
                           growth_rate_per_mutation_I = numeric(),
                           growth_rate_per_mutation_N = numeric(),
                           decrease_rate = numeric(),
                           decrease_rate_per_mutation = numeric(),
                           stringsAsFactors = FALSE) 



for (i in (1:length(R))) {
  for (Rep in (1:n_rep)) {
    File_R = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_R_1e-5_m_0005_scaling_5 R ",R[i]," S ","-0.001"," I ","0.05"," H ","0.0", " S_INV ", "0.05" ," Mu ","1.0e-05"," Rep ",Rep," RecombinationRateSensitivityAnalysis"," 2"," ReverseMutations F"," .csv", sep ="")
    
    aux_inversions_R <- read.csv(File_R, stringsAsFactors = FALSE)
    
    n_I_R = length(aux_inversions_R$number_mutations_I)
    n_N_R = length(aux_inversions_R$number_mutations_N)
    n_inv_R = length(aux_inversions_R$freq_I)
    
    aux_inversions_R <- aux_inversions_R %>%
      mutate(Rep = Rep,
             R = as.numeric(R[i]),
             number_mutations_no_NANs_I = as.numeric(replace_NAN(number_mutations_I)),
             growth_rate_I = c((number_mutations_no_NANs_I[2:n_I_R] - number_mutations_no_NANs_I[1:(n_I_R-1)])/h,NA),
             growth_rate_N = c((number_mutations_N[2:n_N_R] - number_mutations_N[1:(n_N_R-1)])/h,NA),
             growth_rate_per_mutation_I = growth_rate_I/number_mutations_no_NANs_I,
             growth_rate_per_mutation_N = growth_rate_N/number_mutations_N,
             decrease_rate = c(abs(freq_I[2:n_inv_R] - freq_I[1:(n_inv_R-1)])/h,NA),
             decrease_rate_per_mutation = decrease_rate/freq_I)
    
    inversions_R <- inversions_R %>%
      rbind(aux_inversions_R)
    
  }}


n_R = 11  ## first n after inversion

summary_by_Rep_R <- inversions_R %>% 
  group_by(R,Rep) %>% 
  summarize(mean_growth_I_by_Rep = mean(growth_rate_I[n_R:(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_N_by_Rep = mean(growth_rate_N[n_R:(length(growth_rate_N)-1)], na.rm = T),
            mean_growth_rates_per_mutation_I_by_Rep = mean(growth_rate_per_mutation_I[(n_R+1):(length(growth_rate_I)-1)], na.rm = T),
            mean_growth_rates_per_mutation_N_by_Rep = mean(growth_rate_per_mutation_N[(n_R+1):(length(growth_rate_N)-1)], na.rm = T),
            overall_growth_I_by_Rep = (number_mutations_no_NANs_I[length(number_mutations_no_NANs_I)-1] - number_mutations_no_NANs_I[n_R])/(h*((length(number_mutations_no_NANs_I)-1)-n_R)),
            overall_growth_N_by_Rep = (number_mutations_N[length(number_mutations_N)-1] - number_mutations_N[n_R])/(h*(length(number_mutations_N)-1-n_R)),
            mean_decrease_rate_by_Rep = mean(decrease_rate[(n_R+1):(length(decrease_rate)-1)],na.rm = T),
            mean_decrease_rate_per_mutation_by_Rep = mean(decrease_rate_per_mutation[(n_R+1):(length(decrease_rate)-1)],na.rm = T),
            extinction_time = last(extinction_time),
            extincted = !is.na(extinction_time),
            proportion_mutation_in_haplosome_I = mean(as.numeric(number_mutations_I[n_R:(length(number_mutations_I))])/inv_length, na.rm = T),
            proportion_mutation_in_haplosome_N = mean(as.numeric(number_mutations_N[n_R:(length(number_mutations_N))])/inv_length, na.rm = T)) %>% 
  arrange(R) %>% 
  mutate(R = as.factor(R))


summary_by_Rep_I_R <- summary_by_Rep_R[c("R","Rep","mean_growth_I_by_Rep","mean_growth_rates_per_mutation_I_by_Rep","overall_growth_I_by_Rep", "extincted","proportion_mutation_in_haplosome_I")]
summary_by_Rep_N_R <- summary_by_Rep_R[c("R","Rep","mean_growth_N_by_Rep","mean_growth_rates_per_mutation_N_by_Rep","overall_growth_N_by_Rep","extincted","proportion_mutation_in_haplosome_N")]

summary_by_Rep_I_R  <- summary_by_Rep_I_R %>% 
  rename(mean_growth_by_Rep = mean_growth_I_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_I_by_Rep, 
         overall_growth_by_Rep = overall_growth_I_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_I) %>% 
  mutate(arrangement = "I")

summary_by_Rep_N_R <- summary_by_Rep_N_R %>% 
  rename(mean_growth_by_Rep = mean_growth_N_by_Rep, 
         mean_growth_rates_per_mutation_by_Rep = mean_growth_rates_per_mutation_N_by_Rep, 
         overall_growth_by_Rep = overall_growth_N_by_Rep,
         proportion_mutation_in_haplosome = proportion_mutation_in_haplosome_N) %>% 
  mutate(arrangement = "N")


pivot_summary_by_Rep_R = rbind(summary_by_Rep_I_R,summary_by_Rep_N_R) 
pivot_summary_by_Rep_R <- pivot_summary_by_Rep_R %>% 
  mutate(R = as.factor(R))

summary_R <- summary_by_Rep_R %>% 
  group_by(R) %>% 
  summarize(mean_growth_I = mean(mean_growth_I_by_Rep , na.rm = T),
            mean_growth_N = mean(mean_growth_N_by_Rep , na.rm = T),
            mean_growth_rates_per_mutation_I = mean(mean_growth_rates_per_mutation_I_by_Rep, na.rm = T),
            mean_growth_rates_per_mutation_N = mean(mean_growth_rates_per_mutation_N_by_Rep, na.rm = T),
            overall_growth_I = mean(overall_growth_I_by_Rep,na.rm = T),
            overall_growth_N = mean(overall_growth_N_by_Rep,na.rm = T),
            decrease_rate = mean(mean_decrease_rate_by_Rep,na.rm = T),
            decrease_rate_per_mutation = mean(mean_decrease_rate_per_mutation_by_Rep,na.rm = T),
            extinction_rate = sum(extincted)/length(extincted)) %>% 
  arrange(R)


## Plotting summary statistics

mean_growth_R <- pivot_summary_by_Rep_R %>%  ggplot(aes(x = R,y=mean_growth_by_Rep, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean of mutation growth rates for different recombination rates") +
  labs(fill = "Arrangement", x = "Recombination rate", y = "Mean mutation growth rates") +
  theme_cowplot(12)

proportion_mutation_R <- pivot_summary_by_Rep_R %>%  ggplot(aes(x = R,y=proportion_mutation_in_haplosome, fill = arrangement)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Proportion of mutations in haplosomes for different recombination rates") +
  labs(fill = "Arrangement", x = "Recombination rate", y = "Proportion of mutations in haplosomes") +
  theme_cowplot(12)


# overall_growth_R <- pivot_summary_by_Rep_R %>%  ggplot(aes(x = R,y=overall_growth_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Overall growth rate for different recombination rates") +
#   labs(fill = "Arrangement", x = "Recombination rate", y = "Overall mutation growth rates") +
#   theme_cowplot(12)



# mean_growth_per_mutation_R <- pivot_summary_by_Rep_R %>%  ggplot(aes(x = R,y=mean_growth_rates_per_mutation_by_Rep, fill = arrangement)) +
#   geom_boxplot() +
#   geom_point(position=position_dodge(width=0.75),aes(group=arrangement, col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean of mutation growth rates per mutation for different recombination rates") +
#  labs(fill = "Arrangement", x = "Recombination rate", y = "Mean mutation growth rates per mutation") +
#   theme_cowplot(12)


mean_decrease_R <- summary_by_Rep_R %>% ggplot(aes(x = R,y=mean_decrease_rate_by_Rep)) +
  geom_boxplot() +
  geom_point(aes(col = extincted)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Mean inversion frequency decrease rates for different recombination rates") +
  labs(x = "Recombination rate", y = "Mean inversion frequency decrease rates") +
  theme_cowplot(12)



# mean_decrease_per_mutation_R <- summary_by_Rep_R %>% ggplot(aes(x = R,y=mean_decrease_rate_per_mutation_by_Rep)) +
#   geom_boxplot() +
#   geom_point(aes(col = extincted)) +
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   ggtitle("Mean inversion frequency decrease rates per inversion frequency for different recombination rates") +
#   labs(x = "Recombination rate", y = "Mean inversion frequency decrease rates per inversion frequency") +
#   theme_cowplot(12)


extinction_time_R <- summary_by_Rep_R %>% ggplot(aes(x = R,y=extinction_time)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  ggtitle("Mean extinction time for different recombination rates") +
  labs(x = "Recombination rate", y = "Extinction time") +
  theme_cowplot(12)


grid.arrange(mean_growth_R,proportion_mutation_R,mean_decrease_R,extinction_time_R,nrow = 4)





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
# File = paste("/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H ","0.0"," S_INV ","0.05", " Mu ","5.0e-06 ","Rep ","4"," MutationRateSensitivityAnalysis 1 .csv", sep ="")
File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_5e3_N2_5e3_chrom_length_4e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_1 H 0.0 S_INV 0.05 Mu 5.0e-06 Rep 2 Adele_Comparison 1 ReverseMutations T .csv"
# File = "/home/taubier/Documents/SAUVEGARDES_STAGIAIRES/Kilian/Results/Result_files/freq_ output_test_N1_20e3_N2_20e3_chrom_length_3e3_inv_length_2e3_recrate_1e-5_S_0001_m_0005_freqinv_005_scaling_5 H 0.0 S_INV 0.05 Mu 1.0e-05 Rep 5 MutationDominanceSensitivityAnalysis 2 .csv"
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

three_cols = c("black","violet","blue")
two_cols = c("black","blue")
four_cols = c("grey","pink","green","red")


## Find a way to center the titles

inversions %>%  ggplot(aes(x = tick,y=freq_I)) +
  geom_line(col = "red") +
  ggtitle("Inversion frequency evolution") +
  labs(x = "Generations", y = "Frequency of the inversion") +
  theme_cowplot(12) +
  xlim(0,Tmax*5) +
  ylim(0,1.3)

inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=freq_II, col = "black")) +
  geom_line(aes(y=freq_IN, col = "violet")) +
  geom_line(aes(y=freq_NN, col = "blue")) +
  ggtitle("Evolution of the arrangement haplotypes") +
  labs(x = "Generations", y = "Arrangement haplotype proportions") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,1.05) + 
  scale_colour_manual(name = "Haplotype arrangement", values = three_cols, labels = c("II","IN","NN")) 


inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=as.numeric(mean_fitness_mutation_II), col = "black")) +
  geom_line(aes(y=as.numeric(mean_fitness_mutation_IN), col = "violet")) +
  geom_line(aes(y=as.numeric(mean_fitness_mutation_NN), col = "blue")) +
  ggtitle("Evolution of the mean fitness caused by mutations in arrangements") +
  labs(x = "Generations", y = "Mean fitness caused by mutation") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,1.05) + 
  scale_colour_manual(name = "Haplotype arrangement", values = three_cols, labels = c("II","IN","NN")) 


inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(fitness_load_I),col = "black")) +
  geom_line(aes(y = as.numeric(fitness_load_N),col = "blue")) +
  ggtitle("Evolution of the fitness load of haplosomes") +
  labs(x = "Generations", y = "Fitness load") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 


inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(number_mutations_I),col = "black")) +
  geom_line(aes(y = as.numeric(number_mutations_N),col = "blue")) +
  ggtitle("Evolution of the number of mutations in the inversion area") +
  labs(x = "Generations", y = "Number of mutations in the inversion area") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,max*1.2) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = c(as.numeric(growth_rate_I),NA),col = "black")) +
  geom_line(aes(y = c(as.numeric(growth_rate_N),NA),col = "blue")) +
  ggtitle("Evolution of the mutation growth rates") +
  labs(x = "Generations", y = "Mutation growth rates") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(-0.05,0.25) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(marginal_fitness_inversion_I),col = "black")) +
  geom_line(aes(y = as.numeric(marginal_fitness_inversion_N),col = "blue")) +
  ggtitle("Evolution of the marginal fitness of haplosomes") +
  labs(x = "Generations", y = "Marginal fitness of haplosomes") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(1,max_fitness+0.005) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 


inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=as.numeric(freq_homozygous_mut_II), col = "black")) +
  geom_line(aes(y=as.numeric(freq_homozygous_mut_IN), col = "violet")) +
  geom_line(aes(y=as.numeric(freq_homozygous_mut_NN), col = "blue")) +
  ggtitle("Evolution of the homozygote frequency of mutations") +
  labs(x = "Generations", y = "Homozygote frequency of mutations") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,1) +
  scale_colour_manual(name = "Haplotype arrangement", values = three_cols, labels = c("II","IN","NN")) 


inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=as.numeric(mean_fitness_global_II), col = "black")) +
  geom_line(aes(y=as.numeric(mean_fitness_global_IN), col = "violet")) +
  geom_line(aes(y=as.numeric(mean_fitness_global_NN), col = "blue")) +
  ggtitle("Evolution of the global fitness of arrangements") +
  labs(x = "Generations", y = "Global fitness") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0.1,1.5) +
  scale_colour_manual(name = "Haplotype arrangement", values = three_cols, labels = c("II","IN","NN")) 

inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=as.numeric(mean_fitness_global_II)/as.numeric(mean_fitness_global_NN), col = "black")) +
  geom_line(aes(y=as.numeric(mean_fitness_global_IN)/as.numeric(mean_fitness_global_NN), col = "violet")) +
  geom_line(aes(y=as.numeric(mean_fitness_global_NN)/as.numeric(mean_fitness_global_NN), col = "blue")) +
  ggtitle("Evolution of the normalized global fitness of arrangements") +
  labs(x = "Generations", y = "Normalized lobal fitness") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0.1,1.5) +
  scale_colour_manual(name = "Haplotype arrangement", values = three_cols, labels = c("II","IN","NN")) 

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(marginal_fitness_inversion_I)*(1-as.numeric(fitness_load_I_no_NANs)),col = "black")) +
  geom_line(aes(y = as.numeric(marginal_fitness_inversion_N)*(1-as.numeric(fitness_load_N_no_NANs)),col = "blue")) +
  ggtitle("Evolution of the product of marginal fitnesses") +
  labs(x = "Generations", y = "Product of marginal fitnesses") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0.1,1.5) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(covariance_mutation_inversion_I),col = "black")) +
  geom_line(aes(y = as.numeric(covariance_mutation_inversion_N),col = "blue")) +
  ggtitle("Evolution of the covariance mutation-inversion of haplosomes") +
  labs(x = "Generations", y = "Covariance mutation-inversion") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(-0.02,0.02) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(mean_fitness_global_II)*freq_I + as.numeric(mean_fitness_global_IN)*(1-freq_I),col = "black")) +
  geom_line(aes(y = as.numeric(mean_fitness_global_IN)*freq_I + as.numeric(mean_fitness_global_NN)*(1-freq_I),col = "blue")) +
  ggtitle("Evolution of the marginal fitness") +
  labs(x = "Generations", y = "Marginal fitness") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0.1,1.5) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

# plot(x = inversions$tick, y = inversions$covariance_out_of_marginal_fitness_I, col = "black", type = "l", xlab = "Generations", ylab = "Proportion of fitness explained by covariance", main = "Evolution of the proportion of fitness explained by covariance",xlim = c(0,Tmax),ylim = c(-0.02,0.02))
# lines(x = inversions$tick, y = inversions$covariance_out_of_marginal_fitness_N, col = "blue", type = "l")
# legend("topleft", legend=c("Non-inverted haplosome", "Inverted haplosome"),
#        col=c("blue", "black"), lty=1:2, cex=0.8)

inversions %>% ggplot(aes(x = tick)) +
  geom_line(aes(y = as.numeric(number_mutations_I)/inv_length,col = "black")) +
  geom_line(aes(y = as.numeric(number_mutations_N)/inv_length,col = "blue")) +
  ggtitle("Evolution of the proportion of mutations of haplosomes") +
  labs(x = "Generations", y = "Mutation proportion") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,1) +
  scale_colour_manual(name = "Arrangement", values = two_cols, labels = c("Inverted haplosome","Non-inverted haplosome")) 

non_NA_min <- function(vec) {
  if (length(vec) == 0) {
    return(NA)
  } else {
    return(min(vec,na.rm = T))
  }
}


inversions %>%  ggplot(aes(x = tick,y=c(sapply(c(1 - as.numeric(mean_fitness_global_II)/as.numeric(mean_fitness_global_IN), 1 - as.numeric(mean_fitness_global_NN)/as.numeric(mean_fitness_global_IN)), non_NA_min))),NA) +
  geom_line(col = "cyan") +
  ggtitle("Evolution of the associative overdominance") +
  labs(x = "Generations", y = "AOD (associative overdominance)") +
  theme_cowplot(12) +
  xlim(0,Tmax)

inversions %>%  ggplot(aes(x = tick,y=abs(as.numeric(mean_fitness_global_II)-as.numeric(mean_fitness_global_NN))/as.numeric(mean_fitness_global_IN))) +
  geom_line(col = "cornflowerblue") +
  ggtitle("Evolution of the assymetry between homokaryotypic fitnesses") +
  labs(x = "Generations", y = "Assymetry between homokaryotypic fitnesses") +
  theme_cowplot(12) +
  xlim(0,Tmax)

inversions %>%  ggplot(aes(x = tick,y=as.numeric(effective_selection))) +
  geom_line(col = "orange") +
  ggtitle("Evolution of the effective selection") +
  labs(x = "Generations", y = "Effective selection") +
  theme_cowplot(12) +
  xlim(0,Tmax)

inversions %>%  ggplot(aes(x = tick,y=as.numeric(effective_selection_formula))) +
  geom_line(col = "red") +
  ggtitle("Evolution of the effective selection (obtained with formula)") +
  labs(x = "Generations", y = "Effective selection with formula") +
  theme_cowplot(12) +
  xlim(0,Tmax)

inversions %>%  ggplot(aes(x = tick,y=as.numeric(gamma))) +
  geom_line(col = "brown") +
  ggtitle("Evolution of gamma (overdominance)") +
  labs(x = "Generations", y = "Gamma (overdominance)") +
  theme_cowplot(12) +
  xlim(0,Tmax)

s <- as.numeric(inversions$effective_selection)
m <- 0.005
gamma <- as.numeric(inversions$gamma)

## Equilibria ( when s > 0 )

eq_freq_0 = rep(0,115)
eq_freq_1 = (-2*s-3*gamma+m*gamma+sqrt(gamma^2+m^2*gamma^2+2*m*(2*s^2+gamma*(4+gamma)+s*(2+4*gamma))))/(2*(m-1)*(s + 2*gamma))
eq_freq_2 = -(2*s+3*gamma-m*gamma+sqrt(gamma^2+m^2*gamma^2+2*m*(2*s^2+gamma*(4+gamma)+s*(2+4*gamma))))/(2*(m-1)*(s + 2*gamma))


inversions %>%  ggplot(aes(x = tick)) +
  geom_line(aes(y=freq_I, col = "red")) +
  geom_line(aes(y=eq_freq_0, col = "grey")) +
  geom_line(aes(y=eq_freq_1, col = "pink")) +
  geom_line(aes(y=eq_freq_2, col = "green")) +
  ggtitle("Inversion frequencies evolution (simulation and analytical)") +
  labs(x = "Generations", y = "Frequencies of the inversion") +
  theme_cowplot(12) +
  xlim(0,Tmax) +
  ylim(0,4) +
  # scale_colour_manual(name = "Frequencies", values = c("green","grey","pink","red"), labels = c("Equilibrium n°2","Equilibrium n°0","Equilibrium n°1", "Simulated")) 
  scale_colour_manual(name = "Frequencies", values = four_cols, labels = c("Equilibrium n°0","Equilibrium n°1","Equilibrium n°2", "Simulated")) 



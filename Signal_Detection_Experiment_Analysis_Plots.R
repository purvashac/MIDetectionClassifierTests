names = c("pvalues_super_iter",
          "LogLambda_super_iter",
          "Lambda.hat_iter",
          "score_iter",
          "pvalues_super_boot_iter",
          "pvalues_super_boot_cutoff_iter",
          "pvalues_super_score_boot_iter",
          "pvalues_super_score_boot_cutoff_iter",
          "pvalues_super_permute_iter",
          "pvalues_super_permute_cutoff_iter",
          "pvalues_super_score_permute_iter",
          "pvalues_super_score_permute_cutoff_iter",
          "mean_back_iter",
          "sd_back_iter",
          "pvalues_AUC_iter",
          "pvalues_MC_iter",
          "pvalues_LRT_iter",
          "MisClass_iter",
          "AUCvalues_iter",
          "sd.NWald_iter",
          "LogLambda_iter",
          "pvalues_LRT_boot_iter",
          "pvalues_AUC_boot_iter",
          "pvalues_MC_boot_iter",
          "pvalues_LRT_permute_iter",
          "pvalues_AUC_permute_iter",
          "pvalues_MC_permute_iter",
          "pvalues_LRT_boot_cutoff_iter",
          "pvalues_AUC_boot_cutoff_iter",
          "pvalues_MC_boot_cutoff_iter",
          "pvalues_LRT_permute_cutoff_iter",
          "pvalues_AUC_permute_cutoff_iter",
          "pvalues_MC_permute_cutoff_iter",
          "MisClass2_iter",
          "AUCvalues2_iter",
          "LogLambda2_iter",
          "pvalues_LRT_permute2_iter",
          "pvalues_AUC_permute2_iter",
          "pvalues_MC_permute2_iter",
          "pvalues_LRT_permute2_cutoff_iter",
          "pvalues_AUC_permute2_cutoff_iter",
          "pvalues_MC_permute2_cutoff_iter")
lambda = c(0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.01, 0)
result = vector("list", length(lambda))
super_computed = matrix(NBoot, nrow = B, ncol = length(lambda))
semisuper_computed = matrix(NBoot, nrow = B, ncol = length(lambda))
semisuperin_computed = matrix(NBoot, nrow = B, ncol = length(lambda))

for (iter in 1:length(lambda)) {
  result[[iter]] = as.data.frame(matrix(0, nrow = B, ncol = length(names)))
  colnames(result[[iter]]) = names
  
  for(b in 1:50){
    
    # Supervised Tests
    
    Supervised_file = paste0("Iter", iter, "_Supervised_",b, "_Iter.Rdata")
    
    if (file.exists(Supervised_file)){
      load(Supervised_file)
      result[[iter]][b, 1:12] = c(pvalues_super_iter,
                                  LogLambda_super_iter,
                                  Lambda.hat_iter,
                                  score_iter,
                                  pvalues_super_boot_iter,
                                  pvalues_super_boot_cutoff_iter,
                                  pvalues_super_score_boot_iter,
                                  pvalues_super_score_boot_cutoff_iter,
                                  pvalues_super_permute_iter,
                                  pvalues_super_permute_cutoff_iter,
                                  pvalues_super_score_permute_iter,
                                  pvalues_super_score_permute_cutoff_iter)
      super_computed[b, iter] = 1
    }else{
      result[[iter]][b, 1:12] = rep(NA, 12)
      super_computed[b, iter] = 0
    }
    
    # Semi-Supervised Tests
    
    SemiSupervised_file = paste0("Iter", iter, "_SemiSupervised_",b, "_Iter.Rdata")
    
    if (file.exists(SemiSupervised_file)){
      load(SemiSupervised_file)
      result[[iter]][b, 13:33] = c(mean_back_iter,
                                   sd_back_iter,
                                   pvalues_AUC_iter,
                                   pvalues_MC_iter,
                                   pvalues_LRT_iter,
                                   MisClass_iter,
                                   AUCvalues_iter,
                                   sd.NWald_iter,
                                   LogLambda_iter,
                                   pvalues_LRT_boot_iter,
                                   pvalues_AUC_boot_iter,
                                   pvalues_MC_boot_iter,
                                   pvalues_LRT_permute_iter,
                                   pvalues_AUC_permute_iter,
                                   pvalues_MC_permute_iter,
                                   pvalues_LRT_boot_cutoff_iter,
                                   pvalues_AUC_boot_cutoff_iter,
                                   pvalues_MC_boot_cutoff_iter,
                                   pvalues_LRT_permute_cutoff_iter,
                                   pvalues_AUC_permute_cutoff_iter,
                                   pvalues_MC_permute_cutoff_iter)
                                         
      semisuper_computed[b, iter] = 1
    }else{
      result[[iter]][b, 13:33] = rep(NA, 21)
      semisuper_computed[b, iter] = 0
    }
    
    # Semi-Supervised In-Sample Slow Permutation
    
    SemiSupervisedIn_file = paste0("Iter", iter,"_", setNo, 
                                   "_SemiSupervisedInsample_",b, "_Iter.Rdata")
    
    if (file.exists(SemiSupervisedIn_file)){
      load(SemiSupervisedIn_file)
      result[[iter]][b, 34:42] = c(MisClass2_iter,
                                   AUCvalues2_iter,
                                   LogLambda2_iter,
                                   pvalues_LRT_permute2_iter,
                                   pvalues_AUC_permute2_iter,
                                   pvalues_MC_permute2_iter,
                                   pvalues_LRT_permute2_cutoff_iter,
                                   pvalues_AUC_permute2_cutoff_iter,
                                   pvalues_MC_permute2_cutoff_iter)
                                         
      semisuperin_computed[b, iter] = 1
    }else{
      result[[iter]][b, 34:42] = rep(NA, 9)
      semisuperin_computed[b, iter] = 0
    }
  }
}

pcols = which(str_detect(names, "pvalue") == TRUE &
                str_detect(names, "cutoff") == FALSE)
save(pcols, result, super_computed, semisuper_computed, semisuperin_computed,
     file = "Power_Experiment_BigData.Rdata")

reject.table = matrix(0, nrow = length(pcols), ncol = length(lambda))

for(iter in 1:length(lambda)){
  reject.table[,iter] = colSums(result[[iter]][,pcols] < 0.05, na.rm = T)
  reject.table[1:5,iter] = reject.table[1:5,iter]/sum(super_computed[,iter])
  reject.table[6:14,iter] = reject.table[6:14,iter]/sum(semisuper_computed[,iter])
  reject.table[15:17,iter] = reject.table[15:17,iter]/sum(semisuperin_computed[,iter])
}

colPaperOrder = c(1,2,4,3,5,8,9,12,15,6,10,13,16,7,11,14,17)

reject.table = reject.table[colPaperOrder,]

for(i in 1:length(colPaperOrder)){
  print(names[pcols[colPaperOrder[i]]])
  print(paste(c(rbind(reject.table[i,],  rep("&", 8))),  collapse = " "))
}



b_index = sapply(c(1:8), function(i){return(which(semisuperin_computed[,i] > 0))})
save(b_index, file = "Power_b.RData")


# Emperical Distribution of pvalues Plots

library(data.table)
result_final = vector("list", 7)
for(iter in 1:7){
  result_final[[iter]] = result[[iter+1]][b_select[,iter+1],
                                          pcols[colPaperOrder]]
}
plot_pvalues = as.data.frame(rbindlist(result_final, idcol=T))
plot_pvalues = cbind(plot_pvalues[,-1], plot_pvalues[,1])
colnames(plot_pvalues)[ncol(plot_pvalues)] = ".id"



######## pvalue Nice Plots ###############
lambda = c(0.15, 0.1, 0.07, 0.05, 0.03, 0.01, 0)
df = data.frame(lambda = rep(gl(n = length(lambda), k = 50, 
                                labels = lambda), 5),
                Variables = gl(n = 5, k = nrow(plot_pvalues), 
                               labels = c("Asymptotic LRT",
                                          "Bootstrap LRT",
                                          "Permutation LRT",
                                          "Bootstrap Score",
                                          "Permutation Score")),
                Value = c(as.matrix(plot_pvalues[,1:5])))

df$Method = sapply(strsplit(as.character(df$Variables), split = " "),"[" , 1)
df$Method <- factor(df$Method , levels=c("Asymptotic" , "Bootstrap", 
                                         "Permutation"))

df$Statistic = sapply(strsplit(as.character(df$Variables), split = " "),"[" , 2)
df$Statistic <- factor(df$Statistic , levels=c("LRT" , "Score"))

# A data frame with labels for each facet
f_labels <- data.frame(Method = gl(n = 3, k = 2, 
                                   labels = c("Asymptotic" , "Bootstrap", 
                                              "Permutation")),
                       Statistic = rep(c("LRT" , "Score"), 3),
                       label = c("", "N/A", "", "", "", ""),
                       Value = rep(0.5,6))

ann_text <- data.frame(label = "Text",
                       Value = 0.5,
                       Method = factor("Permutation",
                                       levels = c("Asymptotic" , "Bootstrap", 
                                                  "Permutation")),
                       Statistic = factor("Score",
                                          levels = c("LRT" , "Score")))



cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_super = ggplot(df, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Empirical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method, drop = T) 
#geom_text(aes(label = label), data = f_labels)
#geom_text(data = ann_text, aes(label = label))


pdf("Supervised_Higgs_Power_pvalues_BigData.pdf",
    width = 9, height = 5, onefile=F)
plot_super
dev.off()

plot_super2 = ggplot(df, aes(x = Value, lty = Method, color = Method)) +
  stat_ecdf(size = 0.9)+ ylab("Empirical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(lambda),  nrow = 2)
  facet_grid(Statistic ~ lambda) 
#geom_text(aes(label = label), data = f_labels)
#geom_text(data = ann_text, aes(label = label))


pdf("Supervised_Higgs_Power_pvalues_BigDatav2.pdf",
    width = 12, height = 5, onefile=F)
plot_super2
dev.off()



############### Semi-Supervised Plots ##############################

df_semi = data.frame(lambda = rep(gl(n = length(lambda), k = 50, 
                                     labels = lambda), 12),
                     Variables = gl(n = 12, k = nrow(plot_pvalues), 
                                    labels = c("Asymptotic LRT",
                                               "Bootstrap LRT",
                                               "Permutation LRT",
                                               "Slow LRT",
                                               "Asymptotic AUC",
                                               "Bootstrap AUC",
                                               "Permutation AUC",
                                               "Slow AUC",
                                               "Asymptotic MCE",
                                               "Bootstrap MCE",
                                               "Permutation MCE",
                                               "Slow MCE")),
                     Value = c(as.matrix(plot_pvalues[,6:17])))
df_semi$Method = sapply(strsplit(as.character(df_semi$Variables), 
                                 split = " "),"[" , 1)
df_semi$Method[which(df_semi$Method == "Slow")] = "Slow Permutation"
df_semi$Method <- factor(df_semi$Method , levels=c("Asymptotic" , "Bootstrap", 
                                                   "Permutation", "Slow Permutation"))

df_semi$Statistic = sapply(strsplit(as.character(df_semi$Variables), 
                                    split = " "),"[" , 2)
df_semi$Statistic <- factor(df_semi$Statistic , levels=c("LRT" ,"AUC", "MCE"))


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_semisuper = ggplot(df_semi, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Empirical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method)


pdf("SemiSupervised_Higgs_Power_pvalues_BigData.pdf",
    width = 9, height = 6.5, onefile=F)
plot_semisuper
dev.off()

df_semi = df_semi[-which(df_semi$Variables == "Slow LRT"),]

plot_semisuper_noLRT = ggplot(df_semi, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Empirical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method)

pdf("SemiSupervised_Higgs_Power_pvalues_BigData_NoLRT.pdf",
    width = 9, height = 6.5, onefile=F)
plot_semisuper_noLRT
dev.off()

plot_semisuper_noLRT2 = ggplot(df_semi, aes(x = Value, lty = Method, color = Method)) +
  stat_ecdf(size = 0.9)+ ylab("Empirical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ lambda)

pdf("SemiSupervised_Higgs_Power_pvalues_BigData_NoLRTv2.pdf",
    width = 12, height = 6.5, onefile=F)
plot_semisuper_noLRT2
dev.off()



# Mis-specified Power

load("Power_Experiment_BigData_Misspec.Rdata")

reject.table = matrix(0, nrow = length(pcols), ncol = 8)

for(iter in 1:8){
  reject.table[,iter] = colSums(result[[iter]][b_select[,iter],pcols] < 0.05, 
                                na.rm = T)
  reject.table[1:5,iter] = reject.table[1:5,iter]
  reject.table[6:14,iter] = reject.table[6:14,iter]
  reject.table[15:17,iter] = reject.table[15:17,iter]
}

colPaperOrder = c(1,2,4,3,5,8,9,12,15,6,10,13,16,7,11,14,17)


reject.table = reject.table[colPaperOrder,]

for(i in 1:length(colPaperOrder)){
  print(names[pcols[colPaperOrder[i]]])
  print(paste(c(rbind(reject.table[i,],  rep("&", 8))),  collapse = " "))
}

result_final = vector("list", 7)
for(iter in 1:7){
  result_final[[iter]] = result[[iter+1]][b_select[,iter+1],
                                          pcols[colPaperOrder]]
}
plot_pvalues = as.data.frame(rbindlist(result_final, idcol=T))
plot_pvalues = cbind(plot_pvalues[,-1], plot_pvalues[,1])
colnames(plot_pvalues)[ncol(plot_pvalues)] = ".id"



######## pvalue Nice Plots ###############
df = data.frame(lambda = rep(gl(n = length(lambda), k = 50, 
                                labels = lambda), 5),
                Variables = gl(n = 5, k = nrow(plot_pvalues), 
                               labels = c("Asymptotic LRT",
                                          "Bootstrap LRT",
                                          "Permutation LRT",
                                          "Bootstrap Score",
                                          "Permutation Score")),
                Value = c(as.matrix(plot_pvalues[,1:5])))

df$Method = sapply(strsplit(as.character(df$Variables), split = " "),"[" , 1)
df$Method <- factor(df$Method , levels=c("Asymptotic" , "Bootstrap", 
                                         "Permutation"))

df$Statistic = sapply(strsplit(as.character(df$Variables), split = " "),"[" , 2)
df$Statistic <- factor(df$Statistic , levels=c("LRT" , "Score"))

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_super = ggplot(df, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Emperical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method)


pdf("Supervised_Higgs_Power_pvalues_BigData_Misspec.pdf",
    width = 9, height = 5, onefile=F)
plot_super
dev.off()


############### Semi-Supervised Plots ##############################

df_semi = data.frame(lambda = rep(gl(n = length(lambda), k = 50, 
                                     labels = lambda), 12),
                     Variables = gl(n = 12, k = nrow(plot_pvalues), 
                                    labels = c("Asymptotic LRT",
                                               "Bootstrap LRT",
                                               "Permutation LRT",
                                               "Slow LRT",
                                               "Asymptotic AUC",
                                               "Bootstrap AUC",
                                               "Permutation AUC",
                                               "Slow AUC",
                                               "Asymptotic MCE",
                                               "Bootstrap MCE",
                                               "Permutation MCE",
                                               "Slow MCE")),
                     Value = c(as.matrix(plot_pvalues[,6:17])))
df_semi$Method = sapply(strsplit(as.character(df_semi$Variables), 
                                 split = " "),"[" , 1)
df_semi$Method[which(df_semi$Method == "Slow")] = "Slow Permutation"
df_semi$Method <- factor(df_semi$Method , levels=c("Asymptotic" , "Bootstrap", 
                                                   "Permutation", "Slow Permutation"))

df_semi$Statistic = sapply(strsplit(as.character(df_semi$Variables), 
                                    split = " "),"[" , 2)
df_semi$Statistic <- factor(df_semi$Statistic , levels=c("LRT" ,"AUC", "MCE"))


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_semisuper = ggplot(df_semi, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Emperical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method)


pdf("SemiSupervised_Higgs_Power_pvalues_BigData_Misspec.pdf",
    width = 9, height = 6.5, onefile=F)
plot_semisuper
dev.off()

df_semi = df_semi[-which(df_semi$Variables == "Slow LRT"),]

plot_semisuper_noLRT = ggplot(df_semi, aes(x = Value, lty = lambda, color = lambda)) +
  stat_ecdf(size = 0.9)+ ylab("Emperical Distribution") +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  #facet_wrap(vars(Variables),  nrow = 2)
  facet_grid(Statistic ~ Method)

pdf("SemiSupervised_Higgs_Power_pvalues_BigData_Misspec_noLRT.pdf",
    width = 9, height = 6.5, onefile=F)
plot_semisuper_noLRT
dev.off()




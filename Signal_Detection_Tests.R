source("Test_Functions.R")
load("InitializedData.Rdata")
numCores = detectCores()
registerDoParallel(numCores)

for(iter in 1:length(lambda)){
  B = 50
  trans = 0 # Correctly specified signal doesn't transform the tau_pt variable
# trans = 0.7 # For misspecified case, where we transform the tau_pt variable
  
  tic(paste0("Time Taken for lambda = ", lambda[iter]))
  
  
  samp = rbinom(1, n, lambda[iter]) # No. of signal samples in experimental training data
  n.train = nrow(logdata_b) - 2*m 
  if(n.train %% 2 == 0){
    n = n.train/2
    n.train.bg = n # No. of Background Training Samples
    n.train.exp = n # No. of Experimental Training Samples
  }else{
    n = (n.train + 1)/2
    n.train.bg = n - 1 # No. of Background Training Samples
    n.train.exp = n # No. of Experimental Training Samples
  }
  
  result = foreach (b=1:B, .combine=rbind) %dopar% {
    
    
    print(paste0("Lambda no: ", iter, "Iter", b))
    
    # Randomly ordering the data indices to form Background, Signal, and Experimental data sets
    samp_bg_ind = sample(c(1:nrow(logdata_b)), nrow(logdata_b),
                         replace = F)
    data_bg_ind.train = samp_bg_ind[1:n.train.bg]
    data_be_ind.train = samp_bg_ind[(n.train.bg + 1):n.train]
    data_bg_ind.test = samp_bg_ind[(n.train + 1):(n.train + m)]
    data_be_ind.test = samp_bg_ind[(n.train + m + 1):nrow(logdata_b)]
    sig.ind = sample(c(1:nrow(logdata_s)), nrow(logdata_s),
                     replace = F)
    sig.ind.train = sig.ind[1:n.train.bg]
    exp.sig.ind.train = sig.ind[(n.train.bg + 1):((n.train.bg + nrow(logdata_s))/2)]
    exp.sig.ind.test = sig.ind[((n.train.bg + nrow(logdata_s)+2)/2):nrow(logdata_s)]
    back.ind.train = data_bg_ind.train
    
    # Forming the Experimental Training Data's indices
    # Background samples in the Experimental Training Data
    exp.back.ind.train = sample(data_be_ind.train, size = n - samp, 
                                replace = F,
                                prob = logdata_b$Weight[data_be_ind.train]/sum(logdata_b$Weight[data_be_ind.train]))
    # Signal samples in the Experimental Training Data
    exp.sig.ind.train2 = sample(exp.sig.ind.train, size = samp, 
                                replace = F,
                                prob = logdata_s$Weight[exp.sig.ind.train]/sum(logdata_s$Weight[exp.sig.ind.train]))
    
    # Signal Training Data for Model-Dependent Methods
    supervised_signal = logdata_s[sig.ind.train,1:15]
    
    # Transform signal data, except supervised signal data (for misspecified case):
    x = logdata_s[-sig.ind.train,1]
    y = x - (trans*(x - min(x)))
    logdata_s[-sig.ind.train,1] = y
    
    # Background and Experimental Training Data
    background_data_train = logdata_b[back.ind.train,1:15] 
    experimental_data_train = rbind(logdata_b[exp.back.ind.train,1:15],
                                    logdata_s[exp.sig.ind.train2,1:15]) 
    
    
    ############# TEST DATA #############
    
    d = ncol(experimental_data_train)
    samp2 = rbinom(1, m, lambda[iter]) # No. of signal samples in experimental test data
    
    # Background Test Data
    background_data_test = logdata_b[data_bg_ind.test,1:15]
    
    # Experimental Test Data
    back.test.ind = sample(data_be_ind.test, size = m-samp2, 
                           replace = F,
                           prob = logdata_b$Weight[data_be_ind.test]/sum(logdata_b$Weight[data_be_ind.test]))
    exp.test.ind = sample(exp.sig.ind.test, size = samp2, 
                          replace = F,
                          prob = logdata_s$Weight[exp.sig.ind.test]/sum(logdata_s$Weight[exp.sig.ind.test]))
    experimental_data_test = rbind(logdata_b[back.test.ind,1:15],
                                   logdata_s[exp.test.ind, 1:15])
    
    # Saving the indices to create the datasets in the future
    save(sig.ind.train, back.ind.train, 
         exp.back.ind.train, exp.sig.ind.train2,
         data_bg_ind.test, back.test.ind, exp.test.ind,
         file = paste0("Iter", iter, "_Random_samples_",b, ".Rdata"))
    
    
    ############# SUPERVISED TEST ###########
    
    # Asymptotic Test + Computing the test statistic
    Supervised_Test = supervised.test(background_data_train, 
                                      supervised_signal, 
                                      rbind(experimental_data_test,
                                            experimental_data_train))
    pvalues_super_iter = Supervised_Test$pvalue
    LogLambda_super_iter = Supervised_Test$LogLambda
    Lambda.hat_iter = Supervised_Test$lambda.hat
    score_iter = Supervised_Test$score
    prop.train = nrow(background_data_train)/nrow(supervised_signal)
    
    # Supervised Bootstrap Test
    Super_boot = replicate(NBoot, permute.supervised.test(Supervised_Test$randomF,
                                                          background_data_test,
                                                          rbind(experimental_data_test,
                                                                experimental_data_train),
                                                          prop.train, 
                                                          withreplacement = T))
    pvalues_super_boot_iter = length(which(Super_boot[1,] >= Supervised_Test$LogLambda))/NBoot
    pvalues_super_boot_cutoff_iter = quantile(Super_boot[1,], 0.95)
    pvalues_super_score_boot_iter = length(which(Super_boot[2,] >= score_iter))/NBoot
    pvalues_super_score_boot_cutoff_iter = quantile(Super_boot[2,], 0.95)
    
    # Supervised Permutation Test
    Super_permute = replicate(NBoot, permute.supervised.test(Supervised_Test$randomF,
                                                             background_data_test,
                                                             rbind(experimental_data_test,
                                                                   experimental_data_train),
                                                             prop.train))
    pvalues_super_permute_iter = length(which(Super_permute[1,] >= Supervised_Test$LogLambda))/NBoot
    pvalues_super_permute_cutoff_iter = quantile(Super_permute[1,], 0.95)
    pvalues_super_score_permute_iter = length(which(Super_permute[2,] >= score_iter))/NBoot
    pvalues_super_score_permute_cutoff_iter = quantile(Super_permute[2,], 0.95)
    
    
    save(pvalues_super_iter,
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
         pvalues_super_score_permute_cutoff_iter,
         file = paste0("Iter", iter, "_Supervised_",b, "_Iter.Rdata"))
    
    ############## SEMI-SUPERVISED TEST ###########
    
    # Semi-Supervised Asymptotic Tests + Computing Test Statistics
    Semi_Supervised_Test = semisupervised.test(background_data_train, 
                                               experimental_data_train, 
                                               background_data_test,
                                               experimental_data_test,
                                               Statistic = "both")
    prop.train = nrow(background_data_train)/nrow(experimental_data_train)
    pvalues_AUC_iter = Semi_Supervised_Test$pvalue_AUC
    pvalues_LRT_iter = Semi_Supervised_Test$pvalue_LRT
    pvalues_MC_iter = Semi_Supervised_Test$pvalue_MC
    AUCvalues_iter = Semi_Supervised_Test$AUC
    MisClass_iter = Semi_Supervised_Test$MisClass
    sd.NWald_iter = Semi_Supervised_Test$sd_AUC
    LogLambda_iter = Semi_Supervised_Test$LogLambda
    mean_back_iter = Semi_Supervised_Test$mean_back
    sd_back_iter = Semi_Supervised_Test$sd_back
    
    # Semi-Supervised Bootstrap Test
    Semisuper_boot = replicate(NBoot, permute.semisupervised.test(Semi_Supervised_Test$randomF,
                                                                  background_data_test,
                                                                  experimental_data_test,
                                                                  prop.train,
                                                                  withreplacement = T))
    pvalues_LRT_boot_iter = length(which(Semisuper_boot[1,] >= Semi_Supervised_Test$LogLambda))/NBoot
    pvalues_LRT_boot_cutoff_iter = quantile(Semisuper_boot[1,], 0.95)
    pvalues_AUC_boot_iter = length(which(Semisuper_boot[2,] >= Semi_Supervised_Test$AUC))/NBoot
    pvalues_AUC_boot_cutoff_iter = quantile(Semisuper_boot[2,], 0.95)
    pvalues_MC_boot_iter = length(which(Semisuper_boot[3,] <= Semi_Supervised_Test$MisClass))/NBoot
    pvalues_MC_boot_cutoff_iter = quantile(Semisuper_boot[3,], 0.05)
    
    # Semi-Supervised Permutation Test
    Semisuper_permute = replicate(NBoot, permute.semisupervised.test(Semi_Supervised_Test$randomF,
                                                                     background_data_test,
                                                                     experimental_data_test,
                                                                     prop.train))
    pvalues_LRT_permute_iter = length(which(Semisuper_permute[1,] >= Semi_Supervised_Test$LogLambda))/NBoot
    pvalues_AUC_permute_iter = length(which(Semisuper_permute[2,] >= Semi_Supervised_Test$AUC))/NBoot
    pvalues_MC_permute_iter = length(which(Semisuper_permute[3,] <= Semi_Supervised_Test$MisClass))/NBoot
    pvalues_LRT_permute_cutoff_iter = quantile(Semisuper_permute[1,], 0.95)
    pvalues_AUC_permute_cutoff_iter = quantile(Semisuper_permute[2,], 0.95)
    pvalues_MC_permute_cutoff_iter = quantile(Semisuper_permute[3,], 0.05)
    
    save(mean_back_iter,
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
         pvalues_MC_permute_cutoff_iter,
         file = paste0("Iter", iter,"_", "_SemiSupervised_",b, "_Iter.Rdata"))
    
    
    # Semi-Supervised Slow Permutation Test with In-sample Test Statistic
    Semi_Supervised_Test_insample = semisupervised.insample.test(rbind(background_data_train,
                                                                       background_data_test),
                                                                 rbind(experimental_data_train,
                                                                       experimental_data_test),
                                                                 Statistic = "both")
    AUCvalues2_iter = Semi_Supervised_Test_insample$AUC
    MisClass2_iter = Semi_Supervised_Test_insample$MisClass
    LogLambda2_iter = Semi_Supervised_Test_insample$LogLambda
    Semisuper_insample_permute = replicate(NBoot, permute.semisupervised.insample.test(rbind(background_data_train,
                                                                                             background_data_test),
                                                                                       rbind(experimental_data_train,
                                                                                             experimental_data_test),
                                                                                       Statistic = "both"))
    pvalues_LRT_permute2_iter = length(which(Semisuper_insample_permute[1,] >= LogLambda2_iter))/NBoot
    pvalues_AUC_permute2_iter = length(which(Semisuper_insample_permute[2,] >= AUCvalues2_iter))/NBoot
    pvalues_MC_permute2_iter = length(which(Semisuper_insample_permute[3,] <= MisClass2_iter))/NBoot
    pvalues_LRT_permute2_cutoff_iter = quantile(Semisuper_insample_permute[1,], 0.95)
    pvalues_AUC_permute2_cutoff_iter = quantile(Semisuper_insample_permute[2,], 0.95)
    pvalues_MC_permute2_cutoff_iter = quantile(Semisuper_insample_permute[3,], 0.05)
    
    
    
    save(MisClass2_iter,
         AUCvalues2_iter,
         LogLambda2_iter,
         pvalues_LRT_permute2_iter,
         pvalues_AUC_permute2_iter,
         pvalues_MC_permute2_iter,
         pvalues_LRT_permute2_cutoff_iter,
         pvalues_AUC_permute2_cutoff_iter,
         pvalues_MC_permute2_cutoff_iter,
         file = paste0("Iter", iter, "_SemiSupervisedInsample_",b, "_Iter.Rdata"))
    
    
    
    toc() 
  }
    
    save( lambda, n, m, B, NBoot, cols, logind,
          trans, k,
          file = paste0("Supervised_Experiment_PSC2", "Iter", iter, ".Rdata"))



}
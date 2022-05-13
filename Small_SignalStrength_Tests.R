source("Test_Functions.R")
load("GMMApproximation.RData")
numCores = detectCores()
registerDoParallel(numCores)


lambda = c(0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.01, 0.00, 0.005, 0.001)
B = 50
trans = 0
calibration = FALSE
use.ranger = FALSE
m = 2 * 10^4 # Different sample sizes for the different bars
n = 2 * 10^4

for(iter in c(5,7,8, 9, 10)){
  for(b in 1:B){
    
    tic(paste0("Lambda no: ", iter, "Iter", b))
    
    # Generating Training Data
    
    samp = rbinom(1, n, lambda[iter])
    background_data_train = rGMM(n = n, GMM_back) 
    experimental_data_train = rbind(rGMM(n = n - samp, GMM_back),
                                    rGMM(n = samp, GMM_signal)) 
    
    
    
    # Generating Test Data
    
    d = ncol(experimental_data_train)
    samp2 = rbinom(1, m, lambda[iter])
    background_data_test = rGMM(n = m, GMM_back)
    experimental_data_test = rbind(rGMM(n = m - samp2, GMM_back),
                                   rGMM(n = samp2, GMM_signal))
    
    
    ############## SEMI-SUPERVISED TEST ###########
    
    
    Semi_Supervised_Test = semisupervised.test(background_data_train, 
                                               experimental_data_train, 
                                               background_data_test,
                                               experimental_data_test,
                                               Statistic = "both",
                                               calibration = calibration,
                                               use.ranger = use.ranger)
    
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
    
    
    save(mean_back_iter,
         sd_back_iter,
         pvalues_AUC_iter,
         pvalues_MC_iter,
         pvalues_LRT_iter,
         MisClass_iter,
         AUCvalues_iter,
         sd.NWald_iter,
         LogLambda_iter,
         file = paste0("Iter_RF", iter, "_SemiSupervised_",b, ".Rdata"))
    
    #}
    toc()
    
  }
   
}




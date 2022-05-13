source("Test_Functions.R")
load("InitializedData.Rdata")

lambda = c(0, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
skewLevel = 0
B = 1000
BootNoise = rep(10^{-2}, 15)
AddNoise = FALSE
Bootreplace = TRUE
Overlap = FALSE

for(iter in 1:length(lambda)){
  # Background Indices
  set.seed(iter)
  n.train = nrow(logdata_b) - 2*m
  if(n.train %% 2 == 0){
    n = n.train/2
    n.train.bg = n
    n.train.exp = n
  }else{
    n = (n.train + 1)/2
    n.train.bg = n - 1
    n.train.exp = n
  }
  samp = rbinom(1, n, lambda[iter])
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
  exp.back.ind.train = sample(data_be_ind.train, size = n - samp, 
                              replace = F,
                              prob = logdata_b$Weight[data_be_ind.train]/sum(logdata_b$Weight[data_be_ind.train]))
  exp.sig.ind.train2 = sample(exp.sig.ind.train, size = samp, 
                              replace = F,
                              prob = logdata_s$Weight[exp.sig.ind.train]/sum(logdata_s$Weight[exp.sig.ind.train]))
  
  experimental_supervised_signal = logdata_s[sig.ind.train,1:15]
  
  ############# TRANSFORMED DATA ##########
  # Transform signal data, except supervised signal data:
  x = logdata_s[-sig.ind.train,1]
  y = x - (skewLevel*(x - min(x)))
  logdata_s[-sig.ind.train,1] = y
  
  experimental_data_back_train = logdata_b[back.ind.train,1:15] 
  experimental_data_train = rbind(logdata_b[exp.back.ind.train,1:15],
                                  logdata_s[exp.sig.ind.train2,1:15]) 
  
  # TEST DATA
  
  d = ncol(experimental_data_train)
  samp2 = rbinom(1, m, lambda[iter])
  experimental_data_back_test = logdata_b[data_bg_ind.test,1:15]
  
  back.test.ind = sample(data_be_ind.test, size = m-samp2, 
                         replace = F,
                         prob = logdata_b$Weight[data_be_ind.test]/sum(logdata_b$Weight[data_be_ind.test]))
  exp.test.ind = sample(exp.sig.ind.test, size = samp2, 
                        replace = F,
                        prob = logdata_s$Weight[exp.sig.ind.test]/sum(logdata_s$Weight[exp.sig.ind.test]))
  experimental_data_test = rbind(logdata_b[back.test.ind,1:15],
                                 logdata_s[exp.test.ind, 1:15])
  
  # Estimates for original data
  
  result_original = semi.super.lambda.estimates(experimental_data_back_train,
                                                experimental_data_train, 
                                                experimental_data_back_test,
                                                experimental_data_test, 
                                                glm.CI = TRUE)
  
  # Data sets for bootstrapping
  
  background.data = rbind(experimental_data_back_train, 
                          experimental_data_back_test)
  experimental.data = rbind(experimental_data_train,
                            experimental_data_test)
  
  save( lambda, n, m, B, k, NBoot, cols, logind,
        background.data, experimental.data, 
        BootNoise, AddNoise, Bootreplace, Overlap,
        n.train.bg, n.train.exp,
        result_original,
        file = paste0("EstimatingLambdaBigDataBoot_",iter, ".Rdata"))
}


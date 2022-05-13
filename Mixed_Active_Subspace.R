library(sparsepca)
library(np)
library(Hmisc)
library(grf)
source("Test_Functions.R")
source("Active_subspace_functions.R")
source("InitializedData.Rdata")
numCores = detectCores()
registerDoParallel(numCores)

lambda = c(lambda, 0.5)
set.seed(7)  #for lambda = 0.15
iter = 2
h = 2
B = 500
names = c("tau_pt","tau_eta", "tau_phi", "lep_pt", "lep_eta", "lep_phi",
          "met", "met_phi", "met_sumet", "lead_pt", "lead_eta",
          "sublead_pt", "sublead_eta", "sublead_phi",
          "all_pt")
BootNoise = rep(10^{-6}, 15)
AddNoise = FALSE
Bootreplace = TRUE
Overlap = FALSE
colnames(logdata_b)[1:15] = names
colnames(logdata_s)[1:15] = names

lower_quantile = function(x){quantile(x, 0.025)}
upper_quantile = function(x){quantile(x, 0.975)}
condition = "Conditional"
if(AddNoise == TRUE){
  Noise = "withNoise"
}else{
  Noise = ""
}


my_theme = theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


model = "mixed"
pvalue = 1000
seperately = FALSE
y.as.logit = TRUE
train.only = TRUE
test.only = FALSE
AddNoise = FALSE
lower_quantile = function(x){quantile(x, 0.025)}
upper_quantile = function(x){quantile(x, 0.975)}
condition = "Conditional"
if(AddNoise == TRUE){
  Noise = "withNoise"
}else{
  Noise = ""
}


while(pvalue > 0.05){ # Find dataset that finds the signal
  samp = rbinom(1, n, lambda[iter])
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
  
  experimental_data_back_train = logdata_b[back.ind.train,1:15] 
  experimental_data_train = rbind(logdata_b[exp.back.ind.train,1:15],
                                  logdata_s[exp.sig.ind.train2,1:15]) 
  
  ############# TEST DATA #############
  
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
  
  
  # Testing significance of signal
  Semi_Supervised_Test = semisupervised.test(experimental_data_back_train, 
                                             experimental_data_train, 
                                             experimental_data_back_test,
                                             experimental_data_test,
                                             Statistic = "both")
  pvalues_AUC_iter = Semi_Supervised_Test$pvalue_AUC  #0.0004018246
  pvalues_LRT_iter = Semi_Supervised_Test$pvalue_LRT  #0.0001714126
  pvalue = pvalues_AUC_iter
}

######## Fitting the Classifier to the Training Data ##################### 
combdata = rbind(experimental_data_back_train, 
                 experimental_data_train)
nb = nrow(experimental_data_back_train)
ne = nrow(experimental_data_train)
combdata$class = as.factor(c(rep("b",nb), 
                             rep("e",ne)))

randomF = Semi_Supervised_Test$randomF
combdata$membership = randomF$votes[,2]
combdata$predictedclass = randomF$predicted

############ Fitting Classifier on test data ################


combdata2 = rbind(experimental_data_back_test, experimental_data_test)
combdata2$class = as.factor(c(rep("b",nrow(experimental_data_back_test)), 
                              rep("e",nrow(experimental_data_test))))
comp.Prediction = predict(randomF,
                          newdata = combdata2,
                          type="prob")[,2]
combdata2$membership = comp.Prediction
combdata2$class2 = as.factor(c(rep("background",(nrow(experimental_data_back_test)+
                                                   length(back.test.ind))), 
                               rep("signal",length(exp.test.ind))))

####### Active Subspace on Actual Data #################

result_original = Active_Subspace(experimental_data_back_train, 
                                  experimental_data_back_test,
                                  experimental_data_train, 
                                  experimental_data_test, 
                                  names, h,
                                  bw = NULL,
                                  randomF = Semi_Supervised_Test$randomF,
                                  fitting.data = model,
                                  y.as.logit = y.as.logit)

# Data sets for bootstrapping
if(train.only){
  background.data = experimental_data_back_train
  experimental.data = experimental_data_train
}else if(test.only){
  
  background.data = experimental_data_back_test
  experimental.data = experimental_data_test
}else{
  background.data = rbind(experimental_data_back_train, 
                          experimental_data_back_test)
  experimental.data = rbind(experimental_data_train,
                            experimental_data_test)
}

n = nrow(experimental_data_train)
m = nrow(experimental_data_test)
B = 100
Normal_eigenvector1 = c()
Sparse_eigenvector1 = c()
MeanGradient = c()
MeanScaledGradient = c()

for(index in 1:5){
  tic("100 iterations")
  result = foreach (b=1:B) %dopar% {
    
    
    # Bootstrapped Data
    
    if(test.only){
      data = bootstrapped.data(n = m, m = 0, background.data, 
                               experimental.data, BootNoise,
                               AddNoise = AddNoise, 
                               Bootreplace = Bootreplace,
                               train.test.overlap = Overlap,
                               train.only = test.only,
                               seperately = seperately)
    }else{
      Trail.no = 0
      pvalue = 1000
      # Bootstrapped Data
      while(pvalue > 0.05){
        data = bootstrapped.data(n, m, background.data, 
                                 experimental.data, BootNoise,
                                 AddNoise = AddNoise, 
                                 Bootreplace = Bootreplace,
                                 train.test.overlap = Overlap,
                                 train.only = train.only,
                                 seperately = seperately,
                                 n1 = n.train.bg, n2 = n.train.exp)
        
        
        if(train.only){
          Semi_Supervised_Test2 = semisupervised.test(data$back.train, 
                                                     data$exp.train, 
                                                     experimental_data_back_test,
                                                     experimental_data_test,
                                                     Statistic = "both")
        }else{
          Semi_Supervised_Test2 = semisupervised.test(data$back.train, 
                                                     data$exp.train, 
                                                     data$back.test,
                                                     data$exp.test,
                                                     Statistic = "both")
        }
        pvalue = Semi_Supervised_Test2$pvalue_AUC
        Trail.no = Trail.no + 1
      }
    }
    
    
    
    if(train.only){
      active_subspace = Active_Subspace(data$back.train, 
                                        experimental_data_back_test,
                                        data$exp.train, 
                                        experimental_data_test, 
                                        names, h,
                                        bw = NULL,
                                        randomF = Semi_Supervised_Test2$randomF,
                                        fitting.data = model,
                                        no.betas = TRUE,
                                        y.as.logit = y.as.logit)
    }else if(test.only){
      active_subspace = Active_Subspace(experimental_data_back_train,
                                        data$back.train, 
                                        experimental_data_train, 
                                        data$exp.train, 
                                        names, h,
                                        bw = NULL,
                                        randomF = Semi_Supervised_Test$randomF,
                                        fitting.data = model,
                                        no.betas = TRUE,
                                        y.as.logit = y.as.logit)
      
    }else{
      active_subspace = Active_Subspace(data$back.train, data$back.test,
                                        data$exp.train, data$exp.test, 
                                        names, h,
                                        bw = NULL,
                                        randomF = Semi_Supervised_Test2$randomF,
                                        fitting.data = model,
                                        no.betas = TRUE,
                                        y.as.logit = y.as.logit)
    }
    
    active_subspace
    
  }
  if(index == 1){
    save( result,
          result_original,
          names,
          combdata,
          combdata2,
          Semi_Supervised_Test,
          file = paste0("Active_Subspace_BootstrapRFMixed_bw", h,
                        "_Lambda_scaled",iter,"_", index, condition,
                        Noise, ".Rdata"))
  }else{
    save( result,
          #result_original,
          #names,
          #combdata,
          #combdata2,
          #Semi_Supervised_Test,
          file = paste0("Active_Subspace_BootstrapRFMixed_bw", h,
                        "_Lambda_scaled",iter,"_", index, condition,
                        Noise, ".Rdata"))
  }
  
  toc()
  
  for(j in 1:length(result)){
    Normal_eigenvector1 = rbind(Normal_eigenvector1, 
                                result[[j]]$PCA_normal$loadings[,1])
    Sparse_eigenvector1 = rbind(Sparse_eigenvector1,
                                result[[j]]$PCA_sparse$loadings[,1])
    #Betas = rbind(Betas, result[[j]]$G)
    #Scaled_Betas = rbind(Scaled_Betas, 
    #                     result[[j]]$G/(result[[j]]$smooth$gerr+10^{-10}))
    MeanGradient = rbind(MeanGradient, result[[j]]$MeanGradient)
    MeanScaledGradient = rbind(MeanScaledGradient, 
                               result[[j]]$MeanScaledGradient)
  }
}

Normal_eigenvector1 = rbind(Normal_eigenvector1, 
                            result_original$PCA_normal$loadings[,1])
Sparse_eigenvector1 = rbind(Sparse_eigenvector1,
                            result_original$PCA_sparse$loadings[,1])
MeanGradient = rbind(MeanGradient, result_original$MeanGradient)
MeanScaledGradient = rbind(MeanScaledGradient, 
                           result_original$MeanScaledGradient)

for(j in 1:nrow(Normal_eigenvector1)){
  if(Normal_eigenvector1[j,10] <=0){
    Normal_eigenvector1[j,] = - Normal_eigenvector1[j,]
  }
  if(Sparse_eigenvector1[j,10] <=0){
    Sparse_eigenvector1[j,] = - Sparse_eigenvector1[j,]
  }
}

if(result_original$PCA_normal$loadings[10,1] <=0){
  result_original$PCA_normal$loadings[,1] = - result_original$PCA_normal$loadings[,1]
}
if(result_original$PCA_sparse$loadings[10,1] <=0){
  result_original$PCA_sparse$loadings[,1] = - result_original$PCA_sparse$loadings[,1]
}
df = data.frame(Variable = gl(n = length(names), 
                              k = 1, 
                              labels = names),
                Normal_ci_l = apply(Normal_eigenvector1, 2, lower_quantile),
                Normal_ci_u = apply(Normal_eigenvector1, 2, upper_quantile),
                Sparse_ci_l = apply(Sparse_eigenvector1, 2, lower_quantile),
                Sparse_ci_u = apply(Sparse_eigenvector1, 2, upper_quantile),
                MeanGradient_ci_l = apply(MeanGradient, 2, lower_quantile),
                MeanGradient_ci_u = apply(MeanGradient, 2, upper_quantile),
                MeanGradientScaled_ci_l = apply(MeanScaledGradient, 2, lower_quantile),
                MeanGradientScaled_ci_u = apply(MeanScaledGradient, 2, upper_quantile),
                Normal = result_original$PCA_normal$loadings[,1],
                Sparse = result_original$PCA_sparse$loadings[,1],
                MeanGradient = result_original$MeanGradient,
                MeanGradientScaled = result_original$MeanScaledGradient,
                Normal_se = apply(Normal_eigenvector1, 2, sd),
                Sparse_se = apply(Sparse_eigenvector1, 2, sd),
                MeanGradient_se = apply(MeanGradient, 2, sd),
                MeanGradientScaled_se = apply(MeanScaledGradient, 2, sd))

save( Normal_eigenvector1,
      Sparse_eigenvector1,
      MeanGradient,
      MeanScaledGradient,
      df,
      file = paste0("Active_Subspace_Plots_Mixed_bw", h,
                    "_Lambda_scaled", iter, condition,
                    Noise, ".Rdata"))


print("Mixed Experimental and Background Done!")
print(dim(Normal_eigenvector1))

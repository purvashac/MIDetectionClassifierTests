#library(mixtools)
#library(mclust)
library(MASS)
#library(RandPro)
#library(matrixStats)
library(tictoc)
#library(ggplot2)
#library(GGally)
library(randomForest)
library(ranger)
#library(ROCR)
library(data.table)
#library(np)
#library(sparsepca)
#library(devtools)
#library(MTSKNN)
#library(FNN)
library(parallel)
library(foreach)
library(doParallel)


# Defining Supervised Likelihood ratio function

likelihood.supervised.ratio = function(lambda, psi){
  
  return(sum(log(1 - lambda + lambda*psi)))
}

# Defining Semi-Supervised Likelihood ratio function

likelihood.semisupervised.ratio = function(psi){
  
  return(mean(log(psi)))
}


# Grid Search to find MLE

grid.search.MLE <- function(f, a, b, n = 10^3) {
  
  # Grid of values to calculate the MLE.
  x0 = seq(from = a, to = b, length.out = n/2) 
  fx0 = sapply(x0, f)
  return(x0[which.max(fx0)])
}


############# SUPERVISED LRT ###################

supervised.test <- function(back.train, sig.train, exp.test, 
                            calibration = FALSE, use.ranger = FALSE){
  
  # Fitting Classifier
  combdata = rbind(back.train, sig.train)
  combdata$class = as.factor(c(rep("b",nrow(back.train)), 
                               rep("s",nrow(sig.train))))
  
  if(use.ranger){
    randomF = ranger(class~., data = combdata, min.node.size = 100, 
                     write.forest = T,
                     probability = T)
    combdata$membership = randomF$predictions[,2]
  }else{
    randomF = randomForest(class~., data = combdata, nodesize = 100, 
                           proximity = F)
    combdata$membership = predict(randomF,
                                  newdata = combdata,
                                  type="prob")[,2]
  }
  
  
  if(calibration){
    logitModel = glm(class ~ membership, data = combdata, family = "binomial")
    
    combdata$membership = logitModel$fitted.values
  }
  
  
  # Evaluating Random Forest on Experimental Data
  combdata2 = exp.test
  
  if(use.ranger){
    comp.Prediction = predict(randomF,
                              data = combdata2)
    combdata2$membership = comp.Prediction$predictions[,2]
  }else{
    comp.Prediction = predict(randomF,
                              newdata = combdata2,
                              type="prob")[,2]
    combdata2$membership = comp.Prediction
  }
  
  if(calibration){
    combdata2$membership = predict(logitModel, newdata = combdata2, type = "response")
  }
  
  h = combdata2$membership
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi = nrow(back.train)*h/(nrow(sig.train)*(1 - h))
  
  # Lambda MLE
  
  lambda.hat = grid.search.MLE(function(lambda){return(likelihood.supervised.ratio(lambda, psi))}, 
                               a = 0, b = 1, n = 10^4)
  LogLambda = likelihood.supervised.ratio(lambda.hat,psi)
  
  pvalue = 1 - 0.5*pchisq(2*LogLambda, 1) - 0.5*ifelse(2*LogLambda > 0, 1, 0)
  
  if(!calibration){
    logitModel = NULL
    }
  return(list(pvalue = pvalue, LogLambda = LogLambda,
              lambda.hat = lambda.hat, psi = psi,
              combdata = combdata, combdata2 = combdata2,
              randomF = randomF, logitModel = logitModel, score = mean(psi)))
  
}

############# BOOTSTRAP SUPERVISED LRT ###################

boot.supervised.test <- function(n, randomF, back.test, prop.train, 
                                 logitModel = NULL, calibration = FALSE){
  
  exp.test = back.test[sample(1:nrow(back.test), 
                              size = n, replace = T),]
  
  # Evaluating Random Forest on the Data
  combdata2 = exp.test
  comp.Prediction = predict(randomF,
                            data = combdata2)
  combdata2$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata2$membership = predict(logitModel, newdata = combdata2, type = "response")
  }
 
  
  h = combdata2$membership
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi = prop.train*h/(1 - h)
  
  
  lambda.hat = grid.search.MLE(function(lambda){return(likelihood.supervised.ratio(lambda, psi))}, 
                               a = 0, b = 1, n = 10^4)
  
  LogLambda = likelihood.supervised.ratio(lambda.hat, psi)
  
  return(c(LogLambda, mean(psi)))
  
}

############# PERMUTATION SUPERVISED LRT ###################

permute.supervised.test <- function(randomF,  back.test, exp.test, prop.train,
                                    withreplacement = F, 
                                    logitModel = NULL, calibration = FALSE){
  
  permute.exp = rbind(back.test, exp.test)[sample(1:(nrow(back.test)+nrow(exp.test)), 
                                                  size = nrow(exp.test), 
                                                  replace = withreplacement),]
  
  # Evaluating Random Forest on the Data
  combdata2 = permute.exp
  comp.Prediction = predict(randomF,
                            data = combdata2)
  combdata2$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata2$membership = predict(logitModel, newdata = combdata2, type = "response")
  }
  
  h = combdata2$membership
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi =  prop.train*h/(1 - h)
  
  lambda.hat = grid.search.MLE(function(lambda){return(likelihood.supervised.ratio(lambda, psi))}, 
                               a = 0, b = 1, n = 10^4)
  
  LogLambda = likelihood.supervised.ratio(lambda.hat, psi)
  
  return(c(LogLambda, mean(psi)))
  
}



############# SEMI-SUPERVISED LRT AND AUC TEST ###################

semisupervised.test <- function(back.train, exp.train, back.test = NULL, 
                                exp.test, 
                                Statistic = "both", 
                                calibration = FALSE,
                                use.ranger = FALSE){
  # Statistic = "AUC" or "LRT" or "both"
  # Fitting Classifier
  combdata = rbind(back.train, exp.train)
  combdata$class = as.factor(c(rep("b",nrow(back.train)), 
                               rep("e",nrow(exp.train))))
  
  if(use.ranger){
    randomF = ranger(class~., data = combdata, 
                     min.node.size = 100, 
                     write.forest = T,
                     probability = T)
    combdata$membership = randomF$predictions[,2]
  }else{
    randomF = randomForest(class~., data = combdata, 
                           nodesize = 100, 
                           proximity = F)
    combdata$membership = predict(randomF,
                                  newdata = combdata,
                                  type="prob")[,2]
  }
  
  
  
  
  if(calibration){
    logitModel = glm(class ~ membership, data = combdata, family = "binomial")
    combdata$membership = logitModel$fitted.values
  }
  
  
  # Evaluating Random Forest on Test Experimental Data
  combdata_exp = exp.test
  
  if(use.ranger){
    comp.Prediction = predict(randomF,
                              data = combdata_exp)
    combdata_exp$membership = comp.Prediction$predictions[,2]
  }else{
    combdata_exp$membership = predict(randomF,
                                      newdata = combdata_exp,
                                      type="prob")[,2]
  }
  
  if(calibration){
    combdata_exp$membership = predict(logitModel, newdata = combdata_exp, type = "response")
  }
  
  h = combdata_exp$membership
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi = nrow(back.train)*h/(nrow(exp.train)*(1 - h))
  
  # Evaluating Random Forest on Test Background Data
  
  combdata_back = back.test
  
  if(use.ranger){
    comp.Prediction = predict(randomF,
                              data = combdata_back)
    combdata_back$membership = comp.Prediction$predictions[,2]
  }else{
    combdata_back$membership = predict(randomF,
                                       newdata = combdata_back,
                                       type="prob")[,2]
  }
  
  if(calibration){
    combdata_back$membership = predict(logitModel, newdata = combdata_back, type = "response")
  }
  
  
  hb = combdata_back$membership
  hb[which(hb == 1)] = 1 - 10^{-6}
  hb[which(hb == 0)] = 10^{-6}
  psib = nrow(back.train)*hb/(nrow(exp.train)*(1 - hb))
  
  LogLambda = NA
  pvalue_LRT = NA
  theta.hat = NA
  pvalue_AUC = NA
  sd_AUC = NA
  
  # LRT
  if(Statistic == "both" | Statistic == "LRT"){
    
    LogLambda = likelihood.semisupervised.ratio(psi)
    
    pvalue_LRT = 1 - pnorm(sqrt(length(psi))*(LogLambda - mean(log(psib)))/(sqrt(2)*sd(log(psib))))
  }
  if(Statistic == "both" | Statistic == "AUC"){
    
    # AUC Test
    
    nb = nrow(back.test)
    ne = nrow(exp.test)
    
    # Misclassification error test statistic:
    
    predicted.b = combdata_back$membership
    predicted.e = combdata_exp$membership
    
    MisClass = (sum(predicted.b >= 0.5)+ 
                  sum(predicted.e < 0.5))
    
    pvalues_MC = pbinom(MisClass, size = (nb+ne), prob = 0.5)
    
    
    
    # AUC Test Statistic
    
    tic("Calculating AUC and Other Things")
    rank.b = frankv(predicted.b, order = c(-1), ties.method = "average")
    rank.e = frankv(predicted.e, order = c(-1), ties.method = "average")
    rank.bande = frankv(c(predicted.b, predicted.e), 
                        order = c(-1), ties.method = "average")
    vpdot = nb - (rank.bande[(nb+1):(nb+ne)] - rank.e)
    updot = nb - vpdot
    vdotp = rank.bande[1:nb] - rank.b
    udotp = ne - vdotp
    
    maxrank.b = frankv(predicted.b, order = c(-1), ties.method = "max")
    minrank.b = frankv(predicted.b, order = c(-1), ties.method = "min")
    maxrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "max")
    minrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "min")
    pbeqe = exp(log(sum(maxrank.bande[1:nb] - minrank.bande[1:nb] 
                        - maxrank.b + minrank.b)) - log(nb) - log(ne))
    pbneqe = 1 - pbeqe
    toc()
    theta.hat = exp(log(sum(vdotp)) - log(nb) - log(ne))
    
    # Variance Estimate using Newcombe’s Wald Method
    
    capN = mean(c(nb, ne))
    Var.auc.NWald = (theta.hat*(1 - theta.hat))*
      (2*capN - 1 - ((3*capN - 3)/((2 - theta.hat)*(1 + theta.hat))))/
      ((nb - 1)*(ne - 1))
    sd_AUC = sqrt(Var.auc.NWald)
    pvalue_AUC = 1 - pnorm((theta.hat - 0.5)/sqrt(Var.auc.NWald))
  }
  
  if(!calibration){
    logitModel = NULL
  }
  return(list(pvalue_LRT = pvalue_LRT, pvalue_AUC = pvalue_AUC,
              pvalue_MC = pvalues_MC, MisClass = MisClass,
              psi = psi, psib = psib, LogLambda = LogLambda, AUC = theta.hat,
              sd_AUC = sd_AUC, randomF = randomF, logitModel = logitModel, 
              combdata = combdata, combdata_exp = combdata_exp, 
              combdata_back = combdata_back,
              mean_back = mean(log(psib)), sd_back = sd(log(psib))))
  
}



############# BOOTSTRAP SEMI-SUPERVISED LRT ###################

boot.semisupervised.test <- function(randomF, back.test, prop.train, 
                                     logitModel = NULL, calibration = FALSE){
  
  exp.test = back.test[sample(1:nrow(back.test), 
                              size = nrow(back.test), replace = T), ]
  
  # Evaluating Random Forest on Data
  combdata_exp = exp.test
  comp.Prediction = predict(randomF,
                            data = combdata_exp)
  combdata_exp$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata_exp$membership = predict(logitModel, newdata = combdata_exp, type = "response")
  }
  
  
  h = combdata_exp$membership
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi =  prop.train*h/(1 - h)
  LogLambda = likelihood.semisupervised.ratio(psi)
  
  return(LogLambda)
  
}

############# BOOTSTRAP SEMI-SUPERVISED AUC ###################

boot.semisupervised.AUC.test <- function(randomF, back.test, 
                                         logitModel = NULL, calibration = FALSE){
  
  nb = nrow(back.test)
  ne = nrow(back.test)
  exp.test = back.test[sample(1:(nb/2), 
                              size = ne, replace = T), ]
  back.test = back.test[sample(((nb/2) + 1):nb, 
                               size = nb, replace = T), ]
  
  # Evaluating Random Forest on Test Experimental Data
  combdata_exp = exp.test
  comp.Prediction = predict(randomF,
                            data = combdata_exp)
  combdata_exp$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata_exp$membership = predict(logitModel, newdata = combdata_exp, type = "response")
  }
  
  
  # Evaluating Random Forest on Test Background Data
  
  combdata_back = back.test
  comp.Prediction = predict(randomF,
                            data = combdata_back)
  combdata_back$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata_back$membership = predict(logitModel, newdata = combdata_back, type = "response")
  }
  
  
  # AUC Test Statistic
  
  predicted.b = combdata_back$membership
  predicted.e = combdata_exp$membership

  
  MisClass = (sum(predicted.b >= 0.5)+ 
                sum(predicted.e < 0.5))
  
  rank.b = frankv(predicted.b, order = c(-1), ties.method = "average")
  rank.e = frankv(predicted.e, order = c(-1), ties.method = "average")
  rank.bande = frankv(c(predicted.b, predicted.e), 
                      order = c(-1), ties.method = "average")
  vpdot = nb - (rank.bande[(nb+1):(nb+ne)] - rank.e)
  updot = nb - vpdot
  vdotp = rank.bande[1:nb] - rank.b
  udotp = ne - vdotp
  
  maxrank.b = frankv(predicted.b, order = c(-1), ties.method = "max")
  minrank.b = frankv(predicted.b, order = c(-1), ties.method = "min")
  maxrank.bande = frankv(c(predicted.b, predicted.e), 
                         order = c(-1), ties.method = "max")
  minrank.bande = frankv(c(predicted.b, predicted.e), 
                         order = c(-1), ties.method = "min")
  pbeqe = sum(maxrank.bande[1:nb] - minrank.bande[1:nb] 
              - maxrank.b + minrank.b)/(nb*ne)
  pbneqe = 1 - pbeqe
  theta.hat = sum(vdotp)/(nb*ne)
  
  return(theta.hat, MisClass)
  
}



############# PERMUTATION SEMI-SUPERVISED LRT AND AUC ###################

permute.semisupervised.test <- function(randomF, back.test, exp.test, prop.train, 
                                        Statistic = "both", withreplacement = F,
                                        logitModel = NULL, calibration = FALSE){
  nb = nrow(back.test)
  ne = nrow(exp.test)
  permute.exp = rbind(back.test, exp.test)[sample(1:(nb + ne), 
                                                  size = nb + ne, 
                                                  replace = withreplacement),]
  
  # Evaluating Random Forest on Data
  combdata_exp = permute.exp
  comp.Prediction = predict(randomF,
                            data = combdata_exp)
  combdata_exp$membership = comp.Prediction$predictions[,2]
  if(calibration){
    combdata_exp$membership = predict(logitModel, newdata = combdata_exp, type = "response")
  }
  
  comp.Prediction = combdata_exp$membership
  LogLambda = 0
  theta.hat = 0 
  
  if(Statistic == "both" | Statistic == "LRT"){
    
    h = comp.Prediction[(nb+1):(nb+ne)]
    h[which(h == 1)] = 1 - 10^{-6}
    h[which(h == 0)] = 10^{-6}
    psi =  prop.train*h/(1 - h)
    LogLambda = likelihood.semisupervised.ratio(psi)
  }
  if(Statistic == "both" | Statistic == "AUC"){
    
    # AUC Test Statistic
    
    predicted.b = comp.Prediction[1:nb]
    predicted.e = comp.Prediction[(nb + 1):(nb+ne)]
    MisClass = (sum(predicted.b >= 0.5)+ 
                  sum(predicted.e < 0.5))
    rank.b = frankv(predicted.b, order = c(-1), ties.method = "average")
    rank.e = frankv(predicted.e, order = c(-1), ties.method = "average")
    rank.bande = frankv(c(predicted.b, predicted.e), 
                        order = c(-1), ties.method = "average")
    vpdot = nb - (rank.bande[(nb+1):(nb+ne)] - rank.e)
    updot = nb - vpdot
    vdotp = rank.bande[1:nb] - rank.b
    udotp = ne - vdotp
    
    maxrank.b = frankv(predicted.b, order = c(-1), ties.method = "max")
    minrank.b = frankv(predicted.b, order = c(-1), ties.method = "min")
    maxrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "max")
    minrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "min")
    pbeqe = sum(maxrank.bande[1:nb] - minrank.bande[1:nb] 
                - maxrank.b + minrank.b)/(nb*ne)
    pbneqe = 1 - pbeqe
    theta.hat = sum(vdotp)/(nb*ne)
  }
  
  return(c(LogLambda, theta.hat, MisClass))
  
}

############# SEMI-SUPERVISED IN-SAMPLE LRT AND AUC TEST ###################

semisupervised.insample.test <- function(back, exp, Statistic = "both", 
                                         calibration = FALSE){
  # AUC.test = "AUC" or "LRT" or "both"
  # Fitting Classifier
  combdata = rbind(back, exp)
  nb = nrow(back)
  ne = nrow(exp)
  combdata$class = as.factor(c(rep("b",nb), 
                               rep("e",ne)))
  
  randomF = ranger(class~., data = combdata, min.node.size = 100, 
                   write.forest = T,
                   probability = T)
  combdata$membership = randomF$predictions[,2]
  if(calibration){
    logitModel = glm(class ~ membership, data = combdata, family = "binomial")
    combdata$membership = logitModel$fitted.values
  }
  
  
  h = combdata$membership[(nb + 1):(nb + ne)]
  h[which(h == 1)] = 1 - 10^{-6}
  h[which(h == 0)] = 10^{-6}
  psi = nb*h/(ne*(1 - h))
  
  LogLambda = NA
  theta.hat = NA
  
  # LRT
  if(Statistic == "both" | Statistic == "LRT"){
    
    LogLambda = likelihood.semisupervised.ratio(psi)
    
  }
  if(Statistic == "both" | Statistic == "AUC"){
    
    # AUC Test Statistic
    
    predicted.b = combdata$membership[1:nb]
    predicted.e = h
    MisClass = (sum(predicted.b >= 0.5)+ 
                  sum(predicted.e < 0.5))
    rank.b = frankv(predicted.b, order = c(-1), ties.method = "average")
    rank.e = frankv(predicted.e, order = c(-1), ties.method = "average")
    rank.bande = frankv(c(predicted.b, predicted.e), 
                        order = c(-1), ties.method = "average")
    vpdot = nb - (rank.bande[(nb+1):(nb+ne)] - rank.e)
    updot = nb - vpdot
    vdotp = rank.bande[1:nb] - rank.b
    udotp = ne - vdotp
    
    maxrank.b = frankv(predicted.b, order = c(-1), ties.method = "max")
    minrank.b = frankv(predicted.b, order = c(-1), ties.method = "min")
    maxrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "max")
    minrank.bande = frankv(c(predicted.b, predicted.e), 
                           order = c(-1), ties.method = "min")
    pbeqe = sum(maxrank.bande[1:nb] - minrank.bande[1:nb] 
                - maxrank.b + minrank.b)/(nb*ne)
    pbneqe = 1 - pbeqe
    theta.hat = sum(vdotp)/(nb*ne)
  }
  
  if(!calibration){
    logitModel = NULL
  }
  
  return(list(psi = psi, LogLambda = LogLambda, AUC = theta.hat,
              randomF = randomF, logitModel = logitModel, 
              combdata = combdata, MisClass = MisClass))
  
}


############# SLOW PERMUTATION WITH RE-CLASSIFICATION SEMI-SUPERVISED LRT ###################

permute.semisupervised.insample.test <- function( back, exp, 
                                                  Statistic = "both"){
  nb = nrow(back)
  ne = nrow(exp)
  permute = rbind(back, exp)[sample(1:(nb + ne),
                                    size = nb + ne, replace = F),]
  back.permute = permute[1:nb,]
  exp.permute = permute[(nb+1):(nb+ne),]
  test.stat = semisupervised.insample.test(back.permute, exp.permute,
                                           Statistic)
  LogLambda = 0
  theta.hat = 0 
  
  if(Statistic == "both" | Statistic == "LRT"){
    LogLambda = test.stat$LogLambda
  }
  if(Statistic == "both" | Statistic == "AUC"){
    theta.hat = test.stat$AUC
    MisClass = test.stat$MisClass
  }
  
  return(c(LogLambda, theta.hat, MisClass))
  
}


######## Bootstrapped Data Generation for lambda estimation ###########################

bootstrapped.data <- function(n, m, background, experimental,
                              BootNoise = rep(10^{-2}, 15),
                              AddNoise = FALSE,
                              Bootreplace = TRUE,
                              train.test.overlap = TRUE,
                              train.only = FALSE,
                              seperately = FALSE,
                              n1 = n - 1, n2 = n){
  
  d = ncol(background)
  if(seperately){
    back.train.ind = sample(c(1:n1), size = n1, replace = T)
    back.test.ind = sample(c((n1 + 1):nrow(background)), 
                           size = m, replace = T)
    exp.train.ind = sample(c(1:n2), size = n2, replace = T)
    exp.test.ind = sample(c((n2+1):nrow(experimental)), 
                          size = m, replace = T)
  }else{
    if(train.test.overlap){
      rand.back.ind = sample(c(1:nrow(background)), 
                             size = nrow(background), replace = Bootreplace)
      rand.exp.ind = sample(c(1:nrow(experimental)), 
                            size = nrow(experimental), replace = Bootreplace)
      back.train.ind = rand.back.ind[1:n1]
      back.test.ind = rand.back.ind[(n1 + 1):nrow(background)]
      exp.train.ind = rand.exp.ind[1:n2]
      exp.test.ind = rand.exp.ind[(n2+1):nrow(experimental)]
      
    }else{
      rand.back.ind = sample(c(1:nrow(background)), 
                             size = nrow(background), replace = F)
      rand.exp.ind = sample(c(1:nrow(experimental)), 
                            size = nrow(experimental), replace = F)
      back.train.ind = sample(rand.back.ind[1:n1], size = n1, 
                              replace = Bootreplace)
      exp.train.ind = sample(rand.exp.ind[1:n2], size = n2, 
                             replace = Bootreplace)
      if(!train.only){
        back.test.ind = sample(rand.back.ind[(n1 + 1):nrow(background)], size = m, 
                               replace = Bootreplace)
        exp.test.ind = sample(rand.exp.ind[(n2+1):nrow(experimental)], size = m, 
                              replace = Bootreplace)
      }
      
    }
  }
  
  # Background Training Data
  experimental_data_back_train = background[back.train.ind,] 
  if(AddNoise){
    experimental_data_back_train = experimental_data_back_train + 
      mvrnorm(n1, mu = rep(0, d), Sigma = diag(BootNoise))
  }
  experimental_data_back_train = as.data.frame(experimental_data_back_train)
  colnames(experimental_data_back_train) = colnames(background)
  
  # Experimental Training Data
  experimental_data_train = experimental[exp.train.ind,]
  if(AddNoise){
    experimental_data_train = experimental_data_train + 
      mvrnorm(n2, mu = rep(0, d), Sigma = diag(BootNoise))
  }
  experimental_data_train = as.data.frame(experimental_data_train)
  colnames(experimental_data_train) = colnames(background)
  
  if(!train.only){
    # Background Test Data
    experimental_data_back_test = background[back.test.ind,] 
    if(AddNoise){
      experimental_data_back_test = experimental_data_back_test + 
        mvrnorm(m, mu = rep(0, d),Sigma = diag(BootNoise))
    }
    experimental_data_back_test = as.data.frame(experimental_data_back_test)
    colnames(experimental_data_back_test) = colnames(background)
    
    # Experimental Test Data
    experimental_data_test = experimental[exp.test.ind,]
    if(AddNoise){
      experimental_data_test = experimental_data_test + 
        mvrnorm(m, mu = rep(0, d), Sigma = diag(BootNoise))
    }
    experimental_data_test = as.data.frame(experimental_data_test)
    colnames(experimental_data_test) = colnames(background)
  }
  
 
 if(train.only){
   return(list(back.train = experimental_data_back_train,
               exp.train = experimental_data_train))
 }else{
   return(list(back.train = experimental_data_back_train,
               back.test = experimental_data_back_test,
               exp.train = experimental_data_train,
               exp.test = experimental_data_test))
 }
  
  
}

######## SIGNAL STRENGTH ESTIMATION IN MI MODE BY BINNING DATA ###########################

semi.super.lambda.estimates <- function(experimental_data_back_train, 
                                        experimental_data_train, 
                                        experimental_data_back_test,
                                        experimental_data_test, 
                                        glm.CI = TRUE){
  
  Semi_Supervised_Test = semisupervised.test(experimental_data_back_train, 
                                             experimental_data_train, 
                                             experimental_data_back_test,
                                             experimental_data_test,
                                             Statistic = "both")
  
  
  # Semi-Supervised Lambda Method 2
  
  nb = length(Semi_Supervised_Test$psib)
  ne = length(Semi_Supervised_Test$psi)
  rho = rep(0, ne)
  for(j in 1:ne){
    rho[j] = length(which(Semi_Supervised_Test$psib >= Semi_Supervised_Test$psi[j]))
    rho[j] = rho[j]/nb
  }
  
  # Kernel Density Estimation
  Lambda_semisuper2 = semi.super.lambda.kernel.estimation(nb, ne, rho)
  
  # Histogram Estimation (Bin size 5)
  Lambda_semisuper3 = 1 - (sum(rho >= 0.5)/(0.5*length(rho)))
  Lambda_semisuper4 = 1 - (sum(rho >= 0.6)/(0.4*length(rho)))
  Lambda_semisuper5 = 1 - (sum(rho >= 0.7)/(0.3*length(rho)))
  Lambda_semisuper6 = 1 - (sum(rho >= 0.8)/(0.2*length(rho)))
  Lambda_semisuper7 = 1 - (sum(rho >= 0.9)/(0.1*length(rho)))
  
  
  # Histogram Estimation with Linear Regression (Breaks 200, >= 0.8)
  
  Lambda_estimate1 = semi.super.lambda.hist.regression(rho, breaks = 200, 
                                                       cutoff = 0.8,
                                                       include.loglinear = F,
                                                       include.glm.CI = glm.CI)
  
  Lambda_semisuper8 = Lambda_estimate1$linear
  Lambda_semisuper14 = Lambda_estimate1$glm
  if(glm.CI){
    lower1 = Lambda_estimate1$glm.CI[1]
    upper1 = Lambda_estimate1$glm.CI[2]
  }
  
  
  # Histogram Estimation with Linear Regression (Breaks 100, >= 0.8)
  
  Lambda_estimate2 = semi.super.lambda.hist.regression(rho, breaks = 100, 
                                                       cutoff = 0.8,
                                                       include.loglinear = T,
                                                       include.glm.CI = glm.CI)
  
  Lambda_semisuper9 = Lambda_estimate2$linear
  Lambda_semisuper15 = Lambda_estimate2$glm
  if(glm.CI){
    lower2 = Lambda_estimate2$glm.CI[1]
    upper2 = Lambda_estimate2$glm.CI[2]
  }
  Lambda_semisuper12 = Lambda_estimate2$loglinear
  
  # Histogram Estimation with Linear Regression (Breaks 200, >= 0.5)
  
  Lambda_estimate3 = semi.super.lambda.hist.regression(rho, breaks = 200, 
                                                       cutoff = 0.5,
                                                       include.loglinear = F,
                                                       include.glm.CI = glm.CI)
  
  Lambda_semisuper10 = Lambda_estimate3$linear
  Lambda_semisuper16 = Lambda_estimate3$glm
  if(glm.CI){
    lower3 = Lambda_estimate3$glm.CI[1]
    upper3 = Lambda_estimate3$glm.CI[2]
  }
  
  
  
  # Histogram Estimation with Linear Regression (Breaks 100, >= 0.5)
  
  Lambda_estimate4 = semi.super.lambda.hist.regression(rho, breaks = 100, 
                                                       cutoff = 0.5,
                                                       include.loglinear = T,
                                                       include.glm.CI = glm.CI)
  
  Lambda_semisuper11 = Lambda_estimate4$linear
  Lambda_semisuper17 = Lambda_estimate4$glm
  if(glm.CI){
    lower4 = Lambda_estimate4$glm.CI[1]
    upper4 = Lambda_estimate4$glm.CI[2]
  }
  Lambda_semisuper13 = Lambda_estimate4$loglinear
  
  if(glm.CI){
    return(list(Lambdas = c(Lambda_semisuper2,
                            Lambda_semisuper3,
                            Lambda_semisuper4,
                            Lambda_semisuper5,
                            Lambda_semisuper6,
                            Lambda_semisuper7,
                            Lambda_semisuper8,
                            Lambda_semisuper9,
                            Lambda_semisuper10,
                            Lambda_semisuper11,
                            Lambda_semisuper12,
                            Lambda_semisuper13,
                            Lambda_semisuper14,
                            Lambda_semisuper15,
                            Lambda_semisuper16,
                            Lambda_semisuper17,
                            lower1, upper1, 
                            lower2, upper2,
                            lower3, upper3,
                            lower4, upper4),
                Semi_Supervised_Test =  Semi_Supervised_Test,
                rho = rho
    ))
  }else{
    return(list(Lambdas = c(Lambda_semisuper2,
                            Lambda_semisuper3,
                            Lambda_semisuper4,
                            Lambda_semisuper5,
                            Lambda_semisuper6,
                            Lambda_semisuper7,
                            Lambda_semisuper8,
                            Lambda_semisuper9,
                            Lambda_semisuper10,
                            Lambda_semisuper11,
                            Lambda_semisuper12,
                            Lambda_semisuper13,
                            Lambda_semisuper14,
                            Lambda_semisuper15,
                            Lambda_semisuper16,
                            Lambda_semisuper17),
                Semi_Supervised_Test =  Semi_Supervised_Test,
                rho = rho
    ))
  }
 
 
  
}

######## FUNCTIONS NOT PRESENTED IN THE PAPER #########

############# PERMUTE NN 2-SAMPLE TEST ###################

permute.KNN <- function(k, NN_index, nb, ne){
  
  samp = sample(1:(nb + ne), size = nb + ne, replace = F)
  inverse.samp = order(samp)
  NN_index_new = matrix(inverse.samp[NN_index[samp,]], ncol = k)
  
  # Evaluating Test Statistic on Data
  
  TestStat_NN = sum(NN_index_new[1:nb,] <= nb) +
    sum(NN_index_new[(nb + 1):(nb+ne),] > nb)
  
  return(TestStat_NN)
  
}

############# PERMUTE NN 2-SAMPLE TEST ###################

Long.permute.KNNvfast <- function(k,  NN_index, nb, ne){
  
  
  samp = sample(1:(nb + ne), size = nb + ne, replace = F)
  inverse.samp = order(samp)
  NN_index_new = matrix(inverse.samp[NN_index[samp,]], ncol = max(k))
  pvalues_NN_TestStat_permute_iter = rep(1000, length(k))
  
  # Evaluating Test Statistic on Data
  
  for(k0 in 1:length(k)){
    
    pvalues_NN_TestStat_permute_iter[k0] = sum(NN_index_new[1:nb,1:k[k0]] <= nb) +
      sum(NN_index_new[(nb + 1):(nb+ne),1:k[k0]] > nb)
    
  }
  
  return(pvalues_NN_TestStat_permute_iter)
  
}


############# PERMUTE NN 2-SAMPLE TEST ###################

Long.permute.KNN <- function(k, back_data, exp_data){
  
  nb = nrow(back_data)
  ne = nrow(exp_data)
  data = rbind(back_data, exp_data)
  samp = sample(1:(nb + ne), size = nb + ne, replace = F)
  back_data = data[samp[1:nb],]
  exp_data = data[samp[(nb+1):(nb+ne)],]
  Finding_NN = get.knn(rbind(back_data, exp_data), k = max(k),
                       algorithm = "kd_tree")
  pvalues_NN_TestStat_permute_iter = rep(1000, length(k))
  
  # Evaluating Test Statistic on Data
  
  for(k0 in 1:length(k)){
    
    pvalues_NN_TestStat_permute_iter[k0] = sum(Finding_NN$nn.index[1:nb,1:k[k0]] <= nb) +
      sum(Finding_NN$nn.index[(nb + 1):(nb+ne),1:k[k0]] > nb)
    
  }
  
  return(pvalues_NN_TestStat_permute_iter)
  
}

####### Epanechnikov density for estimating Lambda #########
EpanechnikovK  = function(x){
  return(ifelse(x <= 1 & x >= -1, 0.75*(1 - x^2), 0))
}


######## SIGNAL STRENGTH ESTIMATION IN MI MODE USING KERNEL REGRESSION #########


semi.super.lambda.kernel.estimation <- function(nb, ne, rho){
  bw = floor(sqrt(ne))*min(1-rho)
  bw = ifelse(bw == 0, floor(sqrt(ne))*min(1-rho[which(rho!=1)]), bw)
  gqder = mean(EpanechnikovK((1 - rho + bw)/bw)) -
    2*mean(EpanechnikovK((1 - rho)/bw))
  k1 = 3/8
  k2 = 9/15
  bw_opt = (2*k2/(ne*(gqder*k1)^2))^(1/3)
  gq = (2/bw_opt)*mean(EpanechnikovK((1 - rho)/bw_opt))
  Lambda = 1-gq
  return(Lambda)
}

semi.super.lambda.hist.regression <- function(rho, breaks = 200,
                                              cutoff = 0.8,
                                              include.linear = T,
                                              include.glm = T,
                                              include.loglinear = T,
                                              include.glm.CI = T,
                                              CI.aplha = 0.05){
  y = hist(rho, breaks = breaks, plot = F)
  Lambda_predict_df = data.frame(Counts = y$counts[-1],
                                 Density = y$density[-1],
                                 Breaks = y$breaks[-c(1:2)])
  if(include.linear){
    Lambda_predict_model = lm(Density ~ Breaks, 
                              data = subset(Lambda_predict_df, Breaks >= cutoff))
    if(coefficients(Lambda_predict_model)[2]>0){
      Lambda_predict = mean(subset(Lambda_predict_df, Breaks >= cutoff)$Density)
      Lambda_linear = 1 - Lambda_predict[1]
    }else{
      Lambda_predict = predict(Lambda_predict_model, 
                               newdata = data.frame(Breaks = c(1)))
      Lambda_linear = 1 - Lambda_predict[1]
    }
  }else{
    Lambda_linear = NA
  }
  
  if(include.glm){
    Lambda_glm_predict_model = glm(Counts ~ Breaks,
                                   family = poisson(),
                                   data = subset(Lambda_predict_df, Breaks >= cutoff))
    family = family(Lambda_glm_predict_model)
    break.l = Lambda_predict_df$Breaks[2] - Lambda_predict_df$Breaks[1]
    Lambda_predict = predict(Lambda_glm_predict_model, 
                             newdata = data.frame(Breaks = c(1)),
                             se.fit = TRUE)
    if(include.glm.CI){
      lower = family$linkinv(Lambda_predict$fit - 
                               qnorm(1 - CI.aplha/2) * Lambda_predict$se.fit)
      upper = family$linkinv(Lambda_predict$fit + 
                               qnorm(1 - CI.aplha/2) * Lambda_predict$se.fit)
      lower1 = 1 - (upper/(length(rho)*break.l))
      if(lower1 < 0){lower1 = 0}
      if(lower1 > 1){lower1 = 1}
      upper1 = 1 - (lower/(length(rho)*break.l))
      if(upper1 > 1){upper1 = 1}
      if(upper1 < 0){upper1 = 0}
    }
    
    if(coefficients(Lambda_predict_model)[2]>0){
      mean_predict = mean(subset(Lambda_predict_df, Breaks >= cutoff)$Density)
      Lambda_glm = 1 - mean_predict
    }else{
      Lambda_glm = 1 - (family$linkinv(Lambda_predict$fit)/
                                  (length(rho)*break.l))
    }
  }else{
    upper1 = NA
    lower1 = NA
    Lambda_glm = NA
  }
  
  if(include.loglinear){
    logdata = subset(Lambda_predict_df, Breaks >= cutoff)
    index = which(logdata$Density > 0)
    logdata = logdata[index,]
    LogLambda_predict_model = lm(log(Density) ~ Breaks, 
                                 data = logdata)
    if(coefficients(Lambda_predict_model)[2]>0){
      Lambda_predict = exp(mean(log(subset(Lambda_predict_df, 
                                           Breaks >= cutoff)$Density)))
      Lambda_loglinear = 1 - Lambda_predict[1]
    }else{
      Lambda_predict = exp(predict(LogLambda_predict_model, 
                                   newdata = data.frame(Breaks = c(1))))
      Lambda_loglinear = 1 - Lambda_predict[1]
    }
    
  }else{
    Lambda_loglinear = NA
  }
  if(include.glm.CI){
    return(list(linear = Lambda_linear, 
                glm = Lambda_glm, 
                glm.CI = c(lower1, upper1),
                loglinear = Lambda_loglinear))
  }else{
    return(list(linear = Lambda_linear, 
                glm = Lambda_glm,
                loglinear = Lambda_loglinear))
  }
  
}

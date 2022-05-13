library(LaplacesDemon)
Active_Subspace <- function(back.train, back.test,
                            exp.train, exp.test, names, h,
                            bw = NULL,
                            randomF = NULL,
                            fitting.data = "exp",
                            no.betas = FALSE,
                            columns = 1:15,
                            y.as.logit = FALSE){
  
  combineddata = as.data.frame(rbind(back.train, exp.train))
  colnames(combineddata) = names
  combineddata$class = as.factor(c(rep("b", nrow(back.train)), 
                         rep("e", nrow(exp.train))))
  
  # Bootstrapped RF
  if(is.null(randomF)){
    randomF_boot = randomForest(class~., data = combineddata, 
                                ntree = 1000,
                                nodesize = 35, 
                                proximity = F)
  }else{
    randomF_boot = randomF
  }
  
  
  # Bootstrapped test data
  if(fitting.data == "exp"){
    combineddata_test = as.data.frame(exp.test)
  }else if(fitting.data == "back"){
    combineddata_test = as.data.frame(back.test)
  }else{
    combineddata_test = as.data.frame(rbind(back.test, exp.test))
  }
  
  if(is.null(bw)){
    if(length(columns) == 1){
      bw1 = sqrt(var(combineddata_test[,columns]))/h
    }else{
      bw1 = sqrt(diag(var(combineddata_test[,columns])))/h
    }
  }
  
  names(combineddata_test) = names

  combineddata_test$membership = predict(randomF_boot,
                                         newdata = combineddata_test,
                                         type="prob")[,2]
  if(y.as.logit){
    h = combineddata_test$membership
    h.0 = which(combineddata_test$membership <= 10^(-10))
    h.1 = which(combineddata_test$membership >= 1 -  10^(-10))
    h[h.0] = 10^(-10)
    h[h.1] = 1 -  10^(-10)
    smooth = npreg(bws = bw1,
                   txdat = combineddata_test[,columns], 
                   tydat = logit(h),
                   gradients = T)
  }else{
    smooth = npreg(bws = bw1,
                   txdat = combineddata_test[,columns], 
                   tydat = combineddata_test$membership, gradients = T)
  }
  
  G = smooth$grad
  MeanGradient = colMeans(G)
  MeanProjection = as.vector(as.matrix(combineddata_test[,columns]) %*% MeanGradient)
  MeanScaledGradient = colMeans(G/(smooth$gerr+10^{-10}))
  
  ############ Eigenvalues and vectors #################################
  
  #G.pick = scale(G[[i]], center = TRUE, scale = TRUE)
  G.pick = scale(G/(smooth$gerr+10^{-10}), center = TRUE, scale = TRUE)
  
  PCA_normal = spca(G.pick, alpha = 0)
  PCA_sparse = spca(G.pick)
  
  if(no.betas){
    return(list(smooth = NULL, 
                G = NULL,
                PCA_normal = PCA_normal,
                PCA_sparse = PCA_sparse,
                MeanProjection = MeanProjection,
                MeanGradient = MeanGradient,
                MeanScaledGradient = MeanScaledGradient))
  }else{
    return(list(smooth = smooth, 
                G = G,
                PCA_normal = PCA_normal,
                PCA_sparse = PCA_sparse,
                MeanProjection = MeanProjection,
                MeanGradient = MeanGradient,
                MeanScaledGradient = MeanScaledGradient))
  }
  
  
}

Active_Subspace_grf <- function(back.train, back.test,
                                exp.train, exp.test, names, h,
                                bw = NULL,
                                randomF = NULL,
                                fitting.data = "exp",
                                no.betas = FALSE,
                                columns = 1:15){
  
  combineddata = as.data.frame(rbind(back.train, exp.train))
  colnames(combineddata) = names
  combineddata$class = c(rep(0, nrow(back.train)), 
                         rep(1, nrow(exp.train)))

# Bootstrapped RF
randomF_boot = boosted_regression_forest(X = as.matrix(combineddata[,columns]),
                                         Y = as.matrix(combineddata$class))


# Bootstrapped test data
if(fitting.data == "exp"){
  combineddata_test = as.data.frame(exp.test)
}else if(fitting.data == "back"){
  combineddata_test = as.data.frame(back.test)
}else if(fitting.data == "mixed"){
  combineddata_test = as.data.frame(rbind(back.test, exp.test))
}else{
  combineddata_test = as.data.frame(seq(-2,2,length=40))
  colnames(combineddata_test) = names
}

if(is.null(bw)){
  if(length(columns) == 1){
    bw1 = sqrt(var(combineddata_test[,columns]))/h
  }else{
    bw1 = sqrt(diag(var(combineddata_test[,columns])))/h
  }
}

names(combineddata_test) = names

combineddata_test$membership = predict(randomF_boot, 
                                       as.matrix(combineddata_test))$predictions

smooth = npreg(bws = bw1,
               txdat = combineddata_test[,columns], 
               tydat = combineddata_test$membership, gradients = T)
G = smooth$grad
MeanGradient = colMeans(G)
MeanProjection = as.vector(as.matrix(combineddata_test[,columns]) %*% MeanGradient)
MeanScaledGradient = colMeans(G/(smooth$gerr+10^{-10}))

############ Eigenvalues and vectors #################################

#G.pick = scale(G[[i]], center = TRUE, scale = TRUE)
G.pick = scale(G/(smooth$gerr+10^{-10}), center = TRUE, scale = TRUE)

PCA_normal = spca(G.pick, alpha = 0)
PCA_sparse = spca(G.pick)

if(no.betas){
  return(list(smooth = NULL, 
              G = NULL,
              PCA_normal = PCA_normal,
              PCA_sparse = PCA_sparse,
              MeanProjection = MeanProjection,
              MeanGradient = MeanGradient,
              MeanScaledGradient = MeanScaledGradient))
}else{
  return(list(smooth = smooth, 
              G = G,
              PCA_normal = PCA_normal,
              PCA_sparse = PCA_sparse,
              MeanProjection = MeanProjection,
              MeanGradient = MeanGradient,
              MeanScaledGradient = MeanScaledGradient))
}


}



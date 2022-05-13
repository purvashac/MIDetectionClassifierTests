source("Test_Functions.R")
# I ran a total of 1000 bootstrap cycles, running 250 at a time over 4 loops.
# So, there are 8 lambdas times 4 iterations (36 iterations)

for(index in 0:35){
  iter = floor(as.numeric(index)/4)+1
  round = as.numeric(index) %% 4
  set.seed(as.numeric(index))
  
  load(paste0("EstimatingLambdaBigDataBoot_",iter, ".Rdata"))
  B = 250
  
  result = matrix(0, nrow = B, ncol = 16)
  for(b in 1:B){
    data = bootstrapped.data(n, m, background.data, 
                             experimental.data, BootNoise,
                             AddNoise = AddNoise, 
                             Bootreplace = Bootreplace,
                             train.test.overlap = Overlap,
                             n1 = n.train.bg, n2 = n.train.exp) 
    
    lambda.estimates = semi.super.lambda.estimates(data$back.train,
                                                   data$exp.train, 
                                                   data$back.test,
                                                   data$exp.test, 
                                                   glm.CI = FALSE)
    result[b,] = lambda.estimates$Lambdas
  }
  
  save( result, 
        file = paste0("EstimatingLambdaBigDataBoot_",iter, "round", round,".Rdata"))
  
}










  


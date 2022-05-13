new_bootstrap <- function(bootstrap_matrix, Estimate){
  
  # Compute δ∗ for each bootstrap sample
  deltastar = bootstrap_matrix - c(rep(1, nrow(bootstrap_matrix)))%*%t(Estimate)
  
  # Find the 0.1 and 0.9 quantile for deltastar
  d = sapply(1:16,
             function(x){quantile(deltastar[,x],c(0.025, 0.975))})
  d = t(d)
  # Calculate the 80% confidence interval for the mean.
  
  cl = Estimate - d[,2]
  cu = Estimate - d[,1]
  
  cl[which(cl < 0)] = 0
  cu[which(cu < 0)] = 0
  cl[which(cl > 1)] = 1
  cu[which(cu > 1)] = 1
  
  return(list(lower = cl, upper = cu))
}


library(gridExtra)
library(ggplot2)


lambdaFinal = c(0, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)


Bootstderrors = as.data.frame(matrix(0, nrow = 16, 
                                     ncol = length(lambdaFinal)))
colnames(Bootstderrors) = lambdaFinal
Bootstderrors2 = as.data.frame(matrix(0, nrow = 16, 
                                      ncol = length(lambdaFinal)))
colnames(Bootstderrors2) = lambdaFinal
BootCI = as.data.frame(matrix(0, nrow = 16*2, 
                              ncol = length(lambdaFinal)))
colnames(BootCI) = lambdaFinal
BootCI2 = as.data.frame(matrix(0, nrow = 16*2, 
                               ncol = length(lambdaFinal)))
colnames(BootCI2) = lambdaFinal

BootBasicCI = as.data.frame(matrix(0, nrow = 16*2,
                                   ncol = length(lambdaFinal)))
colnames(BootBasicCI) = lambdaFinal

MeanEstimate = as.data.frame(matrix(0, nrow = 16, 
                                    ncol = length(lambdaFinal)))
colnames(MeanEstimate) = lambdaFinal

MeanEstimate2 = as.data.frame(matrix(0, nrow = 16, 
                                     ncol = length(lambdaFinal)))
colnames(MeanEstimate2) = lambdaFinal

PointEstimate = as.data.frame(matrix(0, nrow = 16, 
                                     ncol = length(lambdaFinal)))
colnames(PointEstimate) = lambdaFinal

GLMCI = as.data.frame(matrix(0, nrow = 8, 
                             ncol = length(lambdaFinal)))
colnames(PointEstimate) = lambdaFinal


for(i in 1:length(lambdaFinal)){
  load(paste0("EstimatingLambdaBigDataBoot_",i, ".Rdata"))
  if(abs(lambda[iter[i]]-lambdaFinal[i])>10^(-4)){
    print(paste0("Iteration = ", i, " Normal Bootstrap"))
    print("Error! Lambda's don't match!")
  }
  
  results = matrix(0, nrow = B, ncol = 16)
  for(round in 0:3){
    load(paste0("EstimatingLambdaBigDataBoot_",i, 
                               "round", round,".Rdata"))
    results[(round*250+1):((round+1)*250),] = result
  }
  result = results
  
  MeanEstimate[,i] = colMeans(result[,1:16])
  Bootstderrors[,i] = sqrt(diag(var(result[,1:16])))
  BootCI[,(2*i - 1)] = sapply(1:16, 
                              function(x){max(0,quantile(result[,x],0.025))})
  BootCI[,(2*i)] = sapply(1:16, 
                          function(x){min(1,quantile(result[,x],0.975))})
  
  # Point Estimates and GLM CI
  PointEstimate[,i] = result_original$Lambdas[1:16]
  GLMCI[,i] = result_original$Lambdas[17:24]
  
  # Basic Bootstrap
  Interval = new_bootstrap(result[,1:16], PointEstimate[,i])
  BootBasicCI[,(2*i - 1)] = Interval$lower
  BootBasicCI[,(2*i)] = Interval$upper
  
}

rows = c(1,13:16)
for(i in rows){
  m = rbind(sapply(c(1:length(lambdaFinal)), 
                   function(j){
                     return(paste0(round(PointEstimate[i,j],3)))
                   }), 
            rep("&", ncol(MeanEstimate)))
  print(paste(as.vector(t(t(m))), collapse = " "))
  std1 = rbind(sapply(c(1:length(lambdaFinal)), 
                      function(j){
                        return(paste0(round(Bootstderrors[i,j],3)))
                      }), 
               rep("&", ncol(MeanEstimate)))
  print(paste(as.vector(t(t(std1))), collapse = " "))
  std2 = rbind(sapply(c(1:length(lambdaFinal)), 
                      function(j){
                        return(paste0(round(Bootstderrors2[i,j],3)))
                      }), 
               rep("&", ncol(MeanEstimate)))
  print(paste(as.vector(t(t(std2))), collapse = " "))
  BCI1 = rbind(sapply(c(1:length(lambdaFinal)), 
                      function(j){
                        return(paste0("[", round(BootCI[i,(2*j - 1)],3), ", ", 
                                      round(BootCI[i,(2*j)],3), "]"))
                      }), 
               rep("&", ncol(MeanEstimate)))
  print(paste(as.vector(t(t(BCI1))), collapse = " "))
  BCI2 = rbind(sapply(c(1:length(lambdaFinal)), 
                      function(j){
                        return(paste0("[", round(BootCI2[i,(2*j - 1)],3), ", ", 
                                      round(BootCI2[i,(2*j)],3), "]"))
                      }), 
               rep("&", ncol(MeanEstimate)))
  print(paste(as.vector(t(t(BCI2))), collapse = " "))
  if(i > 12){
    o = rbind(sapply(c(1:length(lambdaFinal)), 
                     function(j){
                       return(paste0("[", round(GLMCI[2*(i-12) - 1,j],3), ", ", 
                                     round(GLMCI[2*(i-12),j],3), "]"))
                     }), 
              rep("&", ncol(MeanEstimate)))
    print(paste(as.vector(t(t(o))), collapse = " "))
  }
}



###### PLOTS ###########
rows = c(13:16, 1)
PointEstimate[PointEstimate < 0] = 0
PointEstimate[PointEstimate > 1] = 1


#################### PLOTS #######################################

LambdaEstimates = data.frame(Lambda = rep(lambdaFinal, 5),
                             Estimates = c(t(PointEstimate[rows,])),
                             BootSE = c(t(Bootstderrors[rows,])),
                             BootCIl = c(t(BootCI[rows,(2*c(1:9)) - 1])),
                             BootCIu = c(t(BootCI[rows,(2*c(1:9))])),
                             BootBasicCIl = c(t(BootBasicCI[rows,(2*c(1:9)) - 1])),
                             BootBasicCIu = c(t(BootBasicCI[rows,(2*c(1:9))])),
                             GLMCIl = c(c(t(GLMCI[(2*c(1:4)) - 1,])), rep(NA, 9)),
                             GLMCIu = c(c(t(GLMCI[(2*c(1:4)),])), rep(NA, 9)),
                             Method = gl(n = 5, k = 9, 
                                         labels = c("T = 0.8, Bins = 200",
                                                    "T = 0.8, Bins = 100",
                                                    "T = 0.5, Bins = 200",
                                                    "T = 0.5, Bins = 100",
                                                    "Kernel")))

LambdaEstimates$SECIl = LambdaEstimates$Estimates - 1.96*LambdaEstimates$BootSE
LambdaEstimates$SECIu = LambdaEstimates$Estimates + 1.96*LambdaEstimates$BootSE
LambdaEstimates$SECIl[LambdaEstimates$SECIl < 0] = 0
LambdaEstimates$SECIl[LambdaEstimates$SECIl > 1] = 1
LambdaEstimates$SECIu[LambdaEstimates$SECIu < 0] = 0
LambdaEstimates$SECIu[LambdaEstimates$SECIu > 1] = 1

df = data.frame(Lambda = rep(lambdaFinal, 5*4),
                Estimates = rep(LambdaEstimates$Estimates, 4),
                Method = rep(LambdaEstimates$Method, 4),
                Boot_l = c(LambdaEstimates$BootBasicCIl,
                           LambdaEstimates$BootCIl,
                           LambdaEstimates$SECIl,
                           LambdaEstimates$GLMCIl),
                Boot_u = c(LambdaEstimates$BootBasicCIu,
                           LambdaEstimates$BootCIu,
                           LambdaEstimates$SECIu,
                           LambdaEstimates$GLMCIu),
                Bootstrap = gl(n = 4, k = 9*5, 
                               labels = c("Basic Bootstrap",
                                          "Bootstrapped Quantiles",
                                          "Bootstrapped SE",
                                          "GLM")))

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.02) # move them .05 to the left and right

# The palette with grey:
cbPalette <- c( "#E69F00","#D55E00", "#0072B2", "#009E73", "#F0E442", "#56B4E9",  "#CC79A7", "#999999")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")

plot_new = ggplot(subset(df, Method!="Kernel"),
                  aes(x=Lambda, y=Estimates,
                      colour=Method, group=Method)) + 
  geom_errorbar(aes(ymin=Boot_l, ymax=Boot_u), 
                width=.02, position=pd) +
  #geom_errorbar(aes(ymin=BootCIl, ymax=BootCIu), width=.1) +
  #geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.1) +
  #geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #ggtitle("CIs Using Bootstrapped SE")+
  #theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1))+
  facet_wrap(vars(Bootstrap),  nrow = 2)


ggsave(plot_new, filename = "EstimatingLambda_BigData_Facets.png",
       width = 7, height = 5)
ggsave(plot_new, filename = "EstimatingLambdaSE_BigData_Facets.pdf",
       width = 7, height = 5)

###################################### OLD PLOTS #########################



##### Get the legend of a plot

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



plot0 = ggplot(subset(LambdaEstimates, Method!="Kernel"),
               aes(x=Lambda, y=Estimates,
                   colour=Method, group=Method)) + 
  geom_errorbar(aes(ymin=SECIl, ymax=SECIu), 
                width=.02, position=pd) +
  #geom_errorbar(aes(ymin=BootCIl, ymax=BootCIu), width=.1) +
  #geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.1) +
  geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("CIs Using Bootstrapped SE")+
  theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  #eliminates background, gridlines, and chart border
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1)
  ) 

legend = get_legend(plot0)

plot1 = plot0 + theme(legend.position = "none")

plot2 = ggplot(subset(LambdaEstimates, Method!="Kernel"),
               aes(x=Lambda, y=Estimates, 
                   colour=Method, group=Method)) + 
  #geom_errorbar(aes(ymin=max(Estimates - 1.96*BootSE, 0), 
  #                  ymax=min(Estimates + 1.96*BootSE, 1)), width=.1) +
  geom_errorbar(aes(ymin=BootCIl, ymax=BootCIu), 
                width=.02, position=pd) +
  #geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.1) +
  geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle("Bootstrapped CIs")+
  theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  #eliminates background, gridlines, and chart border
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1),
    legend.position = "none"
  ) 

plot3 = ggplot(subset(LambdaEstimates, Method!= "Kernel"), 
               aes(x=Lambda, y=Estimates, 
                   colour=Method, group=Method)) + 
  #geom_errorbar(aes(ymin=max(Estimates - 1.96*BootSE, 0), 
  #                  ymax=min(Estimates + 1.96*BootSE, 1)), width=.1) +
  #geom_errorbar(aes(ymin=BootCIl, ymax=BootCIu), width=.1) +
  geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.02, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("GLM CIs")+
  theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  #eliminates background, gridlines, and chart border
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1),
    legend.position = "none"
  ) 

plot4 = ggplot(subset(LambdaEstimates, Method!="Kernel"),
               aes(x=Lambda, y=Estimates, 
                   colour=Method, group=Method)) + 
  #geom_errorbar(aes(ymin=max(Estimates - 1.96*BootSE, 0), 
  #                  ymax=min(Estimates + 1.96*BootSE, 1)), width=.1) +
  geom_errorbar(aes(ymin=BootBasicCIl, ymax=BootBasicCIu), 
                width=.02, position=pd) +
  #geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.1) +
  geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle("Bootstrapped CIs")+
  theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  #eliminates background, gridlines, and chart border
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1),
    legend.position = "none"
  ) 



pdf("EstimatingLambdaBiggerData.pdf", width = 7, height = 5)
grid.arrange(plot1, plot2, plot3, plot4, legend, nrow = 2)
dev.off()


plot0 = ggplot(subset(LambdaEstimates, Method!= "Kernel"),
               aes(x=Lambda, y=Estimates,
                   colour=Method, group=Method)) + 
  geom_errorbar(aes(ymin=SECIl, ymax=SECIu), 
                width=.02, position=pd) +
  #geom_errorbar(aes(ymin=BootCIl, ymax=BootCIu), width=.1) +
  #geom_errorbar(aes(ymin=GLMCIl, ymax=GLMCIu), width=.1) +
  geom_line(position=pd) +
  geom_point(aes(shape = Method), position=pd, size=2) +
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("CIs Using Bootstrapped SE")+
  theme_bw() +
  scale_shape_manual(values=c(15:18, 8)) +
  scale_colour_manual(values=cbPalette) +
  #eliminates background, gridlines, and chart border
  theme(
    panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_rect(colour="black",  linetype = 1)
  ) 

plot0 = plot0 + labs(x = "Signal Strength", y = "Estimate",
                     title = "Estimated Signal Strength",
                     subtitle = "Error bars show standard error") 
ggsave(plot0, filename = "EstimatingLambdaSE.png",
       width = 5, height = 3)
ggsave(plot0, filename = "EstimatingLambdaSE.pdf",
       width = 5, height = 3)
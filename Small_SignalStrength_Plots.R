library(ggplot2)
library(PropCIs)
lambda = c(0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.01, 0.00, 0.005, 0.001)
iter.lambda = c(5, 7, 9, 10, 8)
lambda.plot = lambda[iter.lambda]


reject.table.n4 = matrix(c(46, 6, 2, 0, 2, 
                           33, 5, 2, 0, 1,
                           48, 5, 3, 1, 1), ncol = 3)

rownames(reject.table.n4) = lambda.plot
colnames(reject.table.n4) = c("AUC", "MC", "LRT")


reject.table.n5 = matrix(c(50, 15, 9, 1, 1, 
                           50, 12, 5, 2, 1,
                           50, 17, 9, 2, 2), ncol = 3)

rownames(reject.table.n5) = lambda.plot
colnames(reject.table.n5) = c("AUC", "MC", "LRT")


reject.table.n5Double = matrix(c(50, 40, 7, 2, 1, 
                                 50, 24, 6, 1, 1,
                                 50, 43, 8, 1, 1), ncol = 3)

rownames(reject.table.n5Double) = lambda.plot
colnames(reject.table.n5Double) = c("AUC", "MC", "LRT")

reject.table.n6 = matrix(c(50, 50, 26, 2, 2, 
                           50, 50, 18, 1, 3,
                           50, 50, 35, 2, 2), ncol = 3)

rownames(reject.table.n6) = lambda.plot
colnames(reject.table.n6) = c("AUC", "MC", "LRT")


plot.df = data.frame(Power = c(2*c(reject.table.n4), 2*c(reject.table.n5),
                               2*c(reject.table.n5Double), 2*c(reject.table.n6)),
                     SampleSize = c(rep("n", 3*5), rep("10xn", 3*5),
                                    rep("20xn", 3*5), rep("100xn", 3*5)),
                     Lambda = as.factor(rep(lambda.plot, 3*4)),
                     Method = rep(c(rep("AUC", 5), rep("MCE", 5), rep("LRT", 5)),4))
plot.df$SampleSize = factor(plot.df$SampleSize, levels = c("n", "10xn", "20xn", "100xn"))
plot.df$Upper.CI = sapply(plot.df$Power/2, 
                          function(x){return(exactci(x, 50,conf.level=0.95)$conf.int[2])})
plot.df$Lower.CI = sapply(plot.df$Power/2, 
                          function(x){return(exactci(x, 50,conf.level=0.95)$conf.int[1])})

ggplot(plot.df, aes(x=SampleSize, y=Power, fill=Lambda)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  facet_grid(rows = vars(Method)) +
  theme_minimal()
  
plot_power = ggplot(plot.df, aes(x=Lambda, y=Power, fill=SampleSize)) +
  geom_bar(stat="identity", color="black", position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(aes(ymin = 100*Lower.CI, 
                    ymax = 100*Upper.CI), 
                width = 0.2, position = position_dodge(0.8)) +
  facet_grid(rows = vars(Method)) +
  labs(x = "Signal Strength (\u03bb)",
       fill = "Sample Size")
  
ggsave("PowerPlot.jpg", plot_power, width = 10, height = 5)

ggplot(plot.df, aes(x=Lambda, y=Power, color=SampleSize)) +
  geom_point(position = position_dodge(0.3))+
  geom_errorbar(aes(ymin = Power-(Power*(1 - (Power/100))), 
                    ymax = Power+(Power*(1 - (Power/100)))), 
                width = 0.2, position = position_dodge(0.3)) +
  facet_grid(rows = vars(Method)) +
  theme_minimal() 

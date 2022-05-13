source("Test_Functions.R")

higgsBoson = read.csv("atlas-higgs-challenge-2014-v2.csv") # Higgs Boson Simulated Data 
lambda = c(0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.01, 0) # Signal Strengths for detection
data = higgsBoson[which(higgsBoson$PRI_jet_num == 2),] # Only selecting channels with 2 jets
data_b = data[which(data$Label == "b"),] # Background Data
data_s = data[which(data$Label == "s"),] # Signal Data
m = 20000 # Number of Background Training Samples
n = (nrow(data_b) - 2*m)/2 # Number of Background Test Samples
B = 50 # 50 iterations for Power
NBoot = 1000 # Number of bootstrap and permutation iterations
cols = c(15:23,25:31) # only primary variables from the data

# Log transformation for momentums and energies
logind = c(1,4,7,9, 10, 13, 16) 
logdata_b = data_b
logdata_b[,cols[logind]] = log(data_b[,cols[logind]])

logdata_s = data_s
logdata_s[,cols[logind]] = log(data_s[,cols[logind]])


# Selecting the primary variables
logdata_b = logdata_b[,c(cols,32)]
logdata_s = logdata_s[,c(cols,32)]

# Transforming the phi angles to difference b/w the phi angles 
# and the phi angle for leading jet
logdata_b[,c(3,6,8,15)] = ((logdata_b[,c(3,6,8,15)] + pi - 
                              logdata_b[,12])%%(2*pi)) - pi
logdata_b = logdata_b[,-12]
logdata_s[,c(3,6,8,15)] = ((logdata_s[,c(3,6,8,15)] + pi - 
                              logdata_s[,12])%%(2*pi)) - pi
logdata_s = logdata_s[,-12]

k = c(15, 20, 25, 30)

# Saving the data for use in the experiments in Section 6

save( lambda, n, m, B, NBoot, cols, logind,
      logdata_b, logdata_s, k,
      file = paste0("InitializedData.Rdata"))
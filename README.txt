Code files for Section 6 of Model-Independent Detection of New Physics Signals
using Interpretable Semi-Supervised Classifier Tests

Source for all functions: 
Test_Functions.R - contains all the required functions for all the methods in the paper.

Data: It is in folder Higgs_Boson_DataSet

Section 6.2.1 and 6.2.2 experiments: 
Signal_Detection_Tests.R - Actually runs the tests for the experiment settings.
Signal_Detection_Experiment_Analysis_Plots.R - Creates tables and plots from the experiment outputs.

Section 6.2.3 experiments:
GMMApproximation.RData - Data for the simulations. Contains parameters of the GMM models.
Small_SignalStrength_Tests.R - Actually runs the tests for the experiment settings.
Small_SignalStrength_Plots.R - Creates the plot from the experiment outputs.

Section 6.3 Signal Strength Estimation:
Lambda_Original_Estimates.R - Gives the estimates on the data.
Lambda_Bootstrap_Estimates.R - Gives the estimates on Bootstrapped data.
Lambda_Plots.R - Creates the plots.

Section 6.4 Active Subspace Methods:
Active_subspace_functions.R - File with the functions required for the active subspace methods.
Mixed_Active_Subspace.R - Finding the active subspace on the Higgs Boson data.

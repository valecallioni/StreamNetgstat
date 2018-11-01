# Hypotheses: 
# - We have to deal only with numerical variables (not with categorical ones)
# - We have one network at a time

library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/src copia/interface.so")
source("/vagrant/PACSProject/R/get_SSN_model.R")
load("/vagrant/PACSProject/Data/middlefork.RData")

# library(SSN)
# 
# # Importing the Missouri data set
# missourip = importSSN("MissouriHW.ssn")
# 
# # Importing prediction point data sets
# missourip = importPredpts(missourip, "preds", "ssn")
# 
# # Compute the additive function value (to obtain the spatial weights)
# # based on h2oAreaKm2 (the watershed area), which is a column of the dataframe data,
# # slot of the SpatialStreamNetowrk class, if the variable afvArea is not yet present
# if (!is.element("afvArea", colnames(missourip@data))){
#   names(missourip@data)
#   missourip = additive.function(missourip, "h2oAreaKm2", "computed.adf")
#   names(missourip@data)
# }

# Create the matrix N, the matrix containing the distances among observed points, 
# prediction points and among the two groups
# .... call R function

# Missouri arguments: 17, c("STREAM_AUG", "ELEV", "SLOPE"), "afvArea"
# MiddleFork arguments: 1 o 2c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), c("Summer_mn", "ELEV_DEM", "SLOPE"), "afvArea"

modelParam = get_SSN_model(c("Summer_mn", "ELEV_DEM", "SLOPE"), "afvArea", 
                            c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
                            bin_tables, network_data, obs_points, pred_points, obs_data, pred_data)

print("Theta:")
print(modelParam$optTheta)
print("Beta:")
print(modelParam$betaValues)
#print("Covariance matrix:")
#print(modelParam$covMatrix)

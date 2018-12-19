# Hypotheses: 
# - We have to deal only with numerical variables (not with categorical ones)
# - We have one network at a time

library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/StreamNetgstat/src/interface.so")
source("/vagrant/PACSProject/get_SSN_model.R")
source("/vagrant/PACSProject/get_SSN_model_kriging.R")
source("/vagrant/PACSProject/do_SSN_kriging.R")
load("/vagrant/PACSProject/Data/middlefork_singleNet2.RData")
# load("/vagrant/PACSProject/Data/matrices_missouri.RData")

# file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"), 
#           to = tempdir(), recursive = TRUE, copy.mode = FALSE)
# setwd(tempdir())
# 
# mf04p = importSSN("MiddleFork04.ssn", predpts = "pred1km")
# mf04p = importPredpts(mf04p, "Knapp", "ssn")
# mf04p = importPredpts(mf04p, "CapeHorn", "ssn")

# Missouri arguments: 17, c("STREAM_AUG", "ELEV", "SLOPE"), "afvArea"


# result = get_SSN_model(c("STREAM_AUG", "ELEV"), "afvArea",
#                             c("Exponential.tailup", "Exponential.taildown"),
#                             net_num, bin_tables, network_data, obs_points, obs_data, matrices = dist_matrices)
# 
# result = do_SSN_kriging(net_num, bin_tables, network_data,
#                        obs_points, pred_points, obs_data, pred_data, 
#                        c("STREAM_AUG", "ELEV"), "afvArea",
#                        c("Exponential.tailup", "Exponential.taildown"), TRUE,
#                        theta, covMatrix, matrices = dist_matrices)

# result = get_SSN_model_kriging(c("STREAM_AUG", "ELEV"), "afvArea",
#                             c("Exponential.tailup", "Exponential.taildown"),
#                             net_num, bin_tables, network_data, obs_points, pred_points, obs_data, pred_data)

network_data_loc$NetworkID = rep(1,length(network_data_loc$NetworkID))
obs_points_loc$NetworkID = rep(1,length(obs_points_loc$NetworkID))
pred_points_loc$NetworkID = rep(1,length(pred_points_loc$NetworkID))

network_data_loc = data.matrix(network_data_loc)
obs_points_loc = data.matrix(obs_points_loc)
obs_data_loc = data.matrix(obs_data_loc)
pred_points_loc = data.matrix(pred_points_loc)
pred_data_loc = data.matrix(pred_data_loc)

result = .Call("getSSNModelKriging_MultipleNets", c(1), list(bin.table), network_data_loc,
                obs_points_loc, pred_points_loc, obs_data_loc, pred_data_loc, c(varNames, weightVar), CorModels, TRUE, NULL, NULL)


print("Theta:")
print(result$optTheta)
print("Beta:")
print(result$betaValues)
setwd("/vagrant/PACSProject/")
# write(t(result$covMatrix), file = "covMat.txt", ncolumns = 45)
write(t(result$predictions), file = "predData.txt", ncolumns = 2)

# matrices = .Call("createDistanceMatrices", net_num, bin_tables, network_data, obs_points)
# print(dim(matrices$distGeo))
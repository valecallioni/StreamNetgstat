# Hypotheses: 
# - We have to deal only with numerical variables (not with categorical ones)
# - We have one network at a time

library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/StreamNetgstat/src/interface.so")
# load("/vagrant/PACSProject/Data/middlefork_singleNet2.RData")
# load("/vagrant/PACSProject/Data/matrices_missouri.RData")
load("/vagrant/PACSProject/Data/Rdata_model_kriging.RData")


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

# network_data_loc$NetworkID = rep(1,length(network_data_loc$NetworkID))
# obs_points_loc$NetworkID = rep(1,length(obs_points_loc$NetworkID))
# pred_points_loc$NetworkID = rep(1,length(pred_points_loc$NetworkID))
# 
# network_data_loc = data.matrix(network_data_loc)
# obs_points_loc = data.matrix(obs_points_loc)
# obs_data_loc = data.matrix(obs_data_loc)
# pred_points_loc = data.matrix(pred_points_loc)
# pred_data_loc = data.matrix(pred_data_loc)
# 
# result = .Call("getSSNModelKriging_MultipleNets", c(1), list(bin.table), network_data_loc,
#                 obs_points_loc, pred_points_loc, obs_data_loc, pred_data_loc, c(varNames, weightVar), CorModels, TRUE, NULL, NULL)


# print("Theta:")
# print(result$optTheta)
# print("Beta:")
# print(result$betaValues)
# setwd("/vagrant/PACSProject/")
# write(t(result$covMatrix), file = "covMat.txt", ncolumns = 45)
# write(t(result$predictions), file = "predData.txt", ncolumns = 2)

# matrices = .Call("createDistanceMatrices", net_num, bin_tables, network_data, obs_points)
# print(dim(matrices$distGeo))


result = .Call("getSSNModelKriging_MultipleNets", net_num, bin_tables, data.matrix(network_data),
               data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data),
               c(varNames, weightVar), CorModels, useNugget, matrices, bounds, useCholeskyDec)
print(result$predictions)

# result = .Call("doSSNKriging_MultipleNets", net_num, bin_tables, data.matrix(network_data),
#                data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data), 
#                c(varNames, weightVar), CorModels, useNugget, theta, covMat, matrices)
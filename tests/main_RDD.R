library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/StreamNetgstat/src/interface.so")
load("/vagrant/PACSProject/Data/RData_dokriging.RData")

# res_pred = .Call("getSSNModelKriging_SingleNet", bin_tables[[1]], data.matrix(network_data),
#                  data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data), 
#                  c(varNames, weightVar), CorModels, useNugget, matrices, bounds)

# diag(covMat) = diag(covMat) + rep(theta[3], dim(covMat)[1])
# result = .Call("doSSNKriging_SingleNet", bin.table, data.matrix(network_data),
#                data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data),
#                c(varNames, weightVar), CorModels, useNugget, theta, covMat, matrices)
result = .Call("getSSNModelKriging_SingleNet", bin.table, data.matrix(network_data),
               data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data), 
               c(varNames, weightVar), CorModels, useNugget = TRUE, matrices = NULL, bounds = NULL)
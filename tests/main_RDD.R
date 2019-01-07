library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/StreamNetgstat/src/interface.so")
load("/vagrant/PACSProject/Data/RData_for_kriging_pred.RData")

res_pred = .Call("doSSNKriging_MultipleNets", sort(c(netPred,netObs)), bin_tables, data.matrix(network_data[[2]]),
                 data_obs[[closest_net]]$points, data.matrix(pred_points[[2]]), data_obs[[closest_net]]$data, data.matrix(pred_data[[2]][,-1]), 
                 c(varNames, weightVar), CorModels, useNugget, models[[closest_net]]$param, models[[closest_net]]$covMat, NULL)

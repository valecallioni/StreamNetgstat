# Hypotheses: 
# - We have to deal only with numerical variables (not with categorical ones)
# - We have one network at a time

library(Rcpp)
library(RcppEigen)
dyn.load("/vagrant/PACSProject/StreamNetgstat/src/interface.so")
source("/vagrant/PACSProject/StreamNetgstat/tests/get_SSN_model.R")
load("/vagrant/PACSProject/Data/missouriObs.RData")

# file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"), 
#           to = tempdir(), recursive = TRUE, copy.mode = FALSE)
# setwd(tempdir())
# 
# mf04p = importSSN("MiddleFork04.ssn", predpts = "pred1km")
# mf04p = importPredpts(mf04p, "Knapp", "ssn")
# mf04p = importPredpts(mf04p, "CapeHorn", "ssn")

# Missouri arguments: 17, c("STREAM_AUG", "ELEV", "SLOPE"), "afvArea"


result = get_SSN_model(c("STREAM_AUG", "ELEV"), "afvArea", 
                            c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
                            net_num, bin_tables, network_data, obs_points, obs_data)

print("Theta:")
print(result$optTheta)
print("Beta:")
print(result$betaValues)
# setwd("/vagrant/PACSProject/")
# write(t(result$covMatrix), file = "covMat.txt", ncolumns = 1143)
#write(t(result$predictions), file = "predData.txt", ncolumns = 2)
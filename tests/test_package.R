setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat")
library(roxygen2)
roxygenise()

library(devtools)
devtools::install_github("valecallioni/StreamNetgstat")

library(StreamNetgstat)
#load("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data/missouri.RData")
#source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/R/get_SSN_model.R")
source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/tests/funSSNPackage.R")

# file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"),
#           to = tempdir(), recursive = TRUE, copy.mode = FALSE)
# setwd(tempdir())
# ssn = importSSN("MiddleFork04.ssn", predpts = "pred1km")

setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")
ssn = importSSN("MissouriHW.ssn", predpts = "preds")



# Prova:
result.Exponential.cpp = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                                        CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"))

result.Spherical.cpp = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                                      CorModels = c("Spherical.tailup", "Spherical.taildown", "Spherical.Euclid"))

result.LinearSill.cpp = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                                     CorModels = c("LinearSill.tailup", "LinearSill.taildown"))

result.Exp.up.down.R = funSSNPackage(ssn, formula = STREAM_AUG ~ ELEV, predname = NULL)



###### Confronto:

# Model
library(benchr)
benchmark(
  times = 1,
  cppPackage = get_SSN_model_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                             CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), predpts = "preds"),
  RPackage = funSSNPackage(ssn, formula = STREAM_AUG ~ ELEV, predname = "preds")
)


##### Create distance-matrices - Fit model - Do kriging
dist_matrices = get_plots(ssn, "STREAM_AUG", F)
model = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea", CorModels = c("LinearSill.tailup", "Exponential.taildown"), matrices = dist_matrices)
kriging = do_SSN_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea", predpts = "preds", CorModels = c("LinearSill.tailup", "Exponential.taildown"), theta = model$optTheta, covMat = model$covMatrix, matrices = dist_matrices)

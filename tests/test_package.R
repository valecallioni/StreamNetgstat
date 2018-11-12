setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat")
library(roxygen2)
roxygenise()

library(devtools)
devtools::install_github("valecallioni/StreamNetgstat")


library(benchr)

library(StreamNetgstat)
#load("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data/missouri.RData")
source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/R/get_SSN_model.R")
source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/tests/funSSNPackage.R")

library(SSN)
# file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"), 
#           to = tempdir(), recursive = TRUE, copy.mode = FALSE)
# setwd(tempdir())

setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")

ssn = importSSN("Salmon.ssn")
result = get_SSN_model(ssn, varNames = c("Stream_Aug", "ELEV"), weightVar = "afvArea",
                       CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"))

benchmark(
  RPackage = funSSNPackage(ssn, "preds"),
  cppPackage = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                         CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
                         predpts = "preds", doKriging = TRUE)
)

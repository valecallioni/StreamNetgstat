setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat")
library(roxygen2)
roxygenise()

library(devtools)
devtools::install_github("valecallioni/StreamNetgstat")

library(StreamNetgstat)
#load("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data/missouri.RData")
#source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/R/get_SSN_model.R")
source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/tests/funSSNPackage.R")

library(SSN)
# file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"),
#           to = tempdir(), recursive = TRUE, copy.mode = FALSE)
# setwd(tempdir())
# ssn = importSSN("MiddleFork04.ssn", predpts = "pred1km")

setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")
ssn = importSSN("MissouriHW.ssn")



# Prova:
result = get_SSN_model(ssn, varNames = c("Summer_mn", "ELEV_DEM"), weightVar = "afvArea",
                       CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"))




###### Confronto:

# Model
start.time.cpp = Sys.time()
cppPackage = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                           CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"))
end.time.cpp = Sys.time()
time.taken.cpp = end.time.cpp - start.time.cpp

start.time.R = Sys.time()
RPackage = funSSNPackage(ssn, formula = STREAM_AUG ~ ELEV, predname = NULL)
end.time.R = Sys.time()
time.taken.R = end.time.R - start.time.R


library(benchr)
benchmark(
  times = 1,
  RPackage = funSSNPackage(ssn, formula = STREAM_AUG ~ ELEV, predname = NULL),
  cppPackage = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                         CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"))
)

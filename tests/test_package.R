setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat")
library(roxygen2)
roxygenise()

library(devtools)
devtools::install_github("valecallioni/StreamNetgstat")

library(StreamNetgstat)
#load("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data/missouri.RData")
#source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/R/get_SSN_model.R")
source("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/StreamNetgstat/tests/funSSNPackage.R")

file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"),
          to = tempdir(), recursive = TRUE, copy.mode = FALSE)
setwd(tempdir())
ssn = importSSN("MiddleFork04.ssn", predpts = "pred1km")

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
  times = 5,
  cppPackage = get_SSN_model_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                             CorModels = c("Exponential.tailup", "Exponential.taildown"), predpts = "preds"),
  RPackage = funSSNPackage(ssn, formula = STREAM_AUG ~ ELEV, predname = "preds")
)


##### Create distance-matrices - Fit model - Do kriging
dist_matrices = get_plots(ssn, "STREAM_AUG", F)
model = get_SSN_model(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea", CorModels = c("LinearSill.tailup", "Exponential.taildown"), matrices = dist_matrices)
kriging = do_SSN_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea", predpts = "preds", CorModels = c("LinearSill.tailup", "Exponential.taildown"), theta = model$optTheta, covMat = model$covMatrix, matrices = dist_matrices)



#### SIMULATED DATA
set.seed(12)
setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")
simulated.ssn = createSSN(n = c(50,60), obsDesign = binomialDesign(c(70,85)), predDesign = binomialDesign(c(80,90)), importToR = T, path = "SimNetwork_Results.ssn", treeFunction = iterativeTreeLayout)

setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")
simulated.ssn = importSSN("SimNetwork_Results.ssn", predpts = "preds")

setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Report/img")
pdf(file = "Networks.pdf",width=7,height=6,paper="special") 
plot(simulated.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8, lineCol = "blue", pch = 1, 
     xlab = "x-coordinate (m)", ylab= "y-coordinate (m)", cex = 2)
dev.off()

createDistMat(simulated.ssn, "preds", o.write = TRUE, amongpreds = TRUE)
DFobs = getSSNdata.frame(simulated.ssn, "Obs")
DFpred = getSSNdata.frame(simulated.ssn, "preds")

DFobs[, "X1"] = rnorm(length(DFobs[,1]))
DFpred[, "X1"] = rnorm(length(DFpred[,1]))
DFobs[, "X2"] = rnorm(length(DFobs[,1]))
DFpred[, "X2"] = rnorm(length(DFpred[,1]))

set.seed(102)
sim.out = SimulateOnSSN(simulated.ssn, ObsSimDF = DFobs, PredSimDF = DFpred, PredID = "preds",
                        formula = ~ X1 + X2, coefficients = c(2, 4, 0), 
                        CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), 
                        use.nugget = T, CorParms = c(3, 8, 2, 10, 4, 12, 0.1), addfunccol = "addfunccol")
sim.ssn = sim.out$ssn.object
simDFobs = getSSNdata.frame(sim.ssn, "Obs")
simDFpred = getSSNdata.frame(sim.ssn, "preds")
simpreds = simDFpred[, "Sim_Values"]
simDFpred[, "Sim_Values"] = NA
sim.ssn = putSSNdata.frame(simDFpred, sim.ssn, "preds")


pdf(file = "Sim_Values.pdf",width=7,height=6,paper="special") 
brks = plot(sim.ssn, "Sim_Values", lwdLineCol = "addfunccol", lwdLineEx = 15,
            lineCol = "black", xlab = "x-coordinate (m)", ylab= "y-coordinate (m)")
dev.off()

dist_matrices = get_plots(sim.ssn, "Sim_Values", T, nlag_Torg = 4)
model = get_SSN_model(sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol", CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), matrices = dist_matrices)
kriging = do_SSN_kriging(sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol", predpts = "preds", CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), theta = model$optTheta, covMat = model$covMatrix, matrices = dist_matrices)

analysis = get_SSN_model_kriging(ssn = sim.ssn, varNames = c("Sim_Values", "X1", "X2"),
                                 weightVar = "addfunccol", predpts = "preds",
                                 CorModels = c("Exponential.tailup", "Exponential.taildown"))

pdf(file = "Predictions.pdf",width=7,height=6,paper="special") 
plot.predictions(analysis$ssn.object, VariableName = "Sim_Values", predpts = "preds", SEcex.max = 1, SEcex.min = 0.5/3*2, breaktype = "user", brks = brks)
dev.off()

library(benchr)
benchmark(
  times = 5,
  cppPackage = get_SSN_model_kriging(sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol",
                                     CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), predpts = "preds"),
  RPackage = funSSNPackage(sim.ssn, formula = Sim_Values ~ X1 + X2, predname = "preds")
)

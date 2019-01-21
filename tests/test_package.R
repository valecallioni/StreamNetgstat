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

setwd("~/Desktop/OneDrive - Politecnico di Milano/Tesi/RDD/Data")
simulated.ssn = importSSN("SimNetwork_systematicDesign_02_preds.ssn", predpts = "preds")

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
                        use.nugget = T, CorParms = c(5, 8, 4, 10, 2, 9, 0.1), addfunccol = "addfunccol")
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

pdf(file = "Torg_EmpSemivar.pdf",width=10,height=6,paper="special") 
dist_matrices = get_plots(sim.ssn, "Sim_Values", T)
dev.off()

model = get_SSN_model(ssn = sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol", CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), matrices = dist_matrices, bounds = NULL, useCholeskyDec = T)
kriging = do_SSN_kriging(sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol", predpts = "preds", CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), theta = model$modelParam, covMat = model$modelCovariance, matrices = dist_matrices)

sim.obs = sim.ssn
sim.obs@predpoints@SSNPoints[[1]]@network.point.coords = sim.obs@obspoints@SSNPoints[[1]]@network.point.coords
sim.obs@predpoints@SSNPoints[[1]]@point.coords = sim.obs@obspoints@SSNPoints[[1]]@point.coords
sim.obs@predpoints@SSNPoints[[1]]@point.data = sim.obs@obspoints@SSNPoints[[1]]@point.data
sim.obs@predpoints@SSNPoints[[1]]@points.bbox = sim.obs@obspoints@SSNPoints[[1]]@points.bbox

kriging_obs = do_SSN_kriging(sim.obs, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol", predpts = "preds", CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), theta = model$modelParam, covMat = model$modelCovariance, matrices = dist_matrices)
kriging_obs@predpoints@SSNPoints[[1]]@point.data = cbind.data.frame(kriging_obs@predpoints@SSNPoints[[1]]@point.data,
                                                                            abs(sim.obs@obspoints@SSNPoints[[1]]@point.data["Sim_Values"] - kriging_obs@predpoints@SSNPoints[[1]]@point.data["Sim_Values_pred"]))
colnames(kriging_obs@predpoints@SSNPoints[[1]]@point.data)[length(colnames(analysis$ssn.object@predpoints@SSNPoints[[1]]@point.data))] = "Diff_values_pred"
plot.predictions(kriging_obs, "Diff_values", predpts = "preds", VarPlot = "Predictions", nclasses = 1)

pdf(file = "Predictions_new.pdf",width=7,height=6,paper="special") 
plot.predictions(kriging, VariableName = "Sim_Values", predpts = "preds", SEcex.max = 1, SEcex.min = 0.5/3*2, breaktype = "user", brks = brks)
plot.predictions(kriging_obs, "Diff_values", predpts = "preds", VarPlot = "Predictions", nclasses = 1, color.palette = c("black"), add=T, SEcex.min = 0.2, SEcex.max = 0.2)
dev.off()

kriging@predpoints@SSNPoints[[1]]@point.data = cbind.data.frame(kriging@predpoints@SSNPoints[[1]]@point.data,
                                                                abs(simpreds - kriging@predpoints@SSNPoints[[1]]@point.data["Sim_Values_pred"]))
colnames(kriging@predpoints@SSNPoints[[1]]@point.data)[length(colnames(analysis$ssn.object@predpoints@SSNPoints[[1]]@point.data))] = "Diff_values_pred"
pdf(file = "Prediction_error.pdf",width=7,height=6,paper="special") 
plot.predictions(kriging, "Diff_values", predpts = "preds", VarPlot = "Predictions", nclasses = 4)
dev.off()


library(benchr)
benchmark(
  times = 5,
  cppPackage = get_SSN_model_kriging(sim.ssn, varNames = c("Sim_Values", "X1", "X2"), weightVar = "addfunccol",
                                     CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), predpts = "preds", singleNet = 1),
  RPackage = funSSNPackage(sim.ssn, formula = Sim_Values ~ X1 + X2, predname = "preds")
)

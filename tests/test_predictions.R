library(StreamNetgstat)
library(SSN)

file.copy(system.file(file.path("lsndata", "MiddleFork04.ssn"), package = "SSN"),
          to = tempdir(), recursive = TRUE, copy.mode = FALSE)
setwd(tempdir())
ssn = importSSN("MiddleFork04.ssn")


# setwd("~/Desktop/OneDrive - Politecnico di Milano/PACS/Programming/PACSProject/Data")
# ssn = importSSN("MissouriHW.ssn", predpts = "preds")



## Make predictions on some points that actually were observed
ssn@predpoints@SSNPoints = ssn@obspoints@SSNPoints
col = which(colnames(ssn@predpoints@SSNPoints[[1]]@point.data)=="Summer_mn")
ssn@predpoints@SSNPoints[[1]]@point.data = ssn@predpoints@SSNPoints[[1]]@point.data[,-col]
ssn@predpoints@ID = "obs_predict"

# pred = sort(sample(1:45, size = 10))
# ssn@obspoints@SSNPoints[[1]]@point.coords = ssn@obspoints@SSNPoints[[1]]@point.coords[-pred, ]
# ssn@obspoints@SSNPoints[[1]]@point.data = ssn@obspoints@SSNPoints[[1]]@point.data[-pred, ]
# ssn@obspoints@SSNPoints[[1]]@network.point.coords = ssn@obspoints@SSNPoints[[1]]@network.point.coords[-pred, ]
# ssn@predpoints@SSNPoints[[1]]@point.coords = ssn@predpoints@SSNPoints[[1]]@point.coords[pred, ]
# col = which(colnames(ssn@predpoints@SSNPoints[[1]]@point.data)=="Summer_mn")
# ssn@predpoints@SSNPoints[[1]]@point.data = ssn@predpoints@SSNPoints[[1]]@point.data[pred, -col]
# ssn@predpoints@SSNPoints[[1]]@network.point.coords = ssn@predpoints@SSNPoints[[1]]@network.point.coords[pred, ]


## STREAMNETGSTAT
predictions = get_SSN_model_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea",
                                    predpts = "predpts", CorModels = c("Exponential.tailup", "Exponential.taildown"))
pred_new = do_SSN_kriging(ssn, varNames = c("STREAM_AUG", "ELEV"), weightVar = "afvArea", predpts = "preds", 
                          CorModels = c("Exponential.tailup", "Exponential.taildown"), useNugget = TRUE,
                          theta = predictions$optTheta, covMat = predictions$covMatrix)


## SSN
createDistMat(ssn, predpts = "obs_predict", o.write = TRUE, amongpreds = FALSE)
ssn.glmssn = glmssn(Summer_mn ~ ELEV_DEM, ssn, 
                     CorModels = c("Exponential.tailup", "Exponential.taildown"), 
                    addfunccol = "afvArea")
ssn.pred = predict(ssn.glmssn, "obs_predict")

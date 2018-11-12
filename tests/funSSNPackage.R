funSSNPackage = function(ssn, predname){
  
  createDistMat(ssn, predpts = predname, o.write = TRUE, amongpreds = FALSE)
  
  ssn.glmssn1 = glmssn(STREAM_AUG ~ ELEV, ssn, 
                        CorModels = c("Exponential.tailup", "Exponential.taildown", 
                                      "Exponential.Euclid"), addfunccol = "afvArea")
  
  ssn.pred = predict(ssn.glmssn1, predname)
  
  return (ssn.pred)
  
}
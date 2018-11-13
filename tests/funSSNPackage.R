funSSNPackage = function(ssn, formula, predname){
  
  createDistMat(ssn, predpts = predname, o.write = TRUE, amongpreds = FALSE)
  
  ssn.glmssn1 = glmssn(formula, ssn, 
                        CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"), addfunccol = "afvArea")
  
  if (!is.null(predname)){
    ssn.pred = predict(ssn.glmssn1, predname)
    return (ssn.pred)
  }
  
  return (ssn.glmssn1)
  
}
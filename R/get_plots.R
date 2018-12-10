get_plots = function(ssn, ResponseName, Euclidean = FALSE, EmpVarMeth = "MethMoment")
{
  
  if(Euclidean){
   result = emp_semivario_and_torg(ssn = ssn, ResponseName = ResponseName, EmpVarMeth = EmpVarMeth)
   par(mfrow=c(1,2))
   plot(result$Torg, main = "Torgegram") 
   plot(result$EmpSemiVar$distance, result$SemiVar$gamma, main = "Empirical Semivariogram based on Euclidean distance",
        xlab = "Distance", ylab = "Semi-Variance")
  }
  
  else {
    result = torgegram(ssn = ssn, ResponseName = ResponseName, EmpVarMeth = EmpVarMeth)
    plot(result$Torg, main = "Torgegram") 
  }
  
  return(result$distMatrices)
}
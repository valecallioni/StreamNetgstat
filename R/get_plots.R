#' Plot the Torgegram and the Empirical Semivariogram for a Spatial Stream Network object 
#'
#' @param ssn a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param ResponseName the name of the response variable.
#' @param Euclidean If \code{TRUE} the Empirical Semivariogram based on Euclidean distances is also plotted.
#' @param EmpVarMeth a string indicating the name of the method to be used for computing the semivariogram.
#' @return A vector of distance matrices.
#' @description Given a \link[SSN]{SpatialStreamNetwork-class} object, computes the Torgegram and the Empirical Semivariogram and provides the plots.
#' @details This function works on objects of \link[SSN]{SpatialStreamNetwork-class}.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18. DOI: 10.1198/jasa.2009.ap08248. Rejoinder pgs. 22–24.
#' @useDynLib StreamNetgstat
#' @export


get_plots = function(ssn, ResponseName, Euclidean = FALSE, maxlag = NULL, nlag = 6, 
                     inc = 0, nlagcutoff = 15, EmpVarMeth = "MethMoment")
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
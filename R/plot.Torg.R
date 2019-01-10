#' Plotting Method for Torgegram Objects
#'
#' @param x an object of class \link[SSN]{Torgegram}.
#' @param sp.relationship a string or character vector representing the in-stream spatial relationship to be plotted. "fc" specifies plotting of only flow-connected, and "fu" specifies plotting of only flow-unconnected. Default is both.
#' @param min.cex Minimum character expansion size of the plotting symbols. Symbols are scaled according to how many pairs of points went into computing each bin of the semivariogram. The bin with the smallest sample size will be be plotted with this symbol size. The default is 1.5.
#' @param max.cex Maximum character expansion size of the plotting symbols. Symbols are scaled according to how many pairs of points went into computing each bin of the semivariogram. The bin with the largest sample size will be be plotted with this symbol size. The default is 6.
#' @param leg.auto Logical, default to TRUE. Include a legend.
#' @param main Title for plot.
#' @param ylab Label for y-axis.
#' @param xlab Label for x-axis.
#' @return Plot of empirical semivariogram values
#' @description \code{plot}.Torgegram is a generic plot function that has been adapted for Torgegram objects, which are created using the Torgegram function. A Torgegram object stores information used to construct an empirical semivariogram based on hydrologic distance. The \code{plot.Torgegram} function allows the results to be presented separately for flow-connected and flow-unconnected sites.
#' @details The \code{Torgegram} function creates a list of distances and empirical semivariogram values, along with number of pairs of points in each bin, for both flow-connected and flow-unconnected sites. Flow-connected locations lie on the same stream network (share a common downstream junction) and water flows from one location to the other. Flow-unconnected locations also lie on the same stream network, but do not share flow. The output is of class Torgegram. This is the default plotting method for this class.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18. DOI: 10.1198/jasa.2009.ap08248. Rejoinder pgs. 22–24.
#' @useDynLib StreamNetgstat
#' @export

plot.Torg <- function(x, sp.relationship = c("fc", "fu"),
                           min.cex = 1.5, max.cex = 6, leg.auto = TRUE, main = "", ylab = "",
                           xlab = "Stream Distance")
{
  par.orig <- par(no.readonly = TRUE)
  if(is.null(sp.relationship) | any(is.na(sp.relationship)) |
     length(sp.relationship) > 2) return("sp.relationship mis-specified")
  if(length(sp.relationship)  == 2 & any(sp.relationship != c("fc", "fu")))
    return("sp.relationship mis-specified")
  if(length(sp.relationship)  == 1)
    if(sp.relationship != "fc" & sp.relationship != "fu")
      return("sp.relationship mis-specified")
  if(class(x) != "Torgegram") return(
    "Not a Torgegram object")
  ev <- x
  
  if(is.null(ev$call$EmpVarMeth)) ev$call$EmpVarMeth <- "MethMoment"
  if(is.null(ev$call$nlag)) ev$call$nlag <- 6
  if(is.null(ev$call$inc)) ev$call$inc <- 0
  if(is.null(ev$call$nlagcutoff)) ev$call$nlagcutoff <- 6
  if(is.null(as.list(match.call()[-1])$ylab)) {
    if(ev$call$EmpVarMeth == "Covariance") {
      ylab <- "Covariance"} else
        ylab <- "Semivariance"
  }
  if(is.null(as.list(match.call()[-1])$col)) {
    colr = c("blue","green2")
  }
  else {
    if(length(sp.relationship) == 2 & length(as.list(
      match.call()[-1])$col) == 1) {
      colr <- rep(as.character(as.list(match.call()[-1])$col),
                  times = 2)
    } else {
      colr <- as.character(as.list(match.call()[-1])$col)
      if(length(colr)>1) colr <- colr[-1]
    }
  }
  if(is.null(as.list(match.call()[-1])$pch)) {
    plch = c(19,19)
  }
  else {
    if(length(sp.relationship) == 2 & length(as.list(
      match.call()[-1])$pch) == 1) {
      plch <- rep(as.integer(as.character(as.list(match.call()[-1])$pch)),
                  times = 2)
    } else {
      plch <- as.character(as.list(match.call()[-1])$pch)
      if(length(plch)>1) plch <- plch[-1]
      plch <- as.integer(plch)
    }
  }
  if(length(sp.relationship)==2){
    if(is.null(as.list(match.call()[-1])$main)) {
      main <- paste("Estimation Method:", ev$call$EmpVarMeth)
    } else {
      main <- as.list(match.call()[-1])$main
    }
    plot(c(0, max(ev$distance.connect,ev$distance.unconnect)),
         c(min(0, min(ev$gam.connect,ev$gam.unconnect)),
           max(ev$gam.connect,ev$gam.unconnect)),
         type = "n",
         xlab = xlab,
         ylab = ylab,
         main = main)
    maxnp <- max(ev$np.connect,ev$np.unconnect)
    minnp <- min(ev$np.connect,ev$np.unconnect)
    np.range <- maxnp - minnp
    cex.inc <- max.cex - min.cex
    for(i in 1:length(ev$np.connect)) {
      points(ev$distance.connect[i],
             ev$gam.connect[i], pch = plch[1], col = colr[1],
             cex = min.cex + cex.inc*ev$np.connect[i]/maxnp)
    }
    for(i in 1:length(ev$np.unconnect)) {
      points(ev$distance.unconnect[i],
             ev$gam.unconnect[i], pch = plch[2], col = colr[2],
             cex = min.cex + cex.inc*ev$np.unconnect[i]/maxnp)
    }
    if(leg.auto)
      legend(x = (2/5)*max(ev$distance.connect,ev$distance.unconnect),
             y = 2*min(ev$gam.connect,ev$gam.unconnect),
             legend = c("Flow-connected", "Flow-unconnected"), bty = "n",
             pch = plch,
             col = c(colr[1], colr[2]), y.intersp = 0.2)
  } else if(length(sp.relationship)==1){
    if(sp.relationship == "fc") {
      if(is.null(as.list(match.call()[-1])$main)) {
        main <- paste("Flow-connected,  Estimation Method:",
                      ev$call$EmpVarMeth)
      } else {
        main <- as.list(match.call()[-1])$main
      }
      plot(c(0, max(ev$distance.connect)), c(0, max(ev$gam.connect)),
           type = "n", xlab = xlab, ylab = ylab,
           main = main)
      nlag <- length(ev$distance.connect)
      maxnp <- max(ev$np.connect)
      minnp <- min(ev$np.connect)
      np.range <- maxnp - minnp
      cex.inc <- max.cex - min.cex
      for(i in 1:nlag) {
        points(ev$distance.connect[i],
               ev$gam.connect[i], pch = plch[1],
               cex = min.cex + cex.inc*ev$np.connect[i]/maxnp,
               col = colr[1])
      }
    } else if(sp.relationship == "fu") {
      if(is.null(as.list(match.call()[-1])$main)) {
        main <- paste("Flow-unconnected,  Estimation Method:",
                      ev$call$EmpVarMeth)
      } else {
        main <- as.list(match.call()[-1])$main
      }
      plot(c(0, max(ev$distance.unconnect)), c(0, max(ev$gam.unconnect)),
           type = "n", xlab = xlab, ylab = ylab,
           main = main)
      nlag <- length(ev$distance.unconnect)
      maxnp <- max(ev$np.unconnect)
      minnp <- min(ev$np.unconnect)
      np.range <- maxnp - minnp
      cex.inc <- max.cex - min.cex
      for(i in 1:nlag) {
        points(ev$distance.unconnect[i],
               ev$gam.unconnect[i], pch = plch[1],
               cex = min.cex + cex.inc*ev$np.unconnect[i]/maxnp,
               col = colr[1])
      }
    } else return("sp.relationship mis-specified")
  }
}

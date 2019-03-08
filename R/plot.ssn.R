#' Plot the predicted values along a Spatial Stream Network object
#' 
#' @param ssn.object a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param VariableName name of variable to be plotted, it can be also a factor.
#' @param VarPlot a character argument that must be one of "Both", "Predictions", or "Standard Errors". Default is "Both", which colors predictions by their values and makes their size inversely proportional to the prediction standard errors.
#' @param breaktype the method for breaking the predictions (or standard errors) into classes for coloring while plotting. A character argument that must be one of "quantile" (default), "even", or "user".
#' @param brks if breaktype = "user", the break values must be specified here as a vector or matrix using c(...) or cbind(...). The sorted unique values are used as break points (together with the min and max of the variable being plotted if required)
#' @param nclasses the number of classes for coloring the predictions (or standard errors) according to their value.  The default is 10. If brks = c(...) is specified, then nclasses is automatically set to the number of breaks + 1.
#' @param color.palette a color palette for plotting points. The default is rainbow(nclasses, start = .66, end = .99). The number of colors should equal to the number of classes. See \code{\link{palette}} for many ways to create palettes.
#' @param SEcex.min if VarPlot = "both", the minimum cex value when making point sizes is inversely proportional to the prediction standard errors. See \code{\link{par}} for more on cex.  Also see details below. Default is 1.
#' @param SEcex.max if VarPlot = "both", the maximum cex value when making point sizes inversely proportional to the prediction standard errors. See \code{\link{par}} for more on cex.  Also see details below. Default is 3.
#' @param dec.dig the number of decimal places to print in the legend.  Default is 2.
#' @param add Logical value indicating whether the predictions should be added to an existing plot, such as a plot of colored values for observed data. Default is FALSE.
#' @param \dots Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
#' @details The \command{plot.predictions} function creates a map showing color-coded predictions or prediction standard error values. When VarPlot = "Both", predictions values are colored according to breaks.  The size of the points is inversely proportional to the prediction standard errors.
#' @return Maps of stream networks with prediction and prediction standard error values.

#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18.

#' @useDynLib StreamNetgstat

plot.ssn <-
  function(x, VariableName = NULL, color.palette = NULL, nclasses = NULL,
           breaktype = "quantile", brks = NULL, PredPointsID = NULL, add = FALSE,
           addWithLegend = FALSE, lwdLineCol = NULL, lwdLineEx = 1,
           lineCol = "black", ...)
  {
    if(missing(lwdLineEx)) lwdLineEx <- 1
    if(missing(lwdLineCol))
    {
      x@data$lineWidth <- rep(1, nrow(x@data))
      lwdLineCol <- "lineWidth"
    }
    if(is.null(as.list(match.call()[-1])$pch)) {
      plch = 19
    } else plch <- as.list(match.call()[-1])$pch
    if(is.null(as.list(match.call()[-1])$cex)) {
      chex = 1
    } else chex <- as.list(match.call()[-1])$cex
    if(is.null(as.list(match.call()[-1])$col)) {
      colr = "black"
    } else colr <- as.list(match.call()[-1])$col
    par.orig <- par(no.readonly = TRUE)
    if(!is.null(PredPointsID)){
      for(i in 1:length(x@predpoints@ID)){
        if(x@predpoints@ID[i] == PredPointsID){
          if(add == FALSE & addWithLegend == FALSE) {
            plot(x@bbox[1,],x@bbox[2,], type = "n", ...)
            for(j in 1:length(x@lines))
              for(k in 1:length(x@lines[[j]]))
                if(is.null(lwdLineCol))
                  lines((x@lines[[j]]@Lines[[k]]@coords),
                        col = lineCol, ...)
            else
              lines(x@lines[[j]]@Lines[[k]]@coords,
                    lwd = lwdLineEx*x@data[i,lwdLineCol],
                    col = lineCol, ...)
          }
          #}
          if(add == TRUE) {
            par(new = TRUE)
            plot(x@bbox[1,],x@bbox[2,], type = "n", bty = "n",
                 xlab = "", ylab = "",...)
          }
          if(addWithLegend == TRUE) {
            par(new = TRUE)
            layout(matrix(1:2, nrow = 1), widths = c(4,1))
            par(mar = c(5,5,3,0))
            par(mfg = c(1,1))
            plot(x@bbox[1,],x@bbox[2,], type = "n",
                 bty = "n", xlab = "", ylab = "",...)
          }
          points(x@predpoints@SSNPoints[[i]]@point.coords, pch = plch, cex = chex,
                 col = colr)
        }}
      par(par.orig)
    } else
      if(is.null(VariableName)) {
        plot(x@bbox[1,],x@bbox[2,], type = "n", ...)
        for(i in 1:length(x@lines))
          for(j in 1:length(x@lines[[i]]))
            if(is.null(lwdLineCol))
              lines((x@lines[[i]]@Lines[[j]]@coords),
                    col = lineCol, ...)
        else
          lines(x@lines[[i]]@Lines[[j]]@coords,
                lwd = lwdLineEx*x@data[i,lwdLineCol],
                col = lineCol, ...)
        points(x@obspoints@SSNPoints[[1]]@point.coords, pch = plch, cex = chex,
               col = colr)
        par(par.orig)
      } else {
        layout(matrix(1:2, nrow = 1), widths = c(4,1))
        par(mar = c(5,5,3,0))
        plot(x@bbox[1,],x@bbox[2,], type = "n", ...)
        for(i in 1:length(x@lines))
          for(j in 1:length(x@lines[[i]]))
            if(is.null(lwdLineCol))
              lines((x@lines[[i]]@Lines[[j]]@coords),
                    col = lineCol, ...)
        else
          lines(x@lines[[i]]@Lines[[j]]@coords,
                lwd = lwdLineEx*x@data[i,lwdLineCol],
                col = lineCol, ...)
        data <- x@obspoints@SSNPoints[[1]]@point.data
        if(is.null(nclasses)) nclasses <- 10
        lower.breaks <- matrix(0, nrow = nclasses, ncol = 1)
        upper.breaks <- matrix(0, nrow = nclasses, ncol = 1)
        if (!is.factor(data[,VariableName])){
          if(breaktype == "quantile") {
            brks <- quantile(data[,VariableName],
                             probs = (1:(nclasses-1))/nclasses, na.rm = T)
            lower.breaks <- c(min(data[,VariableName], na.rm = T), brks)
            upper.breaks <- c(brks, max(data[,VariableName], na.rm = T))
          }
          if(breaktype == "even") {
            brks <- min(data[,VariableName]) +
              (max(data[,VariableName]) - min(data[,VariableName])) *
              (1:(nclasses-1))/nclasses
            lower.breaks <- c(min(data[,VariableName], na.rm = T), brks)
            upper.breaks <- c(brks, max(data[,VariableName], na.rm = T))
          }
          if(breaktype == "user") {
            if(is.null(brks)) return("Must specify brks if breaktype = user")
            minD <- min(data[,VariableName], na.rm=TRUE)
            maxD <- max(data[,VariableName], na.rm=TRUE)
            brks <- as.vector(unlist(brks))
            if(minD < min(brks)) brks <- c(brks, minD)
            if(maxD > max(brks)) brks <- c(brks, maxD)
            brks <- sort(unique(unlist(brks)))
            nclasses <- length(brks) - 1
            lower.breaks <- brks[1:nclasses]
            upper.breaks <- brks[2:(nclasses+1)]
          } 
        }
        if(length(color.palette) == 0)
          color.palette <- rainbow(nclasses, start = .66, end = .99)
        if (is.factor(data[,VariableName])){
          leglabs = NULL
          for (j in 1:length(levels(data[,VariableName]))){
            indj = which(data[,VariableName] == j)
            points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F],
                   col = color.palette[j], pch = plch, cex = chex)
            leglabs = c(leglabs, paste(VariableName, levels(data[,VariableName])[j]))
          }
          par(mar = c(0,0,0,0))
          plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n",
               xlab = "", ylab ="", bty = "n")
          legend(x = -1, y = 1.1, legend = leglabs, bty = "n",
                 pch = rep(plch, times = length(leglabs)),
                 col = color.palette, cex = .8)
        } else {
          for (j in 1:nclasses){
            jmax <- upper.breaks[j]
            jmin <- lower.breaks[j]
            indj <- data[,VariableName] >= jmin &
              data[,VariableName] <= jmax
            points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F],
                   col = color.palette[j], pch = plch, cex = chex)
          }
          dec.dig <- 2
          left <- as.character(as.numeric(as.integer(lower.breaks*10^dec.dig))/
                                 10^dec.dig)
          rght <- as.character(as.numeric(as.integer(upper.breaks*10^dec.dig))/
                                 10^dec.dig)
          leglabs <- paste(left,"to",rght)
          par(mar = c(0,0,0,0))
          plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n",
               xlab = "", ylab ="", bty = "n")
          legend(x = -1, y = 1.1, legend = leglabs, bty = "n",
                 pch = rep(plch, times = length(leglabs)),
                 col = color.palette, cex = .8)
        }
        par(par.orig)
        return(invisible(data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)))
      }
    
  }

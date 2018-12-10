#' Make prediction on a Spatial Stream Network object
#' 
#' @param ssn a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param varNames a vector of strings, the names of the variables used in the model.
#' @param weightVar a string indicating the name of the variable to compute the spatial weights.
#' @param predpts a string indicating the name of the group of points in which make predictions.
#' @param CorModels a vector of strings, the names of the covariance models.
#' @param useNugget If \code{FALSE} the nugget effect is not included in the model. Default to \code{TRUE}.
#' @param theta vector of the parameters values used for computing the covariance matrix.
#' @param covMat covariance matrix of the observed points.
#' @return A list with the following fields:
#' \item{\code{optTheta}}{ vector of the parameters values of the fitted model. }
#' \item{\code{betaValues}}{ vector of the beta values of the fitted model. }
#' \item{\code{covMatrix}}{ covariance matrix of the fitted model. }
#' \item{\code{predictions}}{ 2-column matrix containing the predicted values and the kriging variance. }
#' @description Given a \link[SSN]{SpatialStreamNetwork-class} object and the covariance parameters values, makes predictions using universal kriging.
#' @details This function calculates prediction values and kriging variance for prediction sites based on the results of a linear model for the \link[SSN]{SpatialStreamNetwork-class} object.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18. DOI: 10.1198/jasa.2009.ap08248. Rejoinder pgs. 22–24.
#' @useDynLib StreamNetgstat
#' @export

do_SSN_kriging = function(ssn, varNames, weightVar, predpts, CorModels, useNugget = TRUE,
                          theta, covMat){

  library(rlist)
  
  # -------------------------------------------------------------
  # Preprocessing of the data
  
  net_num = as.numeric(levels(ssn@network.line.coords[,"NetworkID"]))

  # Create a list for the binaryId tables (one per network)
  bin_tables = list()
  for (k in net_num){
    driver <- RSQLite::SQLite()
    connect.name <- file.path(ssn@path,"binaryID.db")
    connect <- dbConnect(SQLite(), connect.name)
    if (file.exists(file.path(ssn@path, "binaryID.db")) == FALSE){
      stop("binaryID.db is missing from ssn object")
    }
    
    net.name <- paste("net", k, sep = "")
    bin.table <- dbReadTable(connect, net.name)
    dbDisconnect(connect)
    bin.table = bin.table[order(bin.table[,1]),2]
    bin_tables = list.append(bin_tables, bin.table)
  }
  
  # Create a data.frame for the segments attributes
  network_data = ssn@network.line.coords
  indx <- sapply(network_data, is.factor)
  network_data[indx] <- lapply(network_data[indx], function(x) as.numeric(as.character(x)))
  network_data = data.matrix(network_data[order(network_data$NetworkID, network_data$SegmentID),])
  
  # Create a data.frame for the observed points attributes
  obs_points = cbind(ssn@obspoints@SSNPoints[[1]]@network.point.coords,
                     ssn@obspoints@SSNPoints[[1]]@point.coords)
  indx <- sapply(obs_points, is.factor)
  obs_points[indx] <- lapply(obs_points[indx], function(x) as.numeric(as.character(x)))
  obs_points = data.matrix(obs_points[order(obs_points$NetworkID),])
  
  # Create a data.frame for the observed points data
  obs_data = ssn@obspoints@SSNPoints[[1]]@point.data
  indx <- sapply(obs_data, is.factor)
  obs_data[indx] <- lapply(obs_data[indx], function(x) as.numeric(as.character(x)))
  obs_data = data.matrix(obs_data[order(obs_data$netID, obs_data$pid),c(varNames, weightVar)])
  
  # Create a data.frame for the prediction points attributes and data
  pred_points = NULL
  pred_data = NULL
  if (!is.null(predpts)){
    for (p in 1:length(ssn@predpoints@SSNPoints)){
      if (ssn@predpoints@ID[[p]] == predpts){
        tmp = cbind(ssn@predpoints@SSNPoints[[p]]@network.point.coords,
                    ssn@predpoints@SSNPoints[[p]]@point.coords)
        pred_points = rbind(pred_points, tmp)
        pred_data = rbind(pred_data, ssn@predpoints@SSNPoints[[p]]@point.data)
      }
    }
    indx <- sapply(pred_points, is.factor)
    pred_points[indx] <- lapply(pred_points[indx], function(x) as.numeric(as.character(x)))
    pred_points = data.matrix(pred_points[order(pred_points$NetworkID),])
    
    indx <- sapply(pred_data, is.factor)
    pred_data[indx] <- lapply(pred_data[indx], function(x) as.numeric(as.character(x)))
    pred_data = data.matrix(pred_data[order(pred_data$netID, pred_data$pid),c(varNames[-1], weightVar)])
  }
  
  
  
  # Pass bin_table, network_data, obs_points, pred_points, obs_data, pred_data, c(varNames, weightVar), CorModels
  # to the C++ function
  
  result = .Call("doSSNKriging", net_num, bin_tables, network_data,
                 obs_points, pred_points, obs_data, pred_data, c(varNames, weightVar), CorModels, useNugget,
                 theta, covMat)
  
  return (result)
  
  
  # colnames(bin_table) = list of vectors (one per network) of binaryIDs
  # colnames(network_data) = c("NetworkID", SegmentID", "DistanceUpstream")
  # colnames(obs_points) = c("NetworkID", "SegmentID", "DistanceUpstream", "coords.x1", "coords.x2")
  # colnames(pred_points) = c("NetworkID", "SegmentID", "DistanceUpstream", "coords.x1", "coords.x2")
  # colnames(obs_data) = c(varNames, weightVar)
  # colnames(pred_data) = c(varNames[-1], weightVar)
}

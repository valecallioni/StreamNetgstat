#' Make prediction on a Spatial Stream Network object
#' 
#' @param ssn a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param varNames a vector of strings, the names of the variables used in the model.
#' @param weightVar a string indicating the name of the variable to compute the spatial weights.
#' @param predpts a string indicating the name of the group of points in which make predictions.
#' @param CorModels a vector of strings, the names of the covariance models. The individual models should be of different "types" tail-up, tail-down, Euclidean. The tailup models include: "Exponential.tailup", "LinearSill.tailup", "Spherical.tailup", 
#'
#' "Mariah.tailup" "Epanech.tailup"; tail-down models include: "Exponential.taildown", "LinearSill.taildown", "Spherical.taildown", "Mariah.taildown"; Euclidean distance models include: "Spherical.Euclid", "Gaussian.Euclid", "Exponential.Euclid", 
#'
#' "Cauchy.Euclid". The first 4 tailup and taildown models are described in Ver Hoef and Peterson (2010), and the 4 Euclidean distance models are standard spatial covariance models.
#' @param useNugget If \code{FALSE} the nugget effect is not included in the model. Default to \code{TRUE}.
#' @param singleNet an interger, indicating the network ID that is to be analysed. Default to \code{NULL}, so that the analysis is carried on the entire dataset.
#' @param theta vector of the parameters values used for computing the covariance matrix.
#' @param covMat covariance matrix of the observed points.
#' @param matrices a vector of matrices, containing the flow-connection binary matrix, the hydrological distance matrix and, not necessarily, the Euclidean distance matrix, returned by the function \link[StreamNetgstat]{get_plots}. These matrices consider the relationships between observed points.
#' @return A \code{ssn.object}{ (updated \link[SSN]{SpatialStreamNetwork-class} object), with predicted values included }
#' @description Given a \link[SSN]{SpatialStreamNetwork-class} object and the covariance parameters values, makes predictions using universal kriging.
#' @details This function calculates prediction values and kriging variance for prediction sites based on the results of a linear model for the \link[SSN]{SpatialStreamNetwork-class} object.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18.
#' @useDynLib StreamNetgstat
#' @export

do_SSN_kriging = function(ssn, varNames, weightVar, predpts, CorModels, useNugget = TRUE,
                          theta, covMat, singleNet = NULL, matrices = NULL){

  # -------------------------------------------------------------
  # Preprocessing of the data
  
  net_num = as.numeric(levels(ssn@network.line.coords[,"NetworkID"]))
  
  if (!is.null(singleNet)){
    if (!is.element(singleNet, net_num))
      stop("Network selected not valid")
    else
      net_num = singleNet
  }

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
  if (!is.null(singleNet)) network_data = network_data[which(network_data$NetworkID == singleNet),]
  indx <- sapply(network_data, is.factor)
  network_data[indx] <- lapply(network_data[indx], function(x) as.numeric(as.character(x)))
  network_data = network_data[order(network_data$NetworkID, network_data$SegmentID),]
  
  # Create a data.frame for the observed points attributes
  obs_points = cbind(ssn@obspoints@SSNPoints[[1]]@network.point.coords,
                     ssn@obspoints@SSNPoints[[1]]@point.coords)
  if (!is.null(singleNet)) obs_points = obs_points[which(obs_points$NetworkID == singleNet),]
  indx <- sapply(obs_points, is.factor)
  obs_points[indx] <- lapply(obs_points[indx], function(x) as.numeric(as.character(x)))

  # Create a data.frame for the observed points data
  obs_data = ssn@obspoints@SSNPoints[[1]]@point.data
  if (!is.null(singleNet)) obs_data = obs_data[which(obs_data$netID == singleNet),]
  indx <- sapply(obs_data, is.factor)
  obs_data[indx] <- lapply(obs_data[indx], function(x) as.numeric(as.character(x)))
  obs_points = obs_points[order(obs_points$NetworkID, obs_data$pid),]
  obs_data = obs_data[order(obs_data$netID, obs_data$pid),c(varNames, weightVar)]
  
  # Create a data.frame for the prediction points attributes and data
  pred_points = NULL
  pred_data = NULL
  if (!is.null(predpts)){
    for (p in 1:length(ssn@predpoints@SSNPoints)){
      if (ssn@predpoints@ID[[p]] == predpts){
        id_p = p
        tmp = cbind(ssn@predpoints@SSNPoints[[p]]@network.point.coords,
                    ssn@predpoints@SSNPoints[[p]]@point.coords)
        pred_points = rbind(pred_points, tmp)
        pred_data = rbind(pred_data, ssn@predpoints@SSNPoints[[p]]@point.data)
      }
    }
    if (!is.null(singleNet)) pred_points = pred_points[which(pred_points$NetworkID == singleNet),]
    indx <- sapply(pred_points, is.factor)
    pred_points[indx] <- lapply(pred_points[indx], function(x) as.numeric(as.character(x)))

    if (!is.null(singleNet)) pred_data = pred_data[which(pred_data$netID == singleNet),]
    indx <- sapply(pred_data, is.factor)
    pred_data[indx] <- lapply(pred_data[indx], function(x) as.numeric(as.character(x)))
    # locID_cpp = pred_data$locID
    order_data = order(pred_data$netID, pred_data$pid)
    pred_points = pred_points[order(pred_points$NetworkID, pred_data$pid),]
    pred_data = pred_data[order_data,c(varNames[-1], weightVar)]
  }
  
  if (!is.null(singleNet)) {
    result = .Call("doSSNKriging_SingleNet", bin.table, data.matrix(network_data),
                   data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data),
                   c(varNames, weightVar), CorModels, useNugget, theta, covMat, matrices)
  } else { 
    result = .Call("doSSNKriging_MultipleNets", net_num, bin_tables, data.matrix(network_data),
                   data.matrix(obs_points), data.matrix(pred_points), data.matrix(obs_data), data.matrix(pred_data), 
                   c(varNames, weightVar), CorModels, useNugget, theta, covMat, matrices)
  }

  if (!is.null(singleNet)) {
    id_rows = which(ssn@predpoints@SSNPoints[[id_p]]@network.point.coords$NetworkID == singleNet)
  } else {
    id_rows = 1:dim(ssn@predpoints@SSNPoints[[id_p]]@point.data)[1]
  }
  col_names = colnames(ssn@predpoints@SSNPoints[[id_p]]@point.data)
  ssn@predpoints@SSNPoints[[id_p]]@point.data = cbind.data.frame(ssn@predpoints@SSNPoints[[id_p]]@point.data,
                                                                 rep(NA,dim(ssn@predpoints@SSNPoints[[id_p]]@point.data)[1]),
                                                                 rep(NA,dim(ssn@predpoints@SSNPoints[[id_p]]@point.data)[1]))
  colnames(ssn@predpoints@SSNPoints[[id_p]]@point.data) = c(col_names, paste(varNames[1], "pred", sep="_"), paste(varNames[1], "SE", sep="_"))
  ssn@predpoints@SSNPoints[[id_p]]@point.data[as.logical(id_rows), paste(varNames[1], "pred", sep="_")] = result$predictions[match(id_rows, order_data),1]
  ssn@predpoints@SSNPoints[[id_p]]@point.data[as.logical(id_rows), paste(varNames[1], "se", sep="_")] = result$predictions[match(id_rows, order_data),2]
    
  return (ssn.object = ssn)
}

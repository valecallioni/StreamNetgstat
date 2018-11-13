#' do kriging
#' @useDynLib StreamNetgstat
#' @export

do_SSN_kriging = function(ssn, varNames, weightVar, predpts, CorModels, useNugget = TRUE,
                          theta, covMat){

  library(rlist)
  
  # -------------------------------------------------------------
  # Preprocessing of the data
  
  net_num = as.numeric(levels(ssn@network.line.coords[,1]))

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
    pred_points = data.matrix(pred_points[order(pred_points$NetworkID,),])
    
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

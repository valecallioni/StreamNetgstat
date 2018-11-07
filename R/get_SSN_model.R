#' get model
#' @useDynLib StreamNetgstat
#' @export

#get_SSN_model = function(ssn, net.num, varNames, weightVar){
 
get_SSN_model = function(varNames, weightVar, CorModels,
                         bin_tables, network_data, obs_points, pred_points, obs_data, pred_data){
  
  # # Check to see whether distance folder exists...
  # if (!file.exists(file.path(ssn@path, "distance"))) {
  #   dir.create(file.path(ssn@path, "distance"))
  # }
  # 
  # ## And then whether an observation folder exists
  # if (!file.exists(file.path(ssn@path, "distance", "obs"))) {
  #   dir.create(file.path(ssn@path, "distance", "obs"))
  # }
  # 
  # ## And then whether prediction folder exists
  # # if(length(ssn@predpoints@ID) > 0) {
  # #   pred.num <- which(ssn@predpoints@ID == predpts)
  # #   if (length(pred.num) > 1) {
  # #     stop("SSN contains more than one copy of ", predpts)
  # #   }
  # #   if (length(pred.num) == 0){
  # #     stop(predpts, " does not exist in SSN")
  # #   }
  # # }
  # 
  # net_num = levels(ssn@network.line.coords[,1])
  # 
  # bin_tables = list()
  # for (k in net_num){
  #   driver <- RSQLite::SQLite()
  #   connect.name <- file.path(ssn@path,"binaryID.db")
  #   connect <- dbConnect(SQLite(), connect.name)
  #   if (file.exists(file.path(ssn@path, "binaryID.db")) == FALSE){
  #     stop("binaryID.db is missing from ssn object")
  #   }
  # 
  #   net.name <- paste("net", k, sep = "")
  #   bin.table <- dbReadTable(connect, net.name)
  #   dbDisconnect(connect)
  #   bin.table = bin.table[order(bin.table[,1]),2]
  #   bin_tables[[k]] = bin.table
  # }
  # 
  # network_data = ssn@network.line.coords
  # indx <- sapply(network_data, is.factor)
  # network_data[indx] <- lapply(network_data[indx], function(x) as.numeric(as.character(x)))
  # ord = order(network_data$NetworkID)
  # network_data = data.matrix(network_data[ord,])
  # 
  # obs_points = cbind(ssn@obspoints@SSNPoints[[1]]@network.point.coords,
  #                        ssn@obspoints@SSNPoints[[1]]@point.coords)
  # indx <- sapply(obs_points, is.factor)
  # obs_points[indx] <- lapply(obs_points[indx], function(x) as.numeric(as.character(x)))
  # ord = order(obs_points$NetworkID)
  # obs_points = data.matrix(obs_points[ord,])
  # 
  # obs_data = ssn@obspoints@SSNPoints[[1]]@point.data
  # indx <- sapply(obs_data, is.factor)
  # obs_data[indx] <- lapply(obs_data[indx], function(x) as.numeric(as.character(x)))
  # ord = order(obs_data$netID)
  # obs_data = data.matrix(obs_data[ord, c(varNames, weightVar)])
  # 
  # pred_points = NULL
  # pred_data = NULL
  # for (p in 1:length(ssn@predpoints@SSNPoints)){
  #   if (ssn@predpoints@ID[[p]] == "pred1km"){
  #     tmp = cbind(ssn@predpoints@SSNPoints[[p]]@network.point.coords,
  #                 ssn@predpoints@SSNPoints[[p]]@point.coords)
  #     pred_points = rbind(pred_points, tmp)
  #     pred_data = rbind(pred_data, ssn@predpoints@SSNPoints[[p]]@point.data)
  #   }
  # }
  # indx <- sapply(pred_points, is.factor)
  # pred_points[indx] <- lapply(pred_points[indx], function(x) as.numeric(as.character(x)))
  # ord = order(pred_points$NetworkID)
  # pred_points = data.matrix(pred_points[ord, ])
  # 
  # indx <- sapply(pred_data, is.factor)
  # pred_data[indx] <- lapply(pred_data[indx], function(x) as.numeric(as.character(x)))
  # ord = order(pred_data$netID)
  # pred_data = data.matrix(pred_data[ord, c(varNames[-1], weightVar)])

  
  # Pass network_data, bin.Table, obs_points, pred_points, obs_data, pred_data, c(varNames, weightVar)
  # to the c++ function
  
  result = .Call("getSSNM", 2, bin_tables, network_data, obs_points, pred_points, obs_data, pred_data, c(varNames, weightVar), CorModels)
  return (result)
  
  
  # colnames(bin.table) = c("rid", "binaryID")
  # colnames(network_data) = c("SegmentID", "DistanceUpstream")
  # colnames(obs_points) = c("SegmentID", "DistanceUpstream", "coords.x1", "coords.x2")
  # colnames(pred_points) = c("SegmentID", "DistanceUpstream", "coords.x1", "coords.x2")
  # colnames(obs_data) = c(varNames, weightVar)
}

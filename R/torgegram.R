#' Compute the Torgegram for a Spatial Stream Network object 
#'
#' @param ssn a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param ResponseName the name of the response variable.
#' @param singleNet an interger, indicating the network ID that is to be analysed. Default to \code{NULL}, so that the analysis is carried on the entire dataset.
#' @param maxlag the maximum lag distance to consider when binning pairs of locations by the hydrological distance that separates them. The default is the median distance between all pairs of locations.
#' @param nlag the number of lag bins to create for computing the Torgegram. The hydrological distance between endpoints that define a bin will have equal lengths for all bins. The bin sizes are then determined from the minimum lag in the data, and the specification of maxlag_Torg.
#' @param inc the bin hydrological distance between endpoints. It is possible to specify the bin distance rather than nlag_Torg. In this case, the number of bins is determined by the bin distance and the hydrological distance between the mininum and maximum (maxlag_Torg) lag in the data. 
#' @param nlagcutoff the minimum number of pairs needed to estimate the semivariance for a bin in the Torgegram computation. If the sample sizes is less than this value, the semivariance for the bin is not calculated.
#' @param EmpVarMeth method for computing semivariances. The default is "MethMoment", the classical method of moments, which is just the average difference-squared within bin classes. "Covariance" computes covariance rather than semivariance, but may be more biased because it subtracts off the simple mean of the response variable. "RobustMedian" and "RobustMean" are robust estimators proposed by Cressie and Hawkins (1980). If v is a vector of all pairwise square-roots of absolute differences within a bin class, then RobustMedian computes median(v)^4/.457. "RobustMean" computes mean(v)^4/(.457 + .494/length(v)).
#' @return A list.
#' @description \code{torgegram} computes the empirical semivariogram from the data based on hydrological distance. The results are presented separately for flow-connected and flow-unconnected sites.
#' @details The \code{torgegram} function creates a list of hydrological distances and empirical semivariogram values, along with number of pairs of points in each bin, for both flow-connected and flow-unconnected sites. Flow-connected locations lie on the same stream network (share a common downstream junction) and water flows from one location to the other. Flow-unconnected locations also lie on the same stream network, but do not share flow. The output is of class \link[SSN]{Torgegram}.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18.
#' @useDynLib StreamNetgstat

torgegram = function(ssn, ResponseName, singleNet = NULL, maxlag = NULL, nlag = 6, 
                     inc = 0, nlagcutoff = 15, EmpVarMeth = "MethMoment")
{
  
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
  obs_data = obs_data[order(obs_data$netID, obs_data$pid),ResponseName]
  
  if (!is.null(singleNet)) {
    matrices = .Call("createHydroDistanceMatrices_SingleNet", bin.table, data.matrix(network_data), data.matrix(obs_points))
  } else {
    matrices = .Call("createHydroDistanceMatrices_MultipleNets", net_num, bin_tables, data.matrix(network_data), data.matrix(obs_points))
  }
  
  Dif2s = NULL
  sqrtDifs = NULL
  Dists = NULL
  FCons = NULL
  Covp = NULL
  mnz = mean(obs_data, na.rm = TRUE)
  nsofar = 0
  nObs = table(obs_points[,"NetworkID"])
  for (i in 1:length(nObs)){
    Ds = t(matrices$distHydro[(nsofar+1):(nsofar + nObs[i]), (nsofar+1):(nsofar + nObs[i])]) + matrices$distHydro[(nsofar+1):(nsofar + nObs[i]), (nsofar+1):(nsofar + nObs[i])]
    FCs = matrices$flowMat[(nsofar+1):(nsofar + nObs[i]), (nsofar+1):(nsofar + nObs[i])] + diag(1, nObs[i])
    zs = obs_data[(nsofar+1):(nsofar + nObs[i])]
    nsofar = nsofar + nObs[i]
    ind <- !is.na(zs)
    ni <- sum(ind)
    rs = zs - mnz
    CP = rs%o%rs
    diff = abs(zs%o%rep(1, times = ni) - 
                  rep(1, times = ni)%o%zs)
    diff2 = diff^2
    sqrtd = sqrt(diff)
    Dif2s = c(Dif2s, as.vector(diff2[col(diff2) < row(diff2)]))
    sqrtDifs = c(sqrtDifs, as.vector(sqrtd[col(sqrtd) < row(sqrtd)]))
    Dists = c(Dists, as.vector(Ds[col(Ds) < row(Ds)]))
    FCons = c(FCons, as.vector(FCs[col(FCs) < row(FCs)]))
    Covp = c(Covp, as.vector(CP[col(CP) < row(CP)]))
  }
  
  if(is.null(maxlag)) maxlag = max(Dists)
  indmax = Dists <= maxlag
  Dists = Dists[indmax]
  Dif2s = Dif2s[indmax]
  sqrtDifs = sqrtDifs[indmax]
  FCons = FCons[indmax]
  Covp = Covp[indmax]
  
  if(inc <= 0) inc = maxlag/nlag
  
  store.results = matrix(data = NA, ncol = 6, nrow = nlag,
                          dimnames = list(1:nlag, c("distance.connect", "gam.connect", "np.connect",
                                                    "distance.unconnect", "gam.unconnect","np.unconnect")))
  for ( i in 1:nlag){
    ind1 = (Dists >= (i-1)*inc) & (Dists < i*inc) & (FCons == 1)
    ind2 = (Dists >= (i-1)*inc) & (Dists < i*inc) & (FCons == 0)
    nclass1 = sum(ind1)
    nclass2 = sum(ind2)
    if(EmpVarMeth == "MethMoment") {
      cv1 = mean(Dif2s[ind1])
      cv2 = mean(Dif2s[ind2])
    }
    if(EmpVarMeth == "RobustMean") {
      cv1 = ((mean(sqrtDifs[ind1]))^4)/(.457 + .494/sum(ind1))
      cv2 = ((mean(sqrtDifs[ind2]))^4)/(.457 + .494/sum(ind2))
    }
    if(EmpVarMeth == "RobustMedian") {
      cv1 = (median(sqrtDifs[ind1]))^4/.457
      cv2 = (median(sqrtDifs[ind2]))^4/.457
    }
    if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") {
      cv1 = mean(Covp[ind1])
      cv2 = mean(Covp[ind2])
    }
    mean.dis1 = mean(Dists[ind1])
    mean.dis2 = mean(Dists[ind2])
    if(nclass1 > 0 | nclass2 > 0) store.results[i,] =
      c(mean.dis1, cv1, nclass1, mean.dis2, cv2, nclass2)
  }
  store.results[,"gam.connect"] = store.results[,"gam.connect"]/2
  store.results[,"gam.unconnect"] = store.results[,"gam.unconnect"]/2
  indc = store.results[,"np.connect"] >= nlagcutoff 
  store.results[!indc,1:3] = NA	
  indu = store.results[,"np.unconnect"] >= nlagcutoff
  store.results[!indu,4:6] = NA	
  
  ev = as.data.frame(store.results)
  
  ev[is.nan(ev[,"distance.unconnect"]),"distance.unconnect"] = NA
  ev[is.nan(ev[,"gam.unconnect"]),"gam.unconnect"] = NA
  ev[is.nan(ev[,"distance.connect"]),"distance.connect"] = NA
  ev[is.nan(ev[,"gam.connect"]),"gam.connect"] = NA
  
  ev = as.data.frame(ev)
  ev = as.list(ev)
  for(i in 1:length(ev))
    ev[[i]] = ev[[i]][!is.na(ev[[i]])]
  ev$call = as.list(match.call()[-1])
  class(ev) = "Torgegram"
  return (list(Torg = ev, 
               distMatrices = matrices))
}




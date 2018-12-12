#' Compute the Empirical Semivariogram for a Spatial Stream Network object 
#'
#' @param ssn a \link[SSN]{SpatialStreamNetwork-class} object.
#' @param ResponseName the name of the response variable.
#' @param maxlag_Torg the maximum lag distance to consider when binning pairs of locations by the hydrologic distance that separates them. The default is the median distance between all pairs of locations.
#' @param nlag_Torg the number of lag bins to create for computing the Torgegram. The hydrologic distance between endpoints that define a bin will have equal lengths for all bins. The bin sizes are then determined from the minimum lag in the data, and the specification of maxlag_Torg.
#' @param inc_Torg the bin hydrologic distance between endpoints. It is possible to specify the bin distance rather than nlag_Torg. In this case, the number of bins is determined by the bin distance and the hydrologic distance between the mininum and maximum (maxlag_Torg) lag in the data. 
#' @param nlagcutoff_Torg the minimum number of pairs needed to estimate the semivariance for a bin in the Torgegram computation. If the sample sizes is less than this value, the semivariance for the bin is not calculated.
#' @param maxlag_EmpVar the maximum lag distance to consider when binning pairs of locations by the Euclidean distance that separates them. If the specified maxlag is larger than the maximum Euclidean distance among pairs of points, then maxlag_EmpVar is set to the maximum distance among pairs. If inc_EmpVar is greater than 0, then maxlag_EmpVar is disregarded.
#' @param nlag_EmpVar the number of lag bins to create for computing the Empirical Semivariogram, by direction if directions are specified. The distance between endpoints that define a bin will have equal lengths for all bins. The bin sizes are then determined from the minimum lag in the data, and the specification of maxlag_EmpVar. 
#' @param inc_EmpVar the Euclidean distance increment for each bin class. Default is 0, in which case maxlag and nclasses determine the distance increments. 
#' @param nlagcutoff_EmpVar the minimum number of pairs needed to estimate the semivariance for a bin in the Empirical Semivariogram. If the sample size is less than this value, the semivariance for the bin is not calculated.
#' @param directions directions in degrees clockwise from north that allow lag binning to be directional. Default is \code{c(0, 45, 90, 135)}. Values should be between 0 and 180, as there is radial symmetry in orientation between two points.
#' @param tolerance the angle on either side of the directions to determine if a pair of points falls in that direction class. Note, a pair of points may be in more than one lag bin if tolerances for different directions overlap.
#' @param EmpVarMeth method for computing semivariances. The default is "MethMoment", the classical method of moments, which is just the average difference-squared within bin classes. "Covariance" computes covariance rather than semivariance, but may be more biased because it subtracts off the simple mean of the response variable. "RobustMedian" and "RobustMean" are robust estimators proposed by Cressie and Hawkins (1980). If v is a vector of all pairwise square-roots of absolute differences within a bin class, then RobustMedian computes median(v)^4/.457. "RobustMean" computes mean(v)^4/(.457 + .494/length(v)).
#' @return A list with the following fields:
#' \item{\code{EmpSemiVar}}{ object of class \link[SSN]{EmpiricalSemivariogram}. }
#' \item{\code{Torg}}{ object of class \link[SSN]{Torgegram}. }
#' \item{\code{distMatrices}}{ list of distance-matrices. }
#' @description Given a \link[SSN]{SpatialStreamNetwork-class} object, computes the Torgegram (empirical semivariogram from the data based on hydrologic distance) and the empirical semivariogram from the data based on the Euclidean distance.
#' @details Given a \link[SSN]{SpatialStreamNetwork-class} object, the function creates a list of hydrologic distances and Euclidean distances and empirical semivariogram values for both, along with number of pairs of points in each bin, for flow-connected and flow-unconnected sites when considering the hydrologic distance. Flow-connected locations lie on the same stream network (share a common downstream junction) and water flows from one location to the other. Flow-unconnected locations also lie on the same stream network, but do not share flow. The output is a list, containing an object of class \link[SSN]{Torgegram}, an object of class \link[SSN]{EmpiricalSemivariogram} and a list containing the distance-matrices.
#' @references 
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach to geostatistical modeling in stream networks. Ecology 91(3), 644–651.
#' 
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial statistical models of stream networks (with discussion). Journal of the American Statistical Association 105, 6–18. DOI: 10.1198/jasa.2009.ap08248. Rejoinder pgs. 22–24.
#' @useDynLib StreamNetgstat
#' @export

emp_semivario_and_torg = function(ssn, ResponseName, 
                                  maxlag_Torg = NULL, nlag_Torg = 6, inc_Torg = 0, nlagcutoff_Torg = 15, 
                                  maxlag_EmpVar = 1e32, nlag_EmpVar = 20, inc_EmpVar = 0, nlagcutoff_EmpVar = 1, 
                                  directions = c(0,45,90,135), tolerance = 22.5, EmpVarMeth = "MethMoment")
{
  
  # -------------------------------------------------------------
  # Preprocessing of the data
  
  net_num = as.numeric(levels(ssn@network.line.coords[,"NetworkID"]))
  
  # Create a list for the binaryId tables (one per network)
  bin_tables = list()
  for (k in net_num){
    driver = RSQLite::SQLite()
    connect.name = file.path(ssn@path,"binaryID.db")
    connect = dbConnect(SQLite(), connect.name)
    if (file.exists(file.path(ssn@path, "binaryID.db")) == FALSE){
      stop("binaryID.db is missing from ssn object")
    }
    
    net.name = paste("net", k, sep = "")
    bin.table = dbReadTable(connect, net.name)
    dbDisconnect(connect)
    bin.table = bin.table[order(bin.table[,1]),2]
    bin_tables = list.append(bin_tables, bin.table)
  }
  
  # Create a data.frame for the segments attributes
  network_data = ssn@network.line.coords
  indx = sapply(network_data, is.factor)
  network_data[indx] = lapply(network_data[indx], function(x) as.numeric(as.character(x)))
  network_data = data.matrix(network_data[order(network_data$NetworkID, network_data$SegmentID),])
  
  # Create a data.frame for the observed points attributes
  obs_points = cbind(ssn@obspoints@SSNPoints[[1]]@network.point.coords,
                     ssn@obspoints@SSNPoints[[1]]@point.coords)
  indx = sapply(obs_points, is.factor)
  obs_points[indx] = lapply(obs_points[indx], function(x) as.numeric(as.character(x)))
  obs_points = data.matrix(obs_points[order(obs_points$NetworkID),])
  
  # Create a data.frame for the observed points data
  obs_data = ssn@obspoints@SSNPoints[[1]]@point.data
  indx = sapply(obs_data, is.factor)
  obs_data[indx] = lapply(obs_data[indx], function(x) as.numeric(as.character(x)))
  obs_data = data.matrix(obs_data[order(obs_data$netID, obs_data$pid),ResponseName])
  
  
  matrices = .Call("createDistanceMatrices", net_num, bin_tables, network_data, obs_points)

  #----------------------
  # TORGEGRAM
  
  Dif2s = NULL
  sqrtDifs = NULL
  Dists = NULL
  FCons = NULL
  Covp = NULL
  mnz = mean(obs_data, na.rm = TRUE)
  nsofar = 0
  nObs = table(ssn@obspoints@SSNPoints[[1]]@network.point.coords[,"NetworkID"])
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
  
  if(is.null(maxlag_Torg)) maxlag_Torg = max(Dists)
  indmax = Dists <= maxlag_Torg
  Dists = Dists[indmax]
  Dif2s = Dif2s[indmax]
  sqrtDifs = sqrtDifs[indmax]
  FCons = FCons[indmax]
  Covp = Covp[indmax]
  
  if(inc_Torg <= 0) inc_Torg = maxlag_Torg/nlag_Torg
  
  store.results = matrix(data = NA, ncol = 6, nrow = nlag_Torg,
                         dimnames = list(1:nlag_Torg, c("distance.connect", "gam.connect", "np.connect",
                                                   "distance.unconnect", "gam.unconnect","np.unconnect")))
  for ( i in 1:nlag_Torg){
    ind1 = (Dists >= (i-1)*inc_Torg) & (Dists < i*inc_Torg) & (FCons == 1)
    ind2 = (Dists >= (i-1)*inc_Torg) & (Dists < i*inc_Torg) & (FCons == 0)
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
  indc = store.results[,"np.connect"] >= nlagcutoff_Torg 
  store.results[!indc,1:3] = NA	
  indu = store.results[,"np.unconnect"] >= nlagcutoff_Torg
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
  
  
  
  # -------------------------
  # EMPIRICAL SEMIVARIOGRAM
  
  x = obs_points[,4]
  y = obs_points[,5]
  nObs = dim(obs_points)[1]
  distance = matrices$distGeo
  difx <- -(matrix(y,nrow=nObs,ncol=1) %*%
              matrix(rep(1,times=nObs),nrow=1,ncol=nObs) -
              matrix(rep(1,times=nObs),nrow=nObs,ncol=1) %*%
              matrix(y,nrow=1,ncol=nObs))
  signind <- -(matrix(x,nrow=nObs,ncol=1) %*%
                 matrix(rep(1,times=nObs),nrow=1,ncol=nObs) -
                 matrix(rep(1,times=nObs),nrow=nObs,ncol=1) %*%
                 matrix(x,nrow=1,ncol=nObs)) < 0
  distance <- distance*1.0000000001
  theta.deg <- acos(difx/distance)*180/pi
  # matrix of degrees clockwise from north between locations
  theta.deg[signind] <- 360-theta.deg[signind]
  diff2 <- ( matrix(obs_data,nrow=nObs,ncol=1) %*%
               matrix(rep(1,times=nObs),nrow=1,ncol=nObs) -
               matrix(rep(1,times=nObs),nrow=nObs,ncol=1) %*%
               matrix(obs_data,nrow=1,ncol=nObs) )^2
  sqrtdiff <- sqrt(abs( matrix(obs_data,nrow=nObs,ncol=1) %*%
                          matrix(rep(1,times=nObs),nrow=1,ncol=nObs) -
                          matrix(rep(1,times=nObs),nrow=nObs,ncol=1) %*%
                          matrix(obs_data,nrow=1,ncol=nObs) ) )
  if(EmpVarMeth == "CovMean") temp4cov <- obs_data - mean(obs_data)
  else temp4cov <- obs_data
  covprod <- (matrix(temp4cov,nrow=nObs,ncol=1) %*%
                matrix(rep(1,times=nObs),nrow=1,ncol=nObs)) *
    (matrix(rep(1,times=nObs),nrow=nObs,ncol=1) %*%
       matrix(temp4cov,ncol=nObs,nrow=1))
  # convert to vectors
  distance <- matrix(distance, ncol = 1)
  theta.deg <- matrix(theta.deg, ncol = 1)
  diff2 <- matrix(diff2, ncol = 1)
  sqrtdiff <- matrix(sqrtdiff, ncol = 1)
  covprod <- matrix(covprod, ncol = 1)
  # trim off values greater than maxlag
  indmax <- distance <= maxlag_EmpVar
  distance <- distance[indmax,]
  theta.deg <- theta.deg[indmax,]
  diff2 <- diff2[indmax,]
  sqrtdiff <- sqrtdiff[indmax,]
  covprod <- covprod[indmax,]
  
  maxd<-max(distance)
  if( inc_EmpVar <= 0) inc_EmpVar <- maxd/nlag_EmpVar
  ind <- distance==0
  ndir <- length(directions)
  store.results <- matrix(data = NA, ncol = 6,
                          dimnames = list(NULL, c("distance", "gamma", "np", "azimuth", "hx", "hy")))
  for (j in 1:ndir) {
    for ( i in 1:nlag_EmpVar){
      if( (directions[j]-tolerance)<0 && (directions[j]+tolerance)>0 )
        ind1 <- theta.deg >= 360+directions[j]-tolerance |
          theta.deg < directions[j]+tolerance
      else if( (directions[j]+tolerance)>360 && (directions[j]-tolerance)<360 )
        ind1 <- theta.deg < directions[j]+tolerance-360 |
          theta.deg >= directions[j]-tolerance
      else
        ind1 <- theta.deg >= directions[j]-tolerance &
          theta.deg < directions[j]+tolerance
      ind<-distance>(i-1)*inc_EmpVar & distance<=i*inc_EmpVar &
        !is.na(theta.deg) & ind1
      nclass <- sum(ind)
      if(EmpVarMeth == "MethMoment") cv <- mean(diff2[ind])
      if(EmpVarMeth == "RobustMean") cv <- ((mean(sqrtdiff[ind]))^4)/(.457 + .494/sum(ind))
      if(EmpVarMeth == "RobustMedian") cv <- (median(sqrtdiff[ind]))^4/.457
      if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") cv <- mean(covprod[ind])
      mean.dis <- mean(distance[ind])
      if(nclass > 0) store.results<-rbind(store.results,
                                          c(mean.dis,cv,nclass,directions[j],0,0))
    }
  }
  store.results[,"hx"]<-store.results[,"distance"]*sin(store.results[,"azimuth"]*pi/180)
  store.results[,"hy"]<-store.results[,"distance"]*cos(store.results[,"azimuth"]*pi/180)
  store.results[,"gamma"]<-store.results[,"gamma"]/2
  ind <- store.results[,"np"] >= nlagcutoff_EmpVar
  store.results <- store.results[ind,]
  ind <- !is.na(store.results[,"hx"])
  store.results <- store.results[ind,]
  store.results <- as.data.frame(store.results)
  class(store.results) <- "EmpiricalSemivariogramSSN"
  
  
  
  
  
  return  (list(EmpSemiVar = store.results, 
                Torg = ev,
                distMatrices = matrices))
}
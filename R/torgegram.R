Torgegram = function(ssn, ResponseName,
                      maxlag = NULL, nlag = 6, inc = 0, nlagcutoff = 15, 
                      EmpVarMeth = "MethMoment")
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
  obs_data = data.matrix(obs_data[order(obs_data$netID, obs_data$pid),c(varNames, weightVar)])
  
  
  matrices = .Call("createHydroDistanceMatrices", net_num, bin_tables, network_data, obs_points)
  
  maxlag = NULL
  nlag = 6
  inc = 0
  nlagcutoff = 15
  
  Dif2s = NULL
  sqrtDifs = NULL
  Dists = NULL
  FCons = NULL
  Covp = NULL
  mnz = mean(data[,ResponseName], na.rm = TRUE)
  nsofar = 0
  nObs = table(ssn@obspoints@SSNPoints[[1]]@network.point.coords[,"NetworkID"])
  for (i in net_num){
    Ds = t(matrices$distHydro[(nsofar+1):nObs[i], (nsofar+1):nObs[i]])
    FCs = matrices$flowMat[(nsofar+1):nObs[i], (nsofar+1):nObs[i]] + diag(1, nObs[i])
    zs = obs_data[(nsofar+1):nObs[i],1]
    nsofar = nsofar + nObs[i]
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




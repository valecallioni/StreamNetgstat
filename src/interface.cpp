#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <list>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Dataframe.hpp"
#include "Network.hpp"
#include "Optimizer.hpp"
#include "Kriging.hpp"

/*! \file
* Functions to create distance matrices, fit a spatial linear model and perform Universal kriging on a spatial stream network
* object, when the data set has multiple networks or a single one.
*/

extern "C"{

// ---------------------------------------------------------------------------------------------------------------------------
// MULTIPLE NETWORKS FUNCTIONS

/**
* Functions to create distance matrices, fit a spatial linear model and perform Universal kriging on a spatial stream network
* object, which has multiple networks.
*/

// COMPUTE HYDROLOGICAL DISTANCES
/**
* Compute the hydrological distances between observed points, as well as the flow-connection/unconnection binary matrix
* @param net_num vector containing the networkID per each network of the data set
* @param bin_tables list of vectors of strings, one vector per each network, containing the binaryIDs of the stream segments of
* each network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @return A list with the following fields:
*   - 'flowMat' flow-connection/unconnection binary matrix
*   - 'distHydro' hydrological distance matrix
*/
RcppExport SEXP createHydroDistanceMatrices_MultipleNets (SEXP net_num, SEXP bin_tables, SEXP network_data,
  SEXP obs_points){

  BEGIN_RCPP

  // Stream segments storage
  std::vector<int> nets = Rcpp::as<std::vector<int>> (net_num);
  int netNum = nets.size();
  Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
  std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
  std::vector<std::vector<StreamSegment>> segments(netNum);
  std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
  unsigned int j = 0;
  for (unsigned int k=0; k<netNum; k++){
      unsigned int currentNet = nets[k];
      unsigned int i = 0;
      std::vector<StreamSegment> seg;
      while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
        StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
        seg.push_back(s);
        segmentsMaps[k][networkDataTot(j,1)] = binTables.front()[i];
        i++;
        j++;
      }
      segments[k] = seg;
      binTables.pop_front();
  }
  networkDataTot.resize(0,0);
  binTables.resize(0);


  // Observed points storage
  Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (obs_points);
  std::vector<std::vector<Point>> obsPoints(netNum);
  helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
  obsPointsMat.resize(0,0);
  segmentsMaps.clear();


  // Networks creation
  unsigned int nObsTot(0);
  std::vector<Network> networks(netNum);
  for (unsigned int k=0; k<netNum; k++){
    Network net(k, obsPoints[k], segments[k]);
    networks[k] = net;
    obsPoints[k].clear();
    segments[k].clear();
    networks[k].computeDistances(FALSE);
    nObsTot += networks[k].getNObs();
  }
  obsPoints.clear();
  segments.clear();
  Rcpp::Rcout << "Networks stored. \n";

  std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(FALSE, networks, nObsTot));


  Rcpp::List result = Rcpp::List::create(Rcpp::Named("flowMat") = matrices[0],
                                         Rcpp::Named("distHydro") = matrices[1]);


  return Rcpp::wrap(result);

  END_RCPP

}

// COMPUTE HYDROLOGICAL AND EUCLIDEAN DISTANCES
/**
* Compute the hydrological and Euclidean distances between observed points, as well as the flow-connection/unconnection binary matrix
* @param net_num vector containing the networkID per each network of the data set
* @param bin_tables list of vectors of strings, one vector per each network, containing the binaryIDs of the stream segments of
* each network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @return A list with the following fields:
*   - 'flowMat' flow-connection/unconnection binary matrix
*   - 'distHydro' hydrological distance matrix
*   - 'distGeo' Euclidean distance matrix
*/
RcppExport SEXP createDistanceMatrices_MultipleNets (SEXP net_num, SEXP bin_tables, SEXP network_data,
  SEXP obs_points){

  BEGIN_RCPP

  // Stream segments storage
  std::vector<int> nets = Rcpp::as<std::vector<int>> (net_num);
  int netNum = nets.size();
  Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
  std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
  std::vector<std::vector<StreamSegment>> segments(netNum);
  std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
  unsigned int j = 0;
  for (unsigned int k=0; k<netNum; k++){
      unsigned int currentNet = nets[k];
      unsigned int i = 0;
      std::vector<StreamSegment> seg;
      while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
        StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
        seg.push_back(s);
        segmentsMaps[k][networkDataTot(j,1)] = binTables.front()[i];
        i++;
        j++;
      }
      segments[k] = seg;
      binTables.pop_front();
  }
  networkDataTot.resize(0,0);
  binTables.resize(0);


  // Observed points storage
  Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (obs_points);
  std::vector<std::vector<Point>> obsPoints(netNum);
  helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
  obsPointsMat.resize(0,0);
  segmentsMaps.clear();


  // Networks creation
  unsigned int nObsTot(0);
  std::vector<Network> networks(netNum);
  for (unsigned int k=0; k<netNum; k++){
    Network net(k, obsPoints[k], segments[k]);
    networks[k] = net;
    obsPoints[k].clear();
    segments[k].clear();
    networks[k].computeDistances(TRUE);
    nObsTot += networks[k].getNObs();
  }
  obsPoints.clear();
  segments.clear();
  Rcpp::Rcout << "Networks stored. \n";

  std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(TRUE, networks, nObsTot));


  Rcpp::List result = Rcpp::List::create(Rcpp::Named("flowMat") = matrices[0],
                                         Rcpp::Named("distHydro") = matrices[1],
                                         Rcpp::Named("distGeo") = matrices[2]);


  return Rcpp::wrap(result);

  END_RCPP

}

// CREATE MODEL
/**
* Fit a spatial linear model on a spatial stream network
* @param net_num vector containing the networkID per each network of the data set
* @param bin_tables list of vector of strings, one vector per each network, containing the binaryIDs of the stream segments of
* each network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @param model_bounds vector of upper bounds for the parsills of the covariance models. It can be also NULL
* @param use_cholesky boolean indicating if the Cholesky decomposition is to be always preferred for computing the inverse of positive definite matrices
* @return A list with the following fields:
*   - 'optTheta' vector with the optimal values of the covariance parameters found through the model fitting
*   - 'betaValues' vector with the coefficients of the linear model
*   - 'covMatrix' covariance matrix computed using the optimal values of the covariance parameters
*/
RcppExport SEXP getSSNModel_MultipleNets (SEXP net_num, SEXP bin_tables, SEXP network_data,
  SEXP obs_points, SEXP obs_data, SEXP var_names, SEXP model_names, SEXP nugg,
  SEXP dist_matrices, SEXP model_bounds, SEXP use_cholesky) {

    BEGIN_RCPP

    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    int up = 0;
    int down = 0;
    int euclid = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
         tmp_tailUpModel = tailup_fac.create(name);
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
      }
    }


    // Stream segments storage
    std::vector<int> nets = Rcpp::as<std::vector<int>> (net_num);
    int netNum = nets.size();
    Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
    std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
    std::vector<std::vector<StreamSegment>> segments(netNum);
    std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
    unsigned int j = 0;
    for (unsigned int k=0; k<netNum; k++){
        unsigned int currentNet = nets[k];
        unsigned int i = 0;
        std::vector<StreamSegment> seg;
        while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
          StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
          seg.push_back(s);
          segmentsMaps[k][networkDataTot(j,1)] = binTables.front()[i];
          i++;
          j++;
        }
        segments[k] = seg;
        binTables.pop_front();
    }
    networkDataTot.resize(0,0);
    binTables.resize(0);


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<std::vector<Point>> obsPoints(netNum);
    helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
    obsPointsMat.resize(0,0);
    segmentsMaps.clear();

    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Networks creation
    unsigned int nObsTot(0);
    std::vector<Network> networks(netNum);
    for (unsigned int k=0; k<netNum; k++){
      Network net(k, obsPoints[k], segments[k]);
      networks[k] = net;
      obsPoints[k].clear();
      segments[k].clear();
      if(vecMatrices.isNotNull()){
        distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
        std::vector<Eigen::MatrixXd> tmp_vec;
        for (int i=0; i<distMatrices.size(); i++)
          tmp_vec.push_back(distMatrices[i].block(nObsTot,nObsTot,networks[k].getNObs(),networks[k].getNObs()));
        networks[k].setDistPoints(euclid, tmp_vec);
      }
      else
        networks[k].computeDistances(euclid);
      nObsTot += networks[k].getNObs();
    }
    obsPoints.clear();
    segments.clear();
    Rcpp::Rcout << "Networks stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());


    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));
    // Matrices of observed points only, about connection, distances and weight built using block matrices
    Eigen::MatrixXd weightMatOO(dataObs.computeWeightMat(weightVar));
    Eigen::MatrixXd flowMatOO(nObsTot, nObsTot);
    Eigen::MatrixXd distHydroOO(nObsTot, nObsTot);
    Eigen::MatrixXd distGeoOO(0,0);
    if(vecMatrices.isNotNull()){
      flowMatOO = distMatrices[0];
      distHydroOO = distMatrices[1];
      if (euclid)
        distGeoOO = distMatrices[2];
    }
    else {
      if (euclid){
        std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(TRUE, networks, nObsTot));
        flowMatOO = matrices[0];
        distHydroOO = matrices[1];
        distGeoOO.resize(nObsTot, nObsTot);
        distGeoOO = matrices[2];
      }
      else {
        std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(FALSE, networks, nObsTot));
        flowMatOO = matrices[0];
        distHydroOO = matrices[1];
      }
    }
    weightMatOO = weightMatOO.cwiseProduct(flowMatOO);


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // MODEL FITTING
    bool useCholesky = Rcpp::as<bool> (use_cholesky);
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, nuggetEffect, up+down+euclid,
      std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(distHydroOO), std::make_shared<Eigen::MatrixXd>(distGeoOO), std::make_shared<Eigen::MatrixXd>(weightMatOO), std::make_shared<Eigen::MatrixXi>(flowMatOO.cast<int>()), useCholesky);

    Rcpp::Nullable<std::vector<double>> vecBounds(model_bounds);
    std::vector<double> bounds;
    if(vecBounds.isNotNull()){
      bounds = Rcpp::as<std::vector<double>> (model_bounds);
      solver.setBounds(bounds);
    }

    Rcpp::Rcout << "Model fitting \n";
    solver.glmssn();

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat());


    return Rcpp::wrap(result);
    END_RCPP

}

// DO KRIGING
/**
* Perform Universal kriging on a spatial stream network object
* @param net_num vector containing the networkID per each network of the data set
* @param bin_tables list of vector of strings, one vector per each network, containing the binaryIDs of the stream segments of
* each network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param pred_points matrix whose columns correspond to the following attributes of the prediction points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param pred_data matrix whose columns correspond to the values of the responde variable and model covariates of the prediction points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param param vector with the optimal values of the covariance parameters found through the model fitting
* @param cov_mat covariance matrix computed using the optimal values of the covariance parameters
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @return A list with the following fields:
*   - 'predictions' 2-column matrix containing the predicted values and kriging variance of each prediction point
*/
RcppExport SEXP doSSNKriging_MultipleNets (SEXP net_num, SEXP bin_tables, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names, SEXP nugg, SEXP param, SEXP cov_mat, SEXP dist_matrices) {

    BEGIN_RCPP

    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    Eigen::VectorXd optParam = Rcpp::as<Eigen::VectorXd> (param);
    Eigen::MatrixXd covMatrix = Rcpp::as<Eigen::MatrixXd> (cov_mat);

    int up = 0;
    int down = 0;
    int euclid = 0;
    unsigned int j = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
        tmp_tailUpModel = tailup_fac.create(name);
        tmp_tailUpModel->setSigma2(optParam(j));
        j++;
        tmp_tailUpModel->setAlpha(optParam(j));
        j++;
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
        tmp_tailDownModel->setSigma2(optParam(j));
        j++;
        tmp_tailDownModel->setAlpha(optParam(j));
        j++;
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
        tmp_euclidModel->setSigma2(optParam(j));
        j++;
        tmp_euclidModel->setAlpha(optParam(j));
        j++;
      }
    }


    // Stream segments storage
    std::vector<int> nets = Rcpp::as<std::vector<int>> (net_num);
    int netNum = nets.size();
    Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
    std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
    std::vector<std::vector<StreamSegment>> segments(netNum);
    std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
    j = 0;
    for (unsigned int k=0; k<netNum; k++){
        unsigned int currentNet = nets[k];
        unsigned int i = 0;
        std::vector<StreamSegment> seg;
        while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
          StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
          seg.push_back(s);
          segmentsMaps[k][networkDataTot(j,1)] = binTables.front()[i];
          i++;
          j++;
        }
        segments[k] = seg;
        binTables.pop_front();
    }
    networkDataTot.resize(0,0);
    binTables.resize(0);


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<std::vector<Point>> obsPoints(netNum);
    helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
    obsPointsMat.resize(0,0);

    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<std::vector<Point>> predPoints(netNum);
    helpers::pointsStorage(segmentsMaps, predPointsMat, predPoints);
    predPointsMat.resize(0,0);
    segmentsMaps.clear();

    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Networks creation
    unsigned int nObsTot(0);
    unsigned int nPredTot(0);
    std::vector<Network> networks(netNum);
    for (unsigned int k=0; k<netNum; k++){
      Network net(k, obsPoints[k], predPoints[k], segments[k]);
      networks[k] = net;
      obsPoints[k].clear();
      predPoints[k].clear();
      segments[k].clear();
      if(vecMatrices.isNotNull()){
        distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
        std::vector<Eigen::MatrixXd> tmp_vec;
        for (int i=0; i<distMatrices.size(); i++)
          tmp_vec.push_back(distMatrices[i].block(nObsTot,nObsTot,networks[k].getNObs(),networks[k].getNObs()));
        networks[k].setDistPoints(euclid, tmp_vec);
      }
      else
        networks[k].computeDistances(euclid);
      nObsTot += networks[k].getNObs();
      nPredTot += networks[k].getNPred();
    }
    obsPoints.clear();
    predPoints.clear();
    segments.clear();
    Rcpp::Rcout << "Networks stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());
    std::vector<std::string> covNames(varNames.begin()+1, varNames.end());

    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));
    //Dataframe of the prediction points (without the response variable)
    Dataframe dataPred(covNames, Rcpp::as<Eigen::MatrixXd> (pred_data));

    //Matrices about connection, distances and weight between observed and predicion points built using block matrices
    Eigen::MatrixXd weightMatOP(dataObs.computeWeightMat(weightVar, dataPred[weightVar]));
    Eigen::MatrixXd flowMatOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroPO(nPredTot, nObsTot);
    Eigen::MatrixXd distGeoOP(0,0);
    if (euclid){
      std::vector<Eigen::MatrixXd> matrices = helpers::createDistMatricesOP(TRUE, networks, nObsTot, nPredTot);
      flowMatOP = matrices[0];
      distHydroOP = matrices[1];
      distHydroPO = matrices[2].transpose();
      distGeoOP.resize(nObsTot,nPredTot);
      distGeoOP = matrices[3];
    }
    else {
      std::vector<Eigen::MatrixXd> matrices = helpers::createDistMatricesOP(FALSE, networks, nObsTot, nPredTot);
      flowMatOP = matrices[0];
      distHydroOP = matrices[1];
      distHydroPO = matrices[2].transpose();
    }
    weightMatOP = weightMatOP.cwiseProduct(flowMatOP);


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // KRIGING
    Eigen::MatrixXd designMatPred;
    designMatPred.resize(nPredTot, varNames.size()-1);
    designMatPred.fill(1.0);
    for (unsigned int i=1; i<designMatPred.cols(); i++){
      designMatPred.col(i) = dataPred[varNames[i]];
    }


    Kriging universalKriging(std::make_shared<Eigen::MatrixXd>(designMatPred), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(covMatrix), std::make_shared<Eigen::MatrixXd>(distHydroOP), std::make_shared<Eigen::MatrixXd>(distHydroPO), std::make_shared<Eigen::MatrixXd>(distGeoOP),
      std::make_shared<Eigen::MatrixXd>(weightMatOP), std::make_shared<Eigen::MatrixXi>(flowMatOP.cast<int>()), std::make_shared<Eigen::VectorXd>(optParam), std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]),
      tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, up+down+euclid, nuggetEffect);

    Rcpp::Rcout << "Kriging \n";
    universalKriging.predict();



    Rcpp::List result = Rcpp::List::create(Rcpp::Named("predictions") = universalKriging.getPredictions());

    return Rcpp::wrap(result);
    END_RCPP

}

// CREATE MODEL and DO KRIGING
/**
* Fit a spatial linear model on a spatial stream network and perform Universal kriging
* @param net_num vector containing the networkID per each network of the data set
* @param bin_tables list of vector of strings, one vector per each network, containing the binaryIDs of the stream segments of
* each network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param pred_points matrix whose columns correspond to the following attributes of the prediction points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param pred_data matrix whose columns correspond to the values of the responde variable and model covariates of the prediction points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @param model_bounds vector of upper bounds for the parsills of the covariance models. It can be also NULL
* @param use_cholesky boolean indicating if the Cholesky decomposition is to be always preferred for computing the inverse of positive definite matrices
* @return A list with the following fields:
*   - 'optTheta' vector with the optimal values of the covariance parameters found through the model fitting
*   - 'betaValues' vector with the coefficients of the linear model
*   - 'covMatrix' covariance matrix computed using the optimal values of the covariance parameters
*   - 'predictions' 2-column matrix containing the predicted values and kriging variance of each prediction point
*/
RcppExport SEXP getSSNModelKriging_MultipleNets (SEXP net_num, SEXP bin_tables, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names, SEXP nugg, SEXP dist_matrices, SEXP model_bounds, SEXP use_cholesky) {

    BEGIN_RCPP

    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    int up = 0;
    int down = 0;
    int euclid = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
         tmp_tailUpModel = tailup_fac.create(name);
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
      }
    }


    // Stream segments storage
    std::vector<int> nets = Rcpp::as<std::vector<int>> (net_num);
    int netNum = nets.size();
    Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
    std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
    std::vector<std::vector<StreamSegment>> segments(netNum);
    std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
    unsigned int j = 0;
    for (unsigned int k=0; k<netNum; k++){
        unsigned int currentNet = nets[k];
        unsigned int i = 0;
        std::vector<StreamSegment> seg;
        while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
          StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
          seg.push_back(s);
          segmentsMaps[k][networkDataTot(j,1)] = binTables.front()[i];
          i++;
          j++;
        }
        segments[k] = seg;
        binTables.pop_front();
    }

    networkDataTot.resize(0,0);
    binTables.resize(0);


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<std::vector<Point>> obsPoints(netNum);
    helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
    obsPointsMat.resize(0,0);


    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<std::vector<Point>> predPoints(netNum);
    helpers::pointsStorage(segmentsMaps, predPointsMat, predPoints);
    predPointsMat.resize(0,0);
    segmentsMaps.clear();

    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Networks creation
    unsigned int nObsTot(0);
    unsigned int nPredTot(0);
    std::vector<Network> networks(netNum);
    for (unsigned int k=0; k<netNum; k++){
      Network net(k, obsPoints[k], predPoints[k], segments[k]);
      networks[k] = net;
      obsPoints[k].clear();
      predPoints[k].clear();
      segments[k].clear();
      if(vecMatrices.isNotNull()){
        distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
        std::vector<Eigen::MatrixXd> tmp_vec;
        for (int i=0; i<distMatrices.size(); i++)
          tmp_vec.push_back(distMatrices[i].block(nObsTot,nObsTot,networks[k].getNObs(),networks[k].getNObs()));
        networks[k].setDistPoints(euclid, tmp_vec);
      }
      else
        networks[k].computeDistances(euclid);
      nObsTot += networks[k].getNObs();
      nPredTot += networks[k].getNPred();
    }
    obsPoints.clear();
    predPoints.clear();
    segments.clear();
    Rcpp::Rcout << "Networks stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());
    std::vector<std::string> covNames(varNames.begin()+1, varNames.end());


    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));
    // Matrices of observed points only, about connection, distances and weight built using block matrices
    Eigen::MatrixXd weightMatOO(dataObs.computeWeightMat(weightVar));
    Eigen::MatrixXd flowMatOO(nObsTot, nObsTot);
    Eigen::MatrixXd distHydroOO(nObsTot, nObsTot);
    Eigen::MatrixXd distGeoOO(0,0);
    if(vecMatrices.isNotNull()){
      flowMatOO = distMatrices[0];
      distHydroOO = distMatrices[1];
      if (euclid)
        distGeoOO = distMatrices[2];
    }
    else {
      if (euclid){
        std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(TRUE, networks, nObsTot));
        flowMatOO = matrices[0];
        distHydroOO = matrices[1];
        distGeoOO.resize(nObsTot, nObsTot);
        distGeoOO = matrices[2];
      }
      else {
        std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices(FALSE, networks, nObsTot));
        flowMatOO = matrices[0];
        distHydroOO = matrices[1];
      }
    }
    weightMatOO = weightMatOO.cwiseProduct(flowMatOO);


    //Dataframe of the prediction points (without the response variable)
    Dataframe dataPred(covNames, Rcpp::as<Eigen::MatrixXd> (pred_data));


    //Matrices about connection, distances and weight between observed and predicion points built using block matrices
    Eigen::MatrixXd weightMatOP(dataObs.computeWeightMat(weightVar, dataPred[weightVar]));
    Eigen::MatrixXd flowMatOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroPO(nPredTot, nObsTot);
    Eigen::MatrixXd distGeoOP(0,0);
    if (euclid){
      std::vector<Eigen::MatrixXd> matrices = helpers::createDistMatricesOP(TRUE, networks, nObsTot, nPredTot);
      flowMatOP = matrices[0];
      distHydroOP = matrices[1];
      distHydroPO = matrices[2].transpose();
      distGeoOP.resize(nObsTot,nPredTot);
      distGeoOP = matrices[3];
    }
    else {
      std::vector<Eigen::MatrixXd> matrices = helpers::createDistMatricesOP(FALSE, networks, nObsTot, nPredTot);
      flowMatOP = matrices[0];
      distHydroOP = matrices[1];
      distHydroPO = matrices[2].transpose();
    }
    weightMatOP = weightMatOP.cwiseProduct(flowMatOP);


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }



    // -------------------------------------------------------------------------
    // MODEL FITTING
    bool useCholesky = Rcpp::as<bool> (use_cholesky);
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, nuggetEffect, up+down+euclid,
      std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(distHydroOO), std::make_shared<Eigen::MatrixXd>(distGeoOO), std::make_shared<Eigen::MatrixXd>(weightMatOO), std::make_shared<Eigen::MatrixXi>(flowMatOO.cast<int>()), useCholesky);

    Rcpp::Nullable<std::vector<double>> vecBounds(model_bounds);
    std::vector<double> bounds;
    if(vecBounds.isNotNull()){
      bounds = Rcpp::as<std::vector<double>> (model_bounds);
      solver.setBounds(bounds);
    }

    Rcpp::Rcout << "Model fitting \n";
    solver.glmssn();



    // -------------------------------------------------------------------------
    // KRIGING
    Eigen::MatrixXd designMatPred;
    designMatPred.resize(nPredTot, varNames.size()-1);
    designMatPred.fill(1.0);
    for (unsigned int i=1; i<designMatPred.cols(); i++){
      designMatPred.col(i) = dataPred[varNames[i]];
    }


    Kriging universalKriging(std::make_shared<Eigen::MatrixXd>(designMatPred), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(solver.getCovMat()), std::make_shared<Eigen::MatrixXd>(distHydroOP), std::make_shared<Eigen::MatrixXd>(distHydroPO), std::make_shared<Eigen::MatrixXd>(distGeoOP),
      std::make_shared<Eigen::MatrixXd>(weightMatOP), std::make_shared<Eigen::MatrixXi>(flowMatOP.cast<int>()), std::make_shared<Eigen::VectorXd>(solver.getOptimTheta()), std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]),
      solver.getTailUp(), solver.getTailDown(), solver.getEuclid(), up+down+euclid, nuggetEffect);

    Rcpp::Rcout << "Kriging \n";
    universalKriging.predict();



    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat(),
                                           Rcpp::Named("predictions") = universalKriging.getPredictions());

    return Rcpp::wrap(result);
    END_RCPP

}




// ---------------------------------------------------------------------------------------------------------------------------
// SINGLE NETWORK FUNCTIONS

/**
* Functions to create distance matrices, fit a spatial linear model and perform Universal kriging on a spatial stream network
* object, which has one single network.
*/

// COMPUTE HYDROLOGICAL DISTANCES
/**
* Compute the hydrological distances between observed points, as well as the flow-connection/unconnection binary matrix
* @param bin_table vector of strings containing the binaryIDs of the stream segments of the network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @return A list with the following fields:
*   - 'flowMat' flow-connection/unconnection binary matrix
*   - 'distHydro' hydrological distance matrix
*/
RcppExport SEXP createHydroDistanceMatrices_SingleNet (SEXP bin_table, SEXP network_data,
  SEXP obs_points){

  BEGIN_RCPP

  // Stream segments storage
  Eigen::MatrixXd networkData = Rcpp::as<Eigen::MatrixXd> (network_data);
  std::vector<std::string> binTable = Rcpp::as<std::vector<std::string>> (bin_table);
  std::map<unsigned int, std::string> segmentsMap;
  std::vector<StreamSegment> segments;
  for (unsigned int j = 0; j<networkData.rows(); j++){
    StreamSegment s(networkData(j,0), networkData(j,1), networkData(j,2), binTable[j]);
    segments.push_back(s);
    segmentsMap[networkData(j,1)] = binTable[j];
  }
  networkData.resize(0,0);
  binTable.clear();


  // Observed points storage
  Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
  std::vector<Point> obsPoints;
  for (unsigned int j = 0; j<obsPointsMat.rows(); j++){
    Point p;
    p.setPid(j);
    p.setRid(obsPointsMat(j,1));
    p.setDistUpstream(obsPointsMat(j,2));
    p.setCoordinates(obsPointsMat(j,3), obsPointsMat(j,4));
    auto it = segmentsMap.find(obsPointsMat(j,1));
    p.setBinaryID(it->second);
    obsPoints.push_back(p);
  }
  obsPointsMat.resize(0,0);
  segmentsMap.clear();


  // Network creation
  Network network(0, obsPoints, segments);
  obsPoints.clear();
  segments.clear();
  network.computeDistances(FALSE);
  Rcpp::Rcout << "Network stored. \n";


  Rcpp::List result = Rcpp::List::create(Rcpp::Named("flowMat") = network.getFlowMatOO(),
                                         Rcpp::Named("distHydro") = network.getDistHydroOO());


  return Rcpp::wrap(result);

  END_RCPP

}

// COMPUTE HYDROLOGICAL AND EUCLIDEAN DISTANCES
/**
* Compute the hydrological and Euclidean distances between observed points, as well as the flow-connection/unconnection binary matrix
* @param bin_table vector of strings containing the binaryIDs of the stream segments of the network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @return A list with the following fields:
*   - 'flowMat' flow-connection/unconnection binary matrix
*   - 'distHydro' hydrological distance matrix
*   - 'distGeo' Euclidean distance matrix
*/
RcppExport SEXP createDistanceMatrices_SingleNet (SEXP bin_table, SEXP network_data,
  SEXP obs_points){

  BEGIN_RCPP

  // Stream segments storage
  Eigen::MatrixXd networkData = Rcpp::as<Eigen::MatrixXd> (network_data);
  std::vector<std::string> binTable = Rcpp::as<std::vector<std::string>> (bin_table);
  std::map<unsigned int, std::string> segmentsMap;
  std::vector<StreamSegment> segments;
  for (unsigned int j = 0; j<networkData.rows(); j++){
    StreamSegment s(networkData(j,0), networkData(j,1), networkData(j,2), binTable[j]);
    segments.push_back(s);
    segmentsMap[networkData(j,1)] = binTable[j];
  }
  networkData.resize(0,0);
  binTable.clear();


  // Observed points storage
  Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
  std::vector<Point> obsPoints;
  for (unsigned int j = 0; j<obsPointsMat.rows(); j++){
    Point p;
    p.setPid(j);
    p.setRid(obsPointsMat(j,1));
    p.setDistUpstream(obsPointsMat(j,2));
    p.setCoordinates(obsPointsMat(j,3), obsPointsMat(j,4));
    auto it = segmentsMap.find(obsPointsMat(j,1));
    p.setBinaryID(it->second);
    obsPoints.push_back(p);
  }
  obsPointsMat.resize(0,0);
  segmentsMap.clear();


  // Network creation
  Network network(0, obsPoints, segments);
  obsPoints.clear();
  segments.clear();
  network.computeDistances(TRUE);
  Rcpp::Rcout << "Network stored. \n";


  Rcpp::List result = Rcpp::List::create(Rcpp::Named("flowMat") = network.getFlowMatOO(),
                                         Rcpp::Named("distHydro") = network.getDistHydroOO(),
                                         Rcpp::Named("distGeo") = network.getDistGeoOO());


  return Rcpp::wrap(result);

  END_RCPP

}

// CREATE MODEL
/**
* Fit a spatial linear model on a spatial stream network
* @param bin_table vector of strings containing the binaryIDs of the stream segments of the network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @param model_bounds vector of upper bounds for the parsills of the covariance models. It can be also NULL
* @param use_cholesky boolean indicating if the Cholesky decomposition is to be always preferred for computing the inverse of positive definite matrices
* @return A list with the following fields:
*   - 'optTheta' vector with the optimal values of the covariance parameters found through the model fitting
*   - 'betaValues' vector with the coefficients of the linear model
*   - 'covMatrix' covariance matrix computed using the optimal values of the covariance parameters
*/
RcppExport SEXP getSSNModel_SingleNet (SEXP bin_table, SEXP network_data,
  SEXP obs_points, SEXP obs_data, SEXP var_names, SEXP model_names, SEXP nugg,
  SEXP dist_matrices, SEXP model_bounds, SEXP use_cholesky) {

    BEGIN_RCPP

    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    int up = 0;
    int down = 0;
    int euclid = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
         tmp_tailUpModel = tailup_fac.create(name);
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
      }
    }


    // Stream segments storage
    Eigen::MatrixXd networkData = Rcpp::as<Eigen::MatrixXd> (network_data);
    std::vector<std::string> binTable = Rcpp::as<std::vector<std::string>> (bin_table);
    std::map<unsigned int, std::string> segmentsMap;
    std::vector<StreamSegment> segments;
    for (unsigned int j = 0; j<networkData.rows(); j++){
      StreamSegment s(networkData(j,0), networkData(j,1), networkData(j,2), binTable[j]);
      segments.push_back(s);
      segmentsMap[networkData(j,1)] = binTable[j];
    }
    networkData.resize(0,0);
    binTable.clear();


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<Point> obsPoints;
    for (unsigned int j = 0; j<obsPointsMat.rows(); j++){
      Point p;
      p.setPid(j);
      p.setRid(obsPointsMat(j,1));
      p.setDistUpstream(obsPointsMat(j,2));
      p.setCoordinates(obsPointsMat(j,3), obsPointsMat(j,4));
      auto it = segmentsMap.find(obsPointsMat(j,1));
      p.setBinaryID(it->second);
      obsPoints.push_back(p);
    }
    obsPointsMat.resize(0,0);
    segmentsMap.clear();

    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Network creation
    Network network(0, obsPoints, segments);
    obsPoints.clear();
    segments.clear();
    if(vecMatrices.isNotNull()){
      distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
      std::vector<Eigen::MatrixXd> tmp_vec;
      for (int i=0; i<distMatrices.size(); i++)
        tmp_vec.push_back(distMatrices[i]);
      network.setDistPoints(euclid, tmp_vec);
    }
    else
      network.computeDistances(euclid);
    Rcpp::Rcout << "Network stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());


    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));
    // Matrices of observed points only, about connection, distances and weight built using block matrices
    Eigen::MatrixXd weightMatOO(dataObs.computeWeightMat(weightVar));
    weightMatOO = weightMatOO.cwiseProduct(network.getFlowMatOO().cast<double>());


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(network.getNObs(), varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // MODEL FITTING
    bool useCholesky = Rcpp::as<bool> (use_cholesky);
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, nuggetEffect, up+down+euclid,
      std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(network.getDistHydroOO()), std::make_shared<Eigen::MatrixXd>(network.getDistGeoOO()), std::make_shared<Eigen::MatrixXd>(weightMatOO), std::make_shared<Eigen::MatrixXi>(network.getFlowMatOO()), useCholesky);

    Rcpp::Nullable<std::vector<double>> vecBounds(model_bounds);
    std::vector<double> bounds;
    if(vecBounds.isNotNull()){
      bounds = Rcpp::as<std::vector<double>> (model_bounds);
      solver.setBounds(bounds);
    }

    Rcpp::Rcout << "Model fitting \n";
    solver.glmssn();

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat());


    return Rcpp::wrap(result);
    END_RCPP

}

// DO KRIGING
/**
* Perform Universal kriging on a spatial stream network object
* @param bin_table vector of strings containing the binaryIDs of the stream segments of the network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param pred_points matrix whose columns correspond to the following attributes of the prediction points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param pred_data matrix whose columns correspond to the values of the responde variable and model covariates of the prediction points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param param vector with the optimal values of the covariance parameters found through the model fitting
* @param cov_mat covariance matrix computed using the optimal values of the covariance parameters
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @return A list with the following fields:
*   - 'predictions' 2-column matrix containing the predicted values and kriging variance of each prediction point
*/
RcppExport SEXP doSSNKriging_SingleNet (SEXP bin_table, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names, SEXP nugg, SEXP param, SEXP cov_mat, SEXP dist_matrices) {

    BEGIN_RCPP

    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    Eigen::VectorXd optParam = Rcpp::as<Eigen::VectorXd> (param);
    Eigen::MatrixXd covMatrix = Rcpp::as<Eigen::MatrixXd> (cov_mat);

    int up = 0;
    int down = 0;
    int euclid = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
         tmp_tailUpModel = tailup_fac.create(name);
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
      }
    }


    // Stream segments storage
    Eigen::MatrixXd networkData = Rcpp::as<Eigen::MatrixXd> (network_data);
    std::vector<std::string> binTable = Rcpp::as<std::vector<std::string>> (bin_table);
    std::map<unsigned int, std::string> segmentsMap;
    std::vector<StreamSegment> segments;
    for (unsigned int j = 0; j<networkData.rows(); j++){
      StreamSegment s(networkData(j,0), networkData(j,1), networkData(j,2), binTable[j]);
      segments.push_back(s);
      segmentsMap[networkData(j,1)] = binTable[j];
    }
    networkData.resize(0,0);
    binTable.clear();


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<Point> obsPoints;
    for (unsigned int j = 0; j<obsPointsMat.rows(); j++){
      Point p;
      p.setPid(j);
      p.setRid(obsPointsMat(j,1));
      p.setDistUpstream(obsPointsMat(j,2));
      p.setCoordinates(obsPointsMat(j,3), obsPointsMat(j,4));
      auto it = segmentsMap.find(obsPointsMat(j,1));
      p.setBinaryID(it->second);
      obsPoints.push_back(p);
    }
    obsPointsMat.resize(0,0);


    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<Point> predPoints;
    for (unsigned int j = 0; j<predPointsMat.rows(); j++){
      Point p;
      p.setPid(j);
      p.setRid(predPointsMat(j,1));
      p.setDistUpstream(predPointsMat(j,2));
      p.setCoordinates(predPointsMat(j,3), predPointsMat(j,4));
      auto it = segmentsMap.find(predPointsMat(j,1));
      p.setBinaryID(it->second);
      predPoints.push_back(p);
    }
    predPointsMat.resize(0,0);
    segmentsMap.clear();


    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Network creation
    Network network(0, obsPoints, predPoints, segments);
    obsPoints.clear();
    segments.clear();
    if(vecMatrices.isNotNull()){
      distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
      std::vector<Eigen::MatrixXd> tmp_vec;
      for (int i=0; i<distMatrices.size(); i++)
        tmp_vec.push_back(distMatrices[i]);
      network.setDistPoints(euclid, tmp_vec);
    }
    else
      network.computeDistances(euclid);
    Rcpp::Rcout << "Network stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());
    std::vector<std::string> covNames(varNames.begin()+1, varNames.end());

    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));

    //Dataframe of the prediction points (without the response variable)
    Dataframe dataPred(covNames, Rcpp::as<Eigen::MatrixXd> (pred_data));


    //Matrices about connection, distances and weight between observed and predicion points built using block matrices
    Eigen::MatrixXd weightMatOP(dataObs.computeWeightMat(weightVar, dataPred[weightVar]));
    weightMatOP = weightMatOP.cwiseProduct(network.getFlowMatOP().cast<double>());


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(network.getNObs(), varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // KRIGING
    Eigen::MatrixXd designMatPred;
    designMatPred.resize(network.getNPred(), varNames.size()-1);
    designMatPred.fill(1.0);
    for (unsigned int i=1; i<designMatPred.cols(); i++){
      designMatPred.col(i) = dataPred[varNames[i]];
    }

    Kriging universalKriging(std::make_shared<Eigen::MatrixXd>(designMatPred), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(covMatrix), std::make_shared<Eigen::MatrixXd>(network.getDistHydroOP()), std::make_shared<Eigen::MatrixXd>(network.getDistHydroPO()), std::make_shared<Eigen::MatrixXd>(network.getDistGeoOP()),
      std::make_shared<Eigen::MatrixXd>(weightMatOP), std::make_shared<Eigen::MatrixXi>(network.getFlowMatOP()), std::make_shared<Eigen::VectorXd>(optParam), std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]),
      tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, up+down+euclid, nuggetEffect);

    Rcpp::Rcout << "Kriging \n";
    universalKriging.predict();



    Rcpp::List result = Rcpp::List::create(Rcpp::Named("predictions") = universalKriging.getPredictions());


    return Rcpp::wrap(result);
    END_RCPP

}

// CREATE MODEL and DO KRIGING
/**
* Fit a spatial linear model on a spatial stream network and perform Universal kriging
* @param bin_table vector of strings containing the binaryIDs of the stream segments of the network
* @param network_data matrix whose columns correspond to the following attributes of the stream segments: networkID, segmentID, distance upstream
* @param obs_points matrix whose columns correspond to the following attributes of the observed points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param pred_points matrix whose columns correspond to the following attributes of the prediction points: networkID, segmentID, distance upstream,
* first coordinate, second coordinate
* @param obs_data matrix whose columns correspond to the values of the responde variable and model covariates of the observed points
* @param pred_data matrix whose columns correspond to the values of the responde variable and model covariates of the prediction points
* @param var_names vector of strings indicating first the name of the response variable and then the name of the covariates
* @param model_names vector of strings indicating the covariance models selected
* @param nugg boolean indicating if the nugget effect is to be considered in the mixed model
* @param dist_matrices vector of distance matrices between observed points. It can be also NULL
* @param model_bounds vector of upper bounds for the parsills of the covariance models. It can be also NULL
* @param use_cholesky boolean indicating if the Cholesky decomposition is to be always preferred for computing the inverse of positive definite matrices
* @return A list with the following fields:
*   - 'optTheta' vector with the optimal values of the covariance parameters found through the model fitting
*   - 'betaValues' vector with the coefficients of the linear model
*   - 'covMatrix' covariance matrix computed using the optimal values of the covariance parameters
*   - 'predictions' 2-column matrix containing the predicted values and kriging variance of each prediction point
*/
RcppExport SEXP getSSNModelKriging_SingleNet (SEXP bin_table, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names, SEXP nugg, SEXP dist_matrices, SEXP model_bounds, SEXP use_cholesky) {

    BEGIN_RCPP


    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);

    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

    bool nuggetEffect = Rcpp::as<bool> (nugg);

    int up = 0;
    int down = 0;
    int euclid = 0;
    for (auto name: corModels){
      std::size_t found_up = name.find("up");
      std::size_t found_down = name.find("down");
      std::size_t found_euclid = name.find("Euclid");
      if (found_up!=std::string::npos){
        up++;
         tmp_tailUpModel = tailup_fac.create(name);
      }
      if (found_down!=std::string::npos){
        down++;
        tmp_tailDownModel = taildown_fac.create(name);
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        tmp_euclidModel = euclid_fac.create(name);
      }
    }


    // Stream segments storage
    Eigen::MatrixXd networkData = Rcpp::as<Eigen::MatrixXd> (network_data);

    std::vector<std::string> binTable = Rcpp::as<std::vector<std::string>> (bin_table);
    std::map<unsigned int, std::string> segmentsMap;
    std::vector<StreamSegment> segments;
    for (unsigned int j = 0; j<networkData.rows(); j++){
      StreamSegment s(networkData(j,0), networkData(j,1), networkData(j,2), binTable[j]);
      segments.push_back(s);
      segmentsMap[networkData(j,1)] = binTable[j];
    }
    networkData.resize(0,0);
    binTable.clear();


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::MatrixXd> (obs_points);
    std::vector<Point> obsPoints;
    for (unsigned int j = 0; j<obsPointsMat.rows(); j++){
      Point p;
      p.setPid(j);
      p.setRid(obsPointsMat(j,1));
      p.setDistUpstream(obsPointsMat(j,2));
      p.setCoordinates(obsPointsMat(j,3), obsPointsMat(j,4));
      auto it = segmentsMap.find(obsPointsMat(j,1));
      p.setBinaryID(it->second);
      obsPoints.push_back(p);
    }
    obsPointsMat.resize(0,0);


    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<Point> predPoints;
    for (unsigned int j = 0; j<predPointsMat.rows(); j++){
      Point p;
      p.setPid(j);
      p.setRid(predPointsMat(j,1));
      p.setDistUpstream(predPointsMat(j,2));
      p.setCoordinates(predPointsMat(j,3), predPointsMat(j,4));
      auto it = segmentsMap.find(predPointsMat(j,1));
      p.setBinaryID(it->second);
      predPoints.push_back(p);
    }
    predPointsMat.resize(0,0);
    segmentsMap.clear();


    Rcpp::Nullable<std::vector<Eigen::MatrixXd>> vecMatrices(dist_matrices);
    std::vector<Eigen::MatrixXd> distMatrices;
    // Network creation
    Network network(0, obsPoints, predPoints, segments);
    obsPoints.clear();
    segments.clear();
    if(vecMatrices.isNotNull()){
      distMatrices = Rcpp::as<std::vector<Eigen::MatrixXd>> (dist_matrices);
      std::vector<Eigen::MatrixXd> tmp_vec;
      for (int i=0; i<distMatrices.size(); i++)
        tmp_vec.push_back(distMatrices[i]);
      network.setDistPoints(euclid, tmp_vec);
    }
    else
      network.computeDistances(euclid);
    Rcpp::Rcout << "Network stored. \n";


    // Variables (response variable, covariates and weight variable) names
    std::vector<std::string> varNames = Rcpp::as<std::vector<std::string>> (var_names);
    std::string weightVar(varNames.back());
    std::vector<std::string> covNames(varNames.begin()+1, varNames.end());

    // Dataframe for fitting the model, regarding the observed points
    Dataframe dataObs(varNames, Rcpp::as<Eigen::MatrixXd> (obs_data));
    // Matrices of observed points only, about connection, distances and weight built using block matrices
    Eigen::MatrixXd weightMatOO(dataObs.computeWeightMat(weightVar));
    weightMatOO = weightMatOO.cwiseProduct(network.getFlowMatOO().cast<double>());


    //Dataframe of the prediction points (without the response variable)
    Dataframe dataPred(covNames, Rcpp::as<Eigen::MatrixXd> (pred_data));


    //Matrices about connection, distances and weight between observed and predicion points built using block matrices
    Eigen::MatrixXd weightMatOP(dataObs.computeWeightMat(weightVar, dataPred[weightVar]));
    weightMatOP = weightMatOP.cwiseProduct(network.getFlowMatOP().cast<double>());


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(network.getNObs(), varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // MODEL FITTING
    bool useCholesky = Rcpp::as<bool> (use_cholesky);
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, nuggetEffect, up+down+euclid,
      std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(network.getDistHydroOO()), std::make_shared<Eigen::MatrixXd>(network.getDistGeoOO()), std::make_shared<Eigen::MatrixXd>(weightMatOO), std::make_shared<Eigen::MatrixXi>(network.getFlowMatOO()),useCholesky);

    Rcpp::Nullable<std::vector<double>> vecBounds(model_bounds);
    std::vector<double> bounds;
    if(vecBounds.isNotNull()){
      bounds = Rcpp::as<std::vector<double>> (model_bounds);
      solver.setBounds(bounds);
    }

    Rcpp::Rcout << "Model fitting \n";
    solver.glmssn();


    // -------------------------------------------------------------------------
    // KRIGING
    Eigen::MatrixXd designMatPred;
    designMatPred.resize(network.getNPred(), varNames.size()-1);
    designMatPred.fill(1.0);
    for (unsigned int i=1; i<designMatPred.cols(); i++){
      designMatPred.col(i) = dataPred[varNames[i]];
    }

    Kriging universalKriging(std::make_shared<Eigen::MatrixXd>(designMatPred), std::make_shared<Eigen::MatrixXd>(designMat), std::make_shared<Eigen::MatrixXd>(solver.getCovMat()), std::make_shared<Eigen::MatrixXd>(network.getDistHydroOP()), std::make_shared<Eigen::MatrixXd>(network.getDistHydroPO()), std::make_shared<Eigen::MatrixXd>(network.getDistGeoOP()),
      std::make_shared<Eigen::MatrixXd>(weightMatOP), std::make_shared<Eigen::MatrixXi>(network.getFlowMatOP()), std::make_shared<Eigen::VectorXd>(solver.getOptimTheta()), std::make_shared<Eigen::VectorXd>(dataObs[varNames[0]]),
      solver.getTailUp(), solver.getTailDown(), solver.getEuclid(), up+down+euclid, nuggetEffect);

    Rcpp::Rcout << "Kriging \n";
    universalKriging.predict();



    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat(),
                                           Rcpp::Named("predictions") = universalKriging.getPredictions());


    return Rcpp::wrap(result);
    END_RCPP

}

}

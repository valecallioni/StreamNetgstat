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

extern "C"{

// CREATE MODEL
RcppExport SEXP getSSNModel (SEXP net_num, SEXP bin_tables, SEXP network_data,
  SEXP obs_points, SEXP obs_data, SEXP var_names, SEXP model_names) {

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
    Rcpp::Rcout << "Stream segments stored. \n";
    networkDataTot.resize(0,0);
    binTables.resize(0);


    // Observed points storage
    Eigen::MatrixXd obsPointsMat = Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (obs_points);
    std::vector<std::vector<Point>> obsPoints(netNum);
    helpers::pointsStorage(segmentsMaps, obsPointsMat, obsPoints);
    obsPointsMat.resize(0,0);
    segmentsMaps.clear();
    Rcpp::Rcout << "obsPointsMat stored. \n";


    // Networks creation
    unsigned int nObsTot(0);
    std::vector<Network> networks(netNum);
    for (unsigned int k=0; k<netNum; k++){
      Network net(k, obsPoints[k], segments[k]);
      networks[k] = net;
      obsPoints[k].clear();
      segments[k].clear();
      //networks[k].print();
      networks[k].computeDistances();
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
    Eigen::MatrixXd distGeoOO(nObsTot, nObsTot);
    std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices("obs", networks, nObsTot));
    flowMatOO = matrices[0];
    distHydroOO = matrices[1];
    distGeoOO = matrices[2];
    matrices.clear();
    weightMatOO = weightMatOO.cwiseProduct(flowMatOO);
    //Rcpp::Rcout << "Distance matrices completed. \n";
    //Rcpp::Rcout << "weightMatOO[1:20, 1:10]: \n" << weightMatOO.block(0,0,20,10) << "\n";
    //Rcpp::Rcout << "flowMatOO[1:20, 1:10]: \n" << flowMatOO.block(0,0,20,10) << "\n";



    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

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

    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // -------------------------------------------------------------------------
    // MODEL FITTING
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, TRUE, up+down+euclid,
      dataObs[varNames[0]], designMat, distHydroOO, distGeoOO, weightMatOO, flowMatOO.cast<int>());
    Rcpp::Rcout << "Starting model fitting \n";
    solver.glmssn();
    Rcpp::Rcout << "End model fitting \n";

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat());


    return Rcpp::wrap(result);
    END_RCPP

}


// CREATE MODEL and DO KRIGING
RcppExport SEXP getSSNModelKriging (SEXP net_num, SEXP bin_tables, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names) {

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


    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<std::vector<Point>> predPoints(netNum);
    helpers::pointsStorage(segmentsMaps, predPointsMat, predPoints);
    predPointsMat.resize(0,0);
    segmentsMaps.clear();


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
      //networks[k].print();
      networks[k].computeDistances();
      nObsTot += networks[k].getNObs();
      nPredTot += networks[k].getNPred();
    }
    obsPoints.clear();
    predPoints.clear();
    segments.clear();


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
    Eigen::MatrixXd distGeoOO(nObsTot, nObsTot);
    std::vector<Eigen::MatrixXd> matrices(helpers::createDistMatrices("obs", networks, nObsTot));
    flowMatOO = matrices[0];
    distHydroOO = matrices[1];
    distGeoOO = matrices[2];
    matrices.clear();
    weightMatOO = weightMatOO.cwiseProduct(flowMatOO);


    //Dataframe of the prediction points (without the response variable)
    Dataframe dataPred(covNames, Rcpp::as<Eigen::MatrixXd> (pred_data));
    // Matrices of prediction points only, about connection, distances and weight built using block matrices
    Eigen::MatrixXd weightMatPP(dataPred.computeWeightMat(weightVar));
    Eigen::MatrixXd flowMatPP(nPredTot, nPredTot);
    Eigen::MatrixXd distHydroPP(nPredTot, nPredTot);
    Eigen::MatrixXd distGeoPP(nPredTot, nPredTot);
    matrices = helpers::createDistMatrices("pred", networks, nPredTot);
    flowMatPP = matrices[0];
    distHydroPP = matrices[1];
    distGeoPP = matrices[2];
    matrices.clear();
    weightMatPP = weightMatPP.cwiseProduct(flowMatPP);


    //Matrices about connection, distances and weight between observed and predicion points built using block matrices
    Eigen::MatrixXd weightMatOP(dataObs.computeWeightMat(weightVar, dataPred[weightVar]));
    Eigen::MatrixXd flowMatOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroOP(nObsTot, nPredTot);
    Eigen::MatrixXd distHydroPO(nPredTot, nObsTot);
    Eigen::MatrixXd distGeoOP(nObsTot, nPredTot);
    matrices = helpers::createDistMatricesOP(networks, nObsTot, nPredTot);
    flowMatOP = matrices[0];
    distHydroOP = matrices[1];
    distHydroPO = matrices[2].transpose();
    distGeoOP = matrices[3];
    matrices.clear();
    weightMatOP = weightMatOP.cwiseProduct(flowMatOP);


    // Creation of the factories related to the covariance models chosen
    std::vector<std::string> corModels = Rcpp::as<std::vector<std::string>> (model_names);
    tailup_factory::TailUpFactory& tailup_fac (tailup_factory::TailUpFactory::Instance());
    std::unique_ptr<TailUpModel> tmp_tailUpModel;
    taildown_factory::TailDownFactory& taildown_fac (taildown_factory::TailDownFactory::Instance());
    std::unique_ptr<TailDownModel> tmp_tailDownModel;
    euclidean_factory::EuclideanFactory& euclid_fac(euclidean_factory::EuclideanFactory::Instance());
    std::unique_ptr<EuclideanModel> tmp_euclidModel;

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


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }



    // -------------------------------------------------------------------------
    // MODEL FITTING
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, TRUE, up+down+euclid,
      dataObs[varNames[0]], designMat, distHydroOO, distGeoOO, weightMatOO, flowMatOO.cast<int>());
    solver.glmssn();



    // -------------------------------------------------------------------------
    // KRIGING
    Eigen::MatrixXd designMatPred;
    designMatPred.resize(nPredTot, varNames.size()-1);
    designMatPred.fill(1.0);
    for (unsigned int i=1; i<designMatPred.cols(); i++){
      designMatPred.col(i) = dataPred[varNames[i]];
    }

    Kriging universalKriging(designMatPred, designMat, solver.getCovMat(), distHydroOP, distHydroPO, distGeoOP,
      weightMatOP, flowMatOP.cast<int>(), solver.getOptimTheta(), dataObs[varNames[0]],
      solver.getTailUp(), solver.getTailDown(), solver.getEuclid(), up+down+euclid, TRUE);

    universalKriging.predict();






    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat(),
                                           Rcpp::Named("predictions") = universalKriging.getPredictions());


    return Rcpp::wrap(result);
    END_RCPP

}

}

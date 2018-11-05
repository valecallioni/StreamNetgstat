#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <list>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Dataframe.hpp"
#include "Network.hpp"
#include "Factory.hpp"
#include "FactoryHelpers.hpp"
#include "Helpers.hpp"
#include "Proxy.hpp"
#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"
#include "Optimizer.hpp"

extern "C"{

// CREATE MODEL
RcppExport SEXP getSSNM (SEXP net_num, SEXP bin_tables, SEXP network_data, SEXP obs_points, SEXP pred_points,
  SEXP obs_data, SEXP pred_data, SEXP var_names, SEXP model_names) {

    BEGIN_RCPP

    // Stream segments storage
    int netNum = Rcpp::as<int> (net_num);
    Eigen::MatrixXd networkDataTot = Rcpp::as<Eigen::MatrixXd> (network_data);
    Rcpp::Rcout << "Network data OK. \n";
    std::list<std::vector<std::string>> binTables = Rcpp::as<std::list<std::vector<std::string>>> (bin_tables);
    Rcpp::Rcout << "binTables OK. Dimensions:" << binTables.size() << "\n";
    std::vector<std::vector<StreamSegment>> segments(netNum);
    std::vector<std::map<unsigned int, std::string>> segmentsMaps(netNum);
    unsigned int j = 0;
    for (unsigned int k=0; k<netNum; k++){
        //Rcpp::Rcout << "Ciclo for \n";
        unsigned int currentNet = k+1;
        unsigned int i = 0;
        std::vector<StreamSegment> seg;
        while (j<networkDataTot.rows() && networkDataTot(j,0)==currentNet){
          //Rcpp::Rcout << "While loop, iteration j = " << j << "\n";
          StreamSegment s(currentNet, networkDataTot(j,1), networkDataTot(j,2), binTables.front()[i]);
          //Rcpp::Rcout << "Elemento i = " << i << " del bin.table della net n. " << k+1 << "letto \n";
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
    Rcpp::Rcout << "ObsPoints OK. \n";
    obsPointsMat.resize(0,0);


    // Predicted points (points for prediction) storage
    Eigen::MatrixXd predPointsMat = Rcpp::as<Eigen::MatrixXd> (pred_points);
    std::vector<std::vector<Point>> predPoints(netNum);
    helpers::pointsStorage(segmentsMaps, predPointsMat, predPoints);
    Rcpp::Rcout << "PredPoints OK. \n";
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
    Rcpp::Rcout << "Networks OK. Numero Networks = " << networks.size() << "\n";
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
    //Rcpp::Rcout << "Final distGeo: \n" << distGeoOO << "\n";
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
    Eigen::MatrixXd distGeoOP(nObsTot, nPredTot);
    matrices = helpers::createDistMatricesOP(networks, nObsTot, nPredTot);
    flowMatOP = matrices[0];
    distHydroOP = matrices[1];
    distGeoOP = matrices[2];
    matrices.clear();
    weightMatOP = weightMatOP.cwiseProduct(flowMatOP);


    //varNames.pop_back();
    //covNames.pop_back();


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
        if (up>1){
          std::cerr << "No more than one tail up model is allowed. \n";
          //Stop
        }
        else {
          tmp_tailUpModel = tailup_fac.create(name);
        }
      }
      if (found_down!=std::string::npos){
        down++;
        if (down>1){
          std::cerr << "No more than one tail down model is allowed. \n";
          //Stop
        }
        else {
          tmp_tailDownModel = taildown_fac.create(name);
        }
      }
      if (found_euclid!=std::string::npos){
        euclid++;
        if (euclid>1){
          std::cerr << "No more than one tail up model is allowed. \n";
          //Stop
        }
        else {
          tmp_euclidModel = euclid_fac.create(name);
        }
      }
    }
    Rcpp::Rcout << "Modelli creati \n";


    // Design matrix
    Eigen::MatrixXd designMat;
    designMat.resize(nObsTot, varNames.size()-1);
    designMat.fill(1.0);
    for (unsigned int i=1; i<designMat.cols(); i++){
      designMat.col(i) = dataObs[varNames[i]];
    }


    // Optimizer
    Optimizer solver(tmp_tailUpModel, tmp_tailDownModel, tmp_euclidModel, TRUE, up+down+euclid,
      dataObs[varNames[0]], designMat, distHydroOO, distGeoOO, weightMatOO, flowMatOO.cast<int>());

    //Eigen::VectorXd theta0(solver.thetaInit());
    //solver.computeLogL(theta0);
    solver.glmssn();


    Rcpp::List result = Rcpp::List::create(Rcpp::Named("optTheta") = solver.getOptimTheta(),
                                           Rcpp::Named("betaValues") = solver.getBeta(),
                                           Rcpp::Named("covMatrix") = solver.getCovMat());


    return Rcpp::wrap(result);
    END_RCPP

}
}

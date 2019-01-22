#ifndef PACSPROJECT_KRIGING_HPP
#define PACSPROJECT_KRIGING_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <cmath>
#include "Dataframe.hpp"
#include "Network.hpp"
#include "Factory.hpp"
#include "FactoryHelpers.hpp"
#include "Helpers.hpp"
#include "Proxy.hpp"
#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"

/**
* Kriging class, for performing universal kriging.
* It calculates prediction values and kriging variance for prediction sites based on the results of a linear model.
* The matrices, used for computing the covariance, describe relationships between observed and prediction points.
*/

class Kriging {

private:
  std::shared_ptr<Eigen::MatrixXd> Xpred; ///< design matrix of the prediction points
  std::shared_ptr<Eigen::MatrixXd> Xobs; ///< design matrix of the observed points
  std::shared_ptr<Eigen::MatrixXd> invV; ///< inverse of the covariance matrix
  std::shared_ptr<Eigen::MatrixXi> flowMat; ///< flow-connection/unconnection binary matrix
  std::shared_ptr<Eigen::MatrixXd> weightMat; ///< weight matrix for the tail-up model
  std::shared_ptr<Eigen::MatrixXd> distHydroOP; ///< hydrological distance matrix observed-prediction points
  std::shared_ptr<Eigen::MatrixXd> distHydroPO; ///< hydrological distance matrix prediction-observed points
  std::shared_ptr<Eigen::MatrixXd> distGeo; ///< Euclidean distance matrix observed-prediction points

  std::unique_ptr<TailUpModel> tailUpModel; ///< tail-up model
  std::unique_ptr<TailDownModel> tailDownModel; ///< tail-down model
  std::unique_ptr<EuclideanModel> euclidModel; ///< Euclidean model

  double parsill = 0.0; ///< total parsill
  int nModels; ///< number of covariance models
  bool useNugget; ///< boolean, if TRUE the nugget effect is considered in the mixed model
  std::shared_ptr<Eigen::VectorXd> theta; ///< optimal values of the covariance parameters found
  std::shared_ptr<Eigen::VectorXd> z; ///< vector of the response variable values

  Eigen::MatrixXd Vpred; ///< covariance matrix observed-prediction points
  Eigen::MatrixXd predData; ///< matrix to be filled with the prediction values and the kriging variance
  Eigen::MatrixXd XV; ///< matrix result of the product between Xobs and V
  Eigen::MatrixXd invXVX; ///< inverse matrix of Xobs %*% invV %*% Xobs

  unsigned int nObs; ///< number of observed points
  unsigned int nPred; ///< number of prediction points
  unsigned int p; ///< number of covariates used in the model

public:
  /**
  * Default constructor.
  */
  Kriging() =  default;

  /**
  * Constructor.
  * @param dMatPred design matrix of the prediction points
  * @param dMatObs design matrix of the observed points
  * @param V covariance matrix between observed points
  * @param Nop hydrological distance matrix observed-prediction points
  * @param Npo hydrological distance matrix prediction-observed points
  * @param D Euclidean distance matrix
  * @param wMat weight matrix for tail-up model
  * @param connMat flow-connection/unconnection binary matrix
  * @param param
  * @param y vector of the response variable values
  * @param tailup_ptr pointer to the tail-up model selected
  * @param taildown_ptr pointer to the tail-down model selected
  * @param euclid_ptr pointer to the Euclidean model selected
  * @param nMod number of covariance models
  * @param useNugg boolean, indicating if the nuggect effects is to be considered in the mixed model
  */
  Kriging(const Eigen::MatrixXd& dMatPred, const Eigen::MatrixXd& dMatObs, const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Nop, const Eigen::MatrixXd& Npo, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat,
    const Eigen::VectorXd& param, const Eigen::VectorXd& y, std::unique_ptr<TailUpModel>& tailup_ptr,
    std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, int nMod, bool useNugg);

  /**
  * @return the matrix containing the predicted values and the kriging variance
  */
  const Eigen::MatrixXd& getPredictions() const {return predData;};

  /**
  * Perform Universal kriging
  */
  void predict();

};

#endif

#ifndef PACSPROJECT_OPTIMIZER_HPP
#define PACSPROJECT_OPTIMIZER_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>
#include "Factory.hpp"
#include "FactoryHelpers.hpp"
#include "Helpers.hpp"
#include "Proxy.hpp"
#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"

/**
* Optimizer class, for fitting a spatial linear model with spatially autocorrelated errors using normal likelihood methods.
* It uses matrices regarding relationships between observed points.
*/

class Optimizer {

private:
  std::shared_ptr<Eigen::MatrixXd> distHydro; ///< hydrological distance matrix
  std::shared_ptr<Eigen::MatrixXd> distGeo; ///< Euclidean distance matrix
  std::shared_ptr<Eigen::MatrixXd> weightMat; ///< weight matrix for tail-up model
  std::shared_ptr<Eigen::MatrixXi> flowMat; ///< flow-connection/unconnection binary matrix

  std::shared_ptr<Eigen::VectorXd> z; ///< vector of the response variable values
  std::shared_ptr<Eigen::MatrixXd> X; ///< design matrix of the linear model
  unsigned int n; ///< number of observed points
  unsigned int p; ///< number of covariates

  bool useNugget; ///< boolean, if TRUE the nugget effect is considered in the mixed model
  int nModels; ///< number of covariance models considered

  std::unique_ptr<TailUpModel> tailUpModel; ///< tail-up model
  std::unique_ptr<TailDownModel> tailDownModel; ///< tail-down model
  std::unique_ptr<EuclideanModel> euclidModel; ///< Euclidean model
  double bound_up = 1e+4; ///< upper bound for parsill of tail-up model
  double bound_down = 1e+4; ///< upper bound for parsill of tail-down model
  double bound_eu = 1e+4; ///< upper bound for parsill of Euclidean model

  Eigen::VectorXd optimTheta; ///< optimal values of the covariance parameters found
  Eigen::VectorXd betaValues; ///< coefficients of the linear model
  Eigen::MatrixXd covMat; ///< final covariance matrix

  const double a = 1.0; ///< reflection parameter of Nelder Mead
  const double c = 2.0; ///< expansion parameter of Nelder Mead
  const double r = 0.5; ///< contraction parameter of Nelder Mead
  const double s = 0.5; ///< shrinkage parameter of Nelder Mead

  unsigned int iter = 0; ///< number of iterations needed for convergence
  const unsigned int maxIter = 1e3; ///< maximum number of iterations allowed
  unsigned int funEvals = 0; ///< number of function evaluations needed for convergence
  const unsigned int maxFunEvals = 1e6; ///< maximum number of function evaluations allowed
  double tolTheta = 1e-3; ///< tolerance on difference in norm of theta values
  double tolFun = 1e-3; ///< tolerance on difference in function values

  double maxDistHydro = 1e+10; ///< upper bound for range of tail-up and tail-down models
  double maxDistGeo = 1e+10; ///< upper bound for range of Euclidean model

  bool useLDLT = FALSE; ///< boolean, indicating if the user wants to always use the Cholesky decomposition for positive definite matrix inversion

public:
  /**
  * Constructor.
  * @param tailup_ptr pointer to the tail-up model selected
  * @param taildown_ptr pointer to the tail-down model selected
  * @param euclid_ptr pointer to the Euclidean model selected
  * @param useNugg boolean, indicating if the nuggect effects is to be considered in the mixed model
  * @param n_models number of covariance models
  * @param y vector of the response variable values
  * @param designMat design matrix of the linear model
  * @param N hydrological distance matrix
  * @param D Euclidean distance matrix
  * @param wMat weight matrix for tail-up model
  * @param connMat flow-connection/unconnection binary matrix
  * @param cholesky boolean, indicating if the user wants to always use the Cholesky decomposition when inverting positive definite matrices
  */
  Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg, int n_models,
    const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat, bool cholesky);

  /**
  * Constructor.
  * @param tailup_ptr pointer to the tail-up model selected
  * @param taildown_ptr pointer to the tail-down model selected
  * @param euclid_ptr pointer to the Euclidean model selected
  * @param useNugg boolean, indicating if the nuggect effects is to be considered in the mixed model
  * @param bounds vector of upper bounds for the parsills of the covariance models
  * @param n_models umber of covariance models
  * @param y vector of the response variable values
  * @param designMat designMat design matrix of the linear model
  * @param N hydrological distance matrix
  * @param D Euclidean distance matrix
  * @param wMat weight matrix for tail-up model
  * @param connMat flow-connection/unconnection binary matrix
  * @param cholesky boolean, indicating if the user wants to always use the Cholesky decomposition when inverting positive definite matrices
  */
  Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg, const std::vector<double>& bounds,
    int n_models, const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat, bool cholesky);

  /**
  * @return the vector with the optimal values found for the covariance models
  */
  const Eigen::VectorXd& getOptimTheta() const {return optimTheta;};

  /**
  * @return the vector with the coefficients values for the linear model
  */
  const Eigen::VectorXd& getBeta() const {return betaValues;};

  /**
  * @return the covariance matrix corresponding to the values of the parameters found
  */
  const Eigen::MatrixXd& getCovMat() const {return covMat;};

  /**
  * @return a smart pointer to the tail-up model selected
  */
  std::unique_ptr<TailUpModel>& getTailUp() {return tailUpModel;};

  /**
  * @return a smart pointer to the tail-down model selected
  */
  std::unique_ptr<TailDownModel>& getTailDown() {return tailDownModel;};

  /**
  * @return a smart pointer to the Euclidean model selected
  */
  std::unique_ptr<EuclideanModel>& getEuclid() {return euclidModel;};

  /**
  * @param bounds vector of upper bounds for the parsills of the covariance models
  */
  void setBounds(const std::vector<double>& bounds);

  /**
  * Update the values of the covariance paramters during the steps of the Nelder Mead algorithm
  * @param theta vector with the current values of the covariance parameters
  * @return TRUE if the current values of the covariance paramteres are within the bounds, FALSE otherwise
  */
  bool updateParam(const Eigen::VectorXd& theta);

  /**
  * Initialize the values of the covariance parameters to find the first point of the simplex
  * @return the vector with the initial values of the covariance parameters
  */
  Eigen::VectorXd thetaInit();

  /**
  * Compute the log-likelihood given the current values of the covariance parameters
  * @param theta vector with the current values of the covariance parameters
  * @return the value of the log-likelihood
  */
  double computeLogL(const Eigen::VectorXd& theta);

  /**
  * Compute the initial simplex for the Nelder Mead algorithm
  * @param theta0 vector with the current values of the covariance parameters
  * @return TRUE if the current values of the covariance paramteres are within the bounds, FALSE otherwise
  */
  std::vector<std::pair<double, Eigen::VectorXd>> simplexInit(const Eigen::VectorXd& theta0, const double tau);

  /**
  * Implementation of the Nelder Mead algorithm
  */
  void computeTheta();

  /**
  * Function that calls the Nelder Mead algorithm and stores the optimal values of the covariance parameters,
  * the coefficients and the covariance matrix
  * @see computeTheta
  */
  void glmssn();

};

#endif

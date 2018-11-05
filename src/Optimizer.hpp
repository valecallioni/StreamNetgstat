#ifndef PACSPROJECT_OPTIMIZER_HPP
#define PACSPROJECT_OPTIMIZER_HPP

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

class Optimizer {

private:
  std::shared_ptr<Eigen::MatrixXd> distHydro;
  std::shared_ptr<Eigen::MatrixXd> distGeo;
  std::shared_ptr<Eigen::MatrixXd> weightMat;
  std::shared_ptr<Eigen::MatrixXi> flowMat;

  std::shared_ptr<Eigen::VectorXd> z;
  std::shared_ptr<Eigen::MatrixXd> X;
  unsigned int n;
  unsigned int p;

  bool useNugget;
  int nModels;

  std::unique_ptr<TailUpModel> tailUpModel;
  std::unique_ptr<TailDownModel> tailDownModel;
  std::unique_ptr<EuclideanModel> euclidModel;

  Eigen::VectorXd optimTheta;
  Eigen::VectorXd betaValues;
  Eigen::MatrixXd covMat;

  unsigned int iter = 0;
  unsigned int maxIter = 1e4;
  unsigned int funEvals = 0;
  unsigned int maxFunEvals = 1e6;
  unsigned int restart = 0;
  unsigned int maxRestart = 100;
  double tolTheta = 1e-3;
  double tolFun = 1e-3;

  double maxDistHydro = 1e+10;
  double maxDistGeo = 1e+10;

public:
  Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg, int n_models,
    const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat);

  const Eigen::VectorXd& getOptimTheta() const {return optimTheta;};
  const Eigen::VectorXd& getBeta() const {return betaValues;};
  const Eigen::MatrixXd& getCovMat() const {return covMat;};

  Eigen::VectorXd thetaInit();
  double computeLogL(Eigen::VectorXd& theta);
  std::vector<std::pair<double, Eigen::VectorXd>> simplexInit(const Eigen::VectorXd& theta0, const double tau);
  void computeThetaPaper();
  void computeThetaWiki();
  void glmssn();

  void test(Eigen::VectorXd& theta);
};

#endif

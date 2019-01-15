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

class Kriging {

private:
  std::shared_ptr<Eigen::MatrixXd> Xpred;
  std::shared_ptr<Eigen::MatrixXd> Xobs;
  std::shared_ptr<Eigen::MatrixXd> invV;
  std::shared_ptr<Eigen::MatrixXi> flowMat;
  std::shared_ptr<Eigen::MatrixXd> weightMat;
  std::shared_ptr<Eigen::MatrixXd> distHydroOP;
  std::shared_ptr<Eigen::MatrixXd> distHydroPO;
  std::shared_ptr<Eigen::MatrixXd> distGeo;

  std::unique_ptr<TailUpModel> tailUpModel;
  std::unique_ptr<TailDownModel> tailDownModel;
  std::unique_ptr<EuclideanModel> euclidModel;

  double parsill = 0.0;
  int nModels;
  bool useNugget;
  std::shared_ptr<Eigen::VectorXd> theta;
  std::shared_ptr<Eigen::VectorXd> z;

  Eigen::MatrixXd Vpred;
  Eigen::MatrixXd predData;
  Eigen::MatrixXd XV;
  Eigen::MatrixXd invXVX;

  unsigned int nObs;
  unsigned int nPred;
  unsigned int p;

public:
  Kriging() =  default;
  Kriging(const Eigen::MatrixXd& dMatPred, const Eigen::MatrixXd& dMatObs, const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Nop, const Eigen::MatrixXd& Npo, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat,
    const Eigen::VectorXd& param, const Eigen::VectorXd& y, std::unique_ptr<TailUpModel>& tailup_ptr,
    std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, int nMod, bool useNugg);

  const Eigen::MatrixXd& getPredictions() const {return predData;};
  void predict();

};

#endif

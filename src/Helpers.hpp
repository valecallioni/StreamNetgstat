#include <iostream>
#include <vector>
#include <list>
#include <Eigen/Dense>
#include <cmath>
#include "Network.hpp"

class Network;

namespace helpers {

  double geoDist(const double x11, const double x12, const double x21, const double x22);

  double meanVec(const Eigen::VectorXd& v);
  double meanMat(const Eigen::MatrixXd& m);

  bool operandPair(std::pair<double, Eigen::VectorXd> p1, std::pair<double, Eigen::VectorXd> p2);

  Eigen::MatrixXd geoDistBetweenNets(const Network& net1, const Network& net2);
  std::vector<Eigen::MatrixXd> createDistMatricesObs(const std::vector<Network>& net);
  std::vector<Eigen::MatrixXd> createDistMatricesPred(const std::vector<Network>& net);
}

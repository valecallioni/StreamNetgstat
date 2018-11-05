#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>
#include <cmath>
#include "Network.hpp"

class Network;

namespace helpers {

  double geoDist(const double x11, const double x12, const double x21, const double x22);

  bool operandPair(std::pair<double, Eigen::VectorXd> p1, std::pair<double, Eigen::VectorXd> p2);

  void pointsStorage(std::vector<std::map<unsigned int, std::string>>& segmentsMaps, Eigen::MatrixXd& pointsMat, std::vector<std::vector<Point>>& storage);

  std::vector<Eigen::MatrixXd> returnBlockMatrices(const Points& p);
  Eigen::MatrixXd geoDistBetweenNets(const Points& p1, const Points& p2);

  std::vector<Eigen::MatrixXd> createDistMatrices(const std::string& type, const std::vector<Network>& net, unsigned int nTot);
  std::vector<Eigen::MatrixXd> createDistMatricesOP(const std::vector<Network>& net, unsigned int nObs, unsigned int nPred);

}

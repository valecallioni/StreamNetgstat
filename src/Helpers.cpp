#include "Helpers.hpp"

double helpers::geoDist(const double x11, const double x12, const double x21, const double x22) {
  return sqrt((x11-x21)*(x11-x21) + (x12-x22)*(x12-x22));
}

double helpers::meanVec(const Eigen::VectorXd& v){
  return v.sum()/v.size();
}

double helpers::meanMat(const Eigen::MatrixXd& m){
  return m.sum()/(m.rows()*m.cols());
}

bool helpers::operandPair(std::pair<double, Eigen::VectorXd> p1, std::pair<double, Eigen::VectorXd> p2){
  return (p1.first < p2.first);
}

Eigen::MatrixXd helpers::geoDistBetweenNets(const Network& net1, const Network& net2){
  Eigen::MatrixXd res(net1.getNObs(), net2.getNObs());
  for (unsigned int i=0; i<net1.getNObs(); i++){
    for (unsigned int j=0; j<net2.getNObs(); j++){
      res(i,j) = helpers::geoDist(net1.getObsPoints()[i].getX1(), net1.getObsPoints()[i].getX2(), net2.getObsPoints()[j].getX1(), net2.getObsPoints()[j].getX2());
    }
  }
  return res;
}

std::vector<Eigen::MatrixXd> helpers::createDistMatricesObs(const std::vector<Network>& net){
  Eigen::MatrixXd mat;
  unsigned int n(0);
  for (int k=0; k<net.size(); k++)
    n += net[k].getNObs();
  mat.resize(n,n);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(3,mat);
  unsigned int count(0);
  for (int k=0; k<net.size(); k++){
    res[0].block(count, count, net[k].getNObs(), net[k].getNObs()) = net[k].getFlowMatOO().cast<double>();
    res[1].block(count, count, net[k].getNObs(), net[k].getNObs()) = net[k].getNoo();
    res[2].block(count, count, net[k].getNObs(), net[k].getNObs()) = net[k].getDistGeoOO();
    if (k < net.size()-1){
      res[2].block(count + net[k].getNObs(), count, net[k+1].getNObs(), net[k].getNObs()) = helpers::geoDistBetweenNets(net[k], net[k+1]);
      res[2].block(count, count + net[k].getNObs(), net[k].getNObs(), net[k+1].getNObs()) = res[2].block(count + net[k].getNObs(), count, net[k].getNObs(), net[k].getNObs()).transpose();
    }
    count += net[k].getNObs();
  }
  return res;
}

std::vector<Eigen::MatrixXd> helpers::createDistMatricesPred(const std::vector<Network>& net){
  Eigen::MatrixXd mat;
  unsigned int n(0);
  for (int k=0; k<net.size(); k++)
    n += net[k].getNPred();
  mat.resize(n,n);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(3,mat);
  unsigned int count(0);
  for (int k=0; k<net.size(); k++){
    res[0].block(count, count, net[k].getNPred(), net[k].getNPred()) = net[k].getFlowMatPP().cast<double>();
    res[1].block(count, count, net[k].getNPred(), net[k].getNPred()) = net[k].getNpp();
    res[2].block(count, count, net[k].getNPred(), net[k].getNPred()) = net[k].getDistGeoPP();
    if (k < net.size()-1){
      res[2].block(count + net[k].getNPred(), count, net[k+1].getNPred(), net[k].getNPred()) = helpers::geoDistBetweenNets(net[k], net[k+1]);
      res[2].block(count, count + net[k].getNPred(), net[k].getNPred(), net[k+1].getNPred()) = res[2].block(count + net[k].getNPred(), count, net[k].getNPred(), net[k].getNPred()).transpose();
    }
    count += net[k].getNObs();
  }
  return res;
}

std::vector<Eigen::MatrixXd> helpers::createDistMatricesOP(const std::vector<Network>& net){
  Eigen::MatrixXd mat;
  unsigned int nObs(0);
  unsigned int nPred(0);
  for (int k=0; k<net.size(); k++){
    nObs += net[k].getNObs();
    nPred += net[k].getNPred();
  }
  mat.resize(nObs,nPred);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(3,mat);
  unsigned int countRows(0);
  unsigned int countCols(0);
  for (int k=0; k<net.size(); k++){
    res[0].block(countRows, countCols, net[k].getNObs(), net[k].getNPred()) = net[k].getFlowMatOP().cast<double>();
    res[1].block(countRows, countCols, net[k].getNObs(), net[k].getNPred()) = net[k].getNop();
    res[2].block(countRows, countCols, net[k].getNObs(), net[k].getNPred()) = net[k].getDistGeoOP();
    if (k < net.size()-1){
      res[2].block(countRows + net[k].getNObs(), countCols, net[k+1].getNObs(), net[k].getNPred()) = helpers::geoDistBetweenNets(net[k], net[k+1]);
      res[2].block(count, count + net[k].getNPred(), net[k].getNPred(), net[k+1].getNPred()) = res[2].block(count + net[k].getNPred(), count, net[k].getNPred(), net[k].getNPred()).transpose();
    }
    count += net[k].getNObs();
  }
  return res;
}

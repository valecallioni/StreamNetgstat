#include "Helpers.hpp"

double helpers::geoDist(const double x11, const double x12, const double x21, const double x22) {
  return std::sqrt((x11-x21)*(x11-x21) + (x12-x22)*(x12-x22));
}

bool helpers::operandPair(std::pair<double, Eigen::VectorXd> p1, std::pair<double, Eigen::VectorXd> p2){
  return (p1.first < p2.first);
}

void helpers::pointsStorage(std::vector<std::map<unsigned int, std::string>>& segmentsMaps, Eigen::MatrixXd& pointsMat, std::vector<std::vector<Point>>& storage){
  unsigned int j = 0;
  for (unsigned int k=0; k<storage.size(); k++){
      unsigned int currentNet = k+1;
      unsigned int i = 0;
      std::vector<Point> points;
      while (j<pointsMat.rows() && pointsMat(j,0)==currentNet){
        Point p;
        p.setPid(i);
        p.setRid(pointsMat(j,1));
        p.setDistUpstream(pointsMat(j,2));
        p.setCoordinates(pointsMat(j,3), pointsMat(j,4));
        auto it = segmentsMaps[k].find(pointsMat(j,1));
        p.setBinaryID(it->second);
        points.push_back(p);
        i++;
        j++;
      }
      storage[k] = points;
  }
}

std::vector<Eigen::MatrixXd> helpers::returnBlockMatrices(const Points& p){
  Eigen::MatrixXd mat;
  unsigned int n(p.getN());
  mat.resize(n,n);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(3,mat);
  res[0] = p.getFlowMat().cast<double>();
  res[1] = p.getDistHydro();
  res[2] = p.getDistGeo();
  return res;
}

Eigen::MatrixXd helpers::geoDistBetweenNets(const Points& p1, const Points& p2){
  Eigen::MatrixXd res(p1.getN(), p2.getN());
  res.fill(0.0);
  for (unsigned int i=0; i<p1.getN(); i++){
    for (unsigned int j=0; j<p2.getN(); j++){
      res(i,j) = helpers::geoDist(p1.getPoints()[i].getX1(), p1.getPoints()[i].getX2(), p2.getPoints()[j].getX1(), p2.getPoints()[j].getX2());
    }
  }
  return res;
}

std::vector<Eigen::MatrixXd> helpers::createDistMatrices(const std::string& type, const std::vector<Network>& net, unsigned int nTot){
  Eigen::MatrixXd mat;
  mat.resize(nTot,nTot);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(3,mat);
  unsigned int count(0);

  std::vector<Eigen::MatrixXd> tmp;
  Eigen::MatrixXd geoDist;
  Points p1;
  Points p2;
  for (int k=0; k<net.size(); k++){
    // Chiamo una funzione template che mi restituisce i 3 sottoblocchi
    // Qui capisco se i punti che sto studiando sono obs o pred!!!
    if (type == "obs") p1 = net[k].getObsPoints();
    if (type == "pred") p1 = net[k].getPredPoints();;
    tmp = helpers::returnBlockMatrices(p1);
    res[0].block(count, count, tmp[0].rows(), tmp[0].cols()) = tmp[0];
    res[1].block(count, count, tmp[1].rows(), tmp[1].cols()) = tmp[1];
    res[2].block(count, count, tmp[2].rows(), tmp[2].cols()) = tmp[2];
    if (k < net.size()-1){
      if (type == "obs") p2 = net[k+1].getObsPoints();
      if (type == "pred") p2 = net[k+1].getPredPoints();
      geoDist = helpers::geoDistBetweenNets(p1, p2);
      res[2].block(count, count + p1.getN(), geoDist.rows(), geoDist.cols()) = geoDist;
      res[2].block(count + p1.getN(), count, geoDist.cols(), geoDist.rows()) = geoDist.transpose();
    }
    count += p1.getN();
  }
  return res;
}

std::vector<Eigen::MatrixXd> helpers::createDistMatricesOP(const std::vector<Network>& net, unsigned int nObs, unsigned int nPred){
  Eigen::MatrixXd mat;
  mat.resize(nObs,nPred);
  mat.fill(0.0);
  std::vector<Eigen::MatrixXd> res(4,mat);
  unsigned int countObs(0);
  unsigned int countPred(0);
  std::vector<Eigen::MatrixXd> tmp;
  Eigen::MatrixXd geoDist;
  Points p1;
  Points p2;
  for (int k=0; k<net.size(); k++){
    res[0].block(countObs, countPred, net[k].getNObs(), net[k].getNPred()) = net[k].getFlowMatOP().cast<double>();
    res[1].block(countObs, countPred, net[k].getNObs(), net[k].getNPred()) = net[k].getDistHydroOP();
    res[2].block(countObs, countPred, net[k].getNObs(), net[k].getNPred()) = net[k].getDistHydroPO().transpose();
    res[3].block(countObs, countPred, net[k].getNObs(), net[k].getNPred()) = net[k].getDistGeoOP();
    if (k < net.size()-1){
      p1 = net[k].getObsPoints();
      p2 = net[k+1].getPredPoints();
      geoDist = helpers::geoDistBetweenNets(p1, p2);
      res[3].block(countObs, countPred + net[k].getNPred(), geoDist.rows(), geoDist.cols()) = geoDist;
      p1 = net[k+1].getObsPoints();
      p2 = net[k].getPredPoints();
      geoDist = helpers::geoDistBetweenNets(p1, p2);
      res[3].block(countObs + net[k].getNObs(), countPred, geoDist.rows(), geoDist.cols()) = geoDist;
    }
    countObs += net[k].getNObs();
    countPred += net[k].getNPred();
  }
  return res;
}

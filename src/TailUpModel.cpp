#include "TailUpModel.hpp"

Eigen::MatrixXd TailUpModel::computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMat) {
  unsigned int n(weightMat.rows());
  Eigen::MatrixXd res(n,n);
  res.fill(0.0);
  double h;
  for (unsigned int i=0; i<n; i++){
    for (unsigned int j=i; j<n; j++){
      if (i == j){
        res(i,j) = sigma2;
      }
      else if (distMat(i,j) == 0.0 && distMat(j,i) != 0.0){
        h = distMat(j,i);
        res(i,j) = weightMat(i,j) * computeCov(h);
        res(j,i) = res(i,j);
      }
      else if (distMat(i,j) != 0.0 && distMat(j,i) == 0.0){
        h = distMat(i,j);
        res(i,j) = weightMat(i,j) * computeCov(h);
        res(j,i) = res(i,j);
      }
      else if (distMat(i,j) == 0.0 && distMat(j,i) == 0.0 && weightMat(i,j) == 1) {
        res(i,j) = sigma2;
        res(j,i) = sigma2;
      }
    }
  }
  return res;
}

Eigen::MatrixXd TailUpModel::computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMatOP, const Eigen::MatrixXd& distMatPO) {
  unsigned int n(distMatOP.rows());
  unsigned int m(distMatOP.cols());
  Eigen::MatrixXd res(n,m);
  res.fill(0.0);
  double h;
  for (unsigned int i=0; i<n; i++){
    for (unsigned int j=0; j<m; j++){
      if (distMatOP(i,j) == 0.0 && distMatPO(j,i) != 0.0){
        h = distMatPO(j,i);
        res(i,j) = weightMat(i,j) * computeCov(h);
      }
      else if (distMatOP(i,j) != 0.0 && distMatPO(j,i) == 0.0){
        h = distMatOP(i,j);
        res(i,j) = weightMat(i,j) * computeCov(h);
      }
    }
  }
  return res;
}

double LinearWithSillTU::computeCov(double h) {
  if (h <= alpha){
    return sigma2 * (1.0 - h/alpha);
  }
  else {
    return 0.0;
  }
}

double SphericalTU::computeCov(double h) {
  if (h <= alpha){
    return sigma2 * (1.0 - (3.0/2.0) * (h/alpha) + 0.5 * pow(h,3)/pow(alpha,3));
  }
  else {
    return 0.0;
  }
}

double ExponentialTU::computeCov(double h) {
  return sigma2 * exp(-3.0*h/alpha);
}

double MariahTU::computeCov(double h) {
  if (h <= 1e-6){
    return sigma2;
  }
  else {
    return sigma2 * (log(90.0 * h/alpha + 1.0)/(90.0 * h/alpha));
  }
}

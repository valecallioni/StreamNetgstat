#include "TailDownModel.hpp"

Eigen::MatrixXd TailDownModel::computeMatCov(const Eigen::MatrixXi& flowMat, const Eigen::MatrixXd& distMat) {
  unsigned int n(flowMat.rows());
  Eigen::MatrixXd res(n,n);
  res.fill(0.0);
  double h;
  for (unsigned int i=0; i<n; i++){
    for (unsigned int j=i; j<n; j++){
      if (flowMat(i,j) == 1){
        if (distMat(i,j) == 0.0 && distMat(j,i) != 0.0){
          h = distMat(j,i);
          res(i,j) = computeCov(h);
          res(j,i) = res(i,j);
        }
        else if (distMat(i,j) != 0.0 && distMat(j,i) == 0.0){
          h = distMat(i,j);
          res(i,j) = computeCov(h);
          res(j,i) = res(i,j);
        }
        else if (distMat(i,j) == 0.0 && distMat(j,i) == 0.0){
          res(i,j) = sigma2;
          res(j,i) = sigma2;
        }
      }
      else {
        if (i == j){
          res(i,j) = sigma2;
        }
        else {
          double a(0.0);
          double b(0.0);
          if (distMat(i,j) < distMat(j,i)){
            a = distMat(i,j);
            b = distMat(j,i);
          }
          else {
            b = distMat(i,j);
            a = distMat(j,i);
          }
          if (a != 0.0 && b != 0.0){
            res(i,j) = computeCov(a, b);
            res(j,i) = res(i,j);
          }
        }
      }
    }
  }
  return res;
}

Eigen::MatrixXd TailDownModel::computeMatCov(const Eigen::MatrixXi& flowMat, const Eigen::MatrixXd& distMatOP, const Eigen::MatrixXd& distMatPO) {
  unsigned int n(distMatOP.rows());
  unsigned int m(distMatOP.cols());
  Eigen::MatrixXd res(n,m);
  res.fill(0.0);
  double h;
  for (unsigned int i=0; i<n; i++){
    for (unsigned int j=0; j<m; j++){
      if (flowMat(i,j) == 1){
        if (distMatOP(i,j) == 0.0 && distMatPO(j,i) != 0.0){
          h = distMatPO(j,i);
          res(i,j) = computeCov(h);
        }
        else if (distMatOP(i,j) != 0.0 && distMatPO(j,i) == 0.0){
          h = distMatOP(i,j);
          res(i,j) = computeCov(h);
        }
      }
      else {
        double a(0.0);
        double b(0.0);
        if (distMatOP(i,j) < distMatPO(j,i)){
          a = distMatOP(i,j);
          b = distMatPO(j,i);
        }
        else {
          b = distMatOP(i,j);
          a = distMatPO(j,i);
        }
        if (a != 0.0 && b != 0.0){
          res(i,j) = computeCov(a, b);
        }
      }
    }
  }
  return res;
}

double LinearWithSillTD::computeCov(double a, double b) {
  if (b <= alpha){
    return sigma2 * (1.0 - b/alpha);
  }
  else {
    return 0.0;
  }
}

double LinearWithSillTD::computeCov(double h) {
  if (h <= alpha){
    return sigma2 * (1.0 - h/alpha);
  }
  else {
    return 0.0;
  }
}

double SphericalTD::computeCov(double a, double b) {
  if (b <= alpha){
    return sigma2 * (1.0 - (3.0/2.0) * (a/alpha) + 0.5 * b/alpha) * pow(1.0 - b/alpha, 2);
  }
  else {
    return 0.0;
  }
}

double SphericalTD::computeCov(double h) {
  if (h <= alpha){
    return sigma2 * (1.0 - (3.0/2.0) * (h/alpha) + 0.5 * pow(h,3)/pow(alpha,3));
  }
  else {
    return 0.0;
  }
}

double ExponentialTD::computeCov(double a, double b) {
  return sigma2 * exp(-3.0*(a+b)/alpha);
}

double ExponentialTD::computeCov(double h) {
  return sigma2 * exp(-3.0*h/alpha);
}

double MariahTD::computeCov(double a, double b) {
  if (b - a <= 1e-6){
    return sigma2 / (90.0 * a/alpha + 1);
  }
  else {
    return sigma2 * ((log(90.0 * a/alpha + 1.0) - log(90.0 * b/alpha + 1.0))/(90.0 * (a-b)/alpha));
  }
}

double MariahTD::computeCov(double h) {
  if (h <= 1e-6){
    return sigma2;
  }
  else {
    return sigma2 * (log(90.0 * h/alpha + 1.0)/(90.0 * h/alpha));
  }
}

#include "EuclideanModel.hpp"

Eigen::MatrixXd EuclideanModel::computeMatCov(const Eigen::MatrixXd& distGeo) {
  unsigned int n(distGeo.rows());
  unsigned int m(distGeo.cols());
  Eigen::MatrixXd res(n,m);
  res.fill(0.0);
  for (unsigned int i=0; i<n; i++){
    for (unsigned int j=0; j<m; j++){
      if (i == j)
        res(i,j) = sigma2;
      else
        res(i,j) = computeCov(distGeo(i,j));
      if (res(i,j)<0.0)
        std::cout << "Cov Euclid negative between points " << i << " and " << j << std::endl;
    }
  }
  return res;
}

double CauchyEU::computeCov(double d) {
  return sigma2 * pow(1.0 + 4.4*(pow(d/alpha, 2)), -1);
}

double SphericalEU::computeCov(double d) {
  if (d <= alpha){
      return sigma2 * (1.0 - (3.0/2.0)*(d/alpha) + 0.5*pow(d/alpha, 3));
  }
  else {
    return 0.0;
  }
}

double ExponentialEU::computeCov(double d) {
  return sigma2 * exp(-3.0*d/alpha);
}

double GaussianEU::computeCov(double d) {
  return sigma2 * exp(-3.0*pow(d/alpha, 2));
}

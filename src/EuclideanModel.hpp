#ifndef PACSPROJECT_EUCLIDEANMODEL_HPP
#define PACSPROJECT_EUCLIDEANMODEL_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

class EuclideanModel {

protected:
  double sigma2;
  double alpha;

public:
  EuclideanModel() = default;
  EuclideanModel(double s, double a): sigma2(s), alpha(a){};
  virtual double computeCov(double h) = 0;
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& distGeo);

  void setSigma2(double s) {sigma2 = s;};
  void setAlpha(double a) {alpha = a;};
  double getSigma2() const {return sigma2;};
  double getAlpha() const {return alpha;};
};

class CauchyEU : public EuclideanModel{
public:
  double computeCov(double d) override;
};

class SphericalEU : public EuclideanModel{
public:
  double computeCov(double d) override;
};

class ExponentialEU : public EuclideanModel{
public:
  double computeCov(double d) override;
};

class GaussianEU : public EuclideanModel{
public:
  double computeCov(double d) override;
};

#endif

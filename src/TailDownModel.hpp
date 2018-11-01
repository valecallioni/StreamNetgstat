#ifndef PACSPROJECT_TAILDOWNMODEL_HPP
#define PACSPROJECT_TAILDOWNMODEL_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

class TailDownModel {

protected:
  double sigma2;
  double alpha;

public:
  TailDownModel() = default;
  TailDownModel(double s, double a): sigma2(s), alpha(a){};
  virtual double computeCov(double, double) = 0;
  virtual double computeCov(double h) = 0;
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXi& flowMat, const Eigen::MatrixXd& distMat);

  void setSigma2(double s) {sigma2 = s;};
  void setAlpha(double a) {alpha = a;};
  double getSigma2() const {return sigma2;};
  double getAlpha() const {return alpha;};
};

class LinearWithSillTD : public TailDownModel{
public:
  double computeCov(double, double) override;
  double computeCov(double) override;
};

class SphericalTD : public TailDownModel{
public:
  double computeCov(double, double) override;
  double computeCov(double) override;
};

class ExponentialTD : public TailDownModel{
public:
  double computeCov(double, double) override;
  double computeCov(double) override;
};

class MariahTD : public TailDownModel{
public:
  double computeCov(double, double) override;
  double computeCov(double) override;
};

#endif

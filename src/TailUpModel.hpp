#ifndef PACSPROJECT_TAILUPMODEL_HPP
#define PACSPROJECT_TAILUPMODEL_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

class TailUpModel {

protected:
  double sigma2;
  double alpha;

public:
  TailUpModel() = default;
  TailUpModel(double s, double a): sigma2(s), alpha(a){};
  virtual double computeCov(double h) = 0;
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMat);
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMatOP, const Eigen::MatrixXd& distMatPO);

  void setSigma2(double s) {sigma2 = s;};
  void setAlpha(double a) {alpha = a;};
  double getSigma2() const {return sigma2;};
  double getAlpha() const {return alpha;};
};

class LinearWithSillTU : public TailUpModel{
public:
  double computeCov(double h) override;
};

class SphericalTU : public TailUpModel{
public:
  double computeCov(double h) override;
};

class ExponentialTU : public TailUpModel{
public:
  double computeCov(double h) override;
};

class MariahTU : public TailUpModel{
public:
  double computeCov(double h) override;
};

#endif

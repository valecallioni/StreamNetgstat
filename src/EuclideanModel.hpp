#ifndef PACSPROJECT_EUCLIDEANMODEL_HPP
#define PACSPROJECT_EUCLIDEANMODEL_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

/**
* Abstract class for the Euclidean covariance model.
*/

class EuclideanModel {

protected:
  double sigma2; ///< parsill of the Euclidean model
  double alpha; ///< range of the Euclidean model

public:
  /**
  * Default constructor.
  */
  EuclideanModel() = default;

  /**
  * Default destructor.
  */
  virtual ~EuclideanModel() =  default;

  /**
  * Constructor.
  * @param s initial value of the parsill
  * @param a initial value of the range
  */
  EuclideanModel(double s, double a): sigma2(s), alpha(a){};

  /**
  * A pure virtual member function to compute the covariance between two points
  * @param d Euclidean distance between the points
  */
  virtual double computeCov(double d) = 0;

  /**
  * Computes the covariance matrix from a distance matrix
  * @param distGeo distance matrix
  * @return The covariance matrix
  */
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& distGeo);

  /**
  * @param s new value of the parsill
  */
  void setSigma2(double s) {sigma2 = s;};

  /**
  * @param a new value of the range
  */
  void setAlpha(double a) {alpha = a;};

  /**
  * @return the parsill
  */
  double getSigma2() const {return sigma2;};

  /**
  * @return the range
  */
  double getAlpha() const {return alpha;};
};


/**
* Child class for the Euclidean covariance model, implementing the Cauchy covariance function.
*/
class CauchyEU : public EuclideanModel{
public:
  /**
  * Default destructor.
  */
  ~CauchyEU() =  default;

  /**
  * Compute the covariance between two points
  * @param d Euclidean distance between the points
  */
  double computeCov(double d) override;
};

/**
* Child class for the Euclidean covariance model, implementing the spherical covariance function.
*/
class SphericalEU : public EuclideanModel{
public:
  /**
  * Default destructor.
  */
  ~SphericalEU() =  default;

  /**
  * Compute the covariance between two points
  * @param d Euclidean distance between the points
  */
  double computeCov(double d) override;
};

/**
* Child class for the Euclidean covariance model, implementing the exponential covariance function.
*/
class ExponentialEU : public EuclideanModel{
public:
  /**
  * Default destructor.
  */
  ~ExponentialEU() =  default;

  /**
  * Compute the covariance between two points
  * @param d Euclidean distance between the points
  */
  double computeCov(double d) override;
};

/**
* Child class for the Euclidean covariance model, implementing the Gaussian covariance function.
*/
class GaussianEU : public EuclideanModel{
public:
  /**
  * Default destructor.
  */
  ~GaussianEU() =  default;

  /**
  * Compute the covariance between two points
  * @param d Euclidean distance between the points
  */
  double computeCov(double d) override;
};

#endif

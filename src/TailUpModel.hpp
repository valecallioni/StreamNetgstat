#ifndef PACSPROJECT_TAILUPMODEL_HPP
#define PACSPROJECT_TAILUPMODEL_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

/*! \file
* Abstract class for the tail-up covariance model.
*/

class TailUpModel {

protected:
  double sigma2; ///< parsill of the tail-down model
  double alpha; ///< range of the tail-down model

public:
  /**
  * Default constructor.
  */
  TailUpModel() = default;

  /**
  * Default destructor.
  */
  virtual ~TailUpModel() =  default;

  /**
  * Constructor.
  * @param s initial value of the parsill
  * @param a initial value of the range
  */
  TailUpModel(double s, double a): sigma2(s), alpha(a){};

  /**
  * A pure virtual member function to compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  virtual double computeCov(double h) = 0;

  /**
  * Computes the covariance matrix from a distance matrix and a weight matrix
  * @param weightMat matrix whose elements are the weights associated to pair of points
  * @param distMat hydrological distance matrix
  * @return The covariance matrix
  */
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMat);

  /**
  * Computes the covariance matrix from a distance matrix and a weight matrix between two group of points
  * @param weightMat matrix whose elements are the weights associated to pair of points
  * @param distMatOP hydrological distance matrix, where rows represent the first group of points (in general, the observed), and
  * columns represent the second group (in general, the prediction points)
  * @param distMatPO hydrological distance matrix, where rows represent the second group of points (in general, the prediction points), and
  * columns represent the first group (in general, the observed points)
  * @return The covariance matrix
  */
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXd& weightMat, const Eigen::MatrixXd& distMatOP, const Eigen::MatrixXd& distMatPO);

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
* Child class for the tail-up covariance model, implementing the linear-with-sill covariance function.
*/
class LinearWithSillTU : public TailUpModel{
public:
  /**
  * Default destructor.
  */
  virtual ~LinearWithSillTU() =  default;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double h) override;
};

/**
* Child class for the tail-up covariance model, implementing the spherical covariance function.
*/
class SphericalTU : public TailUpModel{
public:
  /**
  * Default destructor.
  */
  virtual ~SphericalTU() =  default;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double h) override;
};

/**
* Child class for the tail-up covariance model, implementing the exponential covariance function.
*/
class ExponentialTU : public TailUpModel{
public:
  /**
  * Default destructor.
  */
  virtual ~ExponentialTU() =  default;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double h) override;
};

/**
* Child class for the tail-up covariance model, implementing the Mariah covariance function.
*/
class MariahTU : public TailUpModel{
public:
  /**
  * Default destructor.
  */
  virtual ~MariahTU() =  default;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double h) override;
};

#endif

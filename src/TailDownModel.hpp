#ifndef PACSPROJECT_TAILDOWNMODEL_HPP
#define PACSPROJECT_TAILDOWNMODEL_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

/**
* Abstract class for the tail-down covariance model.
*/

class TailDownModel {

protected:
  double sigma2; ///< parsill of the tail-down model
  double alpha; ///< range of the tail-down model

public:
  /**
  * Default constructor.
  */
  TailDownModel() = default;

  /**
  * Default destructor.
  */
  virtual ~TailDownModel() =  default;

  /**
  * Constructor.
  * @param s initial value of the parsill
  * @param a initial value of the range
  */
  TailDownModel(double s, double a): sigma2(s), alpha(a){};

  /**
  * A pure virtual member function to compute the covariance between two flow-unconnected points
  * @param a shorter hydrological distance between the points and the common downstream juction
  * @param b greater hydrological distance between the points and the common downstream juction
  */
  virtual double computeCov(double a, double b) = 0;

  /**
  * A pure virtual member function to compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  virtual double computeCov(double h) = 0;

  /**
  * Computes the covariance matrix from a distance matrix and a flow-connection matrix
  * @param flowMat flow-connection/unconnection binary matrix
  * @param distMat hydrological distance matrix
  * @return The covariance matrix
  */
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXi& flowMat, const Eigen::MatrixXd& distMat);

  /**
  * Computes the covariance matrix from a distance matrix and a flow-connection matrix between two group of points
  * @param flowMat flow-connection/unconnection binary matrix
  * @param distMatOP hydrological distance matrix, where rows represent the first group of points (in general, the observed), and
  * columns represent the second group (in general, the prediction points)
  * @param distMatPO hydrological distance matrix, where rows represent the second group of points (in general, the prediction points), and
  * columns represent the first group (in general, the observed points)
  * @return The covariance matrix
  */
  Eigen::MatrixXd computeMatCov(const Eigen::MatrixXi& flowMat, const Eigen::MatrixXd& distMatOP, const Eigen::MatrixXd& distMatPO);

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
* Child class for the tail-down covariance model, implementing the linear-with-sill covariance function.
*/
class LinearWithSillTD : public TailDownModel{
public:
  /**
  * Default destructor.
  */
  ~LinearWithSillTD() =  default;

  /**
  * Compute the covariance between two flow-unconnected points
  * @param a shorter hydrological distance between the points and the common downstream juction
  * @param b greater hydrological distance between the points and the common downstream juction
  */
  double computeCov(double a, double b) override;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double) override;
};

/**
* Child class for the tail-down covariance model, implementing the spherical covariance function.
*/
class SphericalTD : public TailDownModel{
public:
  /**
  * Default destructor.
  */
  ~SphericalTD() =  default;

  /**
  * Compute the covariance between two flow-unconnected points
  * @param a shorter hydrological distance between the points and the common downstream juction
  * @param b greater hydrological distance between the points and the common downstream juction
  */
  double computeCov(double, double) override;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double) override;
};

/**
* Child class for the tail-down covariance model, implementing the exponential covariance function.
*/
class ExponentialTD : public TailDownModel{
public:
  /**
  * Default destructor.
  */
  ~ExponentialTD() =  default;

  /**
  * Compute the covariance between two flow-unconnected points
  * @param a shorter hydrological distance between the points and the common downstream juction
  * @param b greater hydrological distance between the points and the common downstream juction
  */
  double computeCov(double, double) override;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double) override;
};

/**
* Child class for the tail-down covariance model, implementing the Mariah covariance function.
*/
class MariahTD : public TailDownModel{
public:
  /**
  * Default destructor.
  */
  ~MariahTD() =  default;

  /**
  * Compute the covariance between two flow-unconnected points
  * @param a shorter hydrological distance between the points and the common downstream juction
  * @param b greater hydrological distance between the points and the common downstream juction
  */
  double computeCov(double, double) override;

  /**
  * Compute the covariance between two flow-connected points
  * @param h hydrological distance between the points
  */
  double computeCov(double) override;
};

#endif

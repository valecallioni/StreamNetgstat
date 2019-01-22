#ifndef PACSPROJECT_POINTS_HPP
#define PACSPROJECT_POINTS_HPP

#include <iostream>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <stdexcept>
#include "Point.hpp"
#include "StreamSegment.hpp"

/**
* Points class, representing a group of points, observed or prediction ones.
* Each group has three matrices used in computing the covariance matrix: the
* flow-connection/unconnection binary matrix, the hydrological and Euclidean distance matrices.
*/

class Points {

private:
  std::vector<Point> points; ///< vector with the Point-objects of the group
  unsigned int n = 0; ///< number of points
  Eigen::MatrixXi flowMat; ///< flow-connection/unconnection binary matrix
  Eigen::MatrixXd distHydro; ///< hydrological distance matrix
  Eigen::MatrixXd distGeo; ///< Euclidean distance matrix

public:
  /**
  * Default constructor.
  */
  Points() = default;

  /**
  * Constructor.
  * @param p vector containing the Point-objects
  */
  Points(const std::vector<Point>& p);

  /**
  * @return the number of points within the group
  */
  const unsigned int getN() const {return n;};

  /**
  * @return the flow-connection/unconnection binary matrix
  */
  const Eigen::MatrixXi& getFlowMat() const {return flowMat;};

  /**
  * @return the hydrological distance matrix
  */
  const Eigen::MatrixXd& getDistHydro() const {return distHydro;};

  /**
  * @return the Euclidean distance matrix
  */
  const Eigen::MatrixXd& getDistGeo() const {return distGeo;};

  /**
  * @return the vector containing the Point-objects
  */
  const std::vector<Point>& getPoints() const {return points;};

  /**
  * @param p vector containing the Point-objects
  */
  void setPoints(const std::vector<Point>& p);

  /**
  * Computes the distance matrices within the group of points
  * @param geo boolean, if TRUE then also the Euclidean distance matrix is computed
  * @param segments map that associates to a binaryID the corresponding stream segment, used to
  * determine the mutual position of the points along the network
  */
  void computeDistances(bool geo, const std::map<std::string, StreamSegment>& segments);

  /**
  * Sets the distance matrices
  * @param matrices vector of matrices (flow-connection/unconnection, hydrological distance and, not necessarily, Euclidean distance matrix)
  */
  void setDistances(const std::vector<Eigen::MatrixXd>& matrices);

};


#endif

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>
#include <cmath>
#include "Network.hpp"

class Network;

/**
* Helpers functions.
*/

namespace helpers {

  /**
  * Compute the Euclidean distance between two points.
  * @param x11 first coordinate of the first point
  * @param x12 second coordinate of the first point
  * @param x21 first coordinate of the second point
  * @param x22 second coordinate of the second point
  * @return the Euclidean distance between the points
  */
  double geoDist(const double x11, const double x12, const double x21, const double x22);

  /**
  * Compare pairs of <double, Eigen::VectorXd>, looking just at the double
  * @param p1 first pair
  * @param p2 second pair
  * @return TRUE if the first pair is "less than" the second pair
  */
  bool operandPair(std::pair<double, Eigen::VectorXd> p1, std::pair<double, Eigen::VectorXd> p2);

  /**
  * Store a group of points, with all their attributes, belonging to different networks
  * @param segmentsMaps vector of maps, one per each network, that associates to each segmentID the StreamSegment object
  * @param pointsMat matrix containing the attributes of the points
  * @param storage vector of vectors, one per each network, of Point objects, to be filled
  */
  void pointsStorage(std::vector<std::map<unsigned int, std::string>>& segmentsMaps, Eigen::MatrixXd& pointsMat, std::vector<std::vector<Point>>& storage);

  /**
  * Return the distance matrices and flow-connection/unconnection binary matrix of a group of points
  * @param p Points object
  * @return vector of matrices
  */
  std::vector<Eigen::MatrixXd> returnBlockMatrices(const Points& p);

  /**
  * Compute the Euclidean distance between points belonging to different networks
  * @param p1 Points object of the first network
  * @param p2 Points object of the second network
  * @return Euclidean distance matrix
  */
  Eigen::MatrixXd geoDistBetweenNets(const Points& p1, const Points& p2);

  /**
  * Build the block matrices regarding the hydrological and, if necessary, Euclidean distances and flow-connection/unconnection between observed points
  * @param geo boolean indicating if the Euclidean distance matrix has to be computed
  * @param net vector of networks whose general distance matrices have to be computed
  * @param nTot total number of observed points in the data set
  * @return vector of distance matrices
  */
  std::vector<Eigen::MatrixXd> createDistMatrices(bool geo, const std::vector<Network>& net, unsigned int nTot);

  /**
  * Build the block matrices regarding the hydrological and, if necessary, Euclidean distances and flow-connection/unconnection between observed and prediction points
  * @param geo boolean indicating if the Euclidean distance matrix has to be computed
  * @param net vector of networks whose general distance matrices have to be computed
  * @param nObs total number of observed points in the data set
  * @param nPred total number of prediction points in the data set
  * @return vector of distance matrices
  */
  std::vector<Eigen::MatrixXd> createDistMatricesOP(bool geo, const std::vector<Network>& net, unsigned int nObs, unsigned int nPred);
}

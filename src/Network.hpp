#ifndef PACSPROJECT_NETWORK_HPP
#define PACSPROJECT_NETWORK_HPP

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "Points.hpp"
#include "Helpers.hpp"

/**
* Network class, representing a single network of the data set.
*/

class Network{

private:
    unsigned int ID; ///< networkID

    Points obspoints; ///< group of observed points
    Points predpoints; ///< group of prediction points

    std::map<std::string, StreamSegment> segments; ///< map that associates to each stream segment of the network its binaryID

    unsigned int nObs; ///< number of observed points
    unsigned int nPred; ///< number of prediction points

    Eigen::MatrixXi flowMatOP; ///< flow-connection/unconnection binary matrix between observed and prediction points
    Eigen::MatrixXd distHydroOP; ///< hydrological distance matrix between observed and prediction points
    Eigen::MatrixXd distHydroPO; ///< hydrological distance matrix between prediction and observed points
    Eigen::MatrixXd distGeoOP; ///< Euclidean distance matrix between observed and prediction points

public:
    /**
    * Default constructor.
    */
    Network() = default;

    /**
    * Constructor.
    * @param id networkID
    * @param obs vector of Point objects, representing the observed points
    * @param pred vector of Point objects, representing the prediction points
    * @param seg vector of StreamSegment objects
    */
    Network(const unsigned int id, const std::vector<Point>& obs, const std::vector<Point>& pred, const std::vector<StreamSegment>& seg);

    /**
    * Constructor.
    * @param id networkID
    * @param obs vector of Point objects, representing the observed points
    * @param seg vector of StreamSegment objects
    */
    Network(const unsigned int id, const std::vector<Point>& obs, const std::vector<StreamSegment>& seg);

    /**
    * @param id networkID
    */
    void setID(const unsigned int id) {ID = id;};

    /**
    * @return the networkID
    */
    const unsigned int getID() const {return ID;};

    /**
    * @return the number of observed points
    */
    const unsigned int getNObs() const {return nObs;};

    /**
    * @return the number of prediction points
    */
    const unsigned int getNPred() const {return nPred;};

    /**
    * @return the flow-connection/unconnection binary matrix between observed points
    */
    const Eigen::MatrixXi& getFlowMatOO() const {return obspoints.getFlowMat();};

    /**
    * @return the flow-connection/unconnection binary matrix between observed points and prediction points
    */
    const Eigen::MatrixXi& getFlowMatOP() const {return flowMatOP;};

    /**
    * @return the flow-connection/unconnection binary matrix between prediction points
    */
    const Eigen::MatrixXi& getFlowMatPP() const {return predpoints.getFlowMat();};

    /**
    * @return the Euclidean distance matrix between observed points
    */
    const Eigen::MatrixXd& getDistGeoOO() const {return obspoints.getDistGeo();};

    /**
    * @return the Euclidean distance matrix between observed points and prediction points
    */
    const Eigen::MatrixXd& getDistGeoOP() const {return distGeoOP;};

    /**
    * @return the Euclidean distance matrix between prediction points
    */
    const Eigen::MatrixXd& getDistGeoPP() const {return predpoints.getDistGeo();};

    /**
    * @return the hydrological distance matrix between observed points
    */
    const Eigen::MatrixXd& getDistHydroOO() const {return obspoints.getDistHydro();};

    /**
    * @return the hydrological distance matrix between observed points and prediction points
    */
    const Eigen::MatrixXd& getDistHydroOP() const {return distHydroOP;};

    /**
    * @return the hydrological distance matrix between prediction points and observed points
    */
    const Eigen::MatrixXd& getDistHydroPO() const {return distHydroPO;};

    /**
    * @return the hydrological distance matrix between prediction points
    */
    const Eigen::MatrixXd& getDistHydroPP() const {return predpoints.getDistHydro();};

    /**
    * @return the group of observed points
    */
    const Points& getObsPoints() const {return obspoints;};

    /**
    * @return the group of prediction points
    */
    const Points& getPredPoints() const {return predpoints;};

    /**
    * Compute the distance matrices whose rows represent the observed points and whose columns representing the predicion points
    * @param geo boolean, indicating if the Euclidean distances have to be computed
    */
    void computeDistances(bool geo);

    /**
    * Set the distance matrices regarding the relationships between the observed points and
    * then compute the distance matrices whose rows represent the observed points and whose columns representing the predicion points
    * @param geo boolean, indicating if the Euclidean distances have to be computed
    * @param matrices vector of matrices regarding the group of observed points
    */
    void setDistPoints(bool geo, const std::vector<Eigen::MatrixXd>& matrices);

    /**
    * Printing function to visualize the attributes of the stream segment
    */
    void print() const;
};

#endif //PACSPROJECT_NETWORK_HPP

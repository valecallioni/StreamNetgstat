#ifndef PACSPROJECT_NETWORK_HPP
#define PACSPROJECT_NETWORK_HPP

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "Points.hpp"
#include "Helpers.hpp"

class Network{

private:
    unsigned int ID;

    Points obspoints;
    Points predpoints;

    std::map<std::string, StreamSegment> segments;

    unsigned int nObs;
    unsigned int nPred;

    Eigen::MatrixXi flowMatOP;
    Eigen::MatrixXd distHydroOP;
    Eigen::MatrixXd distHydroPO;
    Eigen::MatrixXd distGeoOP;

public:
    Network() = default;
    Network(const unsigned int id, const std::vector<Point>& obs, const std::vector<Point>& pred, const std::vector<StreamSegment>& seg);

    void setID(const unsigned int id) {ID = id;};

    const unsigned int getID() const {return ID;};
    const unsigned int getNObs() const {return nObs;};
    const unsigned int getNPred() const {return nPred;};
    const Eigen::MatrixXi& getFlowMatOO() const {return obspoints.getFlowMat();};
    const Eigen::MatrixXi& getFlowMatOP() const {return flowMatOP;};
    const Eigen::MatrixXi& getFlowMatPP() const {return predpoints.getFlowMat();};
    //const std::vector<Point>& getObsPoints() const {return obspoints.getPoints();};
    //const std::vector<Point>& getPredPoints() const {return predpoints.getPoints();};
    const Eigen::MatrixXd& getDistGeoOO() const {return obspoints.getDistGeo();};
    const Eigen::MatrixXd& getDistGeoOP() const {return distGeoOP;};
    const Eigen::MatrixXd& getDistGeoPP() const {return predpoints.getDistGeo();};
    const Eigen::MatrixXd& getDistHydroOO() const {return obspoints.getDistHydro();};
    const Eigen::MatrixXd& getDistHydroOP() const {return distHydroOP;};
    const Eigen::MatrixXd& getDistHydroPO() const {return distHydroPO;};
    const Eigen::MatrixXd& getDistHydroPP() const {return predpoints.getDistHydro();};
    const Points& getObsPoints() const {return obspoints;};
    const Points& getPredPoints() const {return predpoints;};


    void computeDistances();

    void print() const;
};

#endif //PACSPROJECT_NETWORK_HPP

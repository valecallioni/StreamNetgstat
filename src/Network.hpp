#ifndef PACSPROJECT_NETWORK_HPP
#define PACSPROJECT_NETWORK_HPP

#include <iostream>
#include <vector>
#include <list>
#include <Eigen/Dense>
#include <cmath>
#include "Coordinates.hpp"
#include "Point.hpp"
#include "StreamSegment.hpp"
#include "Helpers.hpp"

class Network{

private:
    unsigned int ID;

    std::vector<Point> obspoints;
    std::vector<Point> predpoints;

    std::map<std::string, StreamSegment> segments;

    unsigned int nObs;
    unsigned int nPred;

    Eigen::MatrixXi flowMatOO;
    Eigen::MatrixXi flowMatOP;
    Eigen::MatrixXi flowMatPP;

    Eigen::MatrixXd Noo;
    Eigen::MatrixXd Nop;
    Eigen::MatrixXd Npo;
    Eigen::MatrixXd Npp;

    Eigen::MatrixXd distGeoOO;
    Eigen::MatrixXd distGeoOP;
    Eigen::MatrixXd distGeoPP;


public:
    Network() = default;
    Network(const unsigned int id, const std::vector<Point>& obs, const std::vector<Point>& pred, const std::vector<StreamSegment>& seg);

    void setID(const unsigned int id) {ID = id;};

    const unsigned int getID() const {return ID;};
    const unsigned int getNObs() const {return nObs;};
    const unsigned int getNPred() const {return nPred;};
    const Eigen::MatrixXi& getFlowMatOO() const {return flowMatOO;};
    const Eigen::MatrixXi& getFlowMatOP() const {return flowMatOP;};
    const Eigen::MatrixXi& getFlowMatPP() const {return flowMatPP;};
    const std::vector<Point>& getObsPoints() const {return obspoints;};
    const std::vector<Point>& getPredPoints() const {return predpoints;};
    const Eigen::MatrixXd& getDistGeoOO() const {return distGeoOO;};
    const Eigen::MatrixXd& getDistGeoOP() const {return distGeoOP;};
    const Eigen::MatrixXd& getDistGeoPP() const {return distGeoPP;};
    const Eigen::MatrixXd& getNoo() const {return Noo;};
    const Eigen::MatrixXd& getNop() const {return Nop;};
    const Eigen::MatrixXd& getNpo() const {return Npo;};
    const Eigen::MatrixXd& getNpp() const {return Npp;};

    void computeDistances();

    void print() const;
};

#endif //PACSPROJECT_NETWORK_HPP

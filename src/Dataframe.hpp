#ifndef PACSPROJECT_DATAFRAME_HPP
#define PACSPROJECT_DATAFRAME_HPP

#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

class Dataframe{

private:
    std::map<std::string, Eigen::VectorXd> dataframe;
    unsigned int n;
    unsigned int p;

public:
    Dataframe() = default;
    Dataframe(const std::vector<std::string>& colnames, const Eigen::MatrixXd& data);

    unsigned int rows() const;
    unsigned int cols() const;

    const Eigen::VectorXd& operator[](const std::string& key) const;
    //const Eigen::VectorXd& operator[](const unsigned int id) const;

    Eigen::MatrixXd computeWeightMat(const std::string& weightVar);
    Eigen::MatrixXd computeWeightMat(const std::string& weightVar, const Eigen::VectorXd& otherPoints);

    void print() const;
};

#endif //PACSPROJECT_DATAFRAME_HPP

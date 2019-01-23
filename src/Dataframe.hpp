#ifndef PACSPROJECT_DATAFRAME_HPP
#define PACSPROJECT_DATAFRAME_HPP

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

/*! \file
* Dataframe class, representing a dataframe as a matrix with labeled columns.
* Labels are the names of the variables.
*/

class Dataframe{

private:
    std::map<std::string, Eigen::VectorXd> dataframe; /**< map that associates to a colum (Eigen::VectorXd) of a data frame the variable name */
    unsigned int n; ///< number of rows, the statistical units
    unsigned int p; ///< numer of variables (size of the map)

public:
    /**
    * Default constructor.
    */
    Dataframe() = default;

    /**
    * Constructor.
    * @param colnames vector containing the names of the variables
    * @param data matrix, whose columns correspond to the values of the statistical units per each variable
    */
    Dataframe(const std::vector<std::string>& colnames, const Eigen::MatrixXd& data);

    /**
    * @return The data frame number of rows
    */
    unsigned int rows() const;

    /**
    * @return The data frame number of columns
    */
    unsigned int cols() const;

    /**
    * Operator to extract a column of the data frame
    * @param key string indicating the name of the variable whose column is to be extracted
    * @return The column of the data frame corresponding to the variable selected
    */
    const Eigen::VectorXd& operator[](const std::string& key) const;

    /**
    * Compute the weight matrix for the tail-up model corresponding to a single group of points
    * @param weightVar constant string indicating the name of the variable to be used to compute the weights
    * @return A matrix whose elements are the weights for the tail-up model
    */
    Eigen::MatrixXd computeWeightMat(const std::string& weightVar);

    /**
    * Compute the weight matrix for the tail-up model corresponding to two groups of points
    * @param weightVar constant string indicating the name of the variable to be used to compute the weights
    * @param otherPoints vector representing the column of the weight variable of another data frame, for instance the one of
    * the prediction points
    * @return A matrix whose elements are the weights for the tail-up model, where rows correspond to the
    * points of the class (in general, the observed points) and columns correspond to the other group of
    * points (in general, the prediction points)
    */
    Eigen::MatrixXd computeWeightMat(const std::string& weightVar, const Eigen::VectorXd& otherPoints);

    /**
    * Printing function to visualize the data frame as a matrix with labeled columns
    */
    void print() const;
};

#endif

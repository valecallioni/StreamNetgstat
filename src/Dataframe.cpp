#include "Dataframe.hpp"

Dataframe::Dataframe(const std::vector<std::string>& colnames, const Eigen::MatrixXd& data){
    if (colnames.size() != data.cols())
        throw std::length_error("The matrix data and the variable names have not the same dimensions");
    n = data.rows();
    p = colnames.size();
    for (int i=0; i<colnames.size(); i++){
        dataframe[colnames[i]] = data.col(i);
    }
}

unsigned int Dataframe::cols() const {
    return p;
}

unsigned int Dataframe::rows() const {
    return n;
}

const Eigen::VectorXd& Dataframe::operator[](const std::string& key) const{
    return dataframe.at(key);
}

Eigen::MatrixXd Dataframe::computeWeightMat(const std::string& weightVar){
    Eigen::VectorXd& addFuncVal(dataframe[weightVar]);
    Eigen::MatrixXd res(n,n);
    res.fill(1.0);
    unsigned int min = 0;
    unsigned int max = 0;
    for (unsigned int i=0; i<n; i++){
        for (unsigned int j=i+1; j<n; j++){
            if (addFuncVal(i)<addFuncVal(j)) {
                min = i;
                max = j;
            }
            else {
                min = j;
                max = i;
            }
            res(i,j) = sqrt(addFuncVal(min)/addFuncVal(max));
            res(j,i) = sqrt(addFuncVal(min)/addFuncVal(max));
        }
    }
    return res;
}

Eigen::MatrixXd Dataframe::computeWeightMat(const std::string& weightVar, const Eigen::VectorXd& otherPoints) {
    Eigen::VectorXd& addFuncVal(dataframe[weightVar]);
    unsigned int m(otherPoints.size());
    Eigen::MatrixXd res(n,m);
    for (unsigned int i=0; i<n; i++){
        for (unsigned int j=0; j<m; j++){
            if (addFuncVal(i)<otherPoints(j)) {
                res(i,j) = sqrt(addFuncVal(i)/otherPoints(j));
            }
            else {
                res(i,j) = sqrt(otherPoints(j)/addFuncVal(i));
            }

        }
    }
    return res;
}

void Dataframe::print() const {
  std::cout << "Dataframe \n";
  for (auto it=dataframe.cbegin(); it!=dataframe.cend(); it++){
    std::cout << it->first;
  }
  for (auto it=dataframe.cbegin(); it!=dataframe.cend(); it++){
    std::cout << it->second;
  }
  std::cout << "\n";
}

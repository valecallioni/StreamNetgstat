#ifndef PACSPROJECT_POINT_HPP
#define PACSPROJECT_POINT_HPP

#include <string>
#include <vector>
#include <iostream>

class Point{

protected:
    double distUpstream;
    unsigned int rid;   // index of the segment
    unsigned int pid;   // index of the point
    std::vector<char> binaryID;
    std::string ID = "obs";
    double x1;
    double x2;

public:
    Point() = default;
    Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
            const double coord1, const double coord2);
    Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
          const double coord1, const double coord2, const std::string& id);

    const double getDistUpstream() const {return distUpstream;};
    const unsigned int getRid() const {return rid;};
    const unsigned int getPid() const {return pid;};
    const std::vector<char>& getBinaryID() const {return binaryID;};
    const std::string getID() const {return ID;};
    const double getX1() const {return x1;};
    const double getX2() const {return x2;};

    void setRid(const unsigned int r) {rid = r;};
    void setPid(const unsigned int p) {pid = p;};
    void setDistUpstream(const double distUp) {distUpstream = distUp;};
    void setBinaryID(const std::string& binID);
    void setID(const std::string id) {ID = id;};
    void setCoordinates(const double coord1, const double coord2);

    void print() const;

};

#endif

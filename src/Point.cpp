#include "Point.hpp"

Point::Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
             const double coord1, const double coord2): distUpstream(dist), rid(r), pid(p),
             x1(coord1), x2(coord2){
    for (auto c: binID)
        binaryID.push_back(c);
}

Point::Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
             const double coord1, const double coord2, const std::string& id): distUpstream(dist), rid(r),
             pid(p), ID(id), x1(coord1), x2(coord2) {
    for (auto c: binID)
        binaryID.push_back(c);
}

void Point::setCoordinates(const double coord1, const double coord2){
  x1 = coord1;
  x2 = coord2;
}


void Point::setBinaryID(const std::string& binID){
    binaryID.clear();
    for (auto c: binID)
        binaryID.push_back(c);
}

void Point::print() const {
    //std::cout << "Point " << ID << ", pid: " << pid << ", rid: ";
    //std::cout << rid << ", distUpstream " << distUpstream << "and binaryID ";
    //for (auto c: binaryID)
    //    std::cout << c;
    std::cout << "Point " << pid << ": (" << x1 << ", " << x2 << ")" << std::endl;
}

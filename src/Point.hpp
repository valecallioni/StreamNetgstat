#ifndef PACSPROJECT_POINT_HPP
#define PACSPROJECT_POINT_HPP

#include <string>
#include <vector>
#include <iostream>

/**
* Point class, representing a single point, observed or for prediction.
* Each point has information about the binaryID and the rid of the segment it lies on,
* its upstream distance and its coordinates.
*/

class Point{

protected:
    double distUpstream;  ///< distance upstream
    unsigned int rid; ///< ID of the segment the point lies on
    unsigned int pid;  ///< ID of the point
    std::vector<char> binaryID;  ///< binaryID of the segment the point lies on, divided into characters
    std::string ID = "obs";
    double x1;  ///< first coordinate
    double x2;  ///< second coordinate

public:
    /**
    * Default constructor.
    */
    Point() = default;

    /**
    * Constructor for an observed point.
    * @param r ID of the segment
    * @param binID binaryID of the segment
    * @param dist distance upstream
    * @param p ID of the point
    * @param coord1 first coordinate
    * @param coord2 second coordinate
    */
    Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
            const double coord1, const double coord2);

    /**
    * Constructor for a prediction point.
    * @param r ID of the segment
    * @param binID binaryID of the segment
    * @param dist distance upstream
    * @param p ID of the point
    * @param coord1 first coordinate
    * @param coord2 second coordinate
    * @param id string indicating the name of the prediction points group the point belongs to
    */
    Point(const unsigned int r, const std::string& binID, const double dist, const unsigned int p,
          const double coord1, const double coord2, const std::string& id);

    /**
    * @return the distance upstream
    */
    const double getDistUpstream() const {return distUpstream;};

    /**
    * @return the ID of the segment
    */
    const unsigned int getRid() const {return rid;};

    /**
    * @return the ID of the point
    */
    const unsigned int getPid() const {return pid;};

    /**
    * @return the vector of characters representing the binaryID of the segment
    */
    const std::vector<char>& getBinaryID() const {return binaryID;};

    /**
    * @return the string ID of the point ("obs" for observed points)
    */
    const std::string getID() const {return ID;};

    /**
    * @return the first coordinate
    */
    const double getX1() const {return x1;};

    /**
    * @return the second coordinate
    */
    const double getX2() const {return x2;};

    /**
    * @param r ID of the segment
    */
    void setRid(const unsigned int r) {rid = r;};

    /**
    * @param p ID of the point
    */
    void setPid(const unsigned int p) {pid = p;};

    /**
    * @param distUp distance upstream
    */
    void setDistUpstream(const double distUp) {distUpstream = distUp;};

    /**
    * @param binID binaryID of the segment, to be decomposed into a vector of characters;
    */
    void setBinaryID(const std::string& binID);

    /**
    * @param id string ID of the point
    */
    void setID(const std::string id) {ID = id;};

    /**
    * @param coord1 first coordinate
    * @param coord2 second coordinate
    */
    void setCoordinates(const double coord1, const double coord2);

    /**
    * Printing function to visualize the attributes of the stream segment
    */
    void print() const;

};

#endif

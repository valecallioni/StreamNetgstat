#ifndef PACSPROJECT_STREAMSEGMENT_HPP
#define PACSPROJECT_STREAMSEGMENT_HPP

#include <string>
#include <iostream>

/*! \file
* StreamSegment class, representing a segment of a stream network as a continuous line.
* Each segment has an ID (rid) and a binaryID representing the position in the network.
*/

class StreamSegment{

private:
    unsigned int networkID; ///< ID of the network the segment belongs to
    unsigned int segmentID; ///< rid
    double distUpstream; ///< distance upstream of the most upstream point lying on the segment
    std::string binaryID; ///< binaryID of the segment

public:
    /**
    * Default constructor
    */
    StreamSegment() = default;

    /**
    * Constructor.
    * @param netID ID of the network
    * @param segID rid, ID of the segment
    * @param distUp distance upstream of the most upstream point lying on the segment
    * @param binID binaryID of the segment
    */
    StreamSegment(const unsigned int netID, const unsigned int segID,
      const double distUp, const std::string binID): networkID(netID),
      segmentID(segID), distUpstream(distUp), binaryID(binID){};

    /**
    * @return the ID of the network
    */
    const unsigned int getNetworkID() const {return networkID;};

    /**
    * @return the ID of the segment
    */
    const unsigned int getSegmentID() const {return segmentID;};

    /**
    * @return the upstream distance of the most upstream point lying on the segment
    */
    const double getDistUpstream() const {return distUpstream;};

    /**
    * @return the binaryID of the segment
    */
    const std::string getBinaryID() const {return binaryID;};

    /**
    * @param netID ID of the network
    */
    void setNetworkID(const unsigned int netID) {networkID = netID;};

    /**
    * @param segID ID of the segment
    */
    void setSegmentID(const unsigned int segID) {segmentID = segID;};

    /**
    * @param distUp distance upstream of the most upstream point lying on the segment
    */
    void setDistUpstream(const double distUp) {distUpstream = distUp;};

    /**
    * @param binID binaryID of the segment
    */
    void setBinaryID(const std::string binID) {binaryID = binID;};

    /**
    * Printing function to visualize the attributes of the stream segment
    */
    void print() const;
};

#endif

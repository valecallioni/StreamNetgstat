#ifndef PACSPROJECT_STREAMSEGMENT_HPP
#define PACSPROJECT_STREAMSEGMENT_HPP

#include <string>
#include <iostream>

class StreamSegment{

private:
    unsigned int networkID;
    unsigned int segmentID;
    double distUpstream;
    std::string binaryID;

public:
    StreamSegment(const unsigned int netID, const unsigned int segID,
      const double distUp, const std::string binID): networkID(netID),
      segmentID(segID), distUpstream(distUp), binaryID(binID){};
    StreamSegment() = default;

    const unsigned int getNetworkID() const {return networkID;};
    const unsigned int getSegmentID() const {return segmentID;};
    const double getDistUpstream() const {return distUpstream;};
    const std::string getBinaryID() const {return binaryID;};

    void setNetworkID(const unsigned int netID) {networkID = netID;};
    void setSegmentID(const unsigned int segID) {segmentID = segID;};
    void setDistUpstream(const double distUp) {distUpstream = distUp;};
    void setBinaryID(const std::string binID) {binaryID = binID;};

    void print() const;
};

#endif

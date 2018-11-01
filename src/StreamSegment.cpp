#include "StreamSegment.hpp"

void StreamSegment::print() const {
    std::cout << "Segment " << segmentID << " on Network n. " << networkID << " with binID " << binaryID << " and distUpstream " << distUpstream << std::endl;
}

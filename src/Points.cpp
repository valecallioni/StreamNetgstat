#include "Points.hpp"

Points::Points(const std::vector<Point>& p): points(p){
  n = p.size();
  flowMat.resize(n,n);
  distHydro.resize(n,n);
  distGeo.resize(n,n);
}

void Points::setPoints(const std::vector<Point>& p){
  points = p;
  n = p.size();
  flowMat.resize(n,n);
  flowMat.fill(0);
  distHydro.resize(n,n);
  distHydro.fill(0.0);
  distGeo.resize(n,n);
  distGeo.fill(0.0);
}

void Points::computeDistances(const std::map<std::string, StreamSegment>& segments) {

  for (unsigned int i = 0; i < 1; i++) {
      std::vector<char> p1 = points[i].getBinaryID();
      unsigned int n1 = p1.size();
      unsigned int min = n1;

      for (unsigned int j = i + 1; j < n; j++) {
          min = n1;
          std::vector<char> p2 = points[j].getBinaryID();
          unsigned int n2 = p2.size();
          if (n2 < n1)
              min = n2;
          unsigned int k = 0;
          while (k < min && p1[k] == p2[k])
              k++;

          if (k == n1 || k == n2) { // flow-connection

              flowMat(i,j) = 1;
              flowMat(j,i) = 1;
              if (n2 == n1){
                  if (points[i].getDistUpstream() > points[j].getDistUpstream())
                      distHydro(i,j) = points[i].getDistUpstream() - points[j].getDistUpstream();
                  else
                      distHydro(j,i) = points[j].getDistUpstream() - points[i].getDistUpstream();
              }
              else if (n2 < n1)
                  distHydro(i,j) = points[i].getDistUpstream() - points[j].getDistUpstream();
              else
                  distHydro(j,i) = points[j].getDistUpstream() - points[i].getDistUpstream();

          } else { //flow-unconnection

              std::string junc("");
              for (unsigned int c = 0; c < k; c++) {
                  junc.append(1, points[i].getBinaryID()[c]);
              }
              distHydro(i, j) = points[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
              distHydro(j, i) = points[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
              if (distHydro(i,j)<0.0 || distHydro(j,i)<0.0){
                points[i].print();
                points[j].print();
                // std::cout << "points[i].getDistUpstream() = " << points[i].getDistUpstream() << ", where i = " << i << ", points[j].getDistUpstream() = " << points[j].getDistUpstream();
                // std::cout << ", where j = " << j << ", and segments.find(junc)->second.getDistUpstream() = " << segments.find(junc)->second.getDistUpstream();
                // std::cout << ", whose segID is " << segments.find(junc)->second.getSegmentID();
                // std::cout << " and has binID = " << segments.find(junc)->second.getBinaryID() << std::endl;
                // std::cout << "k = " << k << ", n1 = " << n1 << ", n2 = " << n2 << std::endl;
                // std::cout << "binID p1 = ";
                // for (auto c: points[i].getBinaryID()) std::cout << c;
                // std::cout << ", binID p2 = ";
                // for (auto c: points[j].getBinaryID()) std::cout << c;
                // std::cout << "\n";
                //throw std::domain_error("distance negative");
              }
          }

          // Compute Euclidean distances
          //distGeo(i,j) = helpers::geoDist(points[i].getX1(), points[i].getX2(), points[j].getX1(), points[j].getX2());
          distGeo(i,j) = std::sqrt((points[i].getX1() - points[j].getX1())*(points[i].getX1() - points[j].getX1()) + (points[i].getX2() - points[j].getX2())*(points[i].getX2() - points[j].getX2()));
          distGeo(j,i) = distGeo(i,j);

      }
      std::cout << "distHydro[,1]: \n" << distHydro.col(0) << std::endl;
  }

}

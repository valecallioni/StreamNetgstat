#include "Points.hpp"

Points::Points(const std::vector<Point>& p): points(p){
  n = p.size();
}

void Points::setPoints(const std::vector<Point>& p){
  points = p;
  n = p.size();
}

void Points::computeDistances(bool geo, const std::map<std::string, StreamSegment>& segments) {

  flowMat.resize(n,n);
  flowMat.fill(0);
  distHydro.resize(n,n);
  distHydro.fill(0.0);

  if (geo){
    distGeo.resize(n,n);
    distGeo.fill(0.0);
  }

  for (unsigned int i = 0; i < n; i++) {
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
        }

        if (geo) {
          // Compute Euclidean distances
          //distGeo(i,j) = helpers::geoDist(points[i].getX1(), points[i].getX2(), points[j].getX1(), points[j].getX2());
          distGeo(i,j) = std::sqrt((points[i].getX1() - points[j].getX1())*(points[i].getX1() - points[j].getX1()) + (points[i].getX2() - points[j].getX2())*(points[i].getX2() - points[j].getX2()));
          distGeo(j,i) = distGeo(i,j);
        }
      }
  }

}

void Points::setDistances(const std::vector<Eigen::MatrixXd>& matrices){
  if (matrices.size() == 3) {
    flowMat = matrices[0].cast<int>();
    distHydro = matrices[1];
    distGeo = matrices[2];
  }
  else {
    flowMat = matrices[0].cast<int>();
    distHydro = matrices[1];
  }
}

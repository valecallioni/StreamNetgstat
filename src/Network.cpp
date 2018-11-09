#include "Network.hpp"

Network::Network(const unsigned int id, const std::vector<Point>& obs,
  const std::vector<Point>& pred, const std::vector<StreamSegment>& seg): ID(id){

    obspoints.setPoints(obs);
    predpoints.setPoints(pred);

    nObs = obspoints.getN();
    nPred = predpoints.getN();

    flowMatOP.resize(nObs, nPred);
    flowMatOP.fill(0);

    distHydroOP.resize(nObs, nPred);
    distHydroOP.fill(0.0);

    distHydroPO.resize(nPred, nObs);
    distHydroPO.fill(0.0);

    distGeoOP.resize(nObs, nPred);
    distGeoOP.fill(0.0);

    for (unsigned int i=0; i<seg.size(); i++){
        segments[seg[i].getBinaryID()] = seg[i];
    }
}

Network::Network(const unsigned int id, const std::vector<Point>& obs,
  const std::vector<StreamSegment>& seg): ID(id){

    obspoints.setPoints(obs);

    nObs = obspoints.getN();
    nPred = 0;

    flowMatOP.resize(0, 0);
    distHydroOP.resize(0, 0);
    distHydroPO.resize(0, 0);
    distGeoOP.resize(0, 0);

    for (unsigned int i=0; i<seg.size(); i++){
        segments[seg[i].getBinaryID()] = seg[i];
    }
}

void Network::print() const {
    std::cout << "Printing observed points: \n";
    for (auto& p: obspoints.getPoints())
        p.print();
    /*std::cout << "Printing predction points: \n";
    for (unsigned int i=0; i<10; i++)
        predpoints.getPoints()[i].print();
    std::cout << "Printing segments: \n";
    for (auto& s: segments)
        (s.second).print();*/
}

void Network::computeDistances() {

    std::vector<Point> obs(obspoints.getPoints());
    std::vector<Point> pred(predpoints.getPoints());

    obspoints.computeDistances(segments);
    predpoints.computeDistances(segments);

    if (nPred > 0){
      for (unsigned int i=0; i<nObs; i++) {

          std::vector<char> p1 = obs[i].getBinaryID();
          unsigned int n1 = p1.size();
          unsigned int min = n1;

          for (unsigned int j = 0; j < nPred; j++) {
              min = n1;
              std::vector<char> p2 = pred[j].getBinaryID();
              unsigned int n2 = p2.size();
              if (n2 < n1)
                  min = n2;
              unsigned int k = 0;
              while (k < min && p1[k] == p2[k])
                  k++;

              if (k == n1 || k == n2){ // flow-connection
                  flowMatOP(i,j) = 1;
                  if (n2 == n1){
                      if (obs[i].getDistUpstream() > pred[j].getDistUpstream())
                          distHydroOP(i,j) = obs[i].getDistUpstream() - pred[j].getDistUpstream();
                      else
                          distHydroPO(j,i) = pred[j].getDistUpstream() - obs[i].getDistUpstream();
                  }
                  else if (n2 < n1){
                      distHydroOP(i,j) = obs[i].getDistUpstream() - pred[j].getDistUpstream();
                  }
                  else
                      distHydroPO(j,i) = pred[j].getDistUpstream() - obs[i].getDistUpstream();
              }

              else { // flow-unconnection
                  std::string junc("");
                  for (unsigned int c=0; c<k; c++){
                      junc.append(1, obs[i].getBinaryID()[c]);
                  }
                  distHydroOP(i,j) = obs[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                  distHydroPO(j,i) = pred[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
              }

              // Compute Euclidean distances
              distGeoOP(i,j) = std::sqrt((obs[i].getX1() - pred[j].getX1())*(obs[i].getX1() - pred[j].getX1()) + (obs[i].getX2() - pred[j].getX2())*(obs[i].getX2() - pred[j].getX2()));

          }

      }
    }

}

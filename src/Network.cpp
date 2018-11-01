#include "Network.hpp"

Network::Network(const unsigned int id, const std::vector<Point>& obs,
  const std::vector<Point>& pred, const std::vector<StreamSegment>& seg):
ID(id), obspoints(obs), predpoints(pred){

    nObs = obspoints.size();
    nPred = predpoints.size();

    flowMatOO.resize(nObs, nObs);
    flowMatOO.fill(0);
    flowMatOP.resize(nObs, nPred);
    flowMatOP.fill(0);
    flowMatPP.resize(nPred, nPred);
    flowMatPP.fill(0);

    Noo.resize(nObs, nObs);
    Noo.fill(0.0);
    Nop.resize(nObs, nPred);
    Nop.fill(0.0);
    Npo.resize(nPred, nObs);
    Npo.fill(0.0);
    Npp.resize(nPred, nPred);
    Npp.fill(0.0);

    distGeoOO.resize(nObs, nObs);
    distGeoOO.fill(0.0);
    distGeoOP.resize(nObs, nPred);
    distGeoOP.fill(0.0);
    distGeoPP.resize(nPred, nPred);
    distGeoPP.fill(0.0);

    for (unsigned int i=0; i<seg.size(); i++){
        segments[seg[i].getBinaryID()] = seg[i];
    }
}

void Network::print() const {
    std::cout << "Printing observed points: \n";
    for (auto& p: obspoints)
        p.print();
    /*std::cout << "Printing predction points: \n";
    for (unsigned int i=0; i<10; i++)
        predpoints[i].print();
    std::cout << "Printing segments: \n";
    for (auto& s: segments)
        (s.second).print();*/
}

void Network::computeDistances() {

    for (unsigned int i=0; i<nObs; i++) {

        std::vector<char> p1 = obspoints[i].getBinaryID();
        unsigned int n1 = p1.size();
        unsigned int min = n1;
        for (unsigned int j = i + 1; j < nObs; j++) {
            min = n1;
            std::vector<char> p2 = obspoints[j].getBinaryID();
            unsigned int n2 = p2.size();
            if (n2 < n1)
                min = n2;
            unsigned int k = 0;
            while (k < min && p1[k] == p2[k])
                k++;

            if (k == n1 || k == n2) { // flow-connection
                /*relationshipsTableOO(i, j) = (-1.0) * double(k);
                relationshipsTableOO(j, i) = (-1.0) * double(k);*/
                flowMatOO(i,j) = 1;
                flowMatOO(j,i) = 1;
                if (n1 == n2){
                    if (obspoints[i].getDistUpstream() > obspoints[j].getDistUpstream())
                        Noo(i,j) = obspoints[i].getDistUpstream() - obspoints[j].getDistUpstream();
                    else
                        Noo(j,i) = obspoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
                }
                else if (n2 < n1)
                    Noo(i,j) = obspoints[i].getDistUpstream() - obspoints[j].getDistUpstream();
                else
                    Noo(j,i) = obspoints[j].getDistUpstream() - obspoints[i].getDistUpstream();


            } else { // flow-unconnection
                /*relationshipsTableOO(i, j) = double(k);
                relationshipsTableOO(j, i) = double(k);*/
                std::string junc("");
                for (unsigned int c=0; c<k; c++){
                    junc.append(1,obspoints[i].getBinaryID()[c]);
                }
                Noo(i,j) = obspoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Noo(j,i) = obspoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

            // Compute Euclidean distances
            distGeoOO(i,j) = helpers::geoDist(obspoints[i].getX1(), obspoints[i].getX2(), obspoints[j].getX1(), obspoints[j].getX2());
            distGeoOO(j,i) = distGeoOO(i,j);
        }

        for (unsigned int j = 0; j < nPred; j++) {
            min = n1;
            std::vector<char> p2 = predpoints[j].getBinaryID();
            unsigned int n2 = p2.size();
            if (n2 < n1)
                min = n2;
            unsigned int k = 0;
            while (k < min && p1[k] == p2[k])
                k++;

            if (k == n1 || k == n2){ // flow-connection
                //relationshipsTableOP(i, j) = (-1.0) * double(k);
                flowMatOP(i,j) = 1;
                if (n2 == n1){
                    if (obspoints[i].getDistUpstream() > predpoints[j].getDistUpstream())
                        Nop(i,j) = obspoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                    else
                        Npo(j,i) = predpoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
                }
                else if (n2 < n1){
                    Nop(i,j) = obspoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                }
                else
                    Npo(j,i) = predpoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
            }

            else { // flow-unconnection
                //relationshipsTableOP(i, j) = double(k);
                std::string junc("");
                for (unsigned int c=0; c<k; c++){
                    junc.append(1,obspoints[i].getBinaryID()[c]);
                }
                Nop(i,j) = obspoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Npo(j,i) = predpoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

            // Compute Euclidean distances
            distGeoOP(i,j) = helpers::geoDist(obspoints[i].getX1(), obspoints[i].getX2(), predpoints[j].getX1(), predpoints[j].getX2());

        }

    }

    for (unsigned int i = 0; i < nPred; i++) {
        std::vector<char> p1 = predpoints[i].getBinaryID();
        unsigned int n1 = p1.size();
        unsigned int min = n1;

        for (unsigned int j = i + 1; j < nPred; j++) {
            min = n1;
            std::vector<char> p2 = predpoints[j].getBinaryID();
            unsigned int n2 = p2.size();
            if (n2 < n1)
                min = n2;
            unsigned int k = 0;
            while (k < min && p1[k] == p2[k])
                k++;

            if (k == n1 || k == n2) { // flow-connection
                /*relationshipsTablePP(i, j) = (-1.0) * double(k);
                relationshipsTablePP(j, i) = (-1.0) * double(k);*/
                flowMatPP(i,j) = 1;
                flowMatPP(j,i) = 1;
                if (n2 == n1){
                    if (predpoints[i].getDistUpstream() > predpoints[j].getDistUpstream())
                        Npp(i,j) = predpoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                    else
                        Npp(j,i) = predpoints[j].getDistUpstream() - predpoints[i].getDistUpstream();
                }
                else if (n2 < n1)
                    Npp(i,j) = predpoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                else
                    Npp(j,i) = predpoints[j].getDistUpstream() - predpoints[i].getDistUpstream();

            } else { //flow-unconnection
                /*relationshipsTablePP(i, j) = double(k);
                relationshipsTablePP(j, i) = double(k);*/
                std::string junc("");
                for (unsigned int c = 0; c < k; c++) {
                    junc.append(1, predpoints[i].getBinaryID()[c]);
                }
                Npp(i, j) = predpoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Npp(j, i) = predpoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

            // Compute Euclidean distances
            distGeoPP(i,j) = helpers::geoDist(predpoints[i].getX1(), predpoints[i].getX2(), predpoints[j].getX1(), predpoints[j].getX2());
            distGeoPP(j,i) = distGeoPP(i,j);

        }
    }

}

/*
void Network::computeDistances() {

    for (unsigned int i=0; i<nObs; i++){

        for (unsigned int j=i+1; j<nObs; j++){

            if (relationshipsTableOO(i,j)<=0){   // flow-connection
                if (obspoints[j].getBinaryID().size() == obspoints[i].getBinaryID().size()){
                    if (obspoints[i].getDistUpstream() > obspoints[j].getDistUpstream())
                        Noo(i,j) = obspoints[i].getDistUpstream() - obspoints[j].getDistUpstream();
                    else
                        Noo(j,i) = obspoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
                }
                else if (obspoints[j].getBinaryID().size() < obspoints[i].getBinaryID().size())
                    Noo(i,j) = obspoints[i].getDistUpstream() - obspoints[j].getDistUpstream();
                else
                    Noo(j,i) = obspoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
            }

            else {  // flow-unconnection
                std::string junc("");
                for (unsigned int c=0; c<relationshipsTableOO(i,j); c++){
                    junc.append(1,obspoints[i].getBinaryID()[c]);
                }
                Noo(i,j) = obspoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Noo(j,i) = obspoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

        }

        for (unsigned int j=0; j<nPred; j++){

            if (relationshipsTableOP(i,j)<=0){   // flow-connection
                if (predpoints[j].getBinaryID().size() == obspoints[i].getBinaryID().size()){
                    if (obspoints[i].getDistUpstream() > predpoints[j].getDistUpstream())
                        Nop(i,j) = obspoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                    else
                        Npo(j,i) = predpoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
                }
                else if (predpoints[j].getBinaryID().size() < obspoints[i].getBinaryID().size()){
                    Nop(i,j) = obspoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                }
                else
                    Npo(j,i) = predpoints[j].getDistUpstream() - obspoints[i].getDistUpstream();
            }

            else {  // flow-unconnection
                std::string junc("");
                for (unsigned int c=0; c<relationshipsTableOP(i,j); c++){
                    junc.append(1,obspoints[i].getBinaryID()[c]);
                }
                Nop(i,j) = obspoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Npo(j,i) = obspoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

        }

    }

    for (unsigned int i=0; i<nPred; i++) {

        for (unsigned int j = i+1; j < nPred; j++) {

            if (relationshipsTablePP(i,j)<=0){   // flow-connection
                if (predpoints[j].getBinaryID().size() == predpoints[i].getBinaryID().size()){
                    if (predpoints[i].getDistUpstream() > predpoints[j].getDistUpstream())
                        Npp(i,j) = predpoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                    else
                        Npp(j,i) = predpoints[j].getDistUpstream() - predpoints[i].getDistUpstream();
                }
                else if (predpoints[j].getBinaryID().size() < predpoints[i].getBinaryID().size())
                    Npp(i,j) = predpoints[i].getDistUpstream() - predpoints[j].getDistUpstream();
                else
                    Npp(j,i) = predpoints[j].getDistUpstream() - predpoints[i].getDistUpstream();
            } else {  // flow-unconnection
                std::string junc("");
                for (unsigned int c = 0; c < relationshipsTablePP(i, j); c++) {
                    junc.append(1, predpoints[i].getBinaryID()[c]);
                }
                Npp(i, j) = predpoints[i].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
                Npp(j, i) = predpoints[j].getDistUpstream() - segments.find(junc)->second.getDistUpstream();
            }

        }
    }
}*/

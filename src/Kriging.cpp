#include "Kriging.hpp"

Kriging::Kriging(const Eigen::MatrixXd& dMatPred, const Eigen::MatrixXd& dMatObs, const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& Nop, const Eigen::MatrixXd& Npo, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat,
  const Eigen::VectorXd& param, const Eigen::VectorXd& y, std::unique_ptr<TailUpModel>& tailup_ptr,
  std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, int nMod, bool useNugg){

    Xpred = std::make_shared<Eigen::MatrixXd>(dMatPred);
    Xobs = std::make_shared<Eigen::MatrixXd>(dMatObs);
    distHydroOP = std::make_shared<Eigen::MatrixXd>(Nop);
    distHydroPO = std::make_shared<Eigen::MatrixXd>(Npo);
    distGeo = std::make_shared<Eigen::MatrixXd>(D);
    weightMat = std::make_shared<Eigen::MatrixXd>(wMat);
    flowMat = std::make_shared<Eigen::MatrixXi>(connMat);

    z = std::make_shared<Eigen::VectorXd>(y);
    theta = std::make_shared<Eigen::VectorXd>(param);
    nModels = nMod;
    useNugget = useNugg;
    tailUpModel = std::move(tailup_ptr);
    tailDownModel = std::move(taildown_ptr);
    euclidModel = std::move(euclid_ptr);

    nObs = dMatObs.rows();
    nPred = dMatPred.rows();
    p = dMatPred.cols() - 1;

    Eigen::LDLT<Eigen::MatrixXd> solver(nObs);
    solver.compute(V);
    Eigen::MatrixXd Id(nObs,nObs);
    Id.setIdentity();
    invV = std::make_shared<Eigen::MatrixXd>(solver.solve(Id));

    Vpred.resize(nObs, nPred);
    Vpred.fill(0.0);

    predData.resize(nPred, 2);
    predData.fill(0.0);

    // Compute Vpred
    int j = 0;
    if (j < nModels*2 && tailUpModel){
      tailUpModel->setSigma2((*theta)(j));
      parsill += (*theta)(j);
      j++;
      tailUpModel->setAlpha((*theta)(j));
      j++;
    }
    if (j < nModels*2 && tailDownModel){
      tailDownModel->setSigma2((*theta)(j));
      parsill += (*theta)(j);
      j++;
      tailDownModel->setAlpha((*theta)(j));
      j++;
    }
    if (j < nModels*2 && euclidModel){
      euclidModel->setSigma2((*theta)(j));
      parsill += (*theta)(j);
      j++;
      euclidModel->setAlpha((*theta)(j));
      j++;
    }

    if (tailUpModel) Vpred += tailUpModel->computeMatCov(*weightMat, *distHydroOP, *distHydroPO);
    if (tailDownModel) Vpred += tailDownModel->computeMatCov(*flowMat, *distHydroOP, *distHydroPO);
    if (euclidModel) Vpred += euclidModel->computeMatCov(*distGeo);

    //Compute XV and invXVX
    XV = Xobs->transpose() * (*invV);
    solver = Eigen::LDLT<Eigen::MatrixXd>(p+1);
    solver.compute(Xobs->transpose()*(*invV)*(*Xobs));
    Id.resize(p+1,p+1);
    Id.setIdentity();
    invXVX = solver.solve(Id);
}

void Kriging::predict(){

  Eigen::VectorXd r(p+1);
  Eigen::VectorXd m(p+1);
  Eigen::VectorXd lambda(nObs);

  for (unsigned int i=0; i<nPred; i++){
    r = (Xpred->row(i)).transpose() - XV*Vpred.col(i);
    m = invXVX*r;

    lambda = (*invV)*(Vpred.col(i) + (*Xobs)*invXVX*r);

    predData(i,0) = lambda.transpose() * (*z);
    predData(i,1) = std::sqrt(parsill - lambda.transpose()*Vpred.col(i) + m.transpose()*(Xpred->row(i)).transpose());
  }

}

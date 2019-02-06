#include "Kriging.hpp"

Kriging::Kriging(std::shared_ptr<Eigen::MatrixXd> dMatPred, std::shared_ptr<Eigen::MatrixXd> dMatObs, std::shared_ptr<Eigen::MatrixXd> V,
  std::shared_ptr<Eigen::MatrixXd> Nop, std::shared_ptr<Eigen::MatrixXd> Npo, std::shared_ptr<Eigen::MatrixXd> D, std::shared_ptr<Eigen::MatrixXd> wMat, std::shared_ptr<Eigen::MatrixXi> connMat,
  std::shared_ptr<Eigen::VectorXd> param, std::shared_ptr<Eigen::VectorXd> y, std::unique_ptr<TailUpModel>& tailup_ptr,
  std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, int nMod, bool useNugg){

    Xpred = dMatPred;
    Xobs = dMatObs;
    distHydroOP = Nop;
    distHydroPO = Npo;
    if (D->rows() > 0) distGeo = D;
    weightMat = wMat;
    flowMat = connMat;

    z = y;
    theta = param;
    nModels = nMod;
    useNugget = useNugg;
    tailUpModel = std::move(tailup_ptr);
    tailDownModel = std::move(taildown_ptr);
    euclidModel = std::move(euclid_ptr);

    nObs = dMatObs->rows();
    nPred = dMatPred->rows();
    p = dMatPred->cols() - 1;

    Eigen::MatrixXd Id(nObs,nObs);
    Id.setIdentity();
    if (V->determinant()>1e-3){
      Eigen::LDLT<Eigen::MatrixXd> solver(nObs);
      solver.compute(*V);
      invV = std::make_shared<Eigen::MatrixXd>(solver.solve(Id));
      solver.setZero();
    }
    else {
      Eigen::HouseholderQR<Eigen::MatrixXd> solver(nObs, nObs);
      solver.compute(*V);
      invV = std::make_shared<Eigen::MatrixXd>(solver.solve(Id));
    }

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

    if (useNugget) parsill += (*theta)(j)

    if (tailUpModel) Vpred += tailUpModel->computeMatCov(*weightMat, *distHydroOP, *distHydroPO);
    if (tailDownModel) Vpred += tailDownModel->computeMatCov(*flowMat, *distHydroOP, *distHydroPO);
    if (euclidModel) Vpred += euclidModel->computeMatCov(*distGeo);
    if (useNugget) Vpred += Eigen::MatrixXd::Identity(nObs,nPred)*(*theta)(theta->size()-1);

    //Compute XV and invXVX
    XV = Xobs->transpose() * (*invV);
    Eigen::MatrixXd XVX(XV * (*Xobs));
    Id.resize(p+1,p+1);
    Id.setIdentity();
    if (XVX.determinant() >= 1e-3){
      Eigen::LDLT<Eigen::MatrixXd> solver(p+1);
      solver.compute(XVX);
      invXVX = solver.solve(Id);
      solver.setZero();
    }
    else {
      Eigen::HouseholderQR<Eigen::MatrixXd> solver2(p+1,p+1);
      solver2.compute(XVX);
      invXVX = solver2.solve(Id);
    }

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

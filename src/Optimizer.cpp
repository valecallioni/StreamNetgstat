#include "Optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg,
  int n_models, const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat, bool cholesky):
useNugget(useNugg), nModels(n_models), useLDLT(cholesky){

  z = std::make_shared<Eigen::VectorXd>(y);
  X = std::make_shared<Eigen::MatrixXd>(designMat);
  distHydro = std::make_shared<Eigen::MatrixXd>(N);
  if (D.rows() > 0) distGeo = std::make_shared<Eigen::MatrixXd>(D);
  weightMat = std::make_shared<Eigen::MatrixXd>(wMat);
  flowMat = std::make_shared<Eigen::MatrixXi>(connMat);

  tailUpModel = std::move(tailup_ptr);
  tailDownModel = std::move(taildown_ptr);
  euclidModel = std::move(euclid_ptr);

  n = X->rows();
  p = X->cols()-1;

  optimTheta.resize(nModels*2 + 1*useNugget);
  betaValues.resize(p+1);
  covMat.resize(n,n);
  covMat.fill(0.0);

  maxDistHydro = 4.0*distHydro->maxCoeff();
  if (euclidModel) maxDistGeo = 4.0*distGeo->maxCoeff();
}

Optimizer::Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg, const std::vector<double>& bounds,
  int n_models, const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat, bool cholesky):
useNugget(useNugg), nModels(n_models), useLDLT(cholesky){

  z = std::make_shared<Eigen::VectorXd>(y);
  X = std::make_shared<Eigen::MatrixXd>(designMat);
  distHydro = std::make_shared<Eigen::MatrixXd>(N);
  if (D.rows() > 0) distGeo = std::make_shared<Eigen::MatrixXd>(D);
  weightMat = std::make_shared<Eigen::MatrixXd>(wMat);
  flowMat = std::make_shared<Eigen::MatrixXi>(connMat);

  tailUpModel = std::move(tailup_ptr);
  tailDownModel = std::move(taildown_ptr);
  euclidModel = std::move(euclid_ptr);

  int j = 0;
  if (j < nModels*2 && tailUpModel){
    bound_up = bounds[j];
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    bound_down = bounds[j];
    j++;
  }
  if (j < nModels*2 && euclidModel){
    bound_eu = bounds[j];
    j++;
  }

  n = X->rows();
  p = X->cols()-1;

  optimTheta.resize(nModels*2 + 1*useNugget);
  betaValues.resize(p+1);
  covMat.resize(n,n);
  covMat.fill(0.0);

  maxDistHydro = 4.0*distHydro->maxCoeff();
  if (euclidModel) maxDistGeo = 4.0*distGeo->maxCoeff();
}

void Optimizer::setBounds(const std::vector<double>& bounds){
  int j = 0;
  if (j < nModels*2 && tailUpModel){
    bound_up = bounds[j];
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    bound_down = bounds[j];
    j++;
  }
  if (j < nModels*2 && euclidModel){
    bound_eu = bounds[j];
    j++;
  }
}

bool Optimizer::updateParam(const Eigen::VectorXd& theta){
  bool control = true;
  int j = 0;
  if (j < nModels*2 && tailUpModel){
    if (std::exp(theta(j)) < 0.0 || std::exp(theta(j)) > bound_up) control = false;
    tailUpModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro || std::exp(theta(j)) < 0.0) control = false;
    tailUpModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    if (std::exp(theta(j)) < 0.0  || std::exp(theta(j)) > bound_down) control = false;
    tailDownModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro || std::exp(theta(j)) < 0.0) control = false;
    tailDownModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && euclidModel){
    if (std::exp(theta(j)) < 0.0  || std::exp(theta(j)) > bound_eu) control = false;
    euclidModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistGeo || std::exp(theta(j)) < 0.0) control = false;
    euclidModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  return control;
}

Eigen::VectorXd Optimizer::thetaInit() {
  Eigen::VectorXd theta(2*nModels + 1*useNugget);

  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> svd(p+1, p+1, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.compute(X->transpose()*(*X));
  double tol(svd.singularValues()(0)*1.5e-8);
  unsigned int nz(0);
  for (unsigned int i=0; i< svd.singularValues().size(); i++){
    if (svd.singularValues()(i)>tol) nz++;
  }
  Eigen::VectorXd sv(nz);
  for (unsigned int i=0; i<nz; i++)
    sv(i) = svd.singularValues()(i);
  Eigen::MatrixXd U(svd.matrixU().leftCols(nz));
  Eigen::MatrixXd V(svd.matrixV().leftCols(nz));

  Eigen::MatrixXd tmp_inv = V * sv.asDiagonal().inverse() * U.transpose();

  Eigen::VectorXd resid(*z - (*X)*tmp_inv*X->transpose()*(*z));
  resid = resid.array().square();
  double varResid(resid.mean());

  int i = 0;
  if (i < 2*nModels && tailUpModel){
    theta(i) = std::log(0.9/double(nModels) * varResid);
    i++;
    theta(i) = std::log((*distHydro+distHydro->transpose()).mean());
    i++;
  }
  if (i < 2*nModels && tailDownModel){
    theta(i) = std::log(0.9/double(nModels) * varResid);
    i++;
    theta(i) = std::log((*distHydro+distHydro->transpose()).mean());
    i++;
  }
  if (i < 2*nModels && euclidModel){
    theta(i) = std::log(0.9/double(nModels) * varResid);
    i++;
    theta(i) = std::log(distGeo->mean());
    i++;
  }

  if (useNugget){
    theta(i) = std::log(0.1 * varResid);
  }

  //log-scale for every parsill (sigma2) and every range (alpha)
  //log-scale also for the nugget
  return theta;
}

double Optimizer::computeLogL(const Eigen::VectorXd& theta){
  bool check = updateParam(theta);
  if (!check) return 1e+32;

  Eigen::MatrixXd V(n,n);
  V.fill(0.0);
  if (tailUpModel) V += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) V += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) V += euclidModel->computeMatCov(*distGeo);
  if (useNugget) V += Eigen::MatrixXd::Identity(n,n)*std::exp(theta(theta.size()-1));

  Eigen::LDLT<Eigen::MatrixXd> solver(n);
  solver.compute(V);

  Eigen::HouseholderQR<Eigen::MatrixXd> qrV(n,n);
  qrV.compute(V);

  if (!solver.isPositive())
    throw std::domain_error("Covariance matrix not positive definite");

  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  Eigen::MatrixXd invV(n,n);
  double det(V.determinant());
  if (useLDLT){
    invV = solver.solve(Id);
  }
  else if (det >= 1e-3){
    invV = solver.solve(Id);
  }
  else {
    invV = qrV.solve(Id);
    Rcpp::warning("Covariance matrix ill-conditioned. QR decomposition needed.\n");
  }

  Eigen::MatrixXd XVX(X->transpose()*invV*(*X));
  Eigen::MatrixXd invXVX(p+1,p+1);
  Id.resize(p+1,p+1);
  Id.setIdentity();
  if (useLDLT){
    Eigen::LDLT<Eigen::MatrixXd> ldlt(p+1);
    ldlt.compute(XVX);
    invXVX = ldlt.solve(Id);
  }
  else if (XVX.determinant() >= 1e-3){
    Eigen::LDLT<Eigen::MatrixXd> ldlt(p+1);
    ldlt.compute(XVX);
    invXVX = ldlt.solve(Id);
  }
  else {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(p+1, p+1);
    qr.compute(XVX);
    invXVX = qr.solve(Id);
  }

  Eigen::VectorXd beta(invXVX*X->transpose()*invV*(*z));
  Eigen::VectorXd r(*z - (*X)*beta);

  return n*log(2*3.14) + r.transpose()*invV*r + qrV.logAbsDeterminant();

}

std::vector<std::pair<double, Eigen::VectorXd>> Optimizer::simplexInit(const Eigen::VectorXd& theta0, const double tau) {
  unsigned int nParam(theta0.size());
  std::vector<std::pair<double,Eigen::VectorXd>> simplex(nParam+1, std::make_pair(0.0, theta0));
  simplex[0].first = computeLogL(simplex[0].second);
  for (unsigned int i=1; i<nParam+1; i++){
    simplex[i].second(i-1) = theta0(i-1) + tau;
    simplex[i].first = computeLogL(simplex[i].second);
  }
  return simplex;
}

void Optimizer::computeTheta(){

  Eigen::VectorXd theta0(thetaInit());
  unsigned int nParam(theta0.size());
  std::vector<std::pair<double,Eigen::VectorXd>> simplex = simplexInit(theta0, 0.05);

  std::sort(simplex.begin(), simplex.end(), helpers::operandPair);

  Eigen::VectorXd thetaR(nParam);
  double fR(0.0);
  Eigen::VectorXd thetaE(nParam);
  double fE(0.0);
  Eigen::VectorXd thetaC(nParam);
  double fC(0.0);

  double max((simplex[1].second-simplex[0].second).lpNorm<Eigen::Infinity>());
  for (unsigned int i=2; i<nParam+1; i++){
    double norm((simplex[i].second-simplex[0].second).lpNorm<Eigen::Infinity>());
    if (norm > max)
      max = norm;
  }

  bool crit1 = (iter <= maxIter);
  bool crit2 = (funEvals <= maxFunEvals);
  bool crit3 = ((simplex[nParam].first - simplex[0].first) > tolFun);
  bool crit4 = (max > tolTheta);

  // Nelder Mead algorithm
  while (crit1 && crit2 && crit3 && crit4){
    //Calcolo degli 8 theta
    //Calcolo Log-likelihood per tutti e 8
    //Passaggi NelderMead

    fR = 0.0;
    fE = 0.0;
    fC = 0.0;

    //Compunting centroid
    theta0.fill(0.0);
    for (unsigned int i=0; i<nParam; i++){
      for (unsigned int j=0; j<nParam; j++){
        theta0(j) += simplex[i].second(j)/double(nParam);
      }
    }

    //Reflection
    thetaR = theta0 + a*(theta0 - simplex[nParam].second);
    // If the reflected point is better than the second worst, but not better
    // than the best then obtain a new simplex by replacing the worst point
    // with the reflected point
    fR = computeLogL(thetaR);
    if (fR < simplex[nParam-1].first && fR >= simplex[0].first){
      simplex[nParam].second = thetaR;
      simplex[nParam].first = fR;
      // go to order
    }

    //Expansion
    else if (fR < simplex[0].first){
      thetaE = theta0 + c*(theta0 - simplex[nParam].second);
      fE = computeLogL(thetaE);
      if (fE < fR){
        simplex[nParam].second = thetaE;
        simplex[nParam].first = fE;
        funEvals++;
        // go to order
      }
      else {
        simplex[nParam].second = thetaR;
        simplex[nParam].first = fR;
        funEvals++;
        // go to order
      }
    }

    //Contraction (case fR >= simplex[nParam-1].first)
    else {
      thetaC = theta0 + r*(simplex[nParam].second - theta0);
      fC = computeLogL(thetaC);

      if (fC < simplex[nParam].first){
        simplex[nParam].second = thetaC;
        simplex[nParam].first = fC;
        funEvals++;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + s*(simplex[i].second - simplex[0].second);
          simplex[i].first = computeLogL(simplex[i].second);
          funEvals++;
        }
        // go to order
      }
    }

    std::sort(simplex.begin(), simplex.end(), helpers::operandPair);

    max = (simplex[1].second-simplex[0].second).lpNorm<Eigen::Infinity>();
    for (unsigned int i=2; i<nParam+1; i++){
      double norm((simplex[i].second-simplex[0].second).lpNorm<Eigen::Infinity>());
      if (norm > max)
        max = norm;
    }

    iter++;
    crit1 = (iter <= maxIter);
    crit2 = (funEvals <= maxFunEvals);
    crit3 = std::abs(simplex[nParam].first - simplex[0].first) > tolFun;
    crit4 = max > tolTheta;

  }

  if (iter > maxIter) std::cerr << "Reached max number of iterations" << std::endl;
  optimTheta = simplex[0].second.array().exp();

}

void Optimizer::glmssn() {
  computeTheta();

  Eigen::VectorXd logTheta = optimTheta.array().log();
  updateParam(logTheta);

  if (tailUpModel) covMat += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) covMat += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) covMat += euclidModel->computeMatCov(*distGeo);
  if (useNugget) covMat += Eigen::MatrixXd::Identity(n,n)*optimTheta(optimTheta.size()-1);

  Eigen::MatrixXd invV(n,n);
  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  if (useLDLT){
    Eigen::LDLT<Eigen::MatrixXd> solver(n);
    solver.compute(covMat);
    invV = solver.solve(Id);
  }
  else if (covMat.determinant()>1e-3){
    Eigen::LDLT<Eigen::MatrixXd> solver(n);
    solver.compute(covMat);
    invV = solver.solve(Id);
  }
  else {
    Eigen::HouseholderQR<Eigen::MatrixXd> solver(n, n);
    solver.compute(covMat);
    invV = solver.solve(Id);
  }

  Eigen::MatrixXd XVX(X->transpose()*invV*(*X));
  Eigen::MatrixXd invXVX(p+1,p+1);
  Id.resize(p+1,p+1);
  Id.setIdentity();
  if (useLDLT){
    Eigen::LDLT<Eigen::MatrixXd> solver(p+1);
    solver.compute(XVX);
    invXVX = solver.solve(Id);
  }
  else if (XVX.determinant() >= 1e-3){
    Eigen::LDLT<Eigen::MatrixXd> solver(p+1);
    solver.compute(XVX);
    invXVX = solver.solve(Id);
  }
  else {
    Eigen::HouseholderQR<Eigen::MatrixXd> solver(p+1,p+1);
    solver.compute(XVX);
    invXVX = solver.solve(Id);
  }

  betaValues = invXVX*X->transpose()*invV*(*z);
}

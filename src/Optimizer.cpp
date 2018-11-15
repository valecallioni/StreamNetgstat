#include "Optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<TailUpModel>& tailup_ptr, std::unique_ptr<TailDownModel>& taildown_ptr, std::unique_ptr<EuclideanModel>& euclid_ptr, bool useNugg,
  int n_models, const Eigen::VectorXd& y, const Eigen::MatrixXd& designMat, const Eigen::MatrixXd& N, const Eigen::MatrixXd& D, const Eigen::MatrixXd& wMat, const Eigen::MatrixXi& connMat):
useNugget(useNugg), nModels(n_models){

  z = std::make_shared<Eigen::VectorXd>(y);
  X = std::make_shared<Eigen::MatrixXd>(designMat);
  distHydro = std::make_shared<Eigen::MatrixXd>(N);
  distGeo = std::make_shared<Eigen::MatrixXd>(D);
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
  maxDistGeo = 4.0*distGeo->maxCoeff();
}

bool Optimizer::updateParam(const Eigen::VectorXd& theta){
  bool control = true;
  int j = 0;
  if (j < nModels*2 && tailUpModel){
    if (std::exp(theta(j)) < 0.0 || std::exp(theta(j)) > 1e+4) control = false;
    tailUpModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro && std::exp(theta(j)) < 0.0) control = false;
    tailUpModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    if (std::exp(theta(j)) < 0.0  || std::exp(theta(j)) > 1e+4) control = false;
    tailDownModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro && std::exp(theta(j)) < 0.0) control = false;
    tailDownModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && euclidModel){
    if (std::exp(theta(j)) < 0.0  || std::exp(theta(j)) > 1e+4) control = false;
    euclidModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistGeo && std::exp(theta(j)) < 0.0) control = false;
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
  double logl(1e+32);
  bool check = updateParam(theta);

  Eigen::MatrixXd V(n,n);
  V.fill(0.0);
  if (tailUpModel) V += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) V += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) V += euclidModel->computeMatCov(*distGeo);
  if (useNugget) V += Eigen::MatrixXd::Identity(n,n)*std::exp(theta(theta.size()-1));


  Eigen::LDLT<Eigen::MatrixXd> solver(n);
  solver.compute(V);

  if (!solver.isPositive())
    throw std::domain_error("Covariance matrix not positive definite");

  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  Eigen::MatrixXd invV(solver.solve(Id));
  std::cout << "invV" << invV.block(0,0,3,3) << std::endl;

  solver = Eigen::LDLT<Eigen::MatrixXd>(p+1);
  solver.compute(X->transpose()*invV*(*X));
  Id.resize(p+1,p+1);
  Id.setIdentity();
  Eigen::MatrixXd invXVX(solver.solve(Id));

  Eigen::VectorXd beta(invXVX*X->transpose()*invV*(*z));
  Eigen::VectorXd r(*z - (*X)*beta);

  Eigen::HouseholderQR<Eigen::MatrixXd> qrV(n,n);
  qrV.solve(V);
  std::cout << "det(V) = " << qrV.absDeterminant() << std::endl;

  if (check) logl = n*log(2*3.14) + r.transpose()*invV*r + qrV.logAbsDeterminant();
  return logl;

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

void Optimizer::computeThetaWiki(){

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

  double a(1.0);
  double c(2.0);
  double r(0.5);
  double s(0.5);

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

    std::cout << "Entered while loop." << std::endl;

    fR = 0.0;
    fE = 0.0;
    fC = 0.0;

    //Order
    //std::sort(simplex.begin(), simplex.end(), helpers::operandPair);

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
    std::cout << "Print after computing fR. fR = " << fR << std::endl;
    if (fR < simplex[nParam-1].first && fR >= simplex[0].first){
      simplex[nParam].second = thetaR;
      simplex[nParam].first = fR;
      std::cout << "Reflexion." << std::endl;
      // go to order
    }

    //Expansion
    else if (fR < simplex[0].first){
      thetaE = theta0 + c*(theta0 - simplex[nParam].second);
      //thetaE = theta0 + c*(theta0 - simplex[nParam].second);
      fE = computeLogL(thetaE);
      std::cout << "Print after computing fE. fE = " << fE << std::endl;
      if (fE < fR){
        simplex[nParam].second = thetaE;
        simplex[nParam].first = fE;
        funEvals++;
        std::cout << "Expansion." << std::endl;
        // go to order
      }
      else {
        simplex[nParam].second = thetaR;
        simplex[nParam].first = fR;
        funEvals++;
        std::cout << "Reflexion." << std::endl;
        // go to order
      }
    }

    //Contraction (case fR >= simplex[nParam-1].first)
    else {
      thetaC = theta0 + r*(simplex[nParam].second - theta0);
      fC = computeLogL(thetaC);
      std::cout << "Print after computing fC. fC = "<< fC << std::endl;

      if (fC < simplex[nParam].first){
        simplex[nParam].second = thetaC;
        simplex[nParam].first = fC;
        funEvals++;
        std::cout << "Contraction." << std::endl;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + s*(simplex[i].second - simplex[0].second);
          simplex[i].first = computeLogL(simplex[i].second);
          funEvals++;
        }
        std::cout << "Shrink." << std::endl;
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
  std::cout << "Minimization terminated with " << iter << " iterations." << std::endl;
  optimTheta = simplex[0].second.array().exp();

}

void Optimizer::computeThetaPaper(){

  Eigen::VectorXd theta0(thetaInit());
  unsigned int nParam(theta0.size());
  std::vector<std::pair<double,Eigen::VectorXd>> simplex = simplexInit(theta0, 0.05);

  // Ordering
  std::sort(simplex.begin(), simplex.end(), helpers::operandPair);

  double a(1.0);
  double b(1.0 + 2.0/double(nParam));
  double c(0.75 - 1.0/(2.0*double(nParam)));
  double d(1.0 - 1.0/double(nParam));

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

  Eigen::VectorXd thetaR(nParam);
  double fR(0.0);
  Eigen::VectorXd thetaE(nParam);
  double fE(0.0);
  Eigen::VectorXd thetaOC(nParam);
  double fOC(0.0);
  Eigen::VectorXd thetaIC(nParam);
  double fIC(0.0);


  // Nelder Mead algorithm
  while (crit1 && crit2 && crit3 && crit4){

    fR = 0.0;
    fE = 0.0;
    fOC = 0.0;
    fIC = 0.0;

    //Calcolo degli 8 theta
    //Calcolo Log-likelihood per tutti e 8
    //Passaggi NelderMead

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
      funEvals++;
      // go to order
    }

    //Expansion
    else if (fR < simplex[0].first){
      thetaE = theta0 + b*(thetaR - theta0);
      fE = computeLogL(thetaE);
      if (fE != 1e+32 && fE < fR){
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

    //Outside Contraction (f(xn)<=fr<f(xn+1))
    else if (simplex[nParam-1].first <= fR && fR < simplex[nParam].first){
      thetaOC = theta0 + c*(thetaR - theta0);
      fOC = computeLogL(thetaOC);

      if (fOC <= fR){
        simplex[nParam].second = thetaOC;
        simplex[nParam].first = fOC;
        funEvals++;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + d*(simplex[i].second - simplex[0].second);
          simplex[i].first = computeLogL(simplex[i].second);
          funEvals++;
        }
        // go to order
      }
    }

    else { // Inside Contraction (case: fR >= f(x(n+1)))
      thetaIC = theta0 - c*(thetaR - theta0);
      fIC = computeLogL(thetaIC);
      if (fIC != 1e+32 && fIC <= simplex[nParam].first){
        simplex[nParam].second = thetaIC;
        simplex[nParam].first = fIC;
        funEvals++;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + d*(simplex[i].second - simplex[0].second);
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
    crit3 = std::abs(simplex[nParam].first - simplex[0].first) > 0.001;
    crit4 = max > tolTheta;

  }

  if (iter > maxIter) std::cerr << "Reached max number of iterations" << std::endl;
  std::cout << "Minimization terminated with " << iter << " iterations." << std::endl;
  optimTheta = simplex[0].second.array().exp();

}

void Optimizer::glmssn() {
  computeThetaWiki();

  Eigen::VectorXd logTheta = optimTheta.array().log();
  updateParam(logTheta);

  if (tailUpModel) covMat += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) covMat += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) covMat += euclidModel->computeMatCov(*distGeo);
  if (useNugget) covMat += Eigen::MatrixXd::Identity(n,n)*optimTheta(optimTheta.size()-1);

  Eigen::LDLT<Eigen::MatrixXd> solver(n);
  solver.compute(covMat);
  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  Eigen::MatrixXd invV(solver.solve(Id));

  solver = Eigen::LDLT<Eigen::MatrixXd>(p);
  solver.compute(X->transpose()*invV*(*X));
  Id.resize(p+1,p+1);
  Id.setIdentity();
  Eigen::MatrixXd invXVX(solver.solve(Id));

  betaValues = invXVX*X->transpose()*invV*(*z);
}

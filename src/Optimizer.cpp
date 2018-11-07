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

Eigen::VectorXd Optimizer::thetaInit() {
  const double eps(std::numeric_limits<double>::epsilon());
  Eigen::VectorXd theta(2*nModels + 1*useNugget);

  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> svd(p+1, p+1, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.compute(X->transpose()*(*X));
  double tol(svd.singularValues()(0)*1.5e-8);
  unsigned int nz(0);
  for (unsigned int i=0; i< svd.singularValues().size(); i++){
    if (svd.singularValues()(i)>tol) nz++;
  }
  std::cout << "nz = " << nz << std::endl;
  Eigen::VectorXd sv(nz);
  for (unsigned int i=0; i<nz; i++)
    sv(i) = svd.singularValues()(i);
  Eigen::MatrixXd U(svd.matrixU().leftCols(nz));
  Eigen::MatrixXd V(svd.matrixV().leftCols(nz));

  Eigen::MatrixXd tmp_inv = V * sv.asDiagonal().inverse() * U.transpose();

  std::cout << "tmp_inv: \n" << tmp_inv << std::endl;
  Eigen::VectorXd resid(*z - (*X)*tmp_inv*X->transpose()*(*z));
  resid = resid.array().square();
  double varResid(resid.mean());
  std::cout << "varResid = " << varResid << std::endl;

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

double Optimizer::computeLogL(Eigen::VectorXd& theta){
  double logl(1e+32);
  bool check = true;
  int j = 0;
  if (j < nModels*2 && tailUpModel){
    tailUpModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro) check = false;
    tailUpModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    tailDownModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistHydro) check = false;
    tailDownModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && euclidModel){
    euclidModel->setSigma2(std::exp(theta(j)));
    j++;
    if (std::exp(theta(j)) > maxDistGeo) check = false;
    euclidModel->setAlpha(std::exp(theta(j)));
    j++;
  }

  Eigen::MatrixXd V(n,n);
  V.fill(0.0);
  if (tailUpModel) V += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) V += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) V += euclidModel->computeMatCov(*distGeo);
  if (useNugget) V += Eigen::MatrixXd::Identity(n,n)*std::exp(theta(j));

  // Da correggere: posso usare solo isPositive() di LDLT
  Eigen::EigenSolver<Eigen::MatrixXd> eig(n);
  eig.compute(V);
  Eigen::VectorXd values(eig.eigenvalues().real());
  for (unsigned int i=0; i<values.size(); i++){
    if (values(i) < 0.0){
      std::cerr << "Covariance matrix not positive definite" << std::endl;
      check = false;
      break;
    }
  }

  Eigen::LDLT<Eigen::MatrixXd> solver(n);
  solver.compute(V);
  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  Eigen::MatrixXd invV(solver.solve(Id));

  solver = Eigen::LDLT<Eigen::MatrixXd>(p+1);
  solver.compute(X->transpose()*invV*(*X));
  Id.resize(p+1,p+1);
  Id.setIdentity();
  Eigen::MatrixXd invXVX(solver.solve(Id));

  Eigen::VectorXd beta(invXVX*X->transpose()*invV*(*z));
  Eigen::VectorXd r(*z - (*X)*beta);

  Eigen::HouseholderQR<Eigen::MatrixXd> qrV(n,n);
  qrV.solve(V);
  if (check) logl = n*log(2*3.14) + r.transpose()*invV*r + std::log(V.determinant());
  return logl;

}

std::vector<std::pair<double, Eigen::VectorXd>> Optimizer::simplexInit(const Eigen::VectorXd& theta0, const double tau) {
  unsigned int nParam(theta0.size());
  std::vector<std::pair<double,Eigen::VectorXd>> simplex(nParam+1, std::make_pair(0.0, theta0));
  simplex[0].first = computeLogL(simplex[0].second);
  std::cout << "Theta: \n" << simplex[0].second.array().exp() << "with logl = " << simplex[0].first << std::endl;
  for (unsigned int i=1; i<nParam+1; i++){
    simplex[i].second(i-1) = theta0(i-1) + tau;
    simplex[i].first = computeLogL(simplex[i].second);
    std::cout << "Theta: \n" << simplex[i].second.array().exp() << "with logl = " << simplex[i].first << std::endl;
  }
  return simplex;
}

void Optimizer::computeThetaWiki(){

  Eigen::VectorXd theta0(thetaInit());
  std::cout << "Initial theta : \n" << theta0.array().exp() << std::endl;
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

  bool crit1 = (iter <= maxIter);
  bool crit2 = (funEvals <= maxFunEvals);

  double max((simplex[1].second-simplex[0].second).lpNorm<Eigen::Infinity>());
  for (unsigned int i=2; i<nParam+1; i++){
    double norm((simplex[i].second-simplex[0].second).lpNorm<Eigen::Infinity>());
    if (norm > max)
      max = norm;
  }

  bool crit3 = ((simplex[nParam].first - simplex[0].first) > tolFun);
  bool crit4 = (max > tolTheta);

  std::cout << "f(x(n+1)) = " << simplex[nParam].first << ", f(x(1)) = " << simplex[0].first << std::endl;
  std::cout << "Difference = " << simplex[nParam].first - simplex[0].first << std::endl;

  // Nelder Mead algorithm
  while (crit1 && crit2 && crit3){
    std::cout << "Entered in the while loop" << std::endl;
    //Calcolo degli 8 theta
    //Calcolo Log-likelihood per tutti e 8
    //Passaggi NelderMead

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
    std::cout << "fR = " << fR << std::endl;
    if (fR < simplex[nParam-1].first && fR >= simplex[0].first){
      simplex[nParam].second = thetaR;
      simplex[nParam].first = fR;
      std::cout << "Reflection." << std::endl;
      // go to order
    }

    //Expansion
    else if (fR < simplex[0].first){
      thetaE = theta0 + c*(theta0 - simplex[nParam].second);
      //thetaE = theta0 + c*(theta0 - simplex[nParam].second);
      fE = computeLogL(thetaE);
      std::cout << "fE = " << fE << std::endl;
      if (fE < fR){
        simplex[nParam].second = thetaE;
        simplex[nParam].first = fE;
        funEvals++;
        std::cout << "Expansion1." << std::endl;
        // go to order
      }
      else {
        simplex[nParam].second = thetaR;
        simplex[nParam].first = fR;
        funEvals++;
        std::cout << "Expansion2." << std::endl;
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
    crit3 = std::abs(simplex[nParam].first - simplex[0].first) > 0.001;
    crit4 = max > tolTheta;
    std::cout << "Best funEval = " << simplex[0].first << std::endl;

  }

  std::cout << "Exited while loop. Iter = " << iter << ", funEvals = " << funEvals << std::endl;
  std::cout << "Diff. simplex points = " << max << ", diff. LogL values = " << simplex[nParam].first - simplex[0].first << std::endl;

  optimTheta = simplex[0].second.array().exp();

}

void Optimizer::computeThetaPaper(){

  Eigen::VectorXd theta0(thetaInit());
  std::cout << "Initial theta : \n" << theta0.array().exp() << std::endl;
  unsigned int nParam(theta0.size());
  std::vector<std::pair<double,Eigen::VectorXd>> simplex = simplexInit(theta0, 0.05);

  // Ordering
  std::sort(simplex.begin(), simplex.end(), helpers::operandPair);

  double a(1.0);
  double b(1.0 + 2.0/double(nParam));
  double c(0.75 - 1.0/(2.0*double(nParam)));
  double d(1.0 - 1.0/double(nParam));

  bool crit1 = (iter <= maxIter);
  bool crit2 = (funEvals <= maxFunEvals);

  double max((simplex[1].second-simplex[0].second).lpNorm<Eigen::Infinity>());
  for (unsigned int i=2; i<nParam+1; i++){
    double norm((simplex[i].second-simplex[0].second).lpNorm<Eigen::Infinity>());
    if (norm > max)
      max = norm;
  }

  bool crit3 = ((simplex[nParam].first - simplex[0].first) > tolFun);
  bool crit4 = (max > tolTheta);
  bool crit5 = (restart <= maxRestart);

  std::cout << "f(x(n+1)) = " << simplex[nParam].first << ", f(x(1)) = " << simplex[0].first << std::endl;
  std::cout << "Difference = " << simplex[nParam].first - simplex[0].first << std::endl;

  Eigen::VectorXd thetaR(nParam);
  double fR(0.0);
  Eigen::VectorXd thetaE(nParam);
  double fE(0.0);
  Eigen::VectorXd thetaOC(nParam);
  double fOC(0.0);
  Eigen::VectorXd thetaIC(nParam);
  double fIC(0.0);

  int modified(1);

  // Nelder Mead algorithm
  while (crit1 && crit2 && crit3 && crit4 && crit5){
    std::cout << "Entered in the while loop" << std::endl;
    modified = 0;

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
    std::cout << "fR = " << fR << std::endl;
    if (fR < simplex[nParam-1].first && fR >= simplex[0].first){
      simplex[nParam].second = thetaR;
      simplex[nParam].first = fR;
      funEvals++;
      std::cout << "Reflection." << std::endl;
      //modified = 1;
      // go to order
    }

    //Expansion
    else if (fR < simplex[0].first){
      thetaE = theta0 + b*(thetaR - theta0);
      fE = computeLogL(thetaE);
      std::cout << "fE = " << fE << std::endl;
      if (fE != 1e+32 && fE < fR){
        simplex[nParam].second = thetaE;
        simplex[nParam].first = fE;
        funEvals++;
        std::cout << "Expansion1." << std::endl;
        //modified = 1;
        // go to order
      }
      else {
        simplex[nParam].second = thetaR;
        simplex[nParam].first = fR;
        funEvals++;
        std::cout << "Expansion2." << std::endl;
        //modified = 1;
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
        std::cout << "Outside Contraction." << std::endl;
        //modified = 1;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + d*(simplex[i].second - simplex[0].second);
          simplex[i].first = computeLogL(simplex[i].second);
          funEvals++;
        }
        std::cout << "Shrink." << std::endl;
        //modified = 1;
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
        std::cout << "Inside Contraction." << std::endl;
        //modified = 1;
        // go to order
      }
      else {
        //Shrink
        for (unsigned int i=1; i<nParam+1; i++){
          simplex[i].second = simplex[0].second + d*(simplex[i].second - simplex[0].second);
          simplex[i].first = computeLogL(simplex[i].second);
          funEvals++;
        }
        std::cout << "Shrink." << std::endl;
        //modified = 1;
        // go to order
      }

    }

    // if (!modified) {
    //   simplex = simplexInit(simplex[0].second, 0.5);
    //   restart++;
    //   modified = 1;
    //   std::cout << "Restart." << std::endl;
    // }

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
    crit5 = (restart <= maxRestart);

    std::cout << "Best funEval = " << simplex[0].first << std::endl;

  }

  std::cout << "Exited while loop. Iter = " << iter << ", funEvals = " << funEvals << ", restarts = " << restart << std::endl;
  std::cout << "Diff. simplex points = " << max << ", diff. LogL values = " << simplex[nParam].first - simplex[0].first << std::endl;

  optimTheta = simplex[0].second.array().exp();

}

void Optimizer::glmssn() {
  computeThetaWiki();

  int j = 0;
  if (j < nModels*2 && tailUpModel){
    tailUpModel->setSigma2(optimTheta(j));
    j++;
    tailUpModel->setAlpha(optimTheta(j));
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    tailDownModel->setSigma2(optimTheta(j));
    j++;
    tailDownModel->setAlpha(optimTheta(j));
    j++;
  }
  if (j < nModels*2 && euclidModel){
    euclidModel->setSigma2(optimTheta(j));
    j++;
    euclidModel->setAlpha(optimTheta(j));
    j++;
  }

  if (tailUpModel) covMat += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) covMat += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) covMat += euclidModel->computeMatCov(*distGeo);
  if (useNugget) covMat += Eigen::MatrixXd::Identity(n,n)*optimTheta(j);


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

void Optimizer::test(const Eigen::VectorXd& theta){
  int j = 0;
  if (j < nModels*2 && tailUpModel){
    tailUpModel->setSigma2(std::exp(theta(j)));
    j++;
    tailUpModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && tailDownModel){
    tailDownModel->setSigma2(std::exp(theta(j)));
    j++;
    tailDownModel->setAlpha(std::exp(theta(j)));
    j++;
  }
  if (j < nModels*2 && euclidModel){
    euclidModel->setSigma2(std::exp(theta(j)));
    j++;
    euclidModel->setAlpha(std::exp(theta(j)));
    j++;
  }

  if (tailUpModel) covMat += tailUpModel->computeMatCov(*weightMat, *distHydro);
  if (tailDownModel) covMat += tailDownModel->computeMatCov(*flowMat, *distHydro);
  if (euclidModel) covMat += euclidModel->computeMatCov(*distGeo);
  if (useNugget) covMat += Eigen::MatrixXd::Identity(n,n)*exp(theta(j));


  Eigen::LDLT<Eigen::MatrixXd> solver(n);
  solver.compute(covMat);
  Eigen::MatrixXd Id(n,n);
  Id.setIdentity();
  Eigen::MatrixXd invV(solver.solve(Id));

  solver = Eigen::LDLT<Eigen::MatrixXd>(p+1);
  solver.compute(X->transpose()*invV*(*X));
  Id.resize(p+1,p+1);
  Id.setIdentity();
  Eigen::MatrixXd invXVX(solver.solve(Id));

  betaValues = invXVX*X->transpose()*invV*(*z);
}

void Optimizer::setTheta(const Eigen::VectorXd& theta){
  optimTheta = theta;
  test(theta.array().log());
}

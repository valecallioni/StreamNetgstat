#include<vector>
#include<utility>
#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"
#include "FactoryHelpers.hpp"

namespace {
  using tailup_factory::TailUpProxy;
  using taildown_factory::TailDownProxy;
  using euclidean_factory::EuclideanProxy;

  TailUpProxy<LinearWithSillTU> linearTU("LinearSill.tailup");
  TailUpProxy<SphericalTU> sphericalTU("Spherical.tailup");
  TailUpProxy<ExponentialTU> exponentialTU("Exponential.tailup");
  TailUpProxy<MariahTU> mariahTU("Mariah.tailup");

  TailDownProxy<LinearWithSillTD> linearTD("LinearSill.taildown");
  TailDownProxy<SphericalTD> sphericalTD("Spherical.taildown");
  TailDownProxy<ExponentialTD> exponentialTD("Exponential.taildown");
  TailDownProxy<MariahTD> mariahTD("Mariah.taildown");

  EuclideanProxy<CauchyEU> cauchyEU("Cauchy.Euclid");
  EuclideanProxy<SphericalEU> sphericalEU("Spherical.Euclid");
  EuclideanProxy<ExponentialEU> exponentialEU("Exponential.Euclid");
  EuclideanProxy<GaussianEU> gaussianEU("Gaussian.Euclid");
}

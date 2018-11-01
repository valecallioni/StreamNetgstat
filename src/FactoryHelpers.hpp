#ifndef PACSPROJECT_HELPERS_FACTORY_HPP
#define PACSPROJECT_HELPERS_FACTORY_HPP

#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"
#include "Factory.hpp"
#include "Proxy.hpp"

namespace tailup_factory {
  typedef generic_factory::Factory<TailUpModel, std::string> TailUpFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using TailUpProxy = generic_factory::Proxy<TailUpFactory,ConcreteProduct>;
}

namespace taildown_factory {
  typedef generic_factory::Factory<TailDownModel, std::string> TailDownFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using TailDownProxy = generic_factory::Proxy<TailDownFactory,ConcreteProduct>;
}

namespace euclidean_factory {
  typedef generic_factory::Factory<EuclideanModel, std::string> EuclideanFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using EuclideanProxy = generic_factory::Proxy<EuclideanFactory,ConcreteProduct>;
}

#endif

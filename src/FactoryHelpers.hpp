#ifndef PACSPROJECT_HELPERS_FACTORY_HPP
#define PACSPROJECT_HELPERS_FACTORY_HPP

#include "TailUpModel.hpp"
#include "TailDownModel.hpp"
#include "EuclideanModel.hpp"
#include "Factory.hpp"
#include "Proxy.hpp"

/**
 * Typedefs for the factory
*/

namespace tailup_factory {
  /**
  * \var typedef generic_factory::Factory<TailUpModel, std::string> TailUpFactory;
  * Factory for the tail-up model
  */
  typedef generic_factory::Factory<TailUpModel, std::string> TailUpFactory;  // Use standard Builder

  /**
  * Proxy for the tail-up model
  */
  template<typename ConcreteProduct>
  using TailUpProxy = generic_factory::Proxy<TailUpFactory,ConcreteProduct>;
}

namespace taildown_factory {
  /**
  * \var typedef generic_factory::Factory<TailDownModel, std::string> TailDownFactory;
  * Factory for the tail-down model
  */
  typedef generic_factory::Factory<TailDownModel, std::string> TailDownFactory;  // Use standard Builder

  /**
  * Proxy for the tail-down model
  */
  template<typename ConcreteProduct>
  using TailDownProxy = generic_factory::Proxy<TailDownFactory,ConcreteProduct>;
}

namespace euclidean_factory {
  /**
  * \var typedef generic_factory::Factory<EuclideanModel, std::string> EuclideanFactory;
  * Factory for the Euclidean model
  */
  typedef generic_factory::Factory<EuclideanModel, std::string> EuclideanFactory;  // Use standard Builder

  /**
  * Proxy for the Euclidean model
  */
  template<typename ConcreteProduct>
  using EuclideanProxy = generic_factory::Proxy<EuclideanFactory,ConcreteProduct>;
}

#endif

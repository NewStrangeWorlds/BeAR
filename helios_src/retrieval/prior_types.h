/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _prior_types_h
#define _prior_types_h

#include <string>
#include <cmath>
#include <vector>

#define BOOST_IF_CONSTEXPR if
#include "../../_deps/boost_math-src/include/boost/math/special_functions/erf.hpp"


namespace helios {


//types of priors
enum class PriorType { uniform, log_uniform, gaussian, delta };

namespace priors{  
  const std::vector<PriorType> prior_types{
    PriorType::uniform, 
    PriorType::log_uniform, 
    PriorType::gaussian, 
    PriorType::delta};
  const std::vector<std::string> prior_type_strings{
    "uniform", 
    "log_uniform", 
    "gaussian", 
    "delta"};
}


//abstract class describing a general prior
class BasicPrior{
  public:
    virtual double parameterValue(const double& hypercube_parameter) = 0; 
    virtual ~BasicPrior() {}
    std::string distributionType() {return distribution_type;}
    std::string parameterName() {return parameter_name;}
    virtual void printInfo() = 0;
  protected:
    std::string distribution_type = "";
    std::string parameter_name = "";
};


//simple uniform prior
//the constructor requires three parameters: the upper and lower domain boundary and the parameter name
class UniformPrior : public BasicPrior {
  public:
    UniformPrior (
      const std::string& parameter,
      const double& upper_bnd,
      const double& lower_bnd) {
        parameter_name = parameter; 
        upper_bound = upper_bnd; 
        lower_bound = lower_bnd; 
        distribution_type = "uniform prior";}
    virtual ~UniformPrior() {}
    virtual double parameterValue (const double& hypercube_value) {
        return lower_bound + hypercube_value * (upper_bound - lower_bound);}
    virtual void printInfo() {
      std::cout << parameter_name << ": " << distribution_type
        << ", lower boundary: " << lower_bound 
        << ", upper boundary: " << upper_bound << "\n";
    }
  private:
    double lower_bound = 0;
    double upper_bound = 0;
};


//uniform prior in log space
//the constructor requires three parameters: the upper and lower domain boundary and the parameter name
class LogUniformPrior : public BasicPrior {
  public:
    LogUniformPrior (
      const std::string& parameter,
      const double& upper_bnd,
      const double& lower_bnd) {
        parameter_name = parameter; 
        log_upper_bound = std::log10(upper_bnd); 
        log_lower_bound = std::log10(lower_bnd); 
        distribution_type = "log-uniform prior";}
    virtual ~LogUniformPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
        return std::pow(10.0, log_lower_bound + hypercube_value * (log_upper_bound - log_lower_bound));}
    virtual void printInfo() {
      std::cout << parameter_name << ":" << distribution_type
        << ", lower boundary: " << log_lower_bound 
        << ", upper boundary: " << log_upper_bound << "\n";
    }
  private:
    double log_lower_bound = 0;
    double log_upper_bound = 0;
};


//Gaussian prior
//the constructor requires three parameters: the mean value of the Gaussian, its standard deviation, and the parameter name
class GaussianPrior : public BasicPrior {
  public:
    GaussianPrior (
      const std::string& parameter,
      const double& mean,
      const double& standard_deviation) {
        parameter_name = parameter; 
        mu = mean; 
        sigma = standard_deviation; 
        distribution_type = "Gaussian prior";}
    virtual ~GaussianPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
      return mu + sigma * std::sqrt(2.0) * boost::math::erf_inv(2.*hypercube_value - 1.);}
    virtual void printInfo() {
      std::cout << parameter_name << ": " << distribution_type
        << ", mean: " << mu 
        << ", standard deviation: " << sigma << "\n";
    }
  private:
    double mu = 0;
    double sigma = 0;
};


//Delta distribution (constant) prior
//the constructor requires one parameters: constant parameter value
class DeltaPrior : public BasicPrior {
  public:
    DeltaPrior (
      const std::string& parameter,
      const double& value) {
        parameter_name = parameter; const_value = value; distribution_type = "Delta prior";}
    virtual ~DeltaPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
      return const_value;}
    virtual void printInfo() {
      std::cout << parameter_name << ": " << distribution_type
        << ", constant valuie: " << const_value << "\n";
    }
  private:
    double const_value = 0;
};


}

#endif

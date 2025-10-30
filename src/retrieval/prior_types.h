/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAR is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAR directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _prior_types_h
#define _prior_types_h

#include <string>
#include <cmath>
#include <vector>

#include "../additional/physical_const.h"
#include "../additional/exceptions.h"

#define BOOST_IF_CONSTEXPR if
#include "../../_deps/boost_math-src/include/boost/math/special_functions/erf.hpp"



namespace bear {

//abstract class describing a general prior
class BasicPrior{
  public:
    virtual double parameterValue(const double& hypercube_parameter) = 0;
    double parameterPhysicalValue(const double& hypercube_parameter) {
      return parameterValue(hypercube_parameter) * unit;};
    virtual double applyParameterUnit(const double& parameter) {
      return parameter * unit;};
    virtual ~BasicPrior() {}
    std::string distributionType() {return distribution_type;}
    std::string parameterName() {return parameter_name;}
    virtual void printInfo() = 0;
  protected:
    std::string distribution_type = "";
    std::string parameter_name = "";

    double unit = 1.0;
    std::string unit_description = "";

    void setUnit() {
      if (unit_description == "") return;

      if (unit_description == "Me" ||  unit_description == "Mearth") {
        unit = constants::mass_earth;
      } else if (unit_description == "Mj" || unit_description == "Mjupiter") {
        unit = constants::mass_jupiter;
      } else if (unit_description == "Ms" || unit_description == "Msun") {
        unit = constants::mass_sun;
      } else if (unit_description == "Re" || unit_description == "Rearth") {
        unit = constants::radius_earth;
      } else if (unit_description == "Rj" || unit_description == "Rjupiter") {
        unit = constants::radius_jupiter;
      } else if (unit_description == "Rs" || unit_description == "Rsun") {
        unit = constants::radius_sun;
      } else if (unit_description == "ly") {
        unit = constants::light_year;
      } else if (unit_description == "pc") {
        unit = constants::parsec;
      } else {
        std::string error_message = "Unit " + unit_description 
                                  + " for parameter " + parameter_name 
                                  + " unknown!\n";
        throw InvalidInput(std::string ("pirors.config"), error_message);
      }
    };
};


//simple uniform prior
//the constructor requires three parameters: the upper and lower domain boundary and the parameter name
//additionally, it can have an optional unit
class UniformPrior : public BasicPrior {
  public:
    UniformPrior (
      const std::string parameter,
      const double upper_bnd,
      const double lower_bnd,
      const std::string unit_desc) {
        parameter_name = parameter;

        upper_bound = upper_bnd; 
        lower_bound = lower_bnd;

        if (upper_bound < lower_bound)
          std::swap(upper_bound, lower_bound); 

        unit_description = unit_desc;
        distribution_type = "uniform prior";
        
        setUnit();
      };
    virtual ~UniformPrior() {}
    virtual double parameterValue (const double& hypercube_value) {
        return lower_bound + hypercube_value * (upper_bound - lower_bound);}
    virtual void printInfo() {
      std::cout << std::setw(15) << std::left << parameter_name << "  "
                << std::setw(20) << std::left << distribution_type
                << "lower boundary: " << lower_bound 
                << ", upper boundary: " << upper_bound;
      if (unit_description != "")
        std::cout << ", unit: " << unit_description;
      std::cout << "\n";};
  private:
    double lower_bound = 0;
    double upper_bound = 0;
};


//uniform prior in log space
//the constructor requires three parameters: the upper and lower domain boundary and the parameter name
//additionally, it can have an optional unit
class LogUniformPrior : public BasicPrior {
  public:
    LogUniformPrior (
      const std::string parameter,
      const double upper_bnd,
      const double lower_bnd,
      const std::string unit_desc) {
        parameter_name = parameter;

        log_upper_bound = std::log10(upper_bnd);
        log_lower_bound = std::log10(lower_bnd);

        if (log_upper_bound < log_lower_bound)
          std::swap(log_upper_bound, log_lower_bound);

        unit_description = unit_desc;
        distribution_type = "log-uniform prior";

        setUnit();
      };
    virtual ~LogUniformPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
        return std::pow(10.0, log_lower_bound + hypercube_value * (log_upper_bound - log_lower_bound));}
    virtual void printInfo() {
      std::cout << std::setw(15) << std::left << parameter_name << "  "
                << std::setw(20) << std::left << distribution_type
                << "lower boundary: " << log_lower_bound 
                << ", upper boundary: " << log_upper_bound;
      if (unit_description != "")
        std::cout << ", unit: " << unit_description;
      std::cout << "\n";};
  private:
    double log_lower_bound = 0;
    double log_upper_bound = 0;
};


//Gaussian prior
//the constructor requires three parameters: the mean value of the Gaussian, its standard deviation, and the parameter name
//additionally, it can have an optional unit
class GaussianPrior : public BasicPrior {
  public:
    GaussianPrior (
      const std::string parameter,
      const double mean,
      const double standard_deviation,
      const std::string unit_desc) {
        parameter_name = parameter; 
        mu = mean; 
        sigma = standard_deviation;
        unit_description = unit_desc;
        distribution_type = "Gaussian prior";

        setUnit();
      };
    virtual ~GaussianPrior() {}
    inline virtual double parameterValue(const double& hypercube_value) {
      return mu + sigma * std::sqrt(2.0) * boost::math::erf_inv(2.*hypercube_value - 1.);}
    virtual void printInfo() {
      std::cout << std::setw(15) << std::left << parameter_name << "  "
                << std::setw(20) << std::left << distribution_type
                << "mean: " << mu 
                << ", standard deviation: " << sigma;
       if (unit_description != "")
        std::cout << ", unit: " << unit_description;
      std::cout << "\n";};
  private:
    double mu = 0;
    double sigma = 0;
};


//Delta distribution (constant) prior
//the constructor requires one parameters: constant parameter value
//additionally, it can have an optional unit
class DeltaPrior : public BasicPrior {
  public:
    DeltaPrior (
      const std::string parameter,
      const double value,
      const std::string unit_desc) {
        parameter_name = parameter; 
        const_value = value; 
        unit_description = unit_desc; 
        distribution_type = "Delta prior";

        setUnit();
      };
    virtual ~DeltaPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
      return const_value;}
    virtual void printInfo() {
      std::cout << std::setw(15) << std::left << parameter_name << "  "
                << std::setw(20) << std::left << distribution_type
                << "constant value: " << const_value;
      if (unit_description != "")
        std::cout << ", unit: " << unit_description;
      std::cout << "\n";};
  private:
    double const_value = 0;
};



//Prior linked to another prior
//the constructor requires one parameters: the linked prior
class LinkedPrior : public BasicPrior {
  public:
    LinkedPrior (
      const std::string parameter,
      BasicPrior* prior)
      : linked_prior(prior){ 
        parameter_name = parameter;
        distribution_type = "Linked prior";};
    virtual ~LinkedPrior() {}
    virtual double parameterValue(const double& hypercube_value) {
      return linked_prior->parameterValue(hypercube_value);}
    virtual double applyParameterUnit(const double& parameter) {
      return linked_prior->applyParameterUnit(parameter);};
    virtual void printInfo() {
      std::cout << std::setw(15) << std::left << parameter_name << "  "
                << std::setw(20) << std::left << distribution_type 
                << "linked to: "; 
      linked_prior->printInfo();
    }
  private:
    BasicPrior* linked_prior;
};


}

#endif

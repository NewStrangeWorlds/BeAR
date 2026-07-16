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


#ifndef _temperature_h
#define _temperature_h

#include <vector>
#include <string>
#include <cstddef>


namespace bear {


//how the free parameters of a control-point temperature profile map to the
//temperatures at those points:
//  absolute - each parameter is the temperature at its control point
//  relative - parameter[0] is the bottom temperature, each further parameter is
//             a multiplicative factor applied to the level below it
enum class Parametrisation { absolute, relative };


class Temperature{
  public:
    virtual ~Temperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature_profile) = 0;
    size_t nbParameters() {return parameter_names.size();}
    const std::vector<std::string>& parameterNames() const {return parameter_names;}
  protected:
    std::vector<std::string> parameter_names;
    Parametrisation parametrisation = Parametrisation::relative;

    //convert the free parameters (index 0 = bottom/deepest level) into the
    //absolute temperatures at the control points, according to `parametrisation`
    std::vector<double> controlTemperatures(
      const std::vector<double>& parameters) const {
        std::vector<double> temperatures(parameters.size(), 0.0);

        if (parameters.empty()) return temperatures;

        temperatures[0] = parameters[0];

        for (size_t i=1; i<parameters.size(); ++i)
          temperatures[i] = (parametrisation == Parametrisation::relative)
            ? temperatures[i-1] * parameters[i]
            : parameters[i];

        return temperatures;
      };

    //generate parameter_names for n control points: the bottom point is always
    //an absolute temperature (temp_t0); the rest are temp_b<i> (relative factors)
    //or temp_t<i> (absolute temperatures)
    void setControlPointNames(const size_t n) {
      parameter_names.clear();

      for (size_t i=0; i<n; ++i)
      {
        if (i == 0 || parametrisation == Parametrisation::absolute)
          parameter_names.push_back("temp_t" + std::to_string(i));
        else
          parameter_names.push_back("temp_b" + std::to_string(i));
      }
    };

    bool checkProfile(std::vector<double>& temperature) {
      for (auto & i : temperature)
        if (i < 50)
        {
          i = 50;
          return true;
        }
      
      return false;
    };
};


}
#endif 

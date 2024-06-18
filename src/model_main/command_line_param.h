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


struct CommandLineParam{
  bool only_post_process = false;
  bool multinest_resume = false;
};


void parseCommandLine(int argc, char *argv[], CommandLineParam& parameters)
{
  
  //check for addtional command line parameters
  if (argc < 3) return;
  
  for (int i=2; i<argc; ++i)
  { 
    std::string input = argv[i];

    if (input == "-p" || input == "-post") {parameters.only_post_process = true; continue;}
    if (input == "-r" || input == "-restart") {parameters.multinest_resume = true; continue;}

    std::cout << "Command line parameter " << input << " unknown\n";
  }

  
}



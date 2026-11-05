// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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



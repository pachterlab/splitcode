#ifndef SPLITCODE_COMMON_H
#define SPLITCODE_COMMON_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <vector>
#include <iostream>

#ifdef _WIN64
typedef unsigned int uint;
#endif

struct ProgramOptions {
  int threads;
  int nfiles;
  bool no_output;
  std::vector<std::string> files;
  
  ProgramOptions() :
    threads(1),
    nfiles(1),
    no_output(false)
  {}
};

#endif // SPLITCODE_COMMON_H

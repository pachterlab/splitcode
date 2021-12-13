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
  std::vector<std::string> files;
  
  ProgramOptions() :
    threads(1),
    nfiles(1)
  {}
};

#endif // SPLITCODE_COMMON_H

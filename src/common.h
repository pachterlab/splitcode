#ifndef SPLITCODE_COMMON_H
#define SPLITCODE_COMMON_H

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
  std::string barcode_str;
  std::string distance_str;
  std::string location_str;
  std::string barcode_identifiers_str;
  std::string max_finds_str;
  std::string min_finds_str;
  
  ProgramOptions() :
    threads(1),
    nfiles(1),
    no_output(false)
  {}
};

#endif // SPLITCODE_COMMON_H

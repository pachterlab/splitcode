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
  bool pipe;
  bool output_fastq_specified;
  bool verbose;
  bool mod_names;
  bool gzip;
  bool discard;
  bool trim_only;
  bool discard_group;
  bool disable_n;
  std::vector<std::string> files;
  std::vector<std::string> output_files;
  std::string outputb_file;
  std::vector<std::string> unassigned_files;
  std::string barcode_str;
  std::string distance_str;
  std::string location_str;
  std::string barcode_identifiers_str;
  std::string group_identifiers_str;
  std::string max_finds_str;
  std::string min_finds_str;
  std::string max_finds_group_str;
  std::string min_finds_group_str;
  std::string exclude_str;
  std::string after_str;
  std::string config_file;
  std::string mapping_file;
  std::string keep_file;
  std::string keep_group_file;
  std::string append_file;
  std::string left_str;
  std::string right_str;
  std::string empty_read_sequence;
  
  ProgramOptions() :
    threads(1),
    nfiles(1),
    no_output(false),
    pipe(false),
    output_fastq_specified(false),
    verbose(false),
    mod_names(false),
    gzip(false),
    discard(false),
    trim_only(false),
    discard_group(false),
    disable_n(false)
  {}
};

#endif // SPLITCODE_COMMON_H

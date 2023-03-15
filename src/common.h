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
  int input_interleaved_nfiles;
  int quality_trimming_threshold;
  int64_t max_num_reads;
  int compress_level;
  size_t hashmap_limit;
  bool extract_no_chain;
  bool output_fasta;
  bool no_output;
  bool no_output_; // When --no-output is specified but we want to write other things (unassigned, --keep, etc.), just not --output
  bool no_output_barcodes;
  bool no_x_out;
  bool empty_remove;
  bool pipe;
  bool output_fastq_specified;
  bool verbose;
  bool mod_names;
  bool com_names;
  bool bc_names;
  bool seq_names;
  bool x_names;
  bool x_only;
  bool gzip;
  bool discard;
  bool trim_only;
  bool discard_group;
  bool disable_n;
  bool quality_trimming_5;
  bool quality_trimming_3;
  bool quality_trimming_pre;
  bool quality_trimming_naive;
  bool phred64;
  bool keep_fastq_comments;
  bool unlimited_hashmap;
  bool remultiplex;
  bool webasm;
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
  std::string before_str;
  std::string partial5_str;
  std::string partial3_str;
  std::string config_file;
  std::string mapping_file;
  std::string keep_file;
  std::string keep_group_file;
  std::string append_file;
  std::string left_str;
  std::string right_str;
  std::string empty_read_sequence;
  std::string trim_5_str;
  std::string trim_3_str;
  std::string filter_length_str;
  std::string extract_str;
  std::string barcode_prefix;
  std::string summary_file;
  std::string subs_str;
  std::string select_output_files_str;
  std::vector<bool> select_output_files;
  std::vector<std::vector<std::string> > sam_tags;
  std::vector<size_t> sub_assign_vec;
  std::vector<std::string> batch_ids;
  
  ProgramOptions() :
    threads(1),
    nfiles(1),
    input_interleaved_nfiles(0),
    quality_trimming_threshold(-1),
    max_num_reads(0),
    compress_level(1),
    hashmap_limit(320000),
    extract_no_chain(false),
    output_fasta(false),
    no_output(false),
    no_output_(false),
    no_output_barcodes(false),
    no_x_out(false),
    empty_remove(false),
    pipe(false),
    output_fastq_specified(false),
    verbose(false),
    mod_names(false),
    com_names(false),
    bc_names(false),
    seq_names(false),
    x_names(false),
    x_only(false),
    gzip(false),
    discard(false),
    trim_only(true),
    discard_group(false),
    disable_n(true),
    quality_trimming_5(false),
    quality_trimming_3(false),
    quality_trimming_pre(false),
    quality_trimming_naive(false),
    phred64(false),
    keep_fastq_comments(false),
    unlimited_hashmap(false),
    remultiplex(false),
    webasm(false)
  {
    const char* sam_tags_default[5] = {"CB:Z:", "RX:Z:", "BI:i:", "SI:i:", "BC:Z:"};
    sam_tags.push_back(std::vector<std::string>(1, std::string(sam_tags_default[0])));
    sam_tags.push_back(std::vector<std::string>(1, std::string(sam_tags_default[1])));
    sam_tags.push_back(std::vector<std::string>(1, std::string(sam_tags_default[2])));
    sam_tags.push_back(std::vector<std::string>(1, std::string(sam_tags_default[3])));
    sam_tags.push_back(std::vector<std::string>(1, std::string(sam_tags_default[4])));
  }
};

#endif // SPLITCODE_COMMON_H

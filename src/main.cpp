#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <algorithm>
#include <limits>

#include <cstdio>

#include <zlib.h>

#include "common.h"
#include "ProcessReads.h"
#include "SplitCode.h"


//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;


int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}


void PrintCite() {
  cout << "[no citation info to display yet]" 
       << endl;
}

void PrintVersion() {
  cout << "splitcode, version " << 	SPLITCODE_VERSION << endl;
}

void usage() {
  cout << "splitcode " << SPLITCODE_VERSION << endl << endl
       << "Usage: splitcode [arguments] fastq-files" << endl << endl
       << "Options (for configuring on the command-line):" << endl
       << "-b, --barcodes   List of barcode sequences (comma-separated)" << endl
       << "-d, --distances  List of error distance (mismatch:indel:total) thresholds (comma-separated)" << endl
       << "-l, --locations  List of locations (file:pos1:pos2) (comma-separated)" << endl
       << "-i, --ids        List of barcode names/identifiers (comma-separated)" << endl
       << "-g, --groups     List of barcode group names (comma-separated)" << endl
       << "-f, --minFinds   List of minimum times a barcode must be found in a read (comma-separated)" << endl
       << "-F, --maxFinds   List of maximum times a barcode can be found in a read (comma-separated)" << endl
       << "-j, --minFindsG  List of minimum times barcodes in a group must be found in a read (comma-separated group_name:min_times)" << endl
       << "-J, --maxFindsG  List of maximum times barcodes in a group can be found in a read (comma-separated group_name:max_times)" << endl
       << "-e, --exclude    List of what to exclude from final barcode (comma-separated; 1 = exclude, 0 = include)" << endl
       << "-L, --left       List of what barcodes to include when trimming from the left (comma-separated; 1 = include, 0 = exclude)" << endl
       << "-R, --right      List of what barcodes to include when trimming from the right (comma-separated; 1 = include, 0 = exclude)" << endl
       << "                 (Note: for --left/--right, can specify an included barcode as 1:x where x = number of extra bp's to trim" << endl
       << "                 from left/right side if the that included barcode is at the leftmost/rightmost position)" << endl
       << "-a, --next       List of what barcode names must come immediately after each barcode (comma-separated)" << endl
       << "-v, --previous   List of what barcode names must come immediately before each barcode (comma-separated)" << endl
       << "                 (Note: for --next/--previous, specify barcode names as {name} and specify barcode group names as {{group}}" << endl
       << "                 Can also specify the number of base pairs that must appear between the current barcode and the next/previous barcode." << endl
       << "                 E.g. {bc}4-12 means the next/previous barcode is 4-12 bases away and has name 'bc')" << endl
       << "-5, --trim-5     Number of base pairs to trim from the 5′-end of reads (comma-separated; one number per each FASTQ file in a run)" << endl
       << "-3, --trim-3     Number of base pairs to trim from the 3′-end of reads (comma-separated; one number per each FASTQ file in a run)" << endl
       << "Options (configurations supplied in a file):" << endl
       << "-c, --config     Configuration file" << endl
       << "Output Options:" << endl
       << "-m, --mapping    Output file where the mapping between final barcode sequences and names will be written" << endl
       << "-o, --output     FASTQ file(s) where output will be written (comma-separated)" << endl
       << "                 Number of output FASTQ files should equal --nFastqs" << endl
       << "-O, --outb       FASTQ file where final barcodes will be written" << endl
       << "                 If not supplied, final barcodes are prepended to reads of first FASTQ file (or as the first read for --pipe)" << endl
       << "-u, --unassigned FASTQ file(s) where output of unassigned reads will be written (comma-separated)" << endl
       << "                 Number of FASTQ files should equal --nFastqs" << endl
       << "-E, --empty      Sequence to fill in empty reads in output FASTQ files (default: no sequence is used to fill in those reads)" << endl
       << "-p, --pipe       Write to standard output (instead of output FASTQ files)" << endl
       << "    --gzip       Output compressed gzip'ed FASTQ files" << endl
       << "    --no-output  Don't output any sequences (output statistics only)" << endl
       << "    --no-outb    Don't output barcode sequences" << endl
       << "    --mod-names  Modify names of outputted sequences to include identified barcodes" << endl
       << "    --com-names  Modify names of outputted sequences to include final barcode sequence ID" << endl
       << "Other Options:" << endl
       << "-N, --nFastqs    Number of FASTQ file(s) per run" << endl
       << "                 (default: 1) (specify 2 for paired-end)" << endl
       << "-n, --numReads   Maximum number of reads to process from supplied input" << endl
       << "-A, --append     An existing mapping file that will be added on to" << endl
       << "-k, --keep       File containing a list of final barcodes to keep" << endl
       << "-r, --remove     File containing a list of final barcodes to remove/discard" << endl
       << "-y, --keep-grp   File containing a list of final barcode groups to keep" << endl
       << "-Y, --remove-grp File containing a list of final barcode groups to remove/discard" << endl
       << "-t, --threads    Number of threads to use" << endl
       << "-T, --trim-only  All reads are assigned and trimmed regardless of barcode identification" << endl
       << "-h, --help       Displays usage information" << endl
       << "    --inleaved   Specifies that input is an interleaved FASTQ file" << endl
       << "    --disable-n  Disables replacing ambiguous bases with pseudorandom bases" << endl
       << "    --version    Prints version number" << endl
       << "    --cite       Prints citation information" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions& opt) {
  int help_flag = 0;
  int version_flag = 0;
  int cite_flag = 0;
  int no_output_flag = 0;
  int no_output_barcodes_flag = 0;
  int gzip_flag = 0;
  int mod_names_flag = 0;
  int com_names_flag = 0;
  int disable_n_flag = 0;
  int interleaved_flag = 0;

  const char *opt_string = "t:N:n:b:d:i:l:f:F:e:c:o:O:u:m:k:r:A:L:R:E:g:y:Y:j:J:a:v:5:3:Tph";
  static struct option long_options[] = {
    // long args
    {"version", no_argument, &version_flag, 1},
    {"cite", no_argument, &cite_flag, 1},
    {"no-output", no_argument, &no_output_flag, 1},
    {"no-outb", no_argument, &no_output_barcodes_flag, 1},
    {"gzip", no_argument, &gzip_flag, 1},
    {"mod-names", no_argument, &mod_names_flag, 1},
    {"com-names", no_argument, &com_names_flag, 1},
    {"disable-n", no_argument, &disable_n_flag, 1},
    {"inleaved", no_argument, &interleaved_flag, 1},
    // short args
    {"help", no_argument, 0, 'h'},
    {"pipe", no_argument, 0, 'p'},
    {"trim-only", no_argument, 0, 'T'},
    {"threads", required_argument, 0, 't'},
    {"nFastqs", required_argument, 0, 'N'},
    {"numReads", required_argument, 0, 'n'},
    {"barcodes", required_argument, 0, 'b'},
    {"distances", required_argument, 0, 'd'},
    {"locations", required_argument, 0, 'l'},
    {"ids", required_argument, 0, 'i'},
    {"groups", required_argument, 0, 'g'},
    {"maxFinds", required_argument, 0, 'F'},
    {"minFinds", required_argument, 0, 'f'},
    {"maxFindsG", required_argument, 0, 'J'},
    {"minFindsG", required_argument, 0, 'j'},
    {"exclude", required_argument, 0, 'e'},
    {"next", required_argument, 0, 'a'},
    {"after", required_argument, 0, 'a'},
    {"previous", required_argument, 0, 'v'},
    {"before", required_argument, 0, 'v'},
    {"config", required_argument, 0, 'c'},
    {"output", required_argument, 0, 'o'},
    {"outb", required_argument, 0, 'O'},
    {"unassigned", required_argument, 0, 'u'},
    {"mapping", required_argument, 0, 'm'},
    {"keep", required_argument, 0, 'k'},
    {"remove", required_argument, 0, 'r'},
    {"append", required_argument, 0, 'A'},
    {"left", required_argument, 0, 'L'},
    {"right", required_argument, 0, 'R'},
    {"empty", required_argument, 0, 'E'},
    {"keep-grp", required_argument, 0, 'y'},
    {"remove-grp", required_argument, 0, 'Y'},
    {"trim-5", required_argument, 0, '5'},
    {"trim-3", required_argument, 0, '3'},
    {0,0,0,0}
  };
  
  int c;
  int option_index = 0;
  int num_opts_supplied = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }
    
    num_opts_supplied++;
    
    switch (c) {
    case 0:
      break;
    case 'h': {
      help_flag = 1;
      break;
    }
    case 'p': {
      opt.pipe = true;
      break;
    }
    case 'T': {
      opt.trim_only = true;
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'N': {
      stringstream(optarg) >> opt.nfiles;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.max_num_reads;
      if (opt.max_num_reads == 0) {
        opt.max_num_reads = -1;
      }
      break;
    }
    case 'b': {
      stringstream(optarg) >> opt.barcode_str;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.distance_str;
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.location_str;
      break;
    }
    case 'i': {
      stringstream(optarg) >> opt.barcode_identifiers_str;
      break;
    }
    case 'g': {
      stringstream(optarg) >> opt.group_identifiers_str;
      break;
    }
    case 'F': {
      stringstream(optarg) >> opt.max_finds_str;
      break;
    }
    case 'f': {
      stringstream(optarg) >> opt.min_finds_str;
      break;
    }
    case 'J': {
      stringstream(optarg) >> opt.max_finds_group_str;
      break;
    }
    case 'j': {
      stringstream(optarg) >> opt.min_finds_group_str;
      break;
    }
    case 'e': {
      stringstream(optarg) >> opt.exclude_str;
      break;
    }
    case 'L': {
      stringstream(optarg) >> opt.left_str;
      break;
    }
    case 'R': {
      stringstream(optarg) >> opt.right_str;
      break;
    }
    case 'a': {
      stringstream(optarg) >> opt.after_str;
      break;
    }
    case 'v': {
      stringstream(optarg) >> opt.before_str;
      break;
    }
    case 'c': {
      stringstream(optarg) >> opt.config_file;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.mapping_file;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.keep_file;
      break;
    }
    case 'r': {
      opt.discard = true;
      stringstream(optarg) >> opt.keep_file;
      break;
    }
    case 'y': {
      stringstream(optarg) >> opt.keep_group_file;
      break;
    }
    case 'Y': {
      opt.discard_group = true;
      stringstream(optarg) >> opt.keep_group_file;
      break;
    }
    case 'o': {
      std::string files;
      stringstream(optarg) >> files;
      std::stringstream ss(files);
      std::string filename;
      while (std::getline(ss, filename, ',')) { 
        opt.output_files.push_back(filename);
      }
      break;
    }
    case 'O': {
      stringstream(optarg) >> opt.outputb_file;
      break;
    }
    case 'A': {
      stringstream(optarg) >> opt.append_file;
      break;
    }
    case 'E': {
      stringstream(optarg) >> opt.empty_read_sequence;
      for (auto& c: opt.empty_read_sequence) {
        c = toupper(c);
      }
      break;
    }
    case 'u': {
      std::string files;
      stringstream(optarg) >> files;
      std::stringstream ss(files);
      std::string filename;
      while (std::getline(ss, filename, ',')) { 
        opt.unassigned_files.push_back(filename);
      }
      break;
    }
    case '5': {
      stringstream(optarg) >> opt.trim_5_str;
      break;
    }
    case '3': {
      stringstream(optarg) >> opt.trim_3_str;
      break;
    }
    default: break;
    }
  }
  
  if (help_flag || num_opts_supplied == 0) {
    usage();
    exit(0);
  }
  if (version_flag) {
    PrintVersion();
    exit(0);
  }
  if (cite_flag) {
    PrintCite();
    exit(0);
  }
  if (no_output_flag) {
    opt.no_output = true;
  }
  if (mod_names_flag) {
    opt.mod_names = true;
  }
  if (com_names_flag) {
    opt.com_names = true;
  }
  if (gzip_flag) {
    opt.gzip = true;
  }
  if (disable_n_flag) {
    opt.disable_n = true;
  }
  if (interleaved_flag) {
    opt.input_interleaved_nfiles = 1;
  }
  if (no_output_barcodes_flag) {
    opt.no_output_barcodes = true;
  }
  
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
}

bool CheckOptions(ProgramOptions& opt, SplitCode& sc) {
  bool ret = true;
  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }    
  }
  if (opt.files.size() == 0) {
    cerr << ERROR_STR << " Missing read files" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    for (auto& fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << ERROR_STR << " file not found " << fn << endl;
        ret = false;
      }
    }
  }
  if (opt.nfiles <= 0) {
    std::cerr << ERROR_STR << " nFastqs must be a non-zero positive number" << std::endl;
    ret = false;
  } else if (opt.input_interleaved_nfiles != 0) {
    if (opt.files.size() != 1) {
      std::cerr << ERROR_STR << " interleaved input cannot consist of more than one input" << std::endl;
      ret = false;
    }
  } else {
    if (opt.files.size() % opt.nfiles != 0) {
      std::cerr << ERROR_STR << " incorrect number of FASTQ file(s)" << std::endl;
      ret = false;
    }
  }
  if (opt.max_num_reads < 0) {
    std::cerr << ERROR_STR << " --numReads must be a positive number" << std::endl;
    ret = false;
  }
  if (opt.mapping_file.empty() && !opt.trim_only) {
    std::cerr << ERROR_STR << " --mapping must be provided" << std::endl;
    ret = false;
  }
  if (opt.no_output_barcodes && !opt.outputb_file.empty()) {
    std::cerr << ERROR_STR << " --no-outb cannot be specified with --outb" << std::endl;
    ret = false;
  }
  
  bool output_files_specified = opt.output_files.size() > 0 || opt.unassigned_files.size() > 0 || !opt.outputb_file.empty();
  if (opt.output_files.size() == 0 && output_files_specified && !opt.pipe) {
    std::cerr << ERROR_STR << " --output not provided" << std::endl;
    ret = false;
  }
  if (opt.no_output) {
    if (output_files_specified || opt.pipe) {
      std::cerr << ERROR_STR << " Cannot specify an output option when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.mod_names || opt.com_names) {
      std::cerr << ERROR_STR << " Cannot use --mod-names/--com-names when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.gzip) {
      std::cerr << ERROR_STR << " Cannot use --gzip when --no-output is specified" << std::endl;
      ret = false;
    }
  } else {
    if (!output_files_specified && !opt.pipe) {
      std::cerr << ERROR_STR << " Must either specify an output option or --no-output" << std::endl;
      ret = false;
    } else if (opt.pipe) {
      if (opt.output_files.size() > 0 || !opt.outputb_file.empty()) { // Still allow --unassigned with --pipe
        std::cerr << ERROR_STR << " Cannot provide output files when --pipe is specified" << std::endl;
        ret = false;
      } else if (opt.unassigned_files.size() != 0 && opt.unassigned_files.size() % opt.nfiles != 0) {
        std::cerr << ERROR_STR << " Incorrect number of --unassigned output files" << std::endl;
        ret = false;
      }
    } else {
      if (opt.output_files.size() % opt.nfiles != 0 || opt.unassigned_files.size() % opt.nfiles != 0) {
        std::cerr << ERROR_STR << " Incorrect number of output files" << std::endl;
        ret = false;
      }
    }
  }
  if (opt.trim_only && opt.no_output) {
    std::cerr << ERROR_STR << " Cannot use --trim-only with --no-output" << std::endl;
    ret = false;
  }
  if (opt.trim_only && opt.unassigned_files.size() != 0) {
    std::cerr << ERROR_STR << " Cannot use --trim-only with --unassigned" << std::endl;
    ret = false;
  }
  if (opt.trim_only && !opt.outputb_file.empty()) {
    std::cerr << ERROR_STR << " Cannot use --trim-only with --outb" << std::endl;
    ret = false;
  }
  if (opt.trim_only && !opt.mapping_file.empty()) {
    std::cerr << ERROR_STR << " Cannot use --trim-only with --mapping" << std::endl;
    ret = false;
  }
  opt.output_fastq_specified = output_files_specified;
  opt.verbose = !opt.pipe;
  
  int num_groups = 0;
  if (!opt.barcode_str.empty() && !opt.config_file.empty()) {
    std::cerr << ERROR_STR << " Cannot specify both --barcodes and --config" << std::endl;
    ret = false;
  } else if (!opt.barcode_str.empty()) {
    stringstream ss1(opt.barcode_str);
    stringstream ss2(opt.distance_str);
    stringstream ss3(opt.barcode_identifiers_str);
    stringstream ss4(opt.location_str);
    stringstream ss5(opt.max_finds_str);
    stringstream ss6(opt.min_finds_str);
    stringstream ss7(opt.exclude_str);
    stringstream ss8(opt.left_str);
    stringstream ss9(opt.right_str);
    stringstream ss10(opt.group_identifiers_str);
    stringstream ss11(opt.after_str);
    stringstream ss12(opt.before_str);
    while (ss1.good()) {
      uint16_t max_finds = 0;
      uint16_t min_finds = 0;
      bool exclude = false;
      string name = "";
      string group = "";
      string location = "";
      string distance = "";
      string left_str = "";
      string right_str = "";
      string after_str = "";
      string before_str = "";
      bool trim_left, trim_right;
      int trim_left_offset, trim_right_offset;
      int16_t file;
      int32_t pos_start;
      int32_t pos_end;
      int mismatch, indel, total_dist;
      string bc;
      getline(ss1, bc, ',');
      if (!opt.distance_str.empty()) {
        auto currpos = ss2.tellg();
        if (!ss2.good()) {
          std::cerr << ERROR_STR << " Number of values in --distances is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss2, distance, ',');
        if (!ss2.good() && ss1.good() && currpos == 0) {
          ss2.clear();
          ss2.str(opt.distance_str);
        }
      }
      if (!SplitCode::parseDistance(distance, mismatch, indel, total_dist)) {
        std::cerr << ERROR_STR << " --distances is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.barcode_identifiers_str.empty()) {
        if (!ss3.good()) {
          std::cerr << ERROR_STR << " Number of values in --ids is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss3, name, ',');
      }
      if (!opt.location_str.empty()) {
        auto currpos = ss4.tellg();
        if (!ss4.good()) {
          std::cerr << ERROR_STR << " Number of values in --locations is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss4, location, ',');
        if (!ss4.good() && ss1.good() && currpos == 0) {
          ss4.clear();
          ss4.str(opt.location_str);
        }
      }
      if (!SplitCode::parseLocation(location, file, pos_start, pos_end, opt.nfiles)) {
        std::cerr << ERROR_STR << " --locations is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.max_finds_str.empty()) {
        auto currpos = ss5.tellg();
        if (!ss5.good()) {
          std::cerr << ERROR_STR << " Number of values in --maxFinds is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss5, f, ',');
        stringstream(f) >> max_finds;
        if (!ss5.good() && ss1.good() && currpos == 0) {
          ss5.clear();
          ss5.str(opt.max_finds_str);
        }
      }
      if (!opt.min_finds_str.empty()) {
        auto currpos = ss6.tellg();
        if (!ss6.good()) {
          std::cerr << ERROR_STR << " Number of values in --minFinds is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss6, f, ',');
        stringstream(f) >> min_finds;
        if (!ss6.good() && ss1.good() && currpos == 0) {
          ss6.clear();
          ss6.str(opt.min_finds_str);
        }
      }
      if (!opt.exclude_str.empty()) {
        if (!ss7.good()) {
          std::cerr << ERROR_STR << " Number of values in --exclude is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss7, f, ',');
        stringstream(f) >> exclude;
      }
      if (!opt.left_str.empty()) {
        auto currpos = ss8.tellg();
        if (!ss8.good()) {
          std::cerr << ERROR_STR << " Number of values in --left is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss8, f, ',');
        stringstream(f) >> left_str;
        if (!ss8.good() && ss1.good() && currpos == 0) {
          ss8.clear();
          ss8.str(opt.left_str);
        }
      }
      if (!SplitCode::parseTrimStr(left_str, trim_left, trim_left_offset)) {
        std::cerr << ERROR_STR << " --left is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.right_str.empty()) {
        auto currpos = ss9.tellg();
        if (!ss9.good()) {
          std::cerr << ERROR_STR << " Number of values in --right is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss9, f, ',');
        stringstream(f) >> right_str;
        if (!ss9.good() && ss1.good() && currpos == 0) {
          ss9.clear();
          ss9.str(opt.right_str);
        }
      }
      if (!SplitCode::parseTrimStr(right_str, trim_right, trim_right_offset)) {
        std::cerr << ERROR_STR << " --right is invalid" << std::endl;
        ret = false;
        break;
      }
      if (trim_left && trim_right) {
        std::cerr << ERROR_STR << " One of the barcodes has both --left and --right trimming specified" << std::endl;
        ret = false;
        break;
      }
      auto trim_dir = trim_left ? sc.left : (trim_right ? sc.right : sc.nodir);
      auto trim_offset = trim_left ? trim_left_offset : (trim_right ? trim_right_offset : 0);
      if (!opt.group_identifiers_str.empty()) {
        if (!ss10.good()) {
          std::cerr << ERROR_STR << " Number of values in --groups is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss10, group, ',');
        num_groups++;
      }
      if (!opt.after_str.empty()) {
        if (!ss11.good()) {
          std::cerr << ERROR_STR << " Number of values in --next is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss11, after_str, ',');
      }
      if (!SplitCode::validateBeforeAfterStr(after_str)) {
        std::cerr << ERROR_STR << " --next is invalid" << std::endl;
        ret = false;
      }
      if (!opt.before_str.empty()) {
        if (!ss12.good()) {
          std::cerr << ERROR_STR << " Number of values in --previous is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss12, before_str, ',');
      }
      if (!SplitCode::validateBeforeAfterStr(before_str)) {
        std::cerr << ERROR_STR << " --previous is invalid" << std::endl;
        ret = false;
      }
      if (!sc.addTag(bc, name.empty() ? bc : name, group, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude, trim_dir, trim_offset, after_str, before_str)) {
        std::cerr << ERROR_STR << " Could not finish processing supplied barcode list" << std::endl;
        ret = false;
        break;
      }
    }
    if (ret && !opt.distance_str.empty() && ss2.good()) {
      std::cerr << ERROR_STR << " Number of values in --distances is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.barcode_identifiers_str.empty() && ss3.good()) {
      std::cerr << ERROR_STR << " Number of values in --ids is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.location_str.empty() && ss4.good()) {
      std::cerr << ERROR_STR << " Number of values in --locations is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.max_finds_str.empty() && ss5.good()) {
      std::cerr << ERROR_STR << " Number of values in --maxFinds is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.min_finds_str.empty() && ss6.good()) {
      std::cerr << ERROR_STR << " Number of values in --minFinds is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.exclude_str.empty() && ss7.good()) {
      std::cerr << ERROR_STR << " Number of values in --exclude is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.left_str.empty() && ss8.good()) {
      std::cerr << ERROR_STR << " Number of values in --left is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.right_str.empty() && ss9.good()) {
      std::cerr << ERROR_STR << " Number of values in --right is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.group_identifiers_str.empty() && ss10.good()) {
      std::cerr << ERROR_STR << " Number of values in --groups is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.after_str.empty() && ss11.good()) {
      std::cerr << ERROR_STR << " Number of values in --next is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (ret && !opt.before_str.empty() && ss12.good()) {
      std::cerr << ERROR_STR << " Number of values in --previous is greater than that in --barcodes" << std::endl;
      ret = false;
    }
  } else if (!opt.distance_str.empty()) {
    std::cerr << ERROR_STR << " --distances cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.barcode_identifiers_str.empty()) {
    std::cerr << ERROR_STR << " --ids cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.location_str.empty()) {
    std::cerr << ERROR_STR << " --locations cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.max_finds_str.empty()) {
    std::cerr << ERROR_STR << " --maxFinds cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.min_finds_str.empty()) {
    std::cerr << ERROR_STR << " --minFinds cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.exclude_str.empty()) {
    std::cerr << ERROR_STR << " --exclude cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.left_str.empty()) {
    std::cerr << ERROR_STR << " --left cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.right_str.empty()) {
    std::cerr << ERROR_STR << " --right cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.group_identifiers_str.empty()) {
    std::cerr << ERROR_STR << " --groups cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.after_str.empty()) {
    std::cerr << ERROR_STR << " --next cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.before_str.empty()) {
    std::cerr << ERROR_STR << " --previous cannot be supplied unless --barcodes is" << std::endl;
    ret = false;
  } else if (!opt.config_file.empty()) {
    ret = ret && sc.addTags(opt.config_file);
  }
  
  if (ret) { // Now process groups arguments (e.g. maxFindsG/minFindsG) rather than individual barcodes arguments
    if (num_groups > 0) {
      try {
        stringstream ss1(opt.max_finds_group_str);
        stringstream ss2(opt.min_finds_group_str);
        string f;
        while (getline(ss1, f, ',')) {
          int i = 0;
          stringstream ss_(f);
          string g_attribute;
          string group_name = "";
          int group_param = 0;
          while (std::getline(ss_, g_attribute, ':')) {
            if (i == 0) {
              group_name = g_attribute;
            } else if (i == 1) {
              group_param = std::stoi(g_attribute);
            }
            i++;
          }
          if (i != 2 || group_name.empty() || group_param < 0) {
            std::cerr << "Error: --maxFindsG string is malformed; unable to parse \"" << f << "\"" << std::endl;
            ret = false;
            break;
          } else if (!sc.addGroupOptions(group_name, group_param, 0)) {
            std::cerr << ERROR_STR << " Could not finish processing supplied groups options in --maxFindsG" << std::endl;
            ret = false;
            break;
          }
        }
        while (getline(ss2, f, ',')) {
          int i = 0;
          stringstream ss_(f);
          string g_attribute;
          string group_name = "";
          int group_param = 0;
          while (std::getline(ss_, g_attribute, ':')) {
            if (i == 0) {
              group_name = g_attribute;
            } else if (i == 1) {
              group_param = std::stoi(g_attribute);
            }
            i++;
          }
          if (i != 2 || group_name.empty() || group_param < 0) {
            std::cerr << "Error: --minFindsG string is malformed; unable to parse \"" << f << "\"" << std::endl;
            ret = false;
            break;
          } else if (!sc.addGroupOptions(group_name, 0, group_param)) {
            std::cerr << ERROR_STR << " Could not finish processing supplied groups options in --minFindsG" << std::endl;
            ret = false;
            break;
          }
        }
      } catch (std::invalid_argument &e) {
        std::cerr << "Error: Invalid number found in --minFindsG or --maxFindsG" << std::endl;
        ret = false;
      }
    } else if (!opt.max_finds_group_str.empty()) {
      std::cerr << ERROR_STR << " --maxFindsG cannot be supplied unless --groups is" << std::endl;
      ret = false;
    } else if (!opt.min_finds_group_str.empty()) {
      std::cerr << ERROR_STR << " --minFindsG cannot be supplied unless --groups is" << std::endl;
      ret = false;
    }
  }
  
  if (ret && !opt.append_file.empty()) {
    ret = ret && sc.addExistingMapping(opt.append_file);
  }

  if (ret && !opt.keep_file.empty()) {
    ret = ret && sc.addFilterList(opt.keep_file, opt.discard);
  }
  
  if (ret && !opt.keep_group_file.empty()) {
    ret = ret && sc.addFilterListGroup(opt.keep_group_file, opt.discard_group);
  }
  
  if (ret && (sc.getNumTags() == 0 || sc.getMapSize() == 0)) {
    std::cerr << ERROR_STR << " No barcodes found" << std::endl;
    ret = false;
  }
  
  return ret;
}


int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  SplitCode sc(opt.nfiles, opt.trim_only, opt.disable_n, opt.trim_5_str, opt.trim_3_str);
  if (!CheckOptions(opt, sc)) {
    usage();
    exit(1);
  }
  if (opt.input_interleaved_nfiles != 0) {
    opt.input_interleaved_nfiles = opt.nfiles;
    opt.nfiles = 1;
  }
  if (!opt.output_files.empty() && !opt.gzip) {
    bool use_gz = true;
    for (std::string f : opt.output_files) {
      if (!(f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)) {
        use_gz = false;
      }
    }
    if (!opt.outputb_file.empty()) {
      std::string f = opt.outputb_file;
      if (!(f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)) {
        use_gz = false;
      }
    }
    if (!opt.unassigned_files.size()) {
      for (std::string f : opt.unassigned_files) {
        if (!(f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)) {
          use_gz = false;
        }
      }
    }
    if (use_gz) {
      std::cerr << "* Forcing --gzip because all output file names end in .gz" << std::endl;
      opt.gzip = true;
    }
  }
  
  if (opt.verbose) {
    std::cerr << "* Using a list of " << sc.getNumTagsOriginallyAdded() << 
      " barcodes (vector size: " << sc.getNumTags() << 
      "; map size: " << pretty_num(sc.getMapSize()) << 
      "; num elements in map: " << pretty_num(sc.getMapSize(false)) << ")" << std::endl;
  }
  MasterProcessor MP(sc, opt);
  int numreads = ProcessReads(MP, opt);
  fflush(stdout);
  if (!opt.mapping_file.empty()) { // output mapping file:
    if (!(opt.mapping_file.size() > 3 && opt.mapping_file.compare(opt.mapping_file.size() - 3, 3, ".gz") == 0)) {
      sc.writeBarcodeMapping(opt.mapping_file); // output plaintext mapping file
    } else {
      gzFile out_gz = gzopen(opt.mapping_file.c_str(), "wb6");
      std::string o;
      while ((o = sc.fetchNextBarcodeMapping()) != "") {
        gzwrite(out_gz, o.c_str(), o.length()); // output gzip mapping file
      }
      gzclose(out_gz);
    }
  }
  
  if (opt.max_num_reads != 0 && numreads < opt.max_num_reads) {
    std::cerr << "Note: Number of reads processed is less than --numReads: " << opt.max_num_reads << ", returning 1" << std::endl;
    return 1;
  }

  return 0;
}

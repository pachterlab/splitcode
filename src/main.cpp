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
       << "-f, --minFinds   List of minimum times a barcode must be found in a read (comma-separated)" << endl
       << "-F, --maxFinds   List of maximum times a barcode can be found in a read (comma-separated)" << endl
       << "-e, --exclude    List of what to exclude from final barcode (comma-separated; 1 = exclude, 0 = include)" << endl
       << "Options (configurations supplied in a file):" << endl
       << "-c, --config     Configuration file" << endl
       << "Other Options:" << endl
       << "-N, --nFastqs    Number of FASTQ file(s) per run" << endl
       << "                 (default: 1) (specify 2 for paired-end)" << endl
       << "-t, --threads    Number of threads to use" << endl
       << "-h, --help       Displays usage information" << endl
       << "    --version    Prints version number" << endl
       << "    --cite       Prints citation information" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions& opt) {
  int help_flag = 0;
  int version_flag = 0;
  int cite_flag = 0;

  const char *opt_string = "t:N:b:d:i:l:f:F:e:c:h";
  static struct option long_options[] = {
    // long args
    {"version", no_argument, &version_flag, 1},
    {"cite", no_argument, &cite_flag, 1},
    // short args
    {"help", no_argument, 0, 'h'},
    {"threads", required_argument, 0, 't'},
    {"nFastqs", required_argument, 0, 'N'},
    {"barcodes", required_argument, 0, 'b'},
    {"distances", required_argument, 0, 'd'},
    {"locations", required_argument, 0, 'l'},
    {"ids", required_argument, 0, 'i'},
    {"maxFinds", required_argument, 0, 'F'},
    {"minFinds", required_argument, 0, 'f'},
    {"exclude", required_argument, 0, 'e'},
    {"config", required_argument, 0, 'c'},
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
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'N': {
      stringstream(optarg) >> opt.nfiles;
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
    case 'F': {
      stringstream(optarg) >> opt.max_finds_str;
      break;
    }
    case 'f': {
      stringstream(optarg) >> opt.min_finds_str;
      break;
    }
    case 'e': {
      stringstream(optarg) >> opt.exclude_str;
      break;
    }
    case 'c': {
      stringstream(optarg) >> opt.config_file;
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
  }
  else {
    if (opt.files.size() % opt.nfiles != 0) {
      std::cerr << ERROR_STR << " incorrect number of FASTQ file(s)" << std::endl;
      ret = false;
    }
  }
  
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
    while (ss1.good()) {
      uint16_t max_finds = 0;
      uint16_t min_finds = 0;
      bool exclude = false;
      string name = "";
      string location = "";
      string distance = "";
      int16_t file;
      int32_t pos_start;
      int32_t pos_end;
      int mismatch, indel, total_dist;
      if (!opt.distance_str.empty()) {
        if (!ss2.good()) {
          std::cerr << ERROR_STR << " Number of values in --distances is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss2, distance, ',');
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
        if (!ss4.good()) {
          std::cerr << ERROR_STR << " Number of values in --locations is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        getline(ss4, location, ',');
      }
      if (!SplitCode::parseLocation(location, file, pos_start, pos_end, opt.nfiles)) {
        std::cerr << ERROR_STR << " --locations is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.max_finds_str.empty()) {
        if (!ss5.good()) {
          std::cerr << ERROR_STR << " Number of values in --maxFinds is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss5, f, ',');
        stringstream(f) >> max_finds;
      }
      if (!opt.min_finds_str.empty()) {
        if (!ss6.good()) {
          std::cerr << ERROR_STR << " Number of values in --minFinds is less than that in --barcodes" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss6, f, ',');
        stringstream(f) >> min_finds;
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
      string bc;
      getline(ss1, bc, ',');
      if (!sc.addTag(bc, name.empty() ? bc : name, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude)) {
        std::cerr << ERROR_STR << " Could not finish processing supplied barcode list" << std::endl;
        ret = false;
      }
    }
    if (!opt.distance_str.empty() && ss2.good()) {
      std::cerr << ERROR_STR << " Number of values in --distances is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (!opt.barcode_identifiers_str.empty() && ss3.good()) {
      std::cerr << ERROR_STR << " Number of values in --ids is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (!opt.location_str.empty() && ss4.good()) {
      std::cerr << ERROR_STR << " Number of values in --locations is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (!opt.max_finds_str.empty() && ss5.good()) {
      std::cerr << ERROR_STR << " Number of values in --maxFinds is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (!opt.min_finds_str.empty() && ss6.good()) {
      std::cerr << ERROR_STR << " Number of values in --minFinds is greater than that in --barcodes" << std::endl;
      ret = false;
    }
    if (!opt.exclude_str.empty() && ss7.good()) {
      std::cerr << ERROR_STR << " Number of values in --exclude is greater than that in --barcodes" << std::endl;
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
  } else if (!opt.config_file.empty()) {
    ret = ret && sc.addTags(opt.config_file);
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
  SplitCode sc;
  if (!CheckOptions(opt, sc)) {
    usage();
    exit(1);
  }
  std::cerr << "* Using a list of " << sc.getNumTags() << " barcodes (map size: " << pretty_num(sc.getMapSize()) << "; num elements: " << pretty_num(sc.getMapSize(false)) << ")" << std::endl;
  MasterProcessor MP(sc, opt);
  ProcessReads(MP, opt);
  fflush(stdout);

  return 0;
}

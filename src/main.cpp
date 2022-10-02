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

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }
  return res;
}

void PrintVersion() {
  cout << "splitcode, version " << 	SPLITCODE_VERSION << endl;
}

void usage() {
  cout << "splitcode " << SPLITCODE_VERSION << endl << endl
       << "Usage: splitcode [arguments] fastq-files" << endl << endl
       << "Sequence identification options (for configuring on the command-line):" << endl
       << "-b, --tags       List of tag sequences (comma-separated)" << endl
       << "-d, --distances  List of error distance (mismatch:indel:total) thresholds (comma-separated)" << endl
       << "-l, --locations  List of locations (file:pos1:pos2) (comma-separated)" << endl
       << "-i, --ids        List of tag names/identifiers (comma-separated)" << endl
       << "-g, --groups     List of tag group names (comma-separated)" << endl
       << "-f, --minFinds   List of minimum times a tag must be found in a read (comma-separated)" << endl
       << "-F, --maxFinds   List of maximum times a tag can be found in a read (comma-separated)" << endl
       << "-j, --minFindsG  List of minimum times tags in a group must be found in a read (comma-separated group_name:min_times)" << endl
       << "-J, --maxFindsG  List of maximum times tags in a group can be found in a read (comma-separated group_name:max_times)" << endl
       << "-e, --exclude    List of what to exclude from final barcode (comma-separated; 1 = exclude, 0 = include)" << endl
       << "-L, --left       List of what tags to include when trimming from the left (comma-separated; 1 = include, 0 = exclude)" << endl
       << "-R, --right      List of what tags to include when trimming from the right (comma-separated; 1 = include, 0 = exclude)" << endl
       << "                 (Note: for --left/--right, can specify an included tag as 1:x where x = number of extra bp's to trim" << endl
       << "                 from left/right side if that included tag is at the leftmost/rightmost position)" << endl
       << "-a, --next       List of what tag names must come immediately after each tag (comma-separated)" << endl
       << "-v, --previous   List of what tag names must come immediately before each tag (comma-separated)" << endl
       << "                 (Note: for --next/--previous, specify tag names as {name} and specify tag group names as {{group}}" << endl
       << "                 Can also specify the number of base pairs that must appear between the current tag and the next/previous tag." << endl
       << "                 E.g. {bc}4-12 means the next/previous tag is 4-12 bases away and has name 'bc')" << endl
       << "-U, --subs       Specifies sequence to substitute tag with when found in read (. = original sequence) (comma-separated)" << endl
       << "-z, --partial5   Specifies tag may be truncated at the 5′ end (comma-separated min_match:mismatch_freq)" << endl
       << "-Z, --partial3   Specifies tag may be truncated at the 3′ end (comma-separated min_match:mismatch_freq)" << endl
       << "Read modification and extraction options (for configuring on the command-line):" << endl
       << "-x, --extract    Pattern(s) describing how to extract UMI and UMI-like sequences from reads" << endl
       << "                 (E.g. {bc}2<umi_1[5]> means extract a 5-bp UMI sequence, called umi_1, 2 base pairs following the tag named 'bc')" << endl
       << "    --no-chain   If an extraction pattern for a UMI/UMI-like sequence is matched multiple times, only extract based on the first match" << endl
       << "-5, --trim-5     Number of base pairs to trim from the 5′-end of reads (comma-separated; one number per each FASTQ file in a run)" << endl
       << "-3, --trim-3     Number of base pairs to trim from the 3′-end of reads (comma-separated; one number per each FASTQ file in a run)" << endl
       << "-w, --filter-len Filter reads based on length (min_length:max_length)" << endl
       << "-q, --qtrim      Quality trimming threshold" << endl
       << "    --qtrim-5    Perform quality trimming from the 5′-end of reads of each FASTQ file" << endl
       << "    --qtrim-3    Perform quality trimming from the 3′-end of reads of each FASTQ file" << endl
       << "    --qtrim-pre  Perform quality trimming before sequence identification operations" << endl
       << "    --qtrim-naive Perform quality trimming using a naive algorithm (i.e. trim until a base that meets the quality threshold is encountered)" << endl
       << "    --phred64    Use phred+64 encoded quality scores" << endl
       << "-P, --prefix     Bases that will prefix each final barcode sequence (useful for merging separate experiments)" << endl
       << "Options (configurations supplied in a file):" << endl
       << "-c, --config     Configuration file" << endl
       << "Output Options:" << endl
       << "-m, --mapping    Output file where the mapping between final barcode sequences and names will be written" << endl
       << "-o, --output     FASTQ file(s) where output will be written (comma-separated)" << endl
       << "                 Number of output FASTQ files should equal --nFastqs (unless --select is provided)" << endl
       << "-O, --outb       FASTQ file where final barcodes will be written" << endl
       << "                 If not supplied, final barcodes are prepended to reads of first FASTQ file (or as the first read for --pipe)" << endl
       << "-u, --unassigned FASTQ file(s) where output of unassigned reads will be written (comma-separated)" << endl
       << "                 Number of FASTQ files should equal --nFastqs (unless --select is provided)" << endl
       << "-E, --empty      Sequence to fill in empty reads in output FASTQ files (default: no sequence is used to fill in those reads)" << endl
       << "    --empty-remove Empty reads are stripped in output FASTQ files (don't even output an empty sequence)" << endl
       << "-p, --pipe       Write to standard output (instead of output FASTQ files)" << endl
       << "-S, --select     Select which FASTQ files to output (comma-separated) (e.g. 0,1,3 = Output files #0, #1, #3)" << endl
       << "    --gzip       Output compressed gzip'ed FASTQ files" << endl
       << "    --out-fasta  Output in FASTA format rather than FASTQ format" << endl
       << "    --keep-com   Preserve the comments of the read names of the input FASTQ file(s)" << endl
       << "    --no-output  Don't output any sequences (output statistics only)" << endl
       << "    --no-outb    Don't output final barcode sequences" << endl
       << "    --no-x-out   Don't output extracted UMI-like sequences (should be used with --x-names)" << endl
       << "    --mod-names  Modify names of outputted sequences to include identified tag names" << endl
       << "    --com-names  Modify names of outputted sequences to include final barcode sequence ID" << endl
       << "    --seq-names  Modify names of outputted sequences to include the sequences of identified tags" << endl
       << "    --x-names    Modify names of outputted sequences to include extracted UMI-like sequences" << endl
       << "    --x-only     Only output extracted UMI-like sequences" << endl
       << "-X, --sub-assign Assign reads to a secondary sequence ID based on a subset of tags present (must be used with --assign)" << endl
       << "                 (e.g. 0,2 = Generate unique ID based the tags present by subsetting those tags to tag #0 and tag #2 only)" << endl
       << "                 The names of the outputted sequences will be modified to include this secondary sequence ID" << endl
       << "-M  --sam-tags   Modify the default SAM tags (default: CB:Z:,RX:Z:,BI:i:,SI:i:,BC:Z:)" << endl
       << "Other Options:" << endl
       << "-N, --nFastqs    Number of FASTQ file(s) per run" << endl
       << "                 (default: 1) (specify 2 for paired-end)" << endl
       << "-n, --numReads   Maximum number of reads to process from supplied input" << endl
       << "-A, --append     An existing mapping file that will be added on to" << endl
       << "-k, --keep       File containing a list of arrangements of tag names to keep" << endl
       << "-r, --remove     File containing a list of arrangements of tag names to remove/discard" << endl
       << "-y, --keep-grp   File containing a list of arrangements of tag groups to keep" << endl
       << "-Y, --remove-grp File containing a list of arrangements of tag groups to remove/discard" << endl
       << "-t, --threads    Number of threads to use" << endl
       << "-s, --summary    File where summary statistics will be written to" << endl
       << "-h, --help       Displays usage information" << endl
       << "    --assign     Assign reads to a final barcode sequence identifier based on tags present" << endl
       << "    --inleaved   Specifies that input is an interleaved FASTQ file" << endl
       << "    --remultiplex  Turn on remultiplexing mode" << endl
       << "    --version    Prints version number" << endl
       << "    --cite       Prints citation information" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions& opt) {
  int help_flag = 0;
  int version_flag = 0;
  int cite_flag = 0;
  int no_chain_flag = 0;
  int output_fasta_flag = 0;
  int no_output_flag = 0;
  int no_output_barcodes_flag = 0;
  int no_output_extracted_flag = 0;
  int gzip_flag = 0;
  int mod_names_flag = 0;
  int bc_names_flag = 0;
  int com_names_flag = 0;
  int seq_names_flag = 0;
  int x_names_flag = 0;
  int x_only_flag = 0;
  int empty_remove_flag = 0;
  int disable_n_flag = 0;
  int interleaved_flag = 0;
  int qtrim_5_flag = 0;
  int qtrim_3_flag = 0;
  int qtrim_pre_flag = 0;
  int qtrim_naive_flag = 0;
  int phred64_flag = 0;
  int assign_flag = 0;
  int keep_com_flag = 0;
  int remultiplex_flag = 0;
  bool trim_only_specified = false;

  const char *opt_string = "t:N:n:b:d:i:l:f:F:e:c:o:O:u:m:k:r:A:L:R:E:g:y:Y:j:J:a:v:z:Z:5:3:w:x:P:q:s:S:M:U:X:Tph";
  static struct option long_options[] = {
    // long args
    {"version", no_argument, &version_flag, 1},
    {"cite", no_argument, &cite_flag, 1},
    {"no-chain", no_argument, &no_chain_flag, 1},
    {"out-fasta", no_argument, &output_fasta_flag, 1},
    {"no-output", no_argument, &no_output_flag, 1},
    {"no-outb", no_argument, &no_output_barcodes_flag, 1},
    {"no-x-out", no_argument, &no_output_extracted_flag, 1},
    {"gzip", no_argument, &gzip_flag, 1},
    {"mod-names", no_argument, &mod_names_flag, 1},
    {"bc-names", no_argument, &bc_names_flag, 1},
    {"com-names", no_argument, &com_names_flag, 1},
    {"seq-names", no_argument, &seq_names_flag, 1},
    {"x-names", no_argument, &x_names_flag, 1},
    {"x-only", no_argument, &x_only_flag, 1},
    {"empty-remove", no_argument, &empty_remove_flag, 1},
    {"disable-n", no_argument, &disable_n_flag, 1},
    {"inleaved", no_argument, &interleaved_flag, 1},
    {"qtrim-5", no_argument, &qtrim_5_flag, 1},
    {"qtrim-3", no_argument, &qtrim_3_flag, 1},
    {"qtrim-pre", no_argument, &qtrim_pre_flag, 1},
    {"qtrim-naive", no_argument, &qtrim_naive_flag, 1},
    {"phred64", no_argument, &phred64_flag, 1},
    {"keep-com", no_argument, &keep_com_flag, 1},
    {"assign", no_argument, &assign_flag, 1},
    {"remultiplex", no_argument, &remultiplex_flag, 1},
    // short args
    {"help", no_argument, 0, 'h'},
    {"pipe", no_argument, 0, 'p'},
    {"trim-only", no_argument, 0, 'T'},
    {"threads", required_argument, 0, 't'},
    {"nFastqs", required_argument, 0, 'N'},
    {"numReads", required_argument, 0, 'n'},
    {"tags", required_argument, 0, 'b'},
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
    {"subs", required_argument, 0, 'U'},
    {"previous", required_argument, 0, 'v'},
    {"partial5", required_argument, 0, 'z'},
    {"partial3", required_argument, 0, 'Z'},
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
    {"filter-len", required_argument, 0, 'w'},
    {"extract", required_argument, 0, 'x'},
    {"prefix", required_argument, 0, 'P'},
    {"summary", required_argument, 0, 's'},
    {"select", required_argument, 0, 'S'},
    {"sam-tags", required_argument, 0, 'M'},
    {"sub-assign", required_argument, 0, 'X'},
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
      trim_only_specified = true;
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
    case 'U': {
      stringstream(optarg) >> opt.subs_str;
      break;
    }
    case 'z': {
      stringstream(optarg) >> opt.partial5_str;
      break;
    }
    case 'Z': {
      stringstream(optarg) >> opt.partial3_str;
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
    case 'w': {
      stringstream(optarg) >> opt.filter_length_str;
      break;
    }
    case 'x': {
      stringstream(optarg) >> opt.extract_str;
      break;
    }
    case 'P': {
      stringstream(optarg) >> opt.barcode_prefix;
      break;
    }
    case 'q': {
      stringstream(optarg) >> opt.quality_trimming_threshold;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.summary_file;
      break;
    }
    case 'S': {
      stringstream(optarg) >> opt.select_output_files_str;
      break;
    }
    case 'X': {
      std::string subset_n;
      stringstream(optarg) >> subset_n;
      std::stringstream ss(subset_n);
      while (std::getline(ss, subset_n, ',')) { 
        try {
          opt.sub_assign_vec.push_back(std::stoi(subset_n));
        } catch (std::exception &e) { }
      }
      std::sort(opt.sub_assign_vec.begin(), opt.sub_assign_vec.end());
      opt.sub_assign_vec.erase(std::unique(opt.sub_assign_vec.begin(), opt.sub_assign_vec.end()), opt.sub_assign_vec.end());
      break;
    }
    case 'M': {
      std::string m;
      stringstream(optarg) >> m;
      m.erase(remove(m.begin(),m.end(),' '),m.end()); // remove spaces from string
      std::stringstream ss(m);
      std::string s;
      int i = 0;
      while (std::getline(ss, s, ',') && i < 5) {
        if (i == 1) { // Allow multiple tags for extraction (default RX:Z:)
          opt.sam_tags[i].clear();
          std::stringstream ss2(s);
          std::string s2;
          while (std::getline(ss2, s2, '/')) { // Multiple tags separated by '/'
            opt.sam_tags[i].push_back(s2);
          }
          if (opt.sam_tags[i].empty()) {
            opt.sam_tags[i][0] = s;
          }
        } else {
          opt.sam_tags[i][0] = s;
        }
        i++;
      }
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
  if (no_chain_flag) {
    opt.extract_no_chain = true;
  }
  if (output_fasta_flag) {
    opt.output_fasta = true;
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
  if (bc_names_flag) {
    opt.bc_names = true;
  }
  if (seq_names_flag) {
    opt.seq_names = true;
  }
  if (x_names_flag) {
    opt.x_names = true;
  }
  if (x_only_flag) {
    opt.x_only = true;
  }
  if (empty_remove_flag) {
    opt.empty_remove = true;
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
  if (no_output_extracted_flag) {
    opt.no_x_out = true;
  }
  if (qtrim_5_flag) {
    opt.quality_trimming_5 = true;
  }
  if (qtrim_3_flag) {
    opt.quality_trimming_3 = true;
  }
  if (qtrim_pre_flag) {
    opt.quality_trimming_pre = true;
  }
  if (qtrim_naive_flag) {
    opt.quality_trimming_naive = true;
  }
  if (phred64_flag) {
    opt.phred64 = true;
  }
  if (keep_com_flag) {
    opt.keep_fastq_comments = true;
  }
  if (assign_flag && !trim_only_specified) {
    opt.trim_only = false;
  }
  if (remultiplex_flag) {
    opt.remultiplex = true;
  }
  
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  opt.select_output_files.resize(opt.nfiles, true);
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
  if (opt.remultiplex && opt.files.size() != 1) {
    cerr << ERROR_STR << " A single batch file must be supplied (for remultiplexing)" << endl;
    ret = false;
  } else if (opt.remultiplex && !opt.trim_only) {
    cerr << ERROR_STR << " --assign cannot be used with remultiplexing" << endl;
    ret = false;
  } else if (opt.remultiplex) {
    struct stat stFileInfo;
    auto intstat = stat(opt.files[0].c_str(), &stFileInfo);
    if (intstat != 0) {
      cerr << ERROR_STR << " batch file not found " << opt.files[0] << endl;
      ret = false;
    } else {
      // open the file, parse and fill the batch_files values
      std::ifstream bfile(opt.files[0]);
      opt.files.clear();
      std::string line;
      std::string id;
      std::vector<std::string> f_vec;
      size_t num_files = 0;
      size_t num_lines = 0;
      while (std::getline(bfile,line)) {
        if (line.size() == 0) {
          continue;
        }
        std::stringstream ss(line);
        ss >> id;
        if (id[0] == '#') {
          continue;
        }
        std::string f;
        size_t i;
        for (i = 0; ss >> f; i++) {
          f_vec.push_back(f);
        }
        if (num_files == 0) {
          num_files = i;
        }
        if (i == 0 || i != num_files) {
          cerr << ERROR_STR << " batch file malformatted" << endl;
          ret = false;
          break;
        }
        opt.batch_ids.push_back(id);
        num_lines++;
      }
      if (opt.nfiles > 1 && num_files != opt.nfiles) {
        cerr << ERROR_STR << " nFastqs inconsistent with number of files in batch file" << endl;
        ret = false;
      } else {
        opt.nfiles = num_files;
        opt.select_output_files.resize(opt.nfiles, true);
        sc.setNFiles(opt.nfiles);
        for (size_t i = 0; i < num_lines; i++) {
          for (size_t j = 0; j < opt.nfiles; j++) {
            opt.files.push_back(f_vec[i*opt.nfiles+j]);
          }
        }
      }
    }
  }
  if (opt.files.size() == 0) {
    cerr << ERROR_STR << " Missing read files" << endl;
    ret = false;
  } else if (!(opt.files.size() == 1 && opt.files[0] == "-") || opt.remultiplex) { // If not reading from stdin via -
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
    if (opt.remultiplex) {
      std::cerr << ERROR_STR << " interleaved input cannot be used with remultiplexing" << std::endl;
      ret = false;
    } else if (opt.files.size() != 1) {
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
  if (opt.x_only && opt.no_x_out) {
    std::cerr << ERROR_STR << " --x-only cannot be specified with --no-x-out" << std::endl;
    ret = false;
  }
  if (!opt.empty_read_sequence.empty() && opt.empty_remove) {
    std::cerr << ERROR_STR << " --empty cannot be specified with --empty-remove" << std::endl;
    ret = false;
  }
  if (!opt.sub_assign_vec.empty() && opt.trim_only) {
    std::cout << "Cannot use --sub-assign unless --assign is specified" << std::endl;
    ret = false;
  }
  int nf = opt.nfiles;
  if (!opt.select_output_files_str.empty()) {
    nf = 0;
    opt.select_output_files.assign(opt.nfiles, false);
    try {
      std::stringstream ss(opt.select_output_files_str);
      std::string s;
      while (getline(ss, s, ',')) {
        int f = std::stoi(s);
        if (f < 0) {
          std::cerr << ERROR_STR << " --select must contain numbers >= 0" << std::endl;
          ret = false;
          break;
        } else if (f >= opt.nfiles) {
          std::cerr << ERROR_STR << " --select must contain numbers less than --nFastqs" << std::endl;
          ret = false;
          break;
        }
        opt.select_output_files[f] = true;
      }
      for (int i = 0; i < opt.select_output_files.size(); i++) {
        if (opt.select_output_files[i]) nf++;
      }
    } catch (std::exception &e) {
      std::cerr << ERROR_STR << " --select must contain numbers >= 0" << std::endl;
      ret = false;
    }
  }
  
  bool output_files_specified = opt.output_files.size() > 0 || opt.unassigned_files.size() > 0 || !opt.outputb_file.empty();
  if (opt.output_files.size() == 0 && output_files_specified && !opt.pipe && !opt.x_only) {
    std::cerr << ERROR_STR << " --output not provided" << std::endl;
    ret = false;
  }
  if (opt.no_output) {
    if (output_files_specified || opt.pipe) {
      std::cerr << ERROR_STR << " Cannot specify an output option when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.mod_names || opt.com_names || opt.x_names || opt.seq_names) {
      std::cerr << ERROR_STR << " Cannot use --mod-names/--com-names/--seq-names/--x-names when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.gzip) {
      std::cerr << ERROR_STR << " Cannot use --gzip when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.x_only) {
      std::cerr << ERROR_STR << " Cannot use --x-only when --no-output is specified" << std::endl;
      ret = false;
    }
    if (opt.output_fasta) {
      std::cerr << ERROR_STR << " Cannot use --out-fasta when --no-output is specified" << std::endl;
      ret = false;
    }
    if (!opt.select_output_files_str.empty()) {
      std::cerr << ERROR_STR << " Cannot use --select when --no-output is specified" << std::endl;
      ret = false;
    }
  } else {
    if (!output_files_specified && !opt.pipe && !opt.x_only) {
      std::cerr << ERROR_STR << " Must either specify an output option or --no-output" << std::endl;
      ret = false;
    } else if (opt.x_only) {
      if (opt.output_files.size() > 0 || opt.unassigned_files.size() > 0) {
        std::cerr << ERROR_STR << " Cannot provide output files when --x-only is specified" << std::endl;
        ret = false;
      }
    } else if (opt.pipe) {
      if (opt.output_files.size() > 0 || !opt.outputb_file.empty()) { // Still allow --unassigned with --pipe
        std::cerr << ERROR_STR << " Cannot provide output files when --pipe is specified" << std::endl;
        ret = false;
      } else if (opt.unassigned_files.size() != 0 && opt.unassigned_files.size() % nf != 0 || opt.unassigned_files.size() > nf) {
        std::cerr << ERROR_STR << " Incorrect number of --unassigned output files" << std::endl;
        ret = false;
      }
    } else {
      if (opt.output_files.size() % nf != 0 || opt.unassigned_files.size() % nf != 0 || opt.output_files.size() > nf || opt.unassigned_files.size() > nf) {
        std::cerr << ERROR_STR << " Incorrect number of output files" << std::endl;
        ret = false;
      }
    }
  }
  if (opt.trim_only && !opt.outputb_file.empty() && !opt.remultiplex) {
    std::cerr << ERROR_STR << " Cannot use --outb unless --assign is specified" << std::endl;
    ret = false;
  }
  if (opt.trim_only && !opt.mapping_file.empty() && !opt.remultiplex) {
    std::cerr << ERROR_STR << " Cannot use --mapping unless --assign is specified" << std::endl;
    ret = false;
  }
  if (opt.trim_only && opt.com_names && !opt.remultiplex) {
    std::cerr << ERROR_STR << " Cannot use --com-names unless --assign is specified" << std::endl;
    ret = false;
  }
  if (opt.trim_only && !opt.append_file.empty()) {
    std::cerr << ERROR_STR << " Cannot use --append unless --assign is specified" << std::endl;
    ret = false;
  }
  opt.output_fastq_specified = output_files_specified;
  opt.verbose = !opt.pipe;
  
  int num_groups = 0;
  if (!opt.barcode_str.empty() && !opt.config_file.empty()) {
    std::cerr << ERROR_STR << " Cannot specify both --tags and --config" << std::endl;
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
    stringstream ss13(opt.partial5_str);
    stringstream ss14(opt.partial3_str);
    stringstream ss15(opt.subs_str);
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
      string partial5_str = "";
      string partial3_str = "";
      string subs_str = "";
      int partial5_min_match, partial3_min_match;
      double partial5_mismatch_freq, partial3_mismatch_freq;
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
          std::cerr << ERROR_STR << " Number of values in --distances is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --ids is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss3, name, ',');
      }
      if (!opt.location_str.empty()) {
        auto currpos = ss4.tellg();
        if (!ss4.good()) {
          std::cerr << ERROR_STR << " Number of values in --locations is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --maxFinds is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --minFinds is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --exclude is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --left is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --right is less than that in --tags" << std::endl;
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
        std::cerr << ERROR_STR << " One of the tags has both --left and --right trimming specified" << std::endl;
        ret = false;
        break;
      }
      auto trim_dir = trim_left ? sc.left : (trim_right ? sc.right : sc.nodir);
      auto trim_offset = trim_left ? trim_left_offset : (trim_right ? trim_right_offset : 0);
      if (!opt.group_identifiers_str.empty()) {
        if (!ss10.good()) {
          std::cerr << ERROR_STR << " Number of values in --groups is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss10, group, ',');
        num_groups++;
      }
      if (!opt.after_str.empty()) {
        if (!ss11.good()) {
          std::cerr << ERROR_STR << " Number of values in --next is less than that in --tags" << std::endl;
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
          std::cerr << ERROR_STR << " Number of values in --previous is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss12, before_str, ',');
      }
      if (!SplitCode::validateBeforeAfterStr(before_str)) {
        std::cerr << ERROR_STR << " --previous is invalid" << std::endl;
        ret = false;
      }
      if (!opt.partial5_str.empty()) {
        if (!ss13.good()) {
          std::cerr << ERROR_STR << " Number of values in --partial5 is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss13, partial5_str, ',');
      }
      if (!SplitCode::parsePartialStr(partial5_str, partial5_min_match, partial5_mismatch_freq)) {
        std::cerr << ERROR_STR << " --partial5 is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.partial3_str.empty()) {
        if (!ss14.good()) {
          std::cerr << ERROR_STR << " Number of values in --partial3 is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss14, partial3_str, ',');
      }
      if (!SplitCode::parsePartialStr(partial3_str, partial3_min_match, partial3_mismatch_freq)) {
        std::cerr << ERROR_STR << " --partial3 is invalid" << std::endl;
        ret = false;
        break;
      }
      if (!opt.subs_str.empty()) {
        if (!ss15.good()) {
          std::cerr << ERROR_STR << " Number of values in --subs is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        getline(ss15, subs_str, ',');
      }
      if (!sc.addTag(bc, name.empty() ? bc : name, group, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude, trim_dir, trim_offset, after_str, before_str, partial5_min_match, partial5_mismatch_freq, partial3_min_match, partial3_mismatch_freq, subs_str)) {
        std::cerr << ERROR_STR << " Could not finish processing supplied tags list" << std::endl;
        ret = false;
        break;
      }
    }
    if (ret && !opt.distance_str.empty() && ss2.good()) {
      std::cerr << ERROR_STR << " Number of values in --distances is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.barcode_identifiers_str.empty() && ss3.good()) {
      std::cerr << ERROR_STR << " Number of values in --ids is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.location_str.empty() && ss4.good()) {
      std::cerr << ERROR_STR << " Number of values in --locations is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.max_finds_str.empty() && ss5.good()) {
      std::cerr << ERROR_STR << " Number of values in --maxFinds is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.min_finds_str.empty() && ss6.good()) {
      std::cerr << ERROR_STR << " Number of values in --minFinds is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.exclude_str.empty() && ss7.good()) {
      std::cerr << ERROR_STR << " Number of values in --exclude is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.left_str.empty() && ss8.good()) {
      std::cerr << ERROR_STR << " Number of values in --left is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.right_str.empty() && ss9.good()) {
      std::cerr << ERROR_STR << " Number of values in --right is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.group_identifiers_str.empty() && ss10.good()) {
      std::cerr << ERROR_STR << " Number of values in --groups is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.after_str.empty() && ss11.good()) {
      std::cerr << ERROR_STR << " Number of values in --next is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.before_str.empty() && ss12.good()) {
      std::cerr << ERROR_STR << " Number of values in --previous is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.partial5_str.empty() && ss13.good()) {
      std::cerr << ERROR_STR << " Number of values in --partial5 is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.partial3_str.empty() && ss14.good()) {
      std::cerr << ERROR_STR << " Number of values in --partial3 is greater than that in --tags" << std::endl;
      ret = false;
    }
    if (ret && !opt.subs_str.empty() && ss15.good()) {
      std::cerr << ERROR_STR << " Number of values in --subs is greater than that in --tags" << std::endl;
      ret = false;
    }
  } else if (!opt.distance_str.empty()) {
    std::cerr << ERROR_STR << " --distances cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.barcode_identifiers_str.empty()) {
    std::cerr << ERROR_STR << " --ids cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.location_str.empty()) {
    std::cerr << ERROR_STR << " --locations cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.max_finds_str.empty()) {
    std::cerr << ERROR_STR << " --maxFinds cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.min_finds_str.empty()) {
    std::cerr << ERROR_STR << " --minFinds cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.exclude_str.empty()) {
    std::cerr << ERROR_STR << " --exclude cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.left_str.empty()) {
    std::cerr << ERROR_STR << " --left cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.right_str.empty()) {
    std::cerr << ERROR_STR << " --right cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.group_identifiers_str.empty()) {
    std::cerr << ERROR_STR << " --groups cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.after_str.empty()) {
    std::cerr << ERROR_STR << " --next cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.before_str.empty()) {
    std::cerr << ERROR_STR << " --previous cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.partial5_str.empty()) {
    std::cerr << ERROR_STR << " --partial5 cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.partial3_str.empty()) {
    std::cerr << ERROR_STR << " --partial3 cannot be supplied unless --tags is" << std::endl;
    ret = false;
  } else if (!opt.subs_str.empty()) {
    std::cerr << ERROR_STR << " --subs cannot be supplied unless --tags is" << std::endl;
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
    /*std::cerr << ERROR_STR << " No tags found" << std::endl;
    ret = false;*/
    sc.checkInit();
  }
  
  return ret;
}


int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  SplitCode sc(opt.nfiles, opt.summary_file, opt.trim_only, opt.disable_n, opt.trim_5_str, opt.trim_3_str, opt.extract_str, opt.extract_no_chain, opt.barcode_prefix, opt.filter_length_str,
               opt.quality_trimming_5, opt.quality_trimming_3, opt.quality_trimming_pre, opt.quality_trimming_naive, opt.quality_trimming_threshold, opt.phred64, opt.sub_assign_vec);
  bool checkopts = CheckOptions(opt, sc);
  if (!checkopts) {
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
      " tags (vector size: " << sc.getNumTags() << 
      "; map size: " << pretty_num(sc.getMapSize()) << 
      "; num elements in map: " << pretty_num(sc.getMapSize(false)) << ")" << std::endl;
  }
  MasterProcessor MP(sc, opt);
  int numreads = ProcessReads(MP, opt);
  fflush(stdout);
  if (!opt.mapping_file.empty()) { // output mapping file:
    if (opt.remultiplex) { // remultiplexing mapping file
      size_t num_remultiplexed_barcodes = opt.files.size() / opt.nfiles;
      if (!(opt.mapping_file.size() > 3 && opt.mapping_file.compare(opt.mapping_file.size() - 3, 3, ".gz") == 0)) { // plaintext
        std::ofstream of;
        of.open(opt.mapping_file);
        if (!of.is_open()) {
          std::cerr << "Error: Couldn't open file: " << opt.mapping_file << std::endl;
        } else {
          for (size_t i = 0; i < opt.batch_ids.size(); i++) {
            std::string o = sc.binaryToString(sc.getID(i), sc.getBarcodeLength()) + "\t" + opt.batch_ids[i] + "\n";
            of << o;
          }
          of.close();
        }
      } else { // gzip
        gzFile out_gz = gzopen(opt.mapping_file.c_str(), "wb6");
        for (size_t i = 0; i < opt.batch_ids.size(); i++) {
          std::string o = sc.binaryToString(sc.getID(i), sc.getBarcodeLength()) + "\t" + opt.batch_ids[i] + "\n";
          gzwrite(out_gz, o.c_str(), o.length());
        }
        gzclose(out_gz);
      }
    } else if (!(opt.mapping_file.size() > 3 && opt.mapping_file.compare(opt.mapping_file.size() - 3, 3, ".gz") == 0)) {
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
  
  sc.writeSummary(argv_to_string(argc, argv));
  
  if (opt.max_num_reads != 0 && numreads < opt.max_num_reads) {
    std::cerr << "Note: Number of reads processed is less than --numReads: " << opt.max_num_reads << ", returning 1" << std::endl;
    return 1;
  }

  return 0;
}

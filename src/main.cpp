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
#include "LiftWorkflow.h"


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

// This is part of the --unmask non-standard splitcode workflow
void runUnmaskingWorkflow(ProgramOptions& opt) {
  if (opt.files.size() != 2 && opt.files.size() != 3) {
    std::cerr << "--unmask requires two FASTA files supplied: a masked and an unmasked file" << std::endl;
    std::cerr << "\n" << "splitcode --unmask file1.unmasked.fasta[.gz] file2.masked.fasta[.gz] [fasta_header]" << "\n" << std::endl;
    exit(1);
  }
  bool header_supplied = (opt.files.size() == 3);
  std::string header;
  if (header_supplied) header = opt.files[2];

  gzFile fp1 = gzopen(opt.files[0].c_str(), "rb");
  gzFile fp2 = gzopen(opt.files[1].c_str(), "rb");

  if (!fp1 || !fp2) {
    std::cerr << "Error opening files." << std::endl;
    if (fp1) gzclose(fp1);
    if (fp2) gzclose(fp2);
    exit(1);
  }

  kseq_t *seq1 = kseq_init(fp1);
  kseq_t *seq2 = kseq_init(fp2);

  int l1, l2;
  int counter = 1;
  while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0) {
        std::string s1 = seq1->seq.s;
        std::string s2 = seq2->seq.s;
        int diff_start = -1;
        size_t min_len = std::min(s1.size(), s2.size());

        for (size_t i = 0; i < min_len; ++i) {
            if (s1[i] != s2[i]) {
                if (diff_start == -1) diff_start = i; // Start of a new difference
            } else {
                if (diff_start != -1) {
                    std::cout << ">" << (header_supplied ? header : std::to_string(counter++)) << std::endl;
                    std::cout << s1.substr(diff_start, i - diff_start) << std::endl;
                    diff_start = -1; // Reset diff_start
                }
            }
        }

        // Handle difference at the end of the sequence
        if (diff_start != -1) {
            std::cout << ">" << (header_supplied ? header : std::to_string(counter++)) << std::endl;
            std::cout << s1.substr(diff_start, min_len - diff_start) << std::endl;
        }
  }

  kseq_destroy(seq1);
  kseq_destroy(seq2);
  gzclose(fp1);
  gzclose(fp2);

  exit(0);

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
       << "    --revcomp    Specifies tag may be reverse complemented" << endl
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
       << "-D, --min-delta  When matching tags error-tolerantly, specifies how much worse the next best match must be than the best match" << endl
       << "    --from-name  Extract sequences from FASTQ header comments. Format: fastq_number,output_file_number,output_position,pattern." << endl
       << "                 (Example: 0,0,0,::;0,0,0,::+ will extract the nucleotides from 1:N:ATCCC+ATCG and put it into the R1 output)" << endl
       << "    --random     Insert a random sequence. Format: output_file_number,output_position,length." << endl
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
       << "    --out-bam    Output a BAM file rather than FASTQ files (enter the output BAM file name to -o or --output)" << endl
       << "    --keep-com   Preserve the comments of the read names of the input FASTQ file(s)" << endl
       << "    --no-output  Don't output any sequences" << endl
       << "    --no-outb    Don't output final barcode sequences" << endl
       << "    --no-x-out   Don't output extracted UMI-like sequences (should be used with --x-names)" << endl
       << "    --mod-names  Modify names of outputted sequences to include identified tag names" << endl
       << "    --com-names  Modify names of outputted sequences to include final barcode sequence ID" << endl
       << "    --seq-names  Modify names of outputted sequences to include the sequences of identified tags" << endl
       << "    --loc-names  Modify names of outputted sequences to include found tag names and locations" << endl
       << "    --x-names    Modify names of outputted sequences to include extracted UMI-like sequences" << endl
       << "    --x-only     Only output extracted UMI-like sequences" << endl
       << "    --bc-names   Modify names of outputted sequences to include final barcode sequence string" << endl
       << "-X, --sub-assign Assign reads to a secondary sequence ID based on a subset of tags present (must be used with --assign)" << endl
       << "                 (e.g. 0,2 = Generate unique ID based the tags present by subsetting those tags to tag #0 and tag #2 only)" << endl
       << "                 The names of the outputted sequences will be modified to include this secondary sequence ID" << endl
       << "-C  --compress   Set the gzip compression level (default: 1) (range: 1-9)" << endl
       << "-M  --sam-tags   Modify the default SAM tags (default: CB:Z:,RX:Z:,BI:i:,SI:i:,BC:Z:,LX:Z:,YM:Z:)" << endl
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
       << "    --barcode-encode Optimize barcode assignment using a sequence of group names (e.g. group1,group2,group3)" << endl
       << "    --bclen      The length of the final barcode sequence identifier (default: 16)" << endl
       << "    --inleaved   Specifies that input is an interleaved FASTQ file" << endl
       << "    --keep-r1-r2 Use R1.fastq, R2.fastq, etc. file name formats when demultiplexing using --keep or --keep-grp" << endl
       << "    --remultiplex  Turn on remultiplexing mode" << endl
       << "    --unmask       Turn on unmasking mode (extract differences from a masked vs. unmasked FASTA)" << endl
       << "    --lift         Turn lift mode (make variant genomes from VCF files)" << endl
       << "    --version    Prints version number" << endl
       << "    --cite       Prints citation information" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions& opt) {
  int help_flag = 0;
  int version_flag = 0;
  int cite_flag = 0;
  int no_chain_flag = 0;
  int output_fasta_flag = 0;
  int output_bam_flag = 0;
  int no_output_flag = 0;
  int no_output_barcodes_flag = 0;
  int no_output_extracted_flag = 0;
  int gzip_flag = 0;
  int mod_names_flag = 0;
  int mod_names_bam_flag = 0;
  int bc_names_flag = 0;
  int com_names_flag = 0;
  int seq_names_flag = 0;
  int loc_names_flag = 0;
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
  int unmask_flag = 0;
  int keep_r1_r2_flag = 0;
  int webasm_flag = 0;
  int show_not_found_flag = 0;
  bool trim_only_specified = false;
  
  // Some --lift specific options
  int lift_flag = 0;
  int lift_diploid = 0;
  int lift_rename = 0;
  int lift_snvonly = 0;
  int lift_filter = 0;
  std::string lift_ref_gtf;
  std::string lift_out_gtf;
  std::string lift_kmer_length;
  std::string lift_kmer_output;
  std::string lift_kmer_header;
  int lift_kmer_header_num = 0;
  int lift_kmer_sj = 0;

  optind=1; // Reset global variable in case we want to call ParseOptions multiple times

  const char *opt_string = "t:N:n:b:B:d:D:i:l:f:F:e:c:o:O:u:m:k:r:A:L:R:E:g:y:Y:j:J:a:v:z:Z:5:3:w:x:P:q:s:S:M:U:X:C:Tph";
  /*static*/ struct option long_options[] = { // No static keyword because we may want to call ParseOptions multiple times
    // long args
    {"version", no_argument, &version_flag, 1},
    {"cite", no_argument, &cite_flag, 1},
    {"no-chain", no_argument, &no_chain_flag, 1},
    {"out-fasta", no_argument, &output_fasta_flag, 1},
    {"out-bam", no_argument, &output_bam_flag, 1},
    {"no-output", no_argument, &no_output_flag, 1},
    {"no-outb", no_argument, &no_output_barcodes_flag, 1},
    {"no-x-out", no_argument, &no_output_extracted_flag, 1},
    {"gzip", no_argument, &gzip_flag, 1},
    {"mod-names", no_argument, &mod_names_flag, 1},
    {"mod-names-bam", no_argument, &mod_names_bam_flag, 1},
    {"bc-names", no_argument, &bc_names_flag, 1},
    {"com-names", no_argument, &com_names_flag, 1},
    {"seq-names", no_argument, &seq_names_flag, 1},
    {"loc-names", no_argument, &loc_names_flag, 1},
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
    {"unmask", no_argument, &unmask_flag, 1},
    {"keep-r1-r2", no_argument, &keep_r1_r2_flag, 1},
    {"show-not-found", no_argument, &show_not_found_flag, 1},
    {"webasm", no_argument, &webasm_flag, 1},
    {"lift", no_argument, &lift_flag, 1},
    {"diploid", no_argument, &lift_diploid, 1},
    {"rename", no_argument, &lift_rename, 1},
    {"snv-only", no_argument, &lift_snvonly, 1},
    {"filter", no_argument, &lift_filter, 1},
    {"kmer-header-num", no_argument, &lift_kmer_header_num, 1},
    {"kmer-sj", no_argument, &lift_kmer_sj, 1},
    // short args
    {"help", no_argument, 0, 'h'},
    {"pipe", no_argument, 0, 'p'},
    {"trim-only", no_argument, 0, 'T'},
    {"threads", required_argument, 0, 't'},
    {"nFastqs", required_argument, 0, 'N'},
    {"numReads", required_argument, 0, 'n'},
    {"tags", required_argument, 0, 'b'},
    {"tag", required_argument, 0, 'b'},
    {"distances", required_argument, 0, 'd'},
    {"distance", required_argument, 0, 'd'},
    {"min-delta", required_argument, 0, 'D'},
    {"locations", required_argument, 0, 'l'},
    {"location", required_argument, 0, 'l'},
    {"ids", required_argument, 0, 'i'},
    {"id", required_argument, 0, 'i'},
    {"groups", required_argument, 0, 'g'},
    {"group", required_argument, 0, 'g'},
    {"maxFinds", required_argument, 0, 'F'},
    {"minFinds", required_argument, 0, 'f'},
    {"maxFindsG", required_argument, 0, 'J'},
    {"minFindsG", required_argument, 0, 'j'},
    {"exclude", required_argument, 0, 'e'},
    {"barcode-encode", required_argument, 0, 'B'},
    {"next", required_argument, 0, 'a'},
    {"after", required_argument, 0, 'a'},
    {"subs", required_argument, 0, 'U'},
    {"sub", required_argument, 0, 'U'},
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
    {"compress", required_argument, 0, 'C'},
    {"bclen", required_argument, 0, '9'},
    {"revcomp", required_argument, 0, 0},
    {"from-name", required_argument, 0, 0},
    {"random", required_argument, 0, 0},
    {"ref-gtf", required_argument, 0, 0},
    {"out-gtf", required_argument, 0, 0},
    {"kmer-length", required_argument, 0, 0},
    {"kmer-output", required_argument, 0, 0},
    {"kmer-header", required_argument, 0, 0},
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
      if (strcmp(long_options[option_index].name,"revcomp") == 0) opt.revcomp_str = optarg;
      if (strcmp(long_options[option_index].name,"from-name") == 0) opt.from_header_str = optarg;
      if (strcmp(long_options[option_index].name,"random") == 0) opt.random_str = optarg;
      if (strcmp(long_options[option_index].name,"ref-gtf") == 0) lift_ref_gtf = optarg;
      if (strcmp(long_options[option_index].name,"out-gtf") == 0) lift_out_gtf = optarg;
      if (strcmp(long_options[option_index].name,"kmer-length") == 0) lift_kmer_length = optarg;
      if (strcmp(long_options[option_index].name,"kmer-output") == 0) lift_kmer_output = optarg;
      if (strcmp(long_options[option_index].name,"kmer-header") == 0) lift_kmer_header = optarg;
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
      opt.barcode_str = optarg;
      break;
    }
    case 'd': {
      opt.distance_str = optarg;
      break;
    }
    case 'D': {
      stringstream(optarg) >> opt.min_delta;
      if (opt.min_delta < 0) {
        opt.min_delta = -1;
      }
      break;
    }
    case 'B': {
      opt.optimize_assignment_str = optarg;
      break;
    }
    case 'l': {
      opt.location_str = optarg;
      break;
    }
    case 'i': {
      opt.barcode_identifiers_str = optarg;
      break;
    }
    case 'g': {
      opt.group_identifiers_str = optarg;
      break;
    }
    case 'F': {
      opt.max_finds_str = optarg;
      break;
    }
    case 'f': {
      opt.min_finds_str = optarg;
      break;
    }
    case 'J': {
      opt.max_finds_group_str = optarg;
      break;
    }
    case 'j': {
      opt.min_finds_group_str = optarg;
      break;
    }
    case 'e': {
      opt.exclude_str = optarg;
      break;
    }
    case 'L': {
      opt.left_str = optarg;
      break;
    }
    case 'R': {
      opt.right_str = optarg;
      break;
    }
    case 'a': {
      opt.after_str = optarg;
      break;
    }
    case 'v': {
      opt.before_str = optarg;
      break;
    }
    case 'U': {
      opt.subs_str = optarg;
      break;
    }
    case 'z': {
      opt.partial5_str = optarg;
      break;
    }
    case 'Z': {
      opt.partial3_str = optarg;
      break;
    }
    case 'c': {
      opt.config_file = optarg;
      break;
    }
    case 'm': {
      opt.mapping_file = optarg;
      break;
    }
    case 'k': {
      opt.keep_file = optarg;
      break;
    }
    case 'r': {
      opt.discard = true;
      opt.keep_file = optarg;
      break;
    }
    case 'y': {
      opt.keep_group_file = optarg;
      break;
    }
    case 'Y': {
      opt.discard_group = true;
      opt.keep_group_file = optarg;
      break;
    }
    case 'o': {
      std::string files = optarg;
      std::stringstream ss(files);
      std::string filename;
      while (std::getline(ss, filename, ',')) { 
        opt.output_files.push_back(filename);
      }
      break;
    }
    case 'O': {
      opt.outputb_file = optarg;
      break;
    }
    case 'A': {
      opt.append_file = optarg;
      break;
    }
    case 'E': {
      opt.empty_read_sequence = optarg;
      for (auto& c: opt.empty_read_sequence) {
        c = toupper(c);
      }
      break;
    }
    case 'u': {
      std::string files = optarg;
      std::stringstream ss(files);
      std::string filename;
      while (std::getline(ss, filename, ',')) { 
        opt.unassigned_files.push_back(filename);
      }
      if (!files.empty() && files.back() == ',') {
        opt.unassigned_files.push_back(""); // Allow an unspecified file (i.e. no file name)
      }
      break;
    }
    case '5': {
      opt.trim_5_str = optarg;
      break;
    }
    case '3': {
      opt.trim_3_str = optarg;
      break;
    }
    case '9': {
      stringstream(optarg) >> opt.bclen;
      break;
    }
    case 'w': {
      opt.filter_length_str = optarg;
      break;
    }
    case 'x': {
      opt.extract_str = optarg;
      break;
    }
    case 'P': {
      opt.barcode_prefix = optarg;
      break;
    }
    case 'q': {
      stringstream(optarg) >> opt.quality_trimming_threshold;
      break;
    }
    case 's': {
      opt.summary_file = optarg;
      break;
    }
    case 'S': {
      opt.select_output_files_str = optarg;
      break;
    }
    case 'C': {
      std::string compress_level;
      stringstream(optarg) >> compress_level;
      try {
        opt.compress_level = std::stoi(compress_level);
      } catch (std::exception &e) { }
      if (opt.compress_level < 1) opt.compress_level = 1;
      if (opt.compress_level > 9) opt.compress_level = 9;
      break;
    }
    case 'X': {
      std::string subset_n = optarg;
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
      std::string m = optarg;
      m.erase(remove(m.begin(),m.end(),' '),m.end()); // remove spaces from string
      std::stringstream ss(m);
      std::string s;
      int i = 0;
      while (std::getline(ss, s, ',') && i <= 6) {
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
  if (output_bam_flag) {
    opt.outbam = true;
    if (opt.pipe) opt.outbampipe = true; // Write BAM to stdout
    opt.pipe = true; // Always pretend there's a pipe (b/c outputting a BAM is similar to having piped output where only one thing is outputted)
  }
  if (no_output_flag) {
    opt.no_output = true;
  }
  if (mod_names_flag) {
    opt.mod_names = true;
  }
  if (mod_names_bam_flag) {
    opt.mod_names_bam = true;
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
  if (loc_names_flag) {
    opt.write_locations = true;
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
  if (keep_r1_r2_flag) {
    opt.keep_r1_r2 = true;
  }
  if (webasm_flag) {
    opt.threads = 1;
    opt.webasm = true;
  }
  if (show_not_found_flag) {
    opt.show_not_found = true;
  }
  
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  opt.select_output_files.resize(opt.nfiles, true);

  // Now for specialized workflows not part of the main "splitcode" workflow
  if (unmask_flag) {
    runUnmaskingWorkflow(opt);
    exit(0);
  } else if (lift_flag) {
    LiftWorkflow lf(opt.files, (bool)lift_diploid, (bool)lift_rename, (bool)lift_snvonly, lift_ref_gtf, lift_out_gtf, (bool)lift_filter, lift_kmer_length.empty() ? 0 : std::atoi(lift_kmer_length.c_str()), lift_kmer_output, lift_kmer_header, (bool)lift_kmer_header_num, (bool)lift_kmer_sj);
    lf.modify_fasta();
    exit(0);
  }
}

bool CheckOptions(ProgramOptions& opt, SplitCode& sc) {
  bool ret = true;

  if (opt.outbam) { // BAM options
#ifdef NO_HTSLIB
    std::cerr << ERROR_STR << " BAM files not supported because splitcode was compiled without BAM file support" << std::endl;
    return false;
#endif
    opt.pipe = true;
    opt.gzip = false;
    if (!opt.outbampipe) {
      if (opt.output_files.size() == 0) {
        std::cerr << ERROR_STR << " Must supply an output BAM file name to --output" << std::endl;
        ret = false;
      } else if (opt.output_files.size() > 1) {
        std::cerr << ERROR_STR << " Cannot supply multiple BAM file names to --output" << std::endl;
        ret = false;
      }
      opt.outbamfile = opt.output_files[0];
      opt.output_files.clear();
    }
    if (opt.phred64) {
      std::cerr << ERROR_STR << " --phred64 incompatible with writing BAM files" << std::endl;
    }
    if (opt.output_fasta) {
      std::cerr << ERROR_STR << " Cannot use --out-fasta when when writing BAM files" << std::endl;
      ret = false;
    }
    if (opt.x_only) {
      std::cerr << ERROR_STR << " Cannot use --x-only when when writing BAM files" << std::endl;
      ret = false;
    }
    if (opt.keep_fastq_comments && opt.mod_names_bam) {
      std::cerr << ERROR_STR << " Cannot use --mod-names-bam with --keep-com" << std::endl;
      ret = false;
    }
    if (!ret) return ret;
  }

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
  if (opt.bclen != 0) {
    if (!opt.barcode_prefix.empty()) {
      cerr << ERROR_STR << " Cannot specify --prefix with --bclen" << endl;
      ret = false;
    }
    if (opt.bclen >= 32 || opt.bclen < 2) {
      cerr << ERROR_STR << " --bclen must have value between 2 and 32 " << endl;
      ret = false;
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
        } /*else if (f >= opt.nfiles) {
          std::cerr << ERROR_STR << " --select must contain numbers less than --nFastqs" << std::endl;
          ret = false;
          break;
        }*/
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
  
  if (opt.no_output && !opt.pipe && opt.output_files.size() == 0 && opt.outputb_file.empty()) {
    // Override --no-output by replacing it with --pipe (but just don't write to stdout)
    opt.pipe = true;
    opt.no_output = false;
    opt.no_output_ = true;
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
      if (opt.output_files.size() > 0/* || opt.unassigned_files.size() > 0*/) {
        std::cerr << ERROR_STR << " Cannot provide output files when --x-only is specified" << std::endl;
        ret = false;
      }
    } else if (opt.pipe) {
      if (opt.output_files.size() > 0 || !opt.outputb_file.empty()) { // Still allow --unassigned with --pipe
        std::cerr << ERROR_STR << " Cannot provide output files when --pipe is specified" << std::endl;
        ret = false;
      } /* else if (opt.unassigned_files.size() != 0 && opt.unassigned_files.size() % nf != 0 || opt.unassigned_files.size() > nf) {
        std::cerr << ERROR_STR << " Incorrect number of --unassigned output files" << std::endl;
        ret = false;
      } */ 
    } else {
      /* if (opt.output_files.size() % nf != 0 || opt.unassigned_files.size() % nf != 0 || opt.output_files.size() > nf || opt.unassigned_files.size() > nf) {
        std::cerr << ERROR_STR << " Incorrect number of output files" << std::endl;
        ret = false;
      }*/ 
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
  opt.verbose = !opt.pipe || (opt.pipe && opt.no_output_);

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
    stringstream ss16(opt.revcomp_str);
    while (ss1.good()) {
      uint16_t max_finds = 0;
      uint16_t min_finds = 0;
      bool exclude = false;
      bool revcomp = false;
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
      if (!opt.revcomp_str.empty()) {
        if (!ss16.good()) {
          std::cerr << ERROR_STR << " Number of values in --revcomp is less than that in --tags" << std::endl;
          ret = false;
          break;
        }
        string f;
        getline(ss16, f, ',');
        stringstream(f) >> revcomp;
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
      if (!sc.addTag(bc, name.empty() ? bc : name, group, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude, trim_dir, trim_offset, after_str, before_str, partial5_min_match, partial5_mismatch_freq, partial3_min_match, partial3_mismatch_freq, subs_str, revcomp)) {
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
    if (ret && !opt.revcomp_str.empty() && ss16.good()) {
      std::cerr << ERROR_STR << " Number of values in --revcomp is greater than that in --tags" << std::endl;
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
  } else if (!opt.revcomp_str.empty()) {
    std::cerr << ERROR_STR << " --revcomp cannot be supplied unless --tags is" << std::endl;
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
    if (!opt.trim_only && !(opt.no_output || opt.no_output_)) {
      std::cerr << "Error: No tags found even though --assign specified" << std::endl;
      ret = false;
    }
  }
  
  if (ret) { // Validate number of output files supplied
    int nf = sc.getNFiles();
    int nf_u = sc.getNFilesUnassigned();
    if (opt.pipe) {
      if (opt.unassigned_files.size() != 0 && opt.unassigned_files.size() % nf_u != 0 || opt.unassigned_files.size() > nf_u) {
        std::cerr << ERROR_STR << " Incorrect number of --unassigned output files" << std::endl;
        ret = false;
      }
    }
    else {
      if (opt.output_files.size() % nf != 0 || opt.unassigned_files.size() % nf_u != 0 || opt.output_files.size() > nf || opt.unassigned_files.size() > nf_u) {
       std::cerr << ERROR_STR << " Incorrect number of output files" << std::endl;
       ret = false;
       }
    }
    opt.select_output_files.resize(nf, true);
  }
  
  return ret;
}


int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  SplitCode sc(opt.nfiles, opt.summary_file, opt.trim_only, opt.disable_n, opt.trim_5_str, opt.trim_3_str, opt.extract_str, opt.extract_no_chain, opt.barcode_prefix, opt.filter_length_str,
               opt.quality_trimming_5, opt.quality_trimming_3, opt.quality_trimming_pre, opt.quality_trimming_naive, opt.quality_trimming_threshold, opt.phred64, opt.write_locations, opt.sub_assign_vec, opt.bclen, opt.min_delta, !opt.summary_file.empty(),
               opt.x_only, opt.no_x_out, opt.outbam, opt.no_output_barcodes, opt.outputb_file, opt.remultiplex, opt.optimize_assignment_str, opt.from_header_str, opt.random_str, opt.show_not_found);
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
    if (use_gz && !opt.outbam) {
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

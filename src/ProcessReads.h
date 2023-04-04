#ifndef SPLITCODE_PROCESSREADS_H
#define SPLITCODE_PROCESSREADS_H

#include "SplitCode.h"

#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

#include "common.h"

#ifdef SPLITCODE_USE_ZLIB_NG
#ifndef WITH_GZFILEOP
#define WITH_GZFILEOP
#endif
#include "zlib-ng.h"
constexpr auto gzopen = zng_gzopen;
constexpr auto gzclose = zng_gzclose;
constexpr auto gzwrite = zng_gzwrite;
constexpr auto gzdopen = zng_gzdopen;
#else
#include <zlib.h>
#endif

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
#ifndef SPLITCODE_USE_ZLIB_NG
KSEQ_INIT(gzFile, gzread)
#else
KSEQ_INIT(gzFile, zng_gzread)
#endif
#endif
  
class MasterProcessor;

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);

class SequenceReader {
public:
  
  SequenceReader(const ProgramOptions& opt) :
  readbatch_id(-1) {};
  SequenceReader() : state(false), readbatch_id(-1) {};
  virtual ~SequenceReader() {}

  virtual bool empty() = 0;
  virtual void reset();
  virtual void reserveNfiles(int n) = 0;
  virtual bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                              std::vector<std::pair<const char*, int>>& names,
                              std::vector<std::pair<const char*, int>>& quals,
                              std::vector<uint32_t>& flags,
                              int &readbatch_id,
                              bool full=false,
                              bool comments=false) = 0;
  
  
public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:
  
  FastqSequenceReader(const ProgramOptions& opt) : SequenceReader(opt),
  current_file(0), files(opt.files) {
    SequenceReader::state = false;
    interleave_nfiles = opt.input_interleaved_nfiles;
    nfiles = opt.nfiles;
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : SequenceReader(), 
  current_file(0) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();
  
  bool empty();  
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      int &readbatch_id,
                      bool full=false,
                      bool comments=false);
  
public:
  int nfiles = 1;
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  std::vector<std::string> files;
  int current_file;
  std::vector<kseq_t*> seq;
  int interleave_nfiles;
};

class MasterProcessor {
public:
  MasterProcessor (SplitCode &sc, const ProgramOptions& opt)
    : sc(sc), opt(opt), numreads(0), bufsize(1ULL<<23), curr_readbatch_id(0) { 

    SR = new FastqSequenceReader(opt);
    verbose = opt.verbose;
    nfiles = opt.input_interleaved_nfiles == 0 ? opt.nfiles : opt.input_interleaved_nfiles;
    const std::string suffix = opt.output_fasta ? ".fasta" : ".fastq";
    const std::string suffix_gz = opt.output_fasta ? ".fasta.gz" : ".fastq.gz";
    std::string compress_level_str = "wb" + std::to_string(opt.compress_level);
    const char* gz_out_str = compress_level_str.c_str();
    for (auto f : opt.output_files) {
      if (opt.gzip) {
        out_gz.push_back(gzopen(f.c_str(), gz_out_str));
      } else {
        out.push_back(fopen(f.c_str(), "wb"));
      }
    }
    for (auto f : opt.unassigned_files) {
      if (opt.gzip) {
        outu_gz.push_back(gzopen(f.c_str(), gz_out_str));
      } else {
        outu.push_back(fopen(f.c_str(), "wb"));
      }
    }
    if (!opt.pipe && !opt.no_x_out && !opt.no_output) {
      for (auto f : sc.umi_names) {
        if (opt.gzip) {
          outumi_gz.push_back(gzopen((f+suffix_gz).c_str(), gz_out_str));
        } else {
          outumi.push_back(fopen((f+suffix).c_str(), "wb"));
        }
      }
    }
    write_barcode_separate_fastq = !opt.outputb_file.empty() && !opt.no_output_barcodes; // We write separate barcode file if this option is nonempty
    bool write_barcode_separate_fastq_keep = write_barcode_separate_fastq || (opt.pipe && !opt.no_output_barcodes); // For the files specified in the --keep options, we'll also write separate barcode file(s) if the user chose to --pipe
    for (auto f : sc.idmapinv_keep) {
      if (f.second.empty()) {
        continue;
      }
      if (opt.gzip) {
        auto it = out_keep_gz.find(f.second);
        if (it == out_keep_gz.end()) {
          for (int i = 0; i <= nfiles; i++) {
            if (i == 0 && !write_barcode_separate_fastq_keep) {
              out_keep_gz[f.second].push_back(nullptr);
              continue;
            }
            if (i != 0 && !opt.select_output_files[i-1]) {
              out_keep_gz[f.second].push_back(nullptr);
              continue;
            }
            out_keep_gz[f.second].push_back(gzopen((f.second + "_" + (i == 0 ? "barcodes" : std::to_string(i-1)) + suffix_gz).c_str(), gz_out_str));
          }
        }
      } else {
        auto it = out_keep.find(f.second);
        if (it == out_keep.end()) {
          for (int i = 0; i <= nfiles; i++) {
            if (i == 0 && !write_barcode_separate_fastq_keep) {
              out_keep[f.second].push_back(nullptr);
              continue;
            }
            if (i != 0 && !opt.select_output_files[i-1]) {
              out_keep[f.second].push_back(nullptr);
              continue;
            }
            out_keep[f.second].push_back(fopen((f.second + "_" + (i == 0 ? "barcodes" : std::to_string(i-1)) + suffix).c_str(), "wb"));
          }
        }
      }
    }
    for (auto f : sc.groupmapinv_keep) {
      if (f.second.empty()) {
        continue;
      }
      if (opt.gzip) {
        auto it = out_keep_gz.find(f.second);
        if (it == out_keep_gz.end()) {
          for (int i = 0; i <= nfiles; i++) {
            if (i == 0 && !write_barcode_separate_fastq_keep) {
              out_keep_gz[f.second].push_back(nullptr); // zeroth index is the barcodes file or nullptr if we're not writing a separate barcodes file
              continue;
            }
            if (i != 0 && !opt.select_output_files[i-1]) {
              out_keep_gz[f.second].push_back(nullptr);
              continue;
            }
            out_keep_gz[f.second].push_back(gzopen((f.second + "_" + (i == 0 ? "barcodes" : std::to_string(i-1)) + suffix_gz).c_str(), gz_out_str));
          }
        }
      } else {
        auto it = out_keep.find(f.second);
        if (it == out_keep.end()) {
          for (int i = 0; i <= nfiles; i++) {
            if (i == 0 && !write_barcode_separate_fastq_keep) {
              out_keep[f.second].push_back(nullptr);
              continue;
            }
            if (i != 0 && !opt.select_output_files[i-1]) {
              out_keep[f.second].push_back(nullptr);
              continue;
            }
            out_keep[f.second].push_back(fopen((f.second + "_" + (i == 0 ? "barcodes" : std::to_string(i-1)) + suffix).c_str(), "wb"));
          }
        }
      }
    }
    write_output_fastq = opt.output_fastq_specified || opt.pipe || opt.x_only;
    write_unassigned_fastq = outu.size() > 0 || outu_gz.size() > 0;
    if (write_barcode_separate_fastq) {
      if (opt.gzip) {
        outb_gz = gzopen(opt.outputb_file.c_str(), gz_out_str);
      } else {
        outb = fopen(opt.outputb_file.c_str(), "wb");
      }
    }
    
    if (opt.remultiplex) {
      std::unordered_map<std::string,int> batch_map;
      batch_id_mapping.resize(opt.batch_ids.size());
      int j = 0;
      for (int i = 0; i < opt.batch_ids.size(); i++) {
        if (batch_map.find(opt.batch_ids[i]) == batch_map.end()) {
          batch_id_mapping[i] = j;
          batch_map.insert(std::make_pair(opt.batch_ids[i], j));
          j++;
        } else {
          batch_id_mapping[i] = batch_map[opt.batch_ids[i]];
        }
      }
    }
  }
  
  ~MasterProcessor() {
    for (auto& of : out) {
      fclose(of);
    }
    for (auto& of : outu) {
      fclose(of);
    }
    for (auto& of : outumi) {
      fclose(of);
    }
    for (auto& of : out_gz) {
      gzclose(of);
    }
    for (auto& of : outu_gz) {
      gzclose(of);
    }
    for (auto& of : outumi_gz) {
      gzclose(of);
    }
    if (write_barcode_separate_fastq) {
      if (opt.gzip) {
        gzclose(outb_gz);
      } else {
        fclose(outb);
      }
    }
    for (auto e : out_keep) {
      for (auto f : e.second) {
        if (f != nullptr) {
          fclose(f);
        }
      }
    }
    for (auto e : out_keep_gz) {
      for (auto f : e.second) {
        if (f != nullptr) {
          gzclose(f);
        }
      }
    }
    delete SR;
  }
  
  std::mutex reader_lock;
  std::vector<std::mutex> parallel_reader_locks;
  bool parallel_read;
  std::mutex writer_lock;
  std::condition_variable cv;
  
  std::vector<FILE*> out;
  std::vector<gzFile> out_gz;
  FILE* outb;
  gzFile outb_gz;
  std::vector<FILE*> outu;
  std::vector<gzFile> outu_gz;
  std::unordered_map<std::string, std::vector<FILE*>> out_keep;
  std::unordered_map<std::string, std::vector<gzFile>> out_keep_gz;
  std::vector<FILE*> outumi;
  std::vector<gzFile> outumi_gz;
  std::vector<int> batch_id_mapping; // minimal perfect mapping of batch ID
  bool write_output_fastq;
  bool write_barcode_separate_fastq;
  bool write_unassigned_fastq;
  bool verbose;
  
  SequenceReader *SR;
  std::vector<FastqSequenceReader> FSRs;
  SplitCode& sc;

  const ProgramOptions& opt;
  int64_t numreads;
  size_t bufsize;
  int nfiles;
  int curr_readbatch_id;

  void processReads();
  void update(int n, std::vector<SplitCode::Results>& rv,
              std::vector<std::pair<const char*, int>>& seqs,
              std::vector<std::pair<const char*, int>>& names,
              std::vector<std::pair<const char*, int>>& quals,
              std::vector<uint32_t>& flags,
              int readbatch_id);  
  void writeOutput(std::vector<SplitCode::Results>& rv,
                   std::vector<std::pair<const char*, int>>& seqs,
                   std::vector<std::pair<const char*, int>>& names,
                   std::vector<std::pair<const char*, int>>& quals,
                   std::vector<uint32_t>& flags);
};

class ReadProcessor {
public:
  ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  
  size_t bufsize;
  MasterProcessor& mp;
  int64_t numreads;
  
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  
  std::vector<SplitCode::Results> rv;
  bool full;
  bool comments;
  
  /*std::vector<std::vector<int>> newIDs;
  std::vector<std::vector<int>> IDs;*/
  
  void operator()();
  void processBuffer();
  void clear();
};

std::string pretty_num(size_t num);

#endif // SPLITCODE_PROCESSREADS_H

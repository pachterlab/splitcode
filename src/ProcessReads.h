#ifndef SPLITCODE_PROCESSREADS_H
#define SPLITCODE_PROCESSREADS_H

#include "SplitCode.h"

#include <zlib.h>
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


#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
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
                              bool full=false) = 0;
  
  
public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:
  
  FastqSequenceReader(const ProgramOptions& opt) : SequenceReader(opt),
  current_file(0), files(opt.files) {
    SequenceReader::state = false;
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
                      bool full=false);
  
public:
  int nfiles = 1;
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  std::vector<std::string> files;
  int current_file;
  std::vector<kseq_t*> seq;
};

class MasterProcessor {
public:
  MasterProcessor (SplitCode &sc, const ProgramOptions& opt)
    : sc(sc), opt(opt), numreads(0), bufsize(1ULL<<23) { 

    SR = new FastqSequenceReader(opt);
    verbose = opt.verbose;
    for (auto f : opt.output_files) {
      if (opt.gzip) {
        out_gz.push_back(gzopen(f.c_str(), "wb1"));
      } else {
        out.push_back(fopen(f.c_str(), "wb"));
      }
    }
    for (auto f : opt.unassigned_files) {
      if (opt.gzip) {
        outu_gz.push_back(gzopen(f.c_str(), "wb1"));
      } else {
        outu.push_back(fopen(f.c_str(), "wb"));
      }
    }
    write_output_fastq = opt.output_fastq_specified || opt.pipe;
    write_unassigned_fastq = outu.size() > 0 || outu_gz.size() > 0;
    write_barcode_separate_fastq = !opt.outputb_file.empty();
    if (write_barcode_separate_fastq) {
      if (opt.gzip) {
        outb_gz = gzopen(opt.outputb_file.c_str(), "wb1");
      } else {
        outb = fopen(opt.outputb_file.c_str(), "wb");
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
    for (auto& of : out_gz) {
      gzclose(of);
    }
    for (auto& of : outu_gz) {
      gzclose(of);
    }
    if (write_barcode_separate_fastq) {
      if (opt.gzip) {
        gzclose(outb_gz);
      } else {
        fclose(outb);
      }
    }
    delete SR;
  }
  
  std::mutex reader_lock;
  std::vector<std::mutex> parallel_reader_locks;
  bool parallel_read;
  std::mutex writer_lock;
  
  std::vector<FILE*> out;
  std::vector<gzFile> out_gz;
  FILE* outb;
  gzFile outb_gz;
  std::vector<FILE*> outu;
  std::vector<gzFile> outu_gz;
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

  void processReads();
  void update(int n, std::vector<SplitCode::Results>& rv,
              std::vector<std::pair<const char*, int>>& seqs,
              std::vector<std::pair<const char*, int>>& names,
              std::vector<std::pair<const char*, int>>& quals);  
  void writeOutput(std::vector<SplitCode::Results>& rv,
                   std::vector<std::pair<const char*, int>>& seqs,
                   std::vector<std::pair<const char*, int>>& names,
                   std::vector<std::pair<const char*, int>>& quals);
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
  
  /*std::vector<std::vector<int>> newIDs;
  std::vector<std::vector<int>> IDs;*/
  
  void operator()();
  void processBuffer();
  void clear();
};

std::string pretty_num(size_t num);

#endif // SPLITCODE_PROCESSREADS_H

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

    }
  
  ~MasterProcessor() {
    delete SR;
  }
  
  std::mutex reader_lock;
  std::vector<std::mutex> parallel_reader_locks;
  bool parallel_read;
  std::mutex writer_lock;
  
  
  SequenceReader *SR;
  std::vector<FastqSequenceReader> FSRs;
  SplitCode& sc;

  const ProgramOptions& opt;
  int64_t numreads;
  size_t bufsize;

  void processReads();
  void update(int n);  
};

class ReadProcessor {
public:
  ReadProcessor(const SplitCode& sc, const ProgramOptions& opt, MasterProcessor& mp);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  
  size_t bufsize;
  const SplitCode& sc;
  MasterProcessor& mp;
  int64_t numreads;
  
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  
  std::vector<std::vector<int>> newIDs;
  std::vector<std::vector<int>> IDs;
  
  void operator()();
  void processBuffer();
  void clear();
};


#endif // SPLITCODE_PROCESSREADS_H

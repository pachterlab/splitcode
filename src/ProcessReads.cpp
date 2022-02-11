#include <fstream>
#include <limits>
#include <iomanip>
#include "ProcessReads.h"
#include "kseq.h"
#include "common.h"
#include <unordered_set>

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");
  
  if (s.size() <= 3) {
    return s;
  }
  
  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }
  
  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }
  
  ret += s.substr(start_pos, 3);
  
  return ret;
}


//methods


int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt) {

  // need to receive an index map
  std::ios_base::sync_with_stdio(false);

  size_t numreads = 0;
  
  if (MP.verbose) {

  for (int i = 0, si=1; i < opt.files.size(); si++) {
    std::cerr << "* will process sample " << si<< ": ";
    for (int j = 0; j < opt.nfiles; j++,i++) {
      if (j>0) {
        std::cerr << "                         ";
      }
      std::cerr << opt.files[i] << std::endl;
    }
  }

  // for each file
  std::cerr << "* processing the reads ..."; std::cerr.flush();
  }


  MP.processReads();
  numreads = MP.numreads;
  if (MP.verbose) {
  std::cerr << std::endl << "done " << std::endl;
  }
  
  int nummapped = MP.sc.getNumMapped();

  if (MP.verbose) {
    std::cerr << "* processed " << pretty_num(numreads) << " reads";
    if (!MP.sc.always_assign) {
      std::cerr << ", " << pretty_num(nummapped) << " reads had identifiable barcodes";
    }
    std::cerr << std::endl;
  }
  
  return numreads;
}


/** -- read processors -- **/

void MasterProcessor::processReads() {
  

  // start worker threads
  
  std::vector<std::thread> workers;
  parallel_read = opt.threads > 4 && opt.files.size() > opt.nfiles;
  if (parallel_read) {
    delete SR;
    SR = nullptr;
    assert(opt.files.size() % opt.nfiles == 0);
    int nbatches = opt.files.size() / opt.nfiles;
    std::vector<std::mutex> mutexes(nbatches);
    parallel_reader_locks.swap(mutexes);
    for (int i = 0; i < nbatches; i++) {
      FastqSequenceReader fSR(opt);
      fSR.files.erase(fSR.files.begin(), fSR.files.begin()+opt.nfiles*i);
      fSR.files.erase(fSR.files.begin()+opt.nfiles, fSR.files.end());
      assert(fSR.files.size() == opt.nfiles);
      FSRs.push_back(std::move(fSR));
    }
  }
  
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(ReadProcessor(opt,*this)));
  }
  
  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }
  
  /*// now handle the modification of the mincollector
  for (int i = 0; i < bus_ecmap.size(); i++) {
    auto &u = bus_ecmap[i];
    int ec = index.ecmapinv.size();
    auto it = bus_ecmapinv.find(u);
    if (it->second != ec) {
      std::cout << "Error" << std::endl;
      exit(1);
    }      
    index.ecmapinv.insert({u,ec});
    index.ecmap.push_back(u);
  }

  busf_out.close();*/
}

void MasterProcessor::update(int n, std::vector<SplitCode::Results>& rv,
                             std::vector<std::pair<const char*, int>>& seqs,
                             std::vector<std::pair<const char*, int>>& names,
                             std::vector<std::pair<const char*, int>>& quals) {
  // acquire the writer lock
  std::lock_guard<std::mutex> lock(this->writer_lock);

  /*for (int i = 0; i < c.size(); i++) {
    tc.counts[i] += c[i]; // add up ec counts
    nummapped += c[i];
  }

  if (!opt.batch_mode || opt.batch_bus) {
    for(auto &u : newEcs) {
      ++newECcount[u];
    }
  }
  nummapped += newEcs.size();

  if (opt.bus_mode || opt.batch_bus) {
    int bus_bc_sum = 0;
    int bus_umi_sum = 0;
    for (int i = 0; i <= 32; i++) {
      bus_bc_sum += bus_bc_len[i];
      bus_umi_sum += bus_umi_len[i];
    }

    if (bus_bc_sum < 10000 or bus_umi_sum < 10000) {
      for (int i = 0; i < 32; i++) {
        bus_bc_len[i] += bc_len[i];
        bus_umi_len[i] += umi_len[i];
      }
    }

    // add new equiv classes to extra format
    int offset = index.ecmapinv.size();
    for (auto &bp : newBP) {
      auto& u = bp.second;
      int ec = -1;
      auto it = bus_ecmapinv.find(u);
      if (it != bus_ecmapinv.end()) {
        ec = it->second;
      } else {
        ec = offset + bus_ecmapinv.size();
        bus_ecmapinv.insert({u,ec});
        bus_ecmap.push_back(u);
      }
      auto &b = bp.first;
      b.ec = ec;
      bv.push_back(b);
    }

    //copy bus mode information, write to disk or queue up
    writeBUSData(busf_out, bv); 
    //for (auto &bp : newBP) {
    //  newB.push_back(std::move(bp));
    //}
  }*/
  
  sc.update(rv);
  
  if (write_output_fastq) {
    writeOutput(rv, seqs, names, quals);
  }

  numreads += n;
  // releases the lock
}

void MasterProcessor::writeOutput(std::vector<SplitCode::Results>& rv,
                                  std::vector<std::pair<const char*, int>>& seqs,
                                  std::vector<std::pair<const char*, int>>& names,
                                  std::vector<std::pair<const char*, int>>& quals) {
  // Write out fastq
  int incf, jmax;
  incf = opt.nfiles-1;
  jmax = opt.nfiles;
  
  std::vector<const char*> s(jmax, nullptr);
  std::vector<const char*> n(jmax, nullptr);
  std::vector<const char*> nl(jmax, nullptr);
  std::vector<const char*> q(jmax, nullptr);
  std::vector<int> l(jmax,0);
  
  int readnum = 0;
  for (int i = 0; i + incf < seqs.size(); i++, readnum++) {
    auto& r = rv[readnum];
    bool assigned = sc.isAssigned(r);
    std::string mod_name = "";
    if (assigned && opt.mod_names) {
      mod_name = "::" + sc.getNameString(r); // Barcode names
    }
    if (assigned && (write_barcode_separate_fastq || opt.pipe) && !sc.always_assign) { // Write out barcode read
      std::stringstream o;
      // Write out barcode read
      o << "@" << std::string(names[i].first, names[i].second) << mod_name << "\n";
      o << sc.binaryToString(r.id, sc.FAKE_BARCODE_LEN) << "\n";
      o << "+" << "\n";
      o << std::string(sc.FAKE_BARCODE_LEN, sc.QUAL) << "\n";
      const std::string& ostr = o.str();
      size_t ostr_len = ostr.length();
      if (opt.gzip && !opt.pipe) {
        gzwrite(outb_gz, ostr.c_str(), ostr_len);
      } else {
        fwrite(ostr.c_str(), 1, ostr_len, opt.pipe ? stdout : outb);
      }
    }
    for (int j = 0; j < jmax; j++) {
      if (!assigned && !write_unassigned_fastq) {
        break;
      }
      std::stringstream o;
      const char* s = seqs[i+j].first;
      int l = seqs[i+j].second;
      const char* n = names[i+j].first;
      int nl = names[i+j].second;
      const char* q = quals[i+j].first;
      // Write out read
      bool embed_final_barcode = assigned && j==0 && !write_barcode_separate_fastq && !opt.pipe && !sc.always_assign;
      o << "@";
      o << std::string(n,nl) << mod_name << "\n";
      if (embed_final_barcode) {
        o << sc.binaryToString(r.id, sc.FAKE_BARCODE_LEN);
      } else if (l == 0 && !opt.empty_read_sequence.empty()) {
        o << opt.empty_read_sequence;
      }
      o << std::string(s,l) << "\n";
      o << "+" << "\n";
      if (embed_final_barcode) {
        o << std::string(sc.FAKE_BARCODE_LEN, sc.QUAL);
      } else if (l == 0 && !opt.empty_read_sequence.empty()) {
        o << std::string(opt.empty_read_sequence.length(), sc.QUAL);
      }
      o << std::string(q,l) << "\n";
      
      const std::string& ostr = o.str();
      size_t ostr_len = ostr.length();
      if (opt.pipe && assigned) {
        fwrite(ostr.c_str(), 1, ostr_len, stdout);
      } else if (!assigned) {
        if (opt.gzip) {
          gzwrite(outu_gz[j], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, outu[j]);
        }
      } else {
        if (opt.gzip) {
          gzwrite(out_gz[j], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, opt.pipe ? stdout : out[j]);
        }
      }
    }
    i += incf;
  }
}

/*
void MasterProcessor::outputFusion(const std::stringstream &o) {
  std::string os = o.str();
  if (!os.empty()) {
    std::lock_guard<std::mutex> lock(this->writer_lock);
    ofusion << os << "\n";
  }
}
*/

ReadProcessor::ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp) : 
  mp(mp), numreads(0) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];
   seqs.reserve(bufsize/50);
   rv.reserve(1000);
   clear();
}

ReadProcessor::ReadProcessor(ReadProcessor && o) :
  bufsize(o.bufsize),
  mp(o.mp),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)) {
    buffer = o.buffer;
    o.buffer = nullptr;
    o.bufsize = 0;
}

ReadProcessor::~ReadProcessor() {
  if (buffer != nullptr) {
      delete[] buffer;
      buffer = nullptr;
  }
}

void ReadProcessor::operator()() {
  uint64_t parallel_read_counter = 0;
  std::unordered_set<int> parallel_read_empty;
  while (true) {
    int readbatch_id;
    std::vector<std::string> umis;
    // grab the reader lock
    if (mp.parallel_read) {
      int nbatches = mp.opt.files.size() / mp.opt.nfiles;
      int i = parallel_read_counter % nbatches;
      if (parallel_read_empty.size() >= nbatches) {
        return;
      }
      parallel_read_counter++;
      std::lock_guard<std::mutex> lock(mp.parallel_reader_locks[i]);
      if (mp.FSRs[i].empty()) {
        parallel_read_empty.emplace(i);
        continue;
      }
      mp.FSRs[i].fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, !mp.opt.no_output);
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, !mp.opt.no_output);
      }
      // release the reader lock
    }
    
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    mp.update(seqs.size() / mp.opt.nfiles, rv, seqs, names, quals);
    clear();
  }
}

void ReadProcessor::processBuffer() {
  // actually process the sequence
  
  int incf, jmax;
  incf = mp.opt.nfiles-1;
  jmax = mp.opt.nfiles;

  std::vector<const char*> s(jmax, nullptr);
  std::vector<int> l(jmax,0);


  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      s[j] = seqs[i+j].first;
      l[j] = seqs[i+j].second;      
    }
    i += incf;
    numreads++;
    
    SplitCode::Results results;
    mp.sc.processRead(s, l, jmax, results);
    rv.push_back(results);
    
    if (mp.sc.isAssigned(results)) { // Only modify/trim the reads stored in seq if assigned
      mp.sc.modifyRead(seqs, quals, i-incf, results);
    }

    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) { 
        numreads = 0; // reset counter
        int nummapped = mp.sc.getNumMapped();

        std::cerr << '\r' << (mp.numreads/1000000) << "M reads processed";
        if (!mp.sc.always_assign) {
          std::cerr << " (" 
            << std::fixed << std::setw( 3 ) << std::setprecision( 1 ) << ((100.0*nummapped)/double(mp.numreads))
            << "% identified)";
        } else {
          std::cerr << " (running in --trim-only mode)";
        }
        std::cerr.flush();
      }
  }
}

void ReadProcessor::clear() {
  memset(buffer,0,bufsize);
  rv.clear();
}

/** -- sequence readers -- **/

void SequenceReader::reset() {
  state = false;
  readbatch_id = -1;
}

FastqSequenceReader::~FastqSequenceReader() {
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
  }

  for (auto &s : seq) {
    kseq_destroy(s);
  }
}


bool FastqSequenceReader::empty() {
  return (!state && current_file >= files.size());
}

void FastqSequenceReader::reset() {
  SequenceReader::reset();
   
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
    f = nullptr;
  }

  for (auto &ll : l) {
    ll = 0;
  }
  for (auto &nll : nl) {
    nll = 0;
  }
  
  current_file = 0;
  for (auto &s : seq) {
    kseq_destroy(s);
    s = nullptr;
  }
}

void FastqSequenceReader::reserveNfiles(int n) {
  fp.resize(nfiles);
  l.resize(nfiles, 0);
  nl.resize(nfiles, 0);
  seq.resize(nfiles, nullptr);
}

// returns true if there is more left to read from the files
bool FastqSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
  std::vector<std::pair<const char *, int> > &names,
  std::vector<std::pair<const char *, int> > &quals,
  std::vector<uint32_t>& flags,
  int& read_id,
  bool full) {
    
  std::string line;
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
  seqs.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
  flags.clear();
  
  int bufpos = 0;
  int pad = nfiles;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current files
        for (auto &f : fp) {
          if (f) {
            gzclose(f);
          }
        }
        
        // open the next one
        for (int i = 0; i < nfiles; i++) {
          fp[i] = gzopen(files[current_file+i].c_str(), "r");
          seq[i] = kseq_init(fp[i]);
          l[i] = kseq_read(seq[i]);
          
        }
        current_file+=nfiles;
        state = true; 
      }
    }
    // the file is open and we have read into seq1 and seq2
    bool all_l = true;
    int bufadd = nfiles;
    for (int i = 0; i < nfiles; i++) {
      all_l = all_l && l[i] >= 0;
      bufadd += l[i]; // includes seq
    }
    if (all_l) {      
      // fits into the buffer
      if (full) {
        for (int i = 0; i < nfiles; i++) {
          nl[i] = seq[i]->name.l;
          bufadd += l[i] + nl[i]; // includes name and qual
        }
        bufadd += 2*pad;
      }

      if (bufpos+bufadd< limit) {

        for (int i = 0; i < nfiles; i++) {
          char *pi = buf + bufpos;
          memcpy(pi, seq[i]->seq.s, l[i]+1);
          bufpos += l[i]+1;
          seqs.emplace_back(pi,l[i]);

          if (full) {
            pi = buf + bufpos;
            memcpy(pi, seq[i]->qual.s,l[i]+1);
            bufpos += l[i]+1;
            quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          }
        }

        numreads++;
        flags.push_back(numreads-1);
      } else {
        return true; // read it next time
      }

      // read for the next one
      for (int i = 0; i < nfiles; i++) {
        l[i] = kseq_read(seq[i]);
      }        
    } else {
      state = false; // haven't opened file yet
    }
  }
}

// move constructor

FastqSequenceReader::FastqSequenceReader(FastqSequenceReader&& o) :
  nfiles(o.nfiles),
  numreads(o.numreads),
  fp(std::move(o.fp)),
  l(std::move(o.l)),
  nl(std::move(o.nl)),
  files(std::move(o.files)),
  current_file(o.current_file),
  seq(std::move(o.seq)) {

  o.fp.resize(nfiles);
  o.l.resize(nfiles, 0);
  o.nl.resize(nfiles, 0);
  o.seq.resize(nfiles, nullptr);
  o.state = false;
}



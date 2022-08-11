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
      std::cerr << ", " << pretty_num(nummapped) << " reads were assigned";
    }
    std::cerr << std::endl;
  }
  
  MP.sc.setNumReads(numreads);
  
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
    //assert(opt.files.size() % opt.nfiles == 0);
    int nbatches = opt.files.size() / opt.nfiles;
    std::vector<std::mutex> mutexes(nbatches);
    parallel_reader_locks.swap(mutexes);
    for (int i = 0; i < nbatches; i++) {
      FastqSequenceReader fSR(opt);
      fSR.files.erase(fSR.files.begin(), fSR.files.begin()+opt.nfiles*i);
      fSR.files.erase(fSR.files.begin()+opt.nfiles, fSR.files.end());
      //assert(fSR.files.size() == opt.nfiles);
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
}

void MasterProcessor::update(int n, std::vector<SplitCode::Results>& rv,
                             std::vector<std::pair<const char*, int>>& seqs,
                             std::vector<std::pair<const char*, int>>& names,
                             std::vector<std::pair<const char*, int>>& quals) {
  // acquire the writer lock
  std::lock_guard<std::mutex> lock(this->writer_lock);
  
  if (opt.max_num_reads != 0 && numreads+n > opt.max_num_reads) {
    n = opt.max_num_reads-numreads;
    rv.resize(n);
  }
  
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
  incf = nfiles-1;
  jmax = nfiles;
  
  std::vector<const char*> s(jmax, nullptr);
  std::vector<const char*> n(jmax, nullptr);
  std::vector<const char*> nl(jmax, nullptr);
  std::vector<const char*> q(jmax, nullptr);
  std::vector<int> l(jmax,0);
  char start_char = opt.output_fasta ? '>' : '@';
  bool include_quals = !opt.output_fasta;
  
  int readnum = 0;
  for (int i = 0; i + incf < seqs.size() && readnum < rv.size(); i++, readnum++) {
    auto& r = rv[readnum];
    if (!r.passes_filter) {
      i += incf;
      continue;
    }
    if (sc.always_assign) {
      r.ofile = ""; // If always-assign, don't allow writing output to alternative files based on identified tags
    }
    auto& umi_vec = r.umi_data;
    bool assigned = sc.isAssigned(r);
    bool use_pipe = opt.pipe && r.ofile.empty(); // Conditions under which we'll write to stdout
    // Conditions under which we'll write to separate barcode file (either write_barcode_separate_fastq specified previously or we need to write reads out to r.ofile even though user specified --pipe):
    bool write_barcode_separate_fastq_ = write_barcode_separate_fastq || (!r.ofile.empty() && (opt.pipe && !opt.no_output_barcodes));
    std::string mod_name = "";
    if ((assigned || r.discard) && opt.mod_names) {
      mod_name = "::" + sc.getNameString(r); // Barcode names
    }
    bool name_modded = false;
    if ((assigned || r.discard) && opt.seq_names && !r.identified_tags_seqs.empty()) {
      mod_name += " CB:Z:" + r.identified_tags_seqs; // Sequences of identified tags stitched together
      name_modded = true;
    }
    if (assigned && opt.com_names && !sc.always_assign) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += "BI:i:" + std::to_string(sc.getID(r.id));
      //mod_name += "\t" + "CB:Z:" + sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength())
      name_modded = true;
    }
    if (assigned && opt.x_names && !sc.umi_names.empty()) {
      std::string mod_name2 = (name_modded ? "\t" : " ");
      mod_name2 += "RX:Z:";
      bool umi_empty = true;
      for (int umi_index = 0; umi_index < sc.umi_names.size(); umi_index++) { // Iterate through vector of all UMI names
        std::string curr_umi = umi_vec[umi_index];
        if (umi_index != 0 && !(opt.empty_remove && curr_umi.empty())) {
          mod_name2 += "-";
        }
        mod_name2 += curr_umi.empty() ? opt.empty_read_sequence : curr_umi;
        if (!curr_umi.empty()) {
          umi_empty = false;
        }
      }
      if (!umi_empty) {
        mod_name += mod_name2;
      }
    }
    if (assigned && (write_barcode_separate_fastq_ || use_pipe) && !sc.always_assign && !opt.no_output_barcodes) { // Write out barcode read
      std::stringstream o;
      // Write out barcode read
      o << start_char << std::string(names[i].first, names[i].second) << mod_name << "\n";
      o << sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength()) << "\n";
      if (include_quals) {
        o << "+" << "\n";
        o << std::string(sc.getBarcodeLength(), sc.QUAL) << "\n";
      }
      const std::string& ostr = o.str();
      size_t ostr_len = ostr.length();
      if (use_pipe && !write_barcode_separate_fastq_) {
        fwrite(ostr.c_str(), 1, ostr_len, stdout);
      } else if (opt.gzip) {
        gzwrite(r.ofile.empty() ? outb_gz : out_keep_gz[r.ofile][0], ostr.c_str(), ostr_len);
      } else {
        fwrite(ostr.c_str(), 1, ostr_len, r.ofile.empty() ? outb : out_keep[r.ofile][0]);
      }
    }
    if (!sc.umi_names.empty() && assigned && !opt.no_x_out) { // Write out extracted UMIs as needed
      for (int umi_index = 0; umi_index < sc.umi_names.size(); umi_index++) { // Iterate through vector of all UMI names
        std::string curr_umi = umi_vec[umi_index];
        if (curr_umi.empty()) {
          if (opt.empty_remove) {
            continue; // Don't write anything
          }
          curr_umi = opt.empty_read_sequence;
        }
        std::stringstream o;
        o << start_char << std::string(names[i].first, names[i].second) << mod_name << "\n";
        o << curr_umi << "\n";
        if (include_quals) {
          o << "+" << "\n";
          o << std::string(curr_umi.length(), sc.QUAL) << "\n";
        }
        const std::string& ostr = o.str();
        size_t ostr_len = ostr.length();
        if (use_pipe) {
          fwrite(ostr.c_str(), 1, ostr_len, stdout);
        } else if (opt.gzip) {
          gzwrite(outumi_gz[umi_index], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, outumi[umi_index]);
        }
      }
    }
    for (int j = 0; j < jmax; j++) {
      if (!assigned && !write_unassigned_fastq) {
        break;
      }
      if (opt.x_only) {
        break;
      }
      std::stringstream o;
      const char* s = seqs[i+j].first;
      int l = seqs[i+j].second;
      const char* n = names[i+j].first;
      int nl = names[i+j].second;
      const char* q = quals[i+j].first;
      // Write out read
      bool embed_final_barcode = assigned && j==0 && !write_barcode_separate_fastq_ && !use_pipe && !sc.always_assign && !opt.no_output_barcodes;
      o << start_char;
      o << std::string(n,nl) << mod_name << "\n";
      if (embed_final_barcode) {
        o << sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength());
      } else if (l == 0 && !opt.empty_read_sequence.empty()) {
        o << opt.empty_read_sequence;
      } else if (l == 0 && opt.empty_remove) {
        continue; // Don't write anything
      }
      o << std::string(s,l) << "\n";
      if (include_quals) {
        o << "+" << "\n";
        if (embed_final_barcode) {
          o << std::string(sc.getBarcodeLength(), sc.QUAL);
        } else if (l == 0 && !opt.empty_read_sequence.empty()) {
          o << std::string(opt.empty_read_sequence.length(), sc.QUAL);
        }
        o << std::string(q,l) << "\n";
      }
      
      const std::string& ostr = o.str();
      size_t ostr_len = ostr.length();
      if (assigned) {
        if (use_pipe) {
          fwrite(ostr.c_str(), 1, ostr_len, stdout);
        } else {
          if (opt.gzip) {
            gzwrite(r.ofile.empty() ? out_gz[j] : out_keep_gz[r.ofile][j+1], ostr.c_str(), ostr_len);
          } else {
            // note: we use j+1 for out_keep and out_keep_gz because the zeroth index is the barcodes file
            fwrite(ostr.c_str(), 1, ostr_len, r.ofile.empty() ? out[j] : out_keep[r.ofile][j+1]);
          }
        }
      } else {
        if (opt.gzip) {
          gzwrite(outu_gz[j], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, outu[j]);
        }
      }
    }
    i += incf;
  }
}

ReadProcessor::ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp) : 
  mp(mp), numreads(0) {
   // initialize buffer
   bufsize = mp.bufsize;
   buffer = new char[bufsize];
   seqs.reserve(bufsize/50);
   full = !mp.opt.no_output;
   if (full) {
     quals.reserve(bufsize/50);
   }
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
  flags(std::move(o.flags)),
  full(o.full) {
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
      // assert(mp.opt.input_interleaved_nfiles == 0);
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
      mp.FSRs[i].fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full);
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full);
      }
      // release the reader lock
    }
    
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    int nfiles = mp.nfiles;
    mp.update(seqs.size() / nfiles, rv, seqs, names, quals);
    clear();
    if (mp.opt.max_num_reads != 0 && mp.numreads >= mp.opt.max_num_reads) {
      return;
    }
  }
}

void ReadProcessor::processBuffer() {
  // actually process the sequence
  
  int incf, jmax, nfiles;
  nfiles = mp.nfiles;
  incf = nfiles-1;
  jmax = nfiles;

  std::vector<const char*> s(jmax, nullptr);
  std::vector<int> l(jmax,0);
  std::vector<const char*> q(full ? jmax : 0, nullptr);


  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      s[j] = seqs[i+j].first;
      l[j] = seqs[i+j].second;      
      if (full) {
        q[j] = quals[i+j].first;
      }
    }
    i += incf;
    numreads++;
    
    SplitCode::Results results;
    mp.sc.processRead(s, l, jmax, results, q);
    if (mp.sc.isAssigned(results)) { // Only modify/trim the reads stored in seq if assigned
      mp.sc.modifyRead(seqs, quals, i-incf, results);
    }
    rv.push_back(results);

    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) { 
        numreads = 0; // reset counter
        int nummapped = mp.sc.getNumMapped();

        std::cerr << '\r' << (mp.numreads/1000000) << "M reads processed";
        if (!mp.sc.always_assign) {
          std::cerr << " (" 
            << std::fixed << std::setw( 3 ) << std::setprecision( 1 ) << ((100.0*nummapped)/double(mp.numreads))
            << "% assigned)";
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
  int count = 0; // for interleaving
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
          fp[i] = files[0] == "-" && nfiles == 1 && files.size() == 1 ? gzdopen(fileno(stdin), "r") : gzopen(files[current_file+i].c_str(), "r");
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
        if (interleave_nfiles != 0) { // Hack to allow interleaving
          // assert(nfiles == 1);
          if (bufpos+bufadd >= limit-262144 && count % interleave_nfiles == 0) {
            return true;
          }
          count++;
        }

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
        if (interleave_nfiles != 0) {
          std::cerr << "Error: There was an error processing interleaved FASTQ input. Exiting..." << std::endl;
          exit(1);
        }
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



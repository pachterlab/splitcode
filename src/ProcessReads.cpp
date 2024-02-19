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
  
  int64_t nummapped = MP.sc.getNumMapped();

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
  parallel_read = false; // We'll just always turn parallel_read off so FASTQs are processed in a deterministic and sequential order (at the expense of speed saturating)
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
  
  if (!opt.webasm) {
    for (int i = 0; i < opt.threads; i++) {
      workers.emplace_back(std::thread(ReadProcessor(opt,*this)));
    }
    
    // let the workers do their thing
    for (int i = 0; i < opt.threads; i++) {
      workers[i].join(); //wait for them to finish
    }
  } else { // webasm doesn't seem to work with threads
    ReadProcessor rp(opt,*this);
    rp();
  }
}

void MasterProcessor::update(int n, std::vector<SplitCode::Results>& rv,
                             std::vector<std::pair<const char*, int>>& seqs,
                             std::vector<std::pair<const char*, int>>& names,
                             std::vector<std::pair<const char*, int>>& quals,
                             std::vector<uint32_t>& flags,
                             int readbatch_id) {
  // acquire the writer lock
  std::unique_lock<std::mutex> lock(this->writer_lock);
  
  if (opt.max_num_reads != 0 && numreads+n > opt.max_num_reads) {
    n = opt.max_num_reads-numreads;
    rv.resize(n);
  }
  
  sc.update(rv);
  
  if (write_output_fastq) {
    while (readbatch_id != curr_readbatch_id && !parallel_read) {
      cv.wait(lock, [this, readbatch_id]{ return readbatch_id == curr_readbatch_id; });
    }
    writeOutput(rv, seqs, names, quals, flags);
  }

  numreads += n;
  curr_readbatch_id++;
  lock.unlock(); // releases the lock
  cv.notify_all(); // Alert all other threads to check their readbatch_id's!
}

void MasterProcessor::writeOutput(std::vector<SplitCode::Results>& rv,
                                  std::vector<std::pair<const char*, int>>& seqs,
                                  std::vector<std::pair<const char*, int>>& names,
                                  std::vector<std::pair<const char*, int>>& quals,
                                  std::vector<uint32_t>& flags) {
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
  
  size_t readnum = 0;
  for (int i = 0; i + incf < seqs.size() && readnum < rv.size(); i++, readnum++) {
    auto& r = rv[readnum];
    if (!r.passes_filter) {
      i += incf;
      continue;
    }
    auto& umi_vec = r.umi_data;
    bool assigned = sc.isAssigned(r); // Note: r.discard and assigned will both true be in the case of sc.always_assign==true but the read doesn't pass our keep/discard filter (if !sc.always_assign, assigned will be false if r.discard is false)
    bool assigned2 = sc.isAssigned(r, true); // Unlike assigned, assigned2 is false if sc.always_assign==true but the read doesn't pass the keep/discard filter; basically, it's equivalent to: (assigned && !r.discard)
    bool use_pipe = opt.pipe && (r.ofile.empty() || (!r.ofile.empty() && !r.ofile_keep)); // Conditions under which we'll write to stdout
    // Conditions under which we'll write to separate barcode file (either write_barcode_separate_fastq specified previously or we need to write reads out to r.ofile even though user specified --pipe):
    bool write_barcode_separate_fastq_ = write_barcode_separate_fastq || (!r.ofile.empty() && (opt.pipe && !opt.no_output_barcodes));
    bool name_modded = false;
    std::string mod_name = "";
    if (opt.keep_fastq_comments) {
      name_modded = true; // In this option, we always start the fastq comment with a tab
    }
    if ((assigned || r.discard) && opt.mod_names) {
      mod_name += (name_modded ? "\t" : "");
      mod_name += "::" + sc.getNameString(r); // Barcode names
    }
    if ((assigned || r.discard) && opt.seq_names && !r.identified_tags_seqs.empty()) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[0][0] + r.identified_tags_seqs; // Sequences of identified tags stitched together
      name_modded = true;
    }
    if (assigned && opt.com_names && !sc.always_assign) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[2][0] + std::to_string(sc.getID(r.id));
      //mod_name += "\t" + opt.sam_tags[0][0] + sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength())
      name_modded = true;
    } else if (assigned && opt.com_names && opt.remultiplex) { // Add remultiplexed ID to read name
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[2][0] + std::to_string(sc.getID(batch_id_mapping[flags[readnum]]));
      name_modded = true;
    }
    if (assigned && opt.bc_names && !sc.always_assign) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[4][0] + sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength());
      name_modded = true;
    } else if (assigned && opt.bc_names && opt.remultiplex) { // Add remultiplexed ID to read name
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[4][0] + sc.binaryToString(sc.getID(batch_id_mapping[flags[readnum]]), sc.getBarcodeLength());
      name_modded = true;
    }
    if (assigned && r.subassign_id != -1) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[3][0] + std::to_string(r.subassign_id);
      name_modded = true;
    }
    if (opt.write_locations && !r.tag_locations.empty()) {
      mod_name += (name_modded ? "\t" : " ");
      mod_name += opt.sam_tags[5][0];
      for (int ri = 0; ri < r.tag_locations.size(); ri++) {
        mod_name +=  + (ri == 0 ? "" : ",") + r.tag_locations[ri];
      }
      name_modded = true;
    }
    if (assigned2 && opt.x_names && !sc.umi_names.empty()) {
      std::string mod_name2 = (name_modded ? "\t" : " ");
      mod_name2 += opt.sam_tags[1][0];
      bool umi_empty = true;
      for (size_t umi_index = 0; umi_index < sc.umi_names.size(); umi_index++) { // Iterate through vector of all UMI names
        std::string curr_umi = umi_vec[umi_index];
        if (umi_index != 0 && !(opt.empty_remove && curr_umi.empty())) {
          if (opt.sam_tags[1].size() >= umi_index+1) {
            mod_name2 += "\t" + opt.sam_tags[1][umi_index];
          } else {
            mod_name2 += "-";
          }
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
    if (mod_name == " " || mod_name == "\t") {
      mod_name = "";
    }
    if (assigned && (write_barcode_separate_fastq_ || use_pipe) && (!sc.always_assign || (opt.remultiplex && assigned2)) && !opt.no_output_barcodes) { // Write out barcode read
      std::stringstream o;
      // Write out barcode read
      o << start_char << std::string(names[i].first, names[i].second) << mod_name << "\n";
      if (!opt.remultiplex) {
        o << sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength()) << "\n";
      } else { // Write out remultiplexing barcode
        o << sc.binaryToString(sc.getID(batch_id_mapping[flags[readnum]]), sc.getBarcodeLength()) << "\n";
      }
      if (include_quals) {
        o << "+" << "\n";
        o << std::string(sc.getBarcodeLength(), sc.QUAL) << "\n";
      }
      const std::string& ostr = o.str();
      size_t ostr_len = ostr.length();
      if (!opt.outbam || !r.ofile.empty()) { // Only write barcode read if we're not writing BAM files (or we're writing BAM files but we need to write the current read into a FASTQ file for the "keep" demultiplexing)
        if (use_pipe && !write_barcode_separate_fastq_) {
          if (!opt.no_output_) fwrite(ostr.c_str(), 1, ostr_len, stdout);
        } else if (opt.gzip) {
          gzwrite(r.ofile.empty() ? outb_gz : out_keep_gz[r.ofile][0], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, r.ofile.empty() ? outb : out_keep[r.ofile][0]);
        }
      }
    }
    if (!sc.umi_names.empty() && assigned2 && !opt.no_x_out && !opt.outbam) { // Write out extracted UMIs as needed
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
          if (!opt.no_output_) fwrite(ostr.c_str(), 1, ostr_len, stdout);
        } else if (opt.gzip && !outumi_gz.empty()) {
          gzwrite(outumi_gz[umi_index], ostr.c_str(), ostr_len);
        } else if (!outumi.empty()) {
          fwrite(ostr.c_str(), 1, ostr_len, outumi[umi_index]);
        }
      }
    }
    int jj = -1;
    bool no_output = (!(assigned2) && !write_unassigned_fastq) || opt.x_only; // When to not output read sequences
    std::vector<const char*> s_(jmax, nullptr);
    std::vector<const char*> q_(jmax, nullptr);
    std::vector<int> l_(jmax,0);
    std::vector<std::pair<std::string,std::string> > edited_s_vec;
    if (!no_output && r.modsubs.empty()) {
      for (int j = 0; j < jmax; j++) {
        s_[j] = seqs[i+j].first;
        q_[j] = quals[i+j].first;
        l_[j] = seqs[i+j].second;
      }
    } else if (!no_output) { // If we need to make a substitution in the read
      edited_s_vec = sc.getEditedRead(seqs, quals, i, jmax, r, include_quals);
      for (int j = 0; j < jmax; j++) {
        s_[j] = edited_s_vec[j].first.c_str(); // sequence
        q_[j] = edited_s_vec[j].second.c_str(); // quality
        l_[j] = edited_s_vec[j].first.length();
        if (edited_s_vec[j].first == " ") { // Substitution resulted in an empty string
          l_[j] = 0;
        } else if (edited_s_vec[j].first.empty()) { // Sequence underwent no substitutions so just extract the original sequence in the buffer
          s_[j] = seqs[i+j].first;
          q_[j] = quals[i+j].first;
          l_[j] = seqs[i+j].second;
        }
      }
    }
    for (int j = 0; j < jmax; j++) {
      if (no_output) {
        break;
      }
      if (!opt.select_output_files[j]) {
        continue;
      }
      jj++;
      std::stringstream o;
      const char* s = s_[j];
      int l = l_[j];
      const char* n = names[i+j].first;
      int nl = names[i+j].second;
      const char* q = q_[j];
      // Write out read
      bool embed_final_barcode = assigned && jj == 0 && !write_barcode_separate_fastq_ && !use_pipe && (!sc.always_assign || opt.remultiplex) && !opt.no_output_barcodes;
      o << start_char;
      o << std::string(n,nl) << mod_name << "\n";
      if (embed_final_barcode) {
        if (!opt.remultiplex) {
          o << sc.binaryToString(sc.getID(r.id), sc.getBarcodeLength());
        } else {
          o << sc.binaryToString(sc.getID(batch_id_mapping[flags[readnum]]), sc.getBarcodeLength());
        }
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
      if (assigned2) {
        if (opt.outbam && r.ofile.empty()) {
          writeBam(ostr, nl, jmax == 2 ? jj+1 : 0);
        } else if (use_pipe && !opt.outbam) {
          if (!opt.no_output_) fwrite(ostr.c_str(), 1, ostr_len, stdout);
        } else {
          if (opt.gzip) {
            gzwrite(r.ofile.empty() ? out_gz[jj] : out_keep_gz[r.ofile][j+1], ostr.c_str(), ostr_len);
          } else {
            // note: we use j+1 for out_keep and out_keep_gz because the zeroth index is the barcodes file
            fwrite(ostr.c_str(), 1, ostr_len, r.ofile.empty() ? out[jj] : out_keep[r.ofile][j+1]);
          }
        }
      } else {
        if (opt.gzip) {
          gzwrite(outu_gz[jj], ostr.c_str(), ostr_len);
        } else {
          fwrite(ostr.c_str(), 1, ostr_len, outu[jj]);
        }
      }
    }
    i += incf;
  }
}

void MasterProcessor::writeBam(const std::string& ostr, int readNameLen, int readPair) {
  std::istringstream iss(ostr);
  std::string read_header;
  std::string read_sequence;
  std::string separator;
  std::string quality_sequence;
  std::getline(iss, read_header, '\n');
  std::getline(iss, read_sequence, '\n');
  if (std::getline(iss, separator, '\n')) {
    std::getline(iss, quality_sequence, '\n');
  }
  const char* name = read_header.c_str()+1;
  int nlen = strlen(name);
  const char* seq = read_sequence.c_str();
  int slen = strlen(seq);
  const char* qual = quality_sequence.c_str();
  int qlen = strlen(qual);
  
  if (opt.keep_fastq_comments) { // We want to preserve FASTQ comments from before (and we'll add them as BAM tags)
    readNameLen = 0;
    for (int i = 0; name[i] != '\0' && name[i] != '\t' && name[i] != ' '; i++) {
      readNameLen++;
    }
  }

  std::vector<std::string> bam_tags;
  std::string rname; // This is where we store the --mod-names
  int rname_original_len = 0; // The length of the original --mod-names in the read header
  int start_bam_tags = readNameLen+1; // Start position of the BAM tags in the read header
  if (nlen > readNameLen) {
    if (name[readNameLen] == ':' && name[readNameLen+1] == ':') { // We have :: (--mod-names)
      int i;
      rname_original_len = 2;
      for (i = readNameLen+2; i < nlen; i++, rname_original_len++) {
        if (name[i] == ' ') break; // Encountered a space so we end it
        if (name[i] != '[' && name[i] != ']') { // Skip these characters since BAM files don't allow these
          rname += name[i];
        } else if (name[i] == ']') {
          rname += ":"; // We'll use a colon as a separator between tag names
        }
      }
      rname.pop_back(); // Remove final ] from read name
      rname = opt.sam_tags[6][0] + rname;
      bam_tags.push_back(rname);
      start_bam_tags = i+1;
    }
  }
  if (start_bam_tags < nlen) { // position number 2 implies a string of length 3 or more (length 2: 0,1)
    std::string read_header_bam_tags(name+start_bam_tags);
    std::stringstream ss(read_header_bam_tags);
    std::string curr_tag;
    while (ss >> curr_tag) {
      bam_tags.push_back(curr_tag);
    }
  }
  nlen = readNameLen;
  if (opt.mod_names_bam) {
    nlen += rname_original_len; // Add the original --mod-names back in if --mod-names-bam is specified
  }

  static char buf1[32768];
  static char buf2[32768];
  bam1_t b1;
  bam1_core_t core1;
  core1.tid = -1;
  core1.pos = -1;
  core1.bin = 4680; // magic bin for unmapped reads
  core1.qual = 0;
  core1.mtid = -1;
  core1.mpos = -1;
  core1.isize = 0;
  core1.n_cigar = 0;
  core1.l_qname = nlen;
  core1.l_qseq = slen;
  core1.flag = BAM_FUNMAP | BAM_FMUNMAP;
  if (readPair == 1) {
    core1.flag |= BAM_FPAIRED | BAM_FREAD1;
  } else if (readPair == 2) {
    core1.flag |= BAM_FPAIRED | BAM_FREAD2;
  }
  b1.core = std::move(core1);
  uint8_t* buf = (uint8_t*)&buf1[0];
  memcpy(buf, name, nlen);
  int p = core1.l_qname;
  // copy the sequence
  int lseq = (slen+1)>>1;
  uint8_t *seqp = (uint8_t *) (buf+p);
  memset(seqp,0,lseq);
  for (int i = 0; i < slen; ++i) {
    seqp[i>>1] |= seq_nt16_table[(int)seq[i]] << ((~i&1)<<2);
  }
  p += lseq;
  // copy qual
  for (int i = 0; i < qlen; i++) {
    buf[p+i] = qual[i] - 33;
  }
  p += qlen;
  b1.l_data = p;
  b1.m_data = core1.l_qname + ((core1.l_qseq+2)>>1) + core1.l_qseq;
  b1.data = buf; // structure: qname-cigar-seq-qual-aux
  for (auto& bam_tag: bam_tags) { // XX:X:XXXX
    if (bam_tag.length() < 5) continue;
    b1.data[b1.l_data] = bam_tag[0];
    b1.data[b1.l_data+1] = bam_tag[1];
    if (bam_tag[2] != ':') continue;
    if (bam_tag[4] != ':') continue;
    b1.data[b1.l_data+2] = bam_tag[3];
    if (bam_tag[3] == 'Z') { // The value of our BAM tag is a string
      memcpy(b1.data + b1.l_data + 3, bam_tag.c_str()+5, bam_tag.length()+1-5);
      b1.l_data += 3+bam_tag.length()+1-5;
      b1.m_data += 3+bam_tag.length()+1-5;
    } else if (bam_tag[3] == 'i') { // The value of our BAM tag is an integer
      int val = std::atoi(bam_tag.c_str()+5);
      memcpy(b1.data + b1.l_data + 3, &val, 4);
      b1.l_data += 7;
      b1.m_data += 7;
    }
  }
  int ret = 0;
  ret = sam_write1(bamfp, hdr, &b1);
  if (ret < 0) {
    std::cerr << "Error writing to BAM file... exiting" << std::endl;
    exit(1);
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
   comments = mp.opt.keep_fastq_comments;
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
  full(o.full),
  comments(o.comments) {
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
      mp.FSRs[i].fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full, comments);
    } else {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (mp.SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        mp.SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full, comments);
      }
      // release the reader lock
    }
    
    // process our sequences
    processBuffer();

    // update the results, MP acquires the lock
    int nfiles = mp.nfiles;
    mp.update(seqs.size() / nfiles, rv, seqs, names, quals, flags, readbatch_id);
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
      mp.sc.modifyRead(seqs, quals, i-incf, results, true);
    }
    rv.push_back(results);

    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) { 
        numreads = 0; // reset counter
        int64_t nummapped = mp.sc.getNumMapped();

        std::cerr << '\r' << (mp.numreads/1000000) << "M reads processed";
        if (!mp.sc.always_assign) {
          std::cerr << " (" 
            << std::fixed << std::setw( 3 ) << std::setprecision( 1 ) << ((100.0*nummapped)/double(mp.numreads))
            << "% assigned)";
        } else {
          std::cerr << "         ";
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
  bool full,
  bool comments) {
    
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
          nl[i] = seq[i]->name.l + (comments && seq[i]->comment.l != 0 ? seq[i]->comment.l+1 : 0);
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

          if (full && !comments) {
            pi = buf + bufpos;
            if (seq[i]->qual.l != 0) memcpy(pi, seq[i]->qual.s,l[i]+1);
            else { std::string placeholder(seq[i]->seq.l, 'K'); memcpy(pi, placeholder.c_str(),l[i]+1); } // FASTA format, no quality (default to K)
            bufpos += l[i]+1;
            quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          } else if (full && comments) {
            pi = buf + bufpos;
            if (seq[i]->qual.l != 0) memcpy(pi, seq[i]->qual.s,l[i]+1);
            else { std::string placeholder(seq[i]->seq.l, 'K'); memcpy(pi, placeholder.c_str(),l[i]+1); } // FASTA format, no quality (default to K)
            bufpos += l[i]+1;
            quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            if (seq[i]->comment.l == 0) {
              memcpy(pi, seq[i]->name.s, nl[i]+1);
              names.emplace_back(pi, nl[i]);
              bufpos += nl[i]+1;
            } else {
              memcpy(pi, seq[i]->name.s, (nl[i]-(seq[i]->comment.l+1)));
              names.emplace_back(pi, nl[i]);
              bufpos += (nl[i]-(seq[i]->comment.l+1));
              pi = buf + bufpos;
              const char* blank_space = " ";
              memcpy(pi, blank_space, 1);
              bufpos += 1;
              pi = buf + bufpos;
              memcpy(pi, seq[i]->comment.s, seq[i]->comment.l+1);
              bufpos += seq[i]->comment.l+1;
            }
          }
        }

        numreads++;
        flags.push_back((current_file-nfiles) / nfiles); // flags.push_back(numreads-1);
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



#ifndef SPLITCODE_LIFTWORKFLOW_H
#define SPLITCODE_LIFTWORKFLOW_H

#include "SplitCode.h"

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

#include <unordered_set>
#include <sstream>

#ifndef NO_HTSLIB
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#endif

#include "common.h"
 
using namespace std; 

class LiftWorkflow; // This workflow modifies FASTA files with variants and can lift over annotations; essentially recapitulating parts of g2gtools's functionalities

// GTFParser class to parse GTF files
class GTFParser {
public:
    struct Record {
        std::string seqname;
        std::string source;
        std::string feature;
        int start;
        int end;
        std::string remaining_fields; // everything after the fifth column
    };

    GTFParser(const std::string& gtf_file) {
        parse(gtf_file);
    }

    const std::vector<Record>& get_records() const {
        return records;
    }

private:
    std::vector<Record> records;

    void parse(const std::string& gtf_file) {
#ifndef NO_HTSLIB
        htsFile* infile = hts_open(gtf_file.c_str(), "r");
        if (!infile) {
            std::cerr << "Error opening GTF file: " << gtf_file << std::endl;
            exit(1);
        }

        kstring_t str = {0, 0, NULL};
        while (hts_getline(infile, KS_SEP_LINE, &str) >= 0) {
            std::string line(str.s);
            if (line.empty() || line[0] == '#') continue; // skip empty lines and comments
            std::istringstream iss(line);
            Record record;

            // read the first five columns
            if (!(iss >> record.seqname >> record.source >> record.feature >> record.start >> record.end)) {
                std::cerr << "Error parsing GTF line: " << line << std::endl;
                continue;
            }

            // get the rest of the line as remaining_fields
            std::string remaining;
            std::getline(iss, remaining);
            // remove leading spaces
            remaining.erase(0, remaining.find_first_not_of(" \t"));
            record.remaining_fields = remaining;

            records.push_back(record);
        }

        hts_close(infile);
        if (str.s) free(str.s);
#endif
    }
};

class LiftWorkflow {
public:
  LiftWorkflow(const std::vector<std::string>& argv_, bool diploid_, bool rename_, std::string ref_gtf_, std::string out_gtf_, bool do_filtering_) {
    diploid = diploid_;
    ref_gtf = ref_gtf_;
    out_gtf = out_gtf_;
    do_filtering = do_filtering_;
    rename = rename_;
    make_temp_fasta_file = false;
    std::vector<std::string> argv;
    argv.push_back("splitcode --lift");
    temp_file_name_prefix = "";
    for (auto x : argv_) {
      argv.push_back(x);
      temp_file_name_prefix += x + ",";
    }

    size_t argc = argv.size();
    if (argc < 4) {
      std::cout << "Usage: " << argv[0] << " <ref_fasta> <vcf_file> <sample> [--diploid] [--filter] [--ref-gtf <ref_gtf>] [--out-gtf <out_gtf>]" << std::endl;
      exit(1);
    }
    ref_fasta = argv[1];
    vcf_file = argv[2];
    temp_file_name_prefix = temp_file_name_prefix + ref_gtf + "," + out_gtf + "," + ref_fasta + "," + vcf_file + "," + std::to_string(diploid);

    for (int i = 3; i < argc; i++) {
      std::string arg = argv[i];
      samples.push_back(arg);
    }

    if (samples.empty()) {
      std::cerr << "Error: one sample must be specified" << std::endl;
      exit(1);
    } else if (samples.size() > 1) {
      std::cerr << "Error: cannot supply more than one sample" << std::endl;
      exit(1);
    }

    if (!ref_gtf.empty() && out_gtf.empty()) {
      std::cerr << "Error: --out-gtf must be specified if --ref-gtf is provided" << std::endl;
      exit(1);
    }
  }
  
  void modify_fasta() {
#ifndef NO_HTSLIB
    
    // Prepare VCF file
    htsFile* vcf_fp = hts_open(vcf_file.c_str(), "rb");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    unordered_set<string> sample_set(samples.begin(), samples.end());
    vector<int> sample_indices; // Indices in the VCF of samples we're interested in
    vector<string> sample_names;  // Names of the samples we're interested in
    unordered_map<int, int> sample_index_map; // Map from VCF sample index to our sample index
    size_t nsmpl = bcf_hdr_nsamples(vcf_hdr);
    for (int i = 0; i < nsmpl; i++) {
      string sample_name = bcf_hdr_int2id(vcf_hdr, BCF_DT_SAMPLE, i);
      if (sample_set.count(sample_name)) {
        sample_index_map[i] = sample_names.size();
        sample_indices.push_back(i);
        sample_names.push_back(sample_name);
      }
    }

    // Check to make sure everything is valid with regards to number of samples and user-input options before proceeding
    bool proceed = true;
    if (sample_indices.size() != 1) {
      proceed = false;
      std::cerr << "Error: The sample supplied must be present in the VCF" << std::endl;
    }
    
    if (!proceed) {
      bcf_hdr_destroy(vcf_hdr);
      hts_close(vcf_fp);
      exit(1);
    }
    
    // Prepare input FASTA file
    faidx_t* fai = nullptr;
    std::string fname;
    if (ref_fasta.substr(ref_fasta.find_last_of(".") + 1) == "gz") {
      temp_file_name_fasta = generate_tmp_file(temp_file_name_prefix, "./");
      make_temp_fasta_file = true;
      FILE* temp_file = std::fopen(temp_file_name_fasta.c_str(), "wb");
      decompressGzipToFile(ref_fasta, temp_file);
      fclose(temp_file);
      fname = temp_file_name_fasta;
    } else {
      fname = ref_fasta;
    }
    std::string faiFilePath = generate_tmp_file(temp_file_name_prefix + "fai", "./");
    std::string faiFilePath_1 = faiFilePath + ".fai";
    std::string faiFilePath_2 = faiFilePath + ".gzi";
    if (FILE *file = fopen(faiFilePath_1.c_str(), "r")) {
      fclose(file);
      // File exists, attempt to delete it
      std::cerr << "fai file already exists; removing it to build a new one" << std::endl;
      std::remove(faiFilePath_1.c_str());
    }
    if (FILE *file = fopen(faiFilePath_2.c_str(), "r")) {
      fclose(file);
      // File exists, attempt to delete it
      std::cerr << "fai file already exists; removing it to build a new one" << std::endl;
      std::remove(faiFilePath_2.c_str());
    }
    //fai = fai_load3(fname.c_str(), faiFilePath_1.c_str(), faiFilePath_2.c_str(), FAI_CREATE|FAI_CACHE);
    fai = fai_load(fname.c_str());
    
    // Initialize input/output variables
    auto& out_fasta = std::cout;
    bcf1_t* record = bcf_init();
    
    int n_seqs = faidx_fetch_nseq(fai); // number of chromosomes in FASTA file
    std::vector<std::string> chrom_vec;
    for (int idx = 0; idx < n_seqs; ++idx) {
      const char* seq_name = faidx_iseq(fai, idx);
      chrom_vec.push_back(std::string(seq_name));
    }

    std::string chrom;
    size_t i = 0;
    std::vector<std::vector<VarLocation>> var_locations;
    var_locations.resize(diploid ? 2 : sample_names.size());
    std::unordered_set<std::string> seen_chromosomes;
    char* ref_seq = nullptr;
    int ref_len = 0;
    std::vector<int> chrom_len; // for storing the new length of a chromosome (including indels)
    chrom_len.resize(diploid ? 2 : sample_names.size(), 0); 
    while ((modify_fasta_helper(sample_names, sample_index_map, nsmpl, vcf_fp, vcf_hdr, record, fai, chrom, var_locations, ref_seq, ref_len, chrom_len, seen_chromosomes, diploid, i++))) {
      prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, chrom_len, sample_names, out_fasta, fasta_line_length);
    }
    
    prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, chrom_len, sample_names, out_fasta, fasta_line_length);
    
    for (int idx = 0; idx < n_seqs; ++idx) {
      const char* seq_name = faidx_iseq(fai, idx);
      if (!seen_chromosomes.count(std::string(seq_name))) {
        std::string chrom = std::string(seq_name);
        ref_seq = fai_fetch(fai, chrom.c_str(), &ref_len);
        for (int i = 0; i < sample_names.size(); i++) {
          // Include final stretch of sequences (i.e. the stuff after the last variant) into
          VarLocation final_loc;
          final_loc.position = 0;
          final_loc.length = ref_len;
          final_loc.variant = "";
          if (diploid) {
             var_locations[0].push_back(final_loc);
             var_locations[1].push_back(final_loc);
             chrom_len[0] = ref_len;
             chrom_len[1] = ref_len;
          } else {
             var_locations[i].push_back(final_loc);
             chrom_len[i] = ref_len;
          }
        }
        prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, chrom_len, sample_names, out_fasta, fasta_line_length);
      }
    }

    
    bcf_destroy(record);
    bcf_hdr_destroy(vcf_hdr);
    hts_close(vcf_fp);
    fai_destroy(fai);
    if (make_temp_fasta_file) std::remove(temp_file_name_fasta.c_str());
    
    if (!ref_gtf.empty()) {
      lift_over_gtf(ref_gtf, out_gtf);
    }
#endif
  }
  
  
  std::string ref_fasta;
  std::string vcf_file;
  std::vector<string> samples;
  std::string ref_gtf;
  std::string out_gtf;
  bool rename;
  std::string temp_file_name_prefix;
  std::string temp_file_name_fasta;
  bool diploid;
  bool do_filtering;
  bool make_temp_fasta_file;
  
private:
  
  struct VarLocation {
    uint32_t position;
    uint32_t length; // length from position to start of variant
    std::string variant; // allele
  };

  struct CoordinateShift {
    uint32_t orig_pos; // original position in reference genome
    int32_t cumulative_shift; // cumulative shift in coordinate up to this position
  };

  std::unordered_map<std::string, std::vector<std::vector<CoordinateShift>>> coordinate_shifts_all;

  static const size_t fasta_line_length = 60;
  
  
    // sample_names: vector of sample names we're interested in
    // sample_index_map: map from VCF sample index to our sample index
#ifndef NO_HTSLIB
  bool modify_fasta_helper(const vector<string>& sample_names, const unordered_map<int, int>& sample_index_map, size_t nsmpl, htsFile* vcf_fp, bcf_hdr_t* vcf_hdr, bcf1_t* record, faidx_t* fai, std::string& chrom, std::vector<std::vector<VarLocation>>& var_locations, char*& ref_seq, int& ref_len, std::vector<int>& chrom_len, std::unordered_set<std::string>& seen_chromosomes, bool diploid, size_t iteration, bool indel = false) {
    chrom_len.assign(diploid ? 2 : sample_names.size(), 0);
    std::vector<int> prev_start;
    prev_start.assign(diploid ? 2 : sample_names.size(), 0);
    std::string prev_chrom = "";
    int prev_start_ = 0;
    bool started_loop = false;
    bool carry_on = false;
    var_locations.clear(); // Clear everything stored
    var_locations.resize(diploid ? 2 : sample_names.size());
    std::vector<int> deletion_start; // used to check if a variant is within an already deleted portion
    deletion_start.assign(diploid ? 2 : sample_names.size(), -1);
    std::vector<int> deletion_end;
    deletion_end.assign(diploid ? 2 : sample_names.size(), -1);
   
    // Initialize coordinate shifts
    std::vector<std::vector<CoordinateShift>> coordinate_shifts;
    coordinate_shifts.resize(diploid ? 2 : sample_names.size());
    std::vector<int32_t> cumulative_shift;
    cumulative_shift.assign(diploid ? 2 : sample_names.size(), 0);

    while ((carry_on = (iteration != 0 && !started_loop)) || bcf_read(vcf_fp, vcf_hdr, record) == 0) { // Note: For carry_on, = is an assignment operator
      if (!carry_on) { // If carrying on from before, no need to do these things again (bcf_read will not be performed if carrying on)
        bcf_unpack(record, BCF_UN_STR);
        bcf_unpack(record, BCF_UN_INFO);
        bcf_unpack(record, BCF_UN_FMT);
        if (do_filtering) bcf_unpack(record, BCF_UN_FLT);
      }
      chrom = bcf_hdr_id2name(vcf_hdr, record->rid);
      if (!started_loop) { // Only load new chromosome if it's our first time going through this loop
        prev_chrom = chrom;
        ref_seq = fai_fetch(fai, chrom.c_str(), &ref_len);
        if (!ref_seq) {
          std::cerr << "Error: Chromosome " << chrom << " exists in VCF but not FASTA" << std::endl;
          exit(1);
        }
      }
      int start = record->pos;
      vector<string> alleles(record->n_allele);
      for (int i = 0; i < record->n_allele; i++) {
        alleles[i] = std::string(record->d.allele[i]);
      }
      if (iteration == 0 && !started_loop) { // Opening up VCF file for the very first time
        seen_chromosomes.insert(chrom);
      } else if (chrom != prev_chrom) { // Encountering a new chromosome in the VCF file
        if (seen_chromosomes.count(chrom)) {
          std::cerr << "Error: VCF is unsorted, encountered chromosome " << chrom << " separately" << std::endl;
          exit(1);
        }
        seen_chromosomes.insert(chrom);
        // Output everything stored so far
        for (int i = 0; i < (diploid ? 2 : sample_names.size()); i++) {
          // Include final stretch of sequences (i.e. the stuff after the last variant) into
          VarLocation final_loc;
          if (var_locations[i].size() != 0) {
            final_loc.position = prev_start[i]; 
            final_loc.length = ref_len - final_loc.position;
            final_loc.variant = "";
            chrom_len[i] += final_loc.length;
          } else {
            final_loc.position = 0;
            final_loc.length = ref_len;
            final_loc.variant = "";
            chrom_len[i] += final_loc.length;
          }
          var_locations[i].push_back(final_loc);
        }
        // Store coordinate shifts
        coordinate_shifts_all[prev_chrom] = coordinate_shifts;
        chrom = prev_chrom;
        return true;
      } else if (prev_start_ > start) {
        std::cerr << "Error: VCF is unsorted at chromosome " << chrom << ", position: " << start << std::endl;
        exit(1);
      }
      started_loop = true;
      int32_t *gt_arr = NULL, ngt_arr = 0;
      char pass_str[] = "PASS";
      if (do_filtering && bcf_has_filter(vcf_hdr, record, pass_str) != 1) continue;
      int ngt = bcf_get_genotypes(vcf_hdr, record, &gt_arr, &ngt_arr);
      if ( ngt<=0 ) continue; // GT not present 
      int max_ploidy = ngt/nsmpl;
      for (int i = 0; i < nsmpl; i++) {
        auto it = sample_index_map.find(i);
        if (it == sample_index_map.end()) continue; // Not a sample we're interested in (technically not needed since we only permit one sample)
        int32_t *ptr = gt_arr + i*max_ploidy;
        std::vector<int> allele_indices;
        for (int j=0; j<max_ploidy; j++) {
          // if true, the sample has smaller ploidy
          if ( ptr[j]==bcf_int32_vector_end ) {
            std::cerr << "ERROR: incorrect ploidy" << std::endl;
            exit(1);
          }
          int allele_index = bcf_gt_allele(ptr[j]);
          allele_indices.push_back(allele_index);
        }
        std::string allele_1, allele_2;
        if (diploid) {
          if (allele_indices.size() != 2) continue; // Don't do anything unless we have exactly two alleles
          if (allele_indices[0] == -1 || allele_indices[1] == -1) continue; // in case GT is ./.
          char s = toupper(ref_seq[start]);
          allele_1 = alleles[allele_indices[0]];
          allele_2 = alleles[allele_indices[1]];
          std::string ref_allele = alleles[0];

          bool make_loc1 = true; // these variables are true by default and will turn false if any conditions are violated
          bool make_loc2 = true;
          if (allele_indices[0] <= 0) make_loc1 = false; // GT for allele is non-zero
          if (allele_indices[1] <= 0) make_loc2 = false;
          if (start >= deletion_start[0] && start <= deletion_end[0]) make_loc1 = false; // inside deleted segment?
          if (start >= deletion_start[1] && start <= deletion_end[1]) make_loc2 = false;

          // insertion! 
          if (allele_1.length() > ref_allele.length() || allele_2.length() > ref_allele.length()) { 
            if (make_loc1 && allele_1.length() > ref_allele.length()) {
              allele_1 = allele_1.substr(1, allele_1.length() - 1);
            }
            if (make_loc2 && allele_2.length() > ref_allele.length()) {
              allele_2 = allele_2.substr(1, allele_2.length() - 1);
            }
            ref_allele = "";
            start = start + 1;
          }

          // deletion!
          else if (allele_1.length() < ref_allele.length() || allele_2.length() < ref_allele.length()) {  
            ref_allele = ref_allele.substr(1, ref_allele.length() - 1);
            start = start + 1;

            if (make_loc1 && allele_1.length() < alleles[0].length()) {
              allele_1 = "";
              deletion_start[0] = start;
              deletion_end[0] = start + ref_allele.length() - 1;
            }
            else {
              allele_1 = ref_allele;
            }
            if (make_loc2 && allele_2.length() < alleles[0].length()) {
              allele_2 = "";
              deletion_start[1] = start;
              deletion_end[1] = start + ref_allele.length() - 1;
            }
            else {
              allele_2 = ref_allele;
            }
          }

          if (start < prev_start[0]) make_loc1 = false;
          if (start < prev_start[1]) make_loc2 = false;
          
          if (make_loc1) { 
            if (allele_1.length() == 1 && allele_indices[0] <= 0 && s != 'A' && s != 'T' && s != 'C' && s != 'G') allele_1[0] = ref_seq[start];
            VarLocation loc_1;
            loc_1.position = prev_start[0];
            loc_1.length = start-prev_start[0];
            loc_1.variant = allele_1;

            var_locations[0].push_back(loc_1);
            prev_start[0] = start + ref_allele.length();
            chrom_len[0] += loc_1.length + allele_1.length();

            // Update coordinate shifts
            int shift = allele_1.length() - ref_allele.length();
            cumulative_shift[0] += shift;
            CoordinateShift coord_shift;
            coord_shift.orig_pos = start + ref_allele.length();
            coord_shift.cumulative_shift = cumulative_shift[0];
            coordinate_shifts[0].push_back(coord_shift);
          }

          if (make_loc2) { 
            if (allele_2.length() == 1 && allele_indices[1] <= 0 && s != 'A' && s != 'T' && s != 'C' && s != 'G') allele_2[0] = ref_seq[start];
            VarLocation loc_2;
            loc_2.position = prev_start[1];
            loc_2.length = start-prev_start[1];
            loc_2.variant = allele_2;

            var_locations[1].push_back(loc_2);
            prev_start[1] = start + ref_allele.length();
            chrom_len[1] += loc_2.length + allele_2.length();

            // Update coordinate shifts
            int shift = allele_2.length() - ref_allele.length();
            cumulative_shift[1] += shift;
            CoordinateShift coord_shift;
            coord_shift.orig_pos = start + ref_allele.length();
            coord_shift.cumulative_shift = cumulative_shift[1];
            coordinate_shifts[1].push_back(coord_shift);
          }
        } else { // Non-diploid mode
          if (allele_indices.size() != 2) continue; // Only process diploid genotypes
          if (allele_indices[0] == -1 || allele_indices[1] == -1) continue; // in case GT is ./.
          if (allele_indices[0] != allele_indices[1]) continue; // Skip heterozygous genotypes
          int allele_idx = allele_indices[0];
          if (allele_idx == 0) continue; // Homozygous reference, no change
          if (allele_idx > 0) {
           // Homozygous non-reference allele
           char s = toupper(ref_seq[start]);
           std::string allele = alleles[allele_idx];
           std::string ref_allele = alleles[0];
           
           bool make_loc = true;
           if (start >= deletion_start[0] && start <= deletion_end[0]) make_loc = false; // inside deleted segment?

           // Insertion
           if (allele.length() > ref_allele.length()) {
            if (make_loc) allele = allele.substr(1); // Remove the first base
            ref_allele = "";
            start = start + 1;
           }
           // Deletion
           else if (allele.length() < ref_allele.length()) {
            ref_allele = ref_allele.substr(1);
            start = start + 1;
            if (make_loc && allele.length() < alleles[0].length()) {
             allele = "";
             deletion_start[0] = start;
             deletion_end[0] = start+ref_allele.length()-1;
            } else {
             allele = ref_allele;
            }
           }

           if (start < prev_start[0]) make_loc = false;

           if (make_loc) {
            if (allele.length() == 1 && allele_idx <= 0 && s != 'A' && s != 'T' && s != 'C' && s != 'G') allele[0] = ref_seq[start];

            VarLocation loc;
            loc.position = prev_start[0];
            loc.length = start - prev_start[0];
            loc.variant = allele;

            var_locations[0].push_back(loc);
            prev_start[0] = start + ref_allele.length();//allele.length();
            chrom_len[0] += loc.length + allele.length();

            // Update coordinate shifts
            int shift = allele.length() - ref_allele.length();
            cumulative_shift[0] += shift;
            CoordinateShift coord_shift;
            coord_shift.orig_pos = start + ref_allele.length();
            coord_shift.cumulative_shift = cumulative_shift[0];
            coordinate_shifts[0].push_back(coord_shift);
           }
          }
        }
      }
      free(gt_arr);
    }
    // Output everything stored so far
    for (int i = 0; i < (diploid ? 2 : sample_names.size()); i++) {
      if (!diploid && sample_names[i].empty()) continue; // Don't care about this sample
      // Include final stretch of sequences (i.e. the stuff after the last variant) into
      VarLocation final_loc;
      if (var_locations[i].size() != 0) {
        final_loc.position = prev_start[i]; 
        final_loc.length = ref_len - final_loc.position;
        final_loc.variant = "";
        chrom_len[i] += final_loc.length;
      } else {
        final_loc.position = 0;
        final_loc.length = ref_len;
        final_loc.variant = "";
        chrom_len[i] += final_loc.length;
      }
      var_locations[i].push_back(final_loc);
    }
    // Store coordinate shifts
    coordinate_shifts_all[chrom] = coordinate_shifts;
    return false;
  }
  
  void prepareFastaAndPrintChromosome(
      char* ref_seq,
      std::vector<std::vector<VarLocation>>& var_locations,
      bool diploid,
      const std::string& prev_chrom,
      int ref_len,
      std::vector<int>& chrom_len, 
      const std::vector<std::string>& sample_names,
      std::ostream& out_fasta,
      int fasta_line_length
  ) {
    for (int i = 0; i < var_locations.size(); i++) {
      std::string prefix = "";
      std::string suffix = "";

      if (!diploid && sample_names[i].empty()) continue;
      
      if (diploid) suffix = std::string("_") + std::string(i == 0 ? "L" : "R") + std::string(" ") + prev_chrom + std::string("_") + std::string(i == 0 ? "L" : "R") + std::string(":1-") + std::to_string(chrom_len[i]) + std::string(" from|") + prev_chrom + std::string("_") + std::string(i == 0 ? "L" : "R") + std::string(":1-") + std::to_string(ref_len);
      else {
       prefix = sample_names[i] + "_";
       if (!rename) prefix = "";
       // Construct header similar to diploid mode
       suffix = " " + prev_chrom + ":1-" + std::to_string(chrom_len[i]) + " from|" + prev_chrom + ":1-" + std::to_string(ref_len);
      }
      
      out_fasta << ">" << prefix << prev_chrom << suffix << std::endl;
      
      // Print out the chromosome, introducing the variants
      size_t currentLineLength = 0;
      const char* currentChar;
      size_t charsToPrint;
      auto lineLength = fasta_line_length;
      
      for (const auto loc : var_locations[i]) { // Loop through each variant
        currentChar = ref_seq + loc.position; // Pointer to the current string
        size_t remainingChars = loc.length;
        bool variant_printed = false;
        
        while (true) { // Print out sequence (ref+variant)
          if (remainingChars == 0 && !variant_printed) { // Now move onto printing the variant
            currentChar = loc.variant.c_str();
            remainingChars = loc.variant.length();
            variant_printed = true;
          }
          
          if (remainingChars == 0) break;
          
          // Calculate how many characters can be printed on this line
          charsToPrint = std::min(lineLength - currentLineLength, remainingChars);
          
          // Print directly from the pointer, spanning the needed characters
          out_fasta.write(currentChar, charsToPrint);
          currentLineLength += charsToPrint;
          remainingChars -= charsToPrint;
          currentChar += charsToPrint; // Move the pointer forward
          
          // If the line is complete, print a newline and reset line length
          if (currentLineLength >= lineLength) {
            out_fasta << "\n";
            currentLineLength = 0;
          }
        }
      }
      
      // Print a newline if there were any characters printed in the last line
      if (currentLineLength > 0) {
        out_fasta << "\n";
      }
      
      out_fasta.flush(); 
    }
    
    if (ref_seq != nullptr) free(ref_seq);
    var_locations.clear(); // Clear everything stored
    var_locations.resize(diploid ? 2 : sample_names.size());
  }

  void decompressGzipToFile(const std::string& gzipFilePath, FILE* outFile) {
    // Open gzip file for reading
    gzFile gz = gzopen(gzipFilePath.c_str(), "rb");
    if (!gz) {
      std::cerr << "Failed to open gzip file: " << gzipFilePath << std::endl;
      exit(1);
    }
    
    // Buffer to hold data read from the gzip file
    char buffer[4096];  // Buffer size can be adjusted as needed
    int bytesRead;

    // Read the compressed data in chunks
    while ((bytesRead = gzread(gz, buffer, sizeof(buffer))) > 0) {
      // Write the decompressed data to the output file
      if (fwrite(buffer, 1, bytesRead, outFile) != static_cast<size_t>(bytesRead)) {
        std::cerr << "Failed to write decompressed data to file." << std::endl;
        gzclose(gz);
        exit(1);
      }
    }
    
    // Close the gzip file
    gzclose(gz);
  }
#endif
  
    void lift_over_gtf(const string& ref_gtf, const string& out_gtf) {
        GTFParser parser(ref_gtf);
        ofstream out_gtf_file(out_gtf);
        if (!out_gtf_file) {
            std::cerr << "Error opening output GTF file: " << out_gtf << std::endl;
            exit(1);
        }
        const auto& records = parser.get_records();
        for (const auto& record : records) {
            const std::string& chrom = record.seqname;
            if (diploid) {
                // Output records for both haplotypes
                for (int hap = 0; hap < 2; ++hap) {
                    std::string hap_suffix = hap == 0 ? "_L" : "_R";
                    std::string chrom_with_suffix = chrom + hap_suffix;

                    int new_start = record.start;
                    int new_end = record.end;

                    if (coordinate_shifts_all.find(chrom) != coordinate_shifts_all.end()) {
                        // Apply coordinate shifts
                        const auto& shifts = coordinate_shifts_all.at(chrom)[hap];
                        new_start = map_coordinate(record.start, shifts);
                        new_end = map_coordinate(record.end, shifts);
                    }

                    out_gtf_file << chrom_with_suffix << "\t" << record.source << "\t" << record.feature << "\t" << new_start << "\t" << new_end << "\t" << record.remaining_fields << std::endl;
                }
            } else {
                // Non-diploid case
                for (int i = 0; i < samples.size(); ++i) {
                    const std::string& sample_name = samples[i];
                    std::string chrom_with_prefix = sample_name + "_" + chrom;
                    int new_start = record.start;
                    int new_end = record.end;

                    if (coordinate_shifts_all.find(chrom) != coordinate_shifts_all.end()) {
                        if (i < coordinate_shifts_all.at(chrom).size()) {
                            const auto& shifts = coordinate_shifts_all.at(chrom)[i];
                            new_start = map_coordinate(record.start, shifts);
                            new_end = map_coordinate(record.end, shifts);
                        }
                    }
                    out_gtf_file << chrom_with_prefix << "\t" << record.source << "\t" << record.feature << "\t" << new_start << "\t" << new_end << "\t" << record.remaining_fields << std::endl;
                }
            }
        }
        out_gtf_file.close();
    }

    int map_coordinate(int old_pos, const std::vector<CoordinateShift>& shifts) {
        // Find the last shift where orig_pos <= old_pos
        if (shifts.empty() || old_pos < shifts[0].orig_pos) {
            // No shift applies
            return old_pos;
        }
        // Binary search to find the appropriate shift
        size_t left = 0;
        size_t right = shifts.size();
        while (left < right) {
            size_t mid = (left + right) / 2;
            if (shifts[mid].orig_pos <= old_pos) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        // The shift to apply is shifts[left - 1]
        int shift = shifts[left - 1].cumulative_shift;
        return old_pos + shift;
    }

  
  int my_mkdir(const char *path, mode_t mode) {
#ifdef _WIN64
    return mkdir(path);
#else
    return mkdir(path,mode);
#endif
  }
  
  std::string generate_tmp_file(std::string seed, std::string tmp_dir) {
    struct stat stFileInfo;
    auto intStat = stat(tmp_dir.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << tmp_dir << " exists and is not a directory" << endl;
        exit(1);
      }
    } else {
      // create directory
      if (my_mkdir(tmp_dir.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << tmp_dir << endl;
        exit(1);
      }
    }
    std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::string tmp_file = "splitcode.lift.";
    srand((unsigned int)std::hash<std::string>{}(seed));
    int pos;
    while (tmp_file.length() < 32) {
      pos = ((rand() % (base.size() - 1)));
      tmp_file += base.substr(pos, 1);
    }
    return tmp_dir + "/" + tmp_file;
  }
};
#endif // SPLITCODE_LIFTWORKFLOW_H

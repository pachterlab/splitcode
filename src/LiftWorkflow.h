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
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#endif

#include "common.h"
 
using namespace std; 

class LiftWorkflow; // This workflow modifies FASTA files with variants and can lift over annotations; essentially recapitulating parts of g2gtools's functionalities

class LiftWorkflow {
public:
  LiftWorkflow(const std::vector<std::string>& argv_, bool diploid_, std::string indel_vcf_file_, std::string ref_gtf_, std::string out_gtf_) {
    indel_vcf_file = indel_vcf_file_;
    diploid = diploid_;
    ref_gtf = ref_gtf_;
    out_gtf = out_gtf_;
    std::vector<std::string> argv;
    argv.push_back("splitcode --lift");
    temp_file_name = "";
    for (auto x : argv_) {
      argv.push_back(x);
      temp_file_name += x + ",";
    }

    size_t argc = argv.size();
    if (argc < 4) {
      std::cout << "Usage: " << argv[0] << " <ref_fasta> <vcf_file> <sample1> [<sample2> ...] [--diploid] [--indel <vcf_file>] [--ref-gtf <ref_gtf>] [--out-gtf <out_gtf>]" << std::endl;
      exit(1);
    }
    temp_file_name = generate_tmp_file(temp_file_name + ref_gtf + "," + out_gtf + "," + ref_fasta + "," + vcf_file + "," + std::to_string(diploid), "./");
    ref_fasta = argv[1];
    vcf_file = argv[2];

    for (int i = 3; i < argc; i++) {
      std::string arg = argv[i];
      samples.push_back(arg);
    }

    if (samples.empty()) {
      std::cerr << "Error: at least one sample must be specified" << std::endl;
      exit(1);
    }

    if (!ref_gtf.empty() && out_gtf.empty()) {
      std::cerr << "Error: --out-gtf must be specified if --ref-gtf is provided" << std::endl;
      exit(1);
    }
  }
  
  void modify_fasta() {
#ifndef NO_HTSLIB
    
    bool indel = !indel_vcf_file.empty();
    
    // Prepare VCF file
    htsFile* vcf_fp = hts_open(vcf_file.c_str(), "rb");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    unordered_set<string> sample_set(samples.begin(), samples.end());
    vector<int> sample_indices;
    vector<string> sample_names;
    size_t nsmpl = bcf_hdr_nsamples(vcf_hdr);
    for (int i = 0; i < nsmpl; i++) {
      string sample_name = bcf_hdr_int2id(vcf_hdr, BCF_DT_SAMPLE, i);
      if (sample_set.count(sample_name)) {
        sample_indices.push_back(i);
        sample_names.resize(i+1);
        sample_names[i] = sample_name;
      }
    }
    
    htsFile* vcf_fp_indel;
    bcf_hdr_t* vcf_hdr_indel;
    size_t nsmpl_indel;
    vector<int> sample_indices_indel;
    vector<string> sample_names_indel;
    if (indel) {
      vcf_fp_indel = hts_open(indel_vcf_file.c_str(), "rb");
      bcf_hdr_t* vcf_hdr_indel = bcf_hdr_read(vcf_fp_indel);
      nsmpl_indel = bcf_hdr_nsamples(vcf_hdr_indel);
      for (int i = 0; i < nsmpl_indel; i++) {
        string sample_name = bcf_hdr_int2id(vcf_hdr_indel, BCF_DT_SAMPLE, i);
        if (sample_set.count(sample_name)) {
          sample_indices_indel.push_back(i);
          sample_names_indel.resize(i+1);
          sample_names_indel[i] = sample_name;
        }
      }
    }
    
    
    // Check to make sure everything is valid with regards to number of samples and user-input options before proceeding
    bool proceed = true;
    if (sample_indices.size() < 1) {
      proceed = false;
      std::cerr << "Error: Must have at least one sample present in the VCF" << std::endl;
    }
    if (diploid && sample_indices.size() != 1) {
      proceed = false;
      std::cerr << "Error: Cannot use --diploid unless only one sample is provided in the VCF" << std::endl;
    }
    if (indel) {
      if (sample_indices_indel.size() < 1) {
        proceed = false;
        std::cerr << "Error: Must have at least one sample present in the VCF" << std::endl;
      }
      if (diploid && sample_indices_indel.size() != 1) {
        proceed = false;
        std::cerr << "Error: Cannot use --diploid unless only one sample is provided in the VCF" << std::endl;
      }
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
      FILE* temp_file = std::fopen(temp_file_name.c_str(), "wb");
      decompressGzipToFile(ref_fasta, temp_file);
      fname = temp_file_name;
    } else {
      fname = ref_fasta;
    }
    std::string faiFilePath = fname + ".fai";
    if (FILE *file = fopen(faiFilePath.c_str(), "r")) {
      fclose(file);
      // File exists, attempt to delete it
      std::cerr << "fai file already exists; removing it to build a new one" << std::endl;
      std::remove(faiFilePath.c_str());
    }
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
    char* ref_seq_indel = nullptr;
    int ref_len = 0;
    int ref_len_indel = 0;
    while ((modify_fasta_helper(sample_names, nsmpl, vcf_fp, vcf_hdr, record, fai, chrom, var_locations, ref_seq, ref_len, seen_chromosomes, diploid, i++))) {
      prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, sample_names, out_fasta, fasta_line_length);
      if (indel) { // Go up to chromosome chrom; maybe a different out_fasta above too; maybe ingest both original and new ref_seq as well
        while ((modify_fasta_helper(sample_names, nsmpl, vcf_fp, vcf_hdr, record, fai, chrom, var_locations, ref_seq_indel, ref_len_indel, seen_chromosomes, diploid, i++, true))) {
          // Better idea: Go up to chromosome chrom and have two final var_locations vecs, and combine at the end?
        }
      }
    }
    if (indel) { // Go up to the remaining untouched chromosomes
      
    }
    prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, sample_names, out_fasta, fasta_line_length);
    // TODO: indel
    
    for (int idx = 0; idx < n_seqs; ++idx) {
      const char* seq_name = faidx_iseq(fai, idx);
      if (!seen_chromosomes.count(std::string(seq_name))) {
        std::string chrom = std::string(seq_name);
        ref_seq = fai_fetch(fai, chrom.c_str(), &ref_len);
        for (int i = 0; i < sample_names.size(); i++) {
          if (sample_names[i].empty()) continue; // Don't care about this sample
          // Include final stretch of sequences (i.e. the stuff after the last variant) into
          VarLocation final_loc;
          final_loc.position = 0;
          final_loc.length = ref_len;
          final_loc.variant = "";
          if (diploid) { // REMOVE:
             var_locations[0].push_back(final_loc);
             var_locations[1].push_back(final_loc);
          } else var_locations[i].push_back(final_loc);
        }
        prepareFastaAndPrintChromosome(ref_seq, var_locations, diploid, chrom, ref_len, sample_names, out_fasta, fasta_line_length);
      }
    }

    
    bcf_destroy(record);
    bcf_hdr_destroy(vcf_hdr);
    hts_close(vcf_fp);
    fai_destroy(fai);
    std::remove(temp_file_name.c_str());
    
    if (!ref_gtf.empty()) {
      //lift_over_gtf(ref_gtf, out_gtf, mod_ranges[ref_gtf]);
    }
#endif
}
  
  
  std::string ref_fasta;
  std::string vcf_file;
  std::string indel_vcf_file;
  std::vector<string> samples;
  std::string ref_gtf;
  std::string out_gtf;
  std::string temp_file_name;
  bool diploid;
  
private:
  
  struct VarLocation {
    uint32_t position;
    uint32_t length; // length from position to start of variant
    std::string variant; // allele
  };
  
  static const size_t fasta_line_length = 60;
  
  // sample_names: maps sample indices to sample names (if a sample name is empty for a given sample index, we don't care about that sample and skip over it)
  // nsmpl: total number of samples in the VCF (i.e. bcf_hdr_nsamples)
#ifndef NO_HTSLIB
  bool modify_fasta_helper(const vector<string>& sample_names, size_t nsmpl, htsFile* vcf_fp, bcf_hdr_t* vcf_hdr, bcf1_t* record, faidx_t* fai, std::string& chrom, std::vector<std::vector<VarLocation>>& var_locations, char*& ref_seq, int& ref_len, std::unordered_set<std::string>& seen_chromosomes, bool diploid, size_t iteration, bool indel = false) {
    std::vector<int> prev_start;
    prev_start.resize(sample_names.size(), 0);
    std::string prev_chrom = "";
    int prev_start_ = 0;
    bool started_loop = false;
    bool carry_on = false;
    var_locations.clear(); // Clear everything stored
    var_locations.resize(diploid ? 2 : sample_names.size());
    prev_start.resize(sample_names.size(), 0);
    prev_start_ = 0;
    while ((carry_on = (iteration != 0 && !started_loop)) || bcf_read(vcf_fp, vcf_hdr, record) == 0) { // Note: For carry_on, = is an assignment operator
      // Note: We use carry_on to determine if we want to carry on from a previous call to this function (i.e. it's not our first time going through a chromosome and haven't started looping yet)
      if (!carry_on) { // If carrying on from before, no need to do these things again (bcf_read will not be performed if carrying on)
        bcf_unpack(record, BCF_UN_STR);
        bcf_unpack(record, BCF_UN_INFO);
        bcf_unpack(record, BCF_UN_FMT);
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
        if (seen_chromosomes.count(chrom)) { // TODO: store somewhere
          std::cerr << "Error: VCF is unsorted, encountered chromosome " << chrom << " separately" << std::endl;
          exit(1);
        }
        seen_chromosomes.insert(chrom);
        // Output everything stored so far
        for (int i = 0; i < (diploid ? 2 : sample_names.size()); i++) {
          if (!diploid && sample_names[i].empty()) continue; // Don't care about this sample 
          // Include final stretch of sequences (i.e. the stuff after the last variant) into
          VarLocation final_loc;
          if (var_locations[i].size() != 0) {
            auto &second_to_final_loc = var_locations[i][var_locations[i].size()-1];
            final_loc.position = second_to_final_loc.position + second_to_final_loc.length + second_to_final_loc.variant.length();
            final_loc.length = ref_len - final_loc.position;
            final_loc.variant = "";
          } else {
            final_loc.position = 0;
            final_loc.length = ref_len;
            final_loc.variant = "";
          }
          var_locations[i].push_back(final_loc);
        }
        chrom = prev_chrom;
        return true;
      } else if (prev_start_ > start) {
        std::cerr << "Error: VCF is unsorted at chromosome " << chrom << ", position: " << start << std::endl;
        exit(1);
      }
      started_loop = true;
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int ngt = bcf_get_genotypes(vcf_hdr, record, &gt_arr, &ngt_arr);
      if ( ngt<=0 ) continue; // GT not present
      int max_ploidy = ngt/nsmpl;
      for (int i = 0; i < nsmpl; i++) {
        if (i >= sample_names.size() || sample_names[i].empty()) continue; // Don't care about this particular sample
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
          // TODO: what about indels? non-diploid?
          // We set . to REF, aka 0 (e.g. ./. = 0/0)
          allele_1 = alleles[allele_indices[0] != -1 ? allele_indices[0]  : 0];
          allele_2 = alleles[allele_indices[1] != -1 ? allele_indices[1] : 0];
          // For SNPs, if FASTA contains a non-ATCG base AND the genotype is 0 (aka REF) or -1 (aka .), just use that non-ATCG base
          char s = toupper(ref_seq[start+i]);
          if (allele_1.length() == 1 && allele_indices[0] <= 0 && s != 'A' && s != 'T' && s != 'C' && s != 'G') allele_1[0] = ref_seq[start+i];
          if (allele_2.length() == 1 && allele_indices[1] <= 0 && s != 'A' && s != 'T' && s != 'C' && s != 'G') allele_2[0] = ref_seq[start+i];
          VarLocation loc_1;
          loc_1.position = prev_start[0];
          loc_1.length = start-prev_start[0];
          loc_1.variant = allele_1;
          VarLocation loc_2;
          loc_2.position = prev_start[1];
          loc_2.length = start-prev_start[1];
          loc_2.variant = allele_2;
          if (!indel && !(allele_1.length() == 1 && allele_2.length() == 1)) continue; // Skip because not a SNP and indel mode not active
          if (allele_1 == allele_2) continue; // Skip because alleles are identical
          var_locations[0].push_back(loc_1);
          var_locations[1].push_back(loc_2);
          prev_start[0] = start + allele_1.length();
          prev_start[1] = start + allele_2.length();
        } else {
          
        }
      }
      free(gt_arr);
    }
    // if started_loop and we actually have stuff (TODO: How to determine if we actually have stuff?)
    // Output everything stored so far
    for (int i = 0; i < (diploid ? 2 : sample_names.size()); i++) {
      if (!diploid && sample_names[i].empty()) continue; // Don't care about this sample
      // Include final stretch of sequences (i.e. the stuff after the last variant) into
      VarLocation final_loc;
      if (var_locations[i].size() != 0) {
        auto &second_to_final_loc = var_locations[i][var_locations[i].size()-1];
        final_loc.position = second_to_final_loc.position + second_to_final_loc.length + second_to_final_loc.variant.length();
        final_loc.length = ref_len - final_loc.position;
        final_loc.variant = "";
      } else {
        final_loc.position = 0;
        final_loc.length = ref_len;
        final_loc.variant = "";
      }
      var_locations[i].push_back(final_loc);
    }
    return false;
  }
  
  void prepareFastaAndPrintChromosome(
      char* ref_seq,
      std::vector<std::vector<VarLocation>>& var_locations,
      bool diploid,
      const std::string& prev_chrom,
      int ref_len,
      const std::vector<std::string>& sample_names,
      std::ostream& out_fasta,
      int fasta_line_length
  ) {
    for (int i = 0; i < var_locations.size(); i++) {
      std::string prefix = "";
      std::string suffix = "";

      if (!diploid && sample_names[i].empty()) continue;
      
      if (diploid) suffix = std::string("_") + std::string(i == 0 ? "L" : "R") + std::string(" ") + prev_chrom + std::string(":1-") + std::to_string(ref_len);
      else prefix = sample_names[i];
      
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
  
  void lift_over_gtf(const string& ref_gtf, const string& out_gtf, const vector<pair<int, int>>& mod_ranges) {
    ifstream ref_gtf_file(ref_gtf);
    ofstream out_gtf_file(out_gtf);
    
    string line;
    while (getline(ref_gtf_file, line)) {
      istringstream iss(line);
      GTFRecord record;
      iss >> record.seqname >> record.source >> record.feature >> record.start >> record.end >> record.score >> record.strand >> record.frame >> record.attributes;
      
      int start = record.start;
      int end = record.end;
      
      for (const auto& range : mod_ranges) {
        if (end <= range.first) {
          continue;
        } else if (start >= range.second) {
          start += range.second - range.first;
          end += range.second - range.first;
        } else {
          if (start < range.first) {
            end += range.second - range.first;
          } else {
            end = range.second + (end - range.first);
          }
          start = range.first;
        }
      }
      
      record.start = start;
      record.end = end;
      
      out_gtf_file << record.seqname << "\t" << record.source << "\t" << record.feature << "\t" << record.start << "\t" << record.end << "\t" << record.score << "\t" << record.strand << "\t" << record.frame << "\t" << record.attributes << endl;
    }
    
    ref_gtf_file.close();
    out_gtf_file.close();
  }
  
  struct GTFRecord {
    std::string seqname;
    std::string source;
    std::string feature;
    int start;
    int end;
    double score;
    char strand;
    std::string frame;
    std::string attributes;
  };
  
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

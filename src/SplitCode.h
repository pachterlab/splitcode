#ifndef SPLITCODE_H
#define SPLITCODE_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "hash.hpp"

struct SplitCode {
  SplitCode() {
    
  }
  
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& IDs, const std::vector<std::vector<int>>& newIDs) {
    
  }
  
  struct SortedVectorHasher {
    size_t operator()(const std::vector<int>& v) const {
      uint64_t r = 0;
      int i=0;
      for (auto x : v) {
        uint64_t t;
        MurmurHash3_x64_64(&x,sizeof(x), 0,&t);
        t = (x>>i) | (x<<(64-i));
        r = r ^ t;
        i = (i+1)%64;
      }
      return r;
    }
  };
  
  struct SplitCodeTag {
    bool initiator;
    bool terminator;
    std::string name;
    int16_t file;
    int32_t pos_start;
    int32_t pos_end;
    bool discard_read_if_not_present;
    bool not_include_in_barcode;
    bool is_original_sequence;
  };

  
  void generate_hamming_mismatches(std::string seq, int dist, std::vector<std::string>& results, std::vector<size_t> pos = std::vector<size_t>()) {
    if (dist == 0) {
      return;
    }
    size_t bc = seq.length();
    for (size_t i = 0; i < bc; ++i) {
      if (std::find(pos.begin(), pos.end(),i)==pos.end()) {
        char bases[] = {'A','T','C','G'};
        for (int d = 0; d < 4; d++) {
          if (seq[i] != bases[d]) {
            std::string y = seq;
            y[i] = bases[d];
            results.push_back(y);
            pos.push_back(i);
            generate_hamming_mismatches(y, dist-1, results, pos);
          }
        }
      }
    }
  }
  
  bool addTag(std::string seq, std::string name, uint8_t mismatch_dist, 
              int16_t file, int32_t pos_start, int32_t pos_end,
              bool discard_read_if_not_present, bool not_include_in_barcode) {
    SplitCodeTag new_tag;
    new_tag.initiator = false;
    new_tag.terminator = false;
    new_tag.is_original_sequence = false;
    
    if (seq.length() > 0 && seq[0] == '*') {
      new_tag.initiator = true;
      seq.erase(0,1);
    }
    if (seq.length() > 0 && seq[seq.size()-1] == '*') {
      new_tag.terminator = true;
      seq.erase(seq.end()-1);
    }
    if (seq.length() == 0) {
      std::cerr << "Error: Sequence: " << name << " is empty" << std::endl;
      return false;
    }
    if (seq.length() > 32) {
      std::cerr << "Error: Sequence: " << name << " cannot be longer than 32 bp's" << std::endl;
      return false;
    }
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    for (int i = 0; i < seq.size(); i++) {
      if (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C' && seq[i] != 'G') {
        std::cerr << "Error: Sequence: " << name << " contains a non-ATCG character" << std::endl;
        return false;
      }
    }

    new_tag.name = name;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.discard_read_if_not_present = discard_read_if_not_present;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    if (tags.find(seq) != tags.end()) {
      std::cerr << "Error: Sequence: " << name << " collides with sequence: " << tags[seq].name << std::endl;
      return false;
    }
    for (auto& it: tags) {
      if (name == it.second.name) {
        std::cerr << "Error: Sequence name: " << name << " is present more than once" << std::endl;
        return false;
      }
    }

    std::vector<std::string> mismatches;
    generate_hamming_mismatches(seq, mismatch_dist, mismatches);
    for (std::string mismatch_seq : mismatches) {
      if (tags.find(mismatch_seq) != tags.end()) {
        if (tags[mismatch_seq].is_original_sequence) {
          std::cerr << "Error: Sequence: " << name << " collides with sequence: " << tags[mismatch_seq].name << std::endl;
          return false;
        } else {
          tags.erase(mismatch_seq);
        }
      }
      tags.insert({mismatch_seq,new_tag});
    }
    
    new_tag.is_original_sequence = true;
    tags.insert({seq,new_tag});
    return true;
  }
  
  void addTags(/*file*/) {
    //for loop: addTag(...)
  }
  
  int getNumTags() {
    std::vector<std::string> names;
    for (auto& it : tags)
      names.push_back(it.second.name);
    std::sort(names.begin(), names.end());
    return std::unique(names.begin(), names.end()) - names.begin();
  }
  
  std::unordered_map<std::string, SplitCodeTag> tags;
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
};


#endif // SPLITCODE_H

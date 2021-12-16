#ifndef SPLITCODE_H
#define SPLITCODE_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

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
    uint64_t seq_hash;
    std::string name;
    uint8_t hamming;
    int16_t file;
    int32_t pos_start;
    int32_t pos_end;
    bool discard_read_if_not_present;
    bool not_include_in_barcode;
  };
  
  uint64_t stringToBinary(const char* s, const size_t len) {
    uint64_t r = 0;
    int numN = 0;
    size_t k = len;
    if (k > 32) {
      k = 32;
    }
    for (size_t i = 0; i < k; ++i) {
      uint64_t x = ((*s) & 4) >> 1;
      if (((*s) & 3) == 2) {
        ++numN;
      }
      r = r << 2;
      r |= (x + ((x ^ (*s & 2)) >>1));
      s++;
    }
    return r;
  }
  
  bool addTag(std::string seq, std::string name, uint8_t hamming, 
              int16_t file, int32_t pos_start, int32_t pos_end,
              bool discard_read_if_not_present, bool not_include_in_barcode) {
    SplitCodeTag new_tag;
    new_tag.initiator = false;
    new_tag.terminator = false;
    
    if (seq[0] == '*') {
      new_tag.terminator = true;
      seq.erase(0,1);
    } else if (seq[seq.size()-1] == '*') {
      new_tag.initiator = true;
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
    
    uint64_t seq_hash = stringToBinary(seq.c_str(), seq.length());
    
    new_tag.name = name;
    new_tag.seq_hash = seq_hash;
    new_tag.hamming = hamming;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.discard_read_if_not_present = discard_read_if_not_present;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    if (tags.find(seq_hash) != tags.end()) {
      std::cerr << "Error: Sequence: " << name << " is a duplicate" << std::endl;
      return false;
    }
    for (auto& it: tags) {
      if (name == it.second.name) {
        std::cerr << "Error: Sequence name: " << name << " is present more than once" << std::endl;
        return false;
      }
    }

    // TODO: Hamming distance computations (and check for collisions)
    tags.insert({seq_hash,new_tag});
    return true;
  }
  
  void addTags(/*file*/) {
    //for loop: addTag(...)
  }
  
  std::unordered_map<uint64_t, SplitCodeTag> tags;
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
};


#endif // SPLITCODE_H

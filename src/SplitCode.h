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
    uint64_t seq_hash;
    std::string name;
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
  
  std::string binaryToString(uint64_t x, size_t len) {
    std::string s(len, 'N');
    size_t sh = len-1;
    for (size_t i = 0; i < len; i++) {
      char c = 'N';
      switch((x >> (2*sh)) & 0x03ULL) {
      case 0x00: c = 'A'; break;
      case 0x01: c = 'C'; break;
      case 0x02: c = 'G'; break;
      case 0x03: c = 'T'; break;
      }
      sh--;
      s.at(i) = c;
    }
    return std::move(s);
  }

  
  void generate_hamming_mismatches(std::string seq, int dist, std::vector<uint64_t>& results, std::vector<size_t> pos = std::vector<size_t>()) {
    if (dist == 0) {
      return;
    }
    const char *s = seq.c_str();
    size_t bc = seq.length();
    size_t sh = bc - 1;
    uint64_t b = stringToBinary(s, bc);
    for (size_t i = 0; i < bc; ++i) {
      if (std::find(pos.begin(), pos.end(),i)==pos.end()) {
        for (uint64_t d = 1; d <= 3; d++) {
          uint64_t y = b ^ (d << (2 * sh));
          results.push_back(y);
          pos.push_back(i);
          generate_hamming_mismatches(binaryToString(y,bc), dist-1, results, pos);
        }
      }
      sh--;
    }
  }
  
  bool addTag(std::string seq, std::string name, uint8_t mismatch_dist, 
              int16_t file, int32_t pos_start, int32_t pos_end,
              bool discard_read_if_not_present, bool not_include_in_barcode) {
    SplitCodeTag new_tag;
    new_tag.initiator = false;
    new_tag.terminator = false;
    
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
    
    uint64_t seq_hash = stringToBinary(seq.c_str(), seq.length());
    
    new_tag.name = name;
    new_tag.seq_hash = seq_hash;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.discard_read_if_not_present = discard_read_if_not_present;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    if (tags.find(seq_hash) != tags.end()) {
      std::cerr << "Error: Sequence: " << name << " collides with sequence: " << tags[seq_hash].name << std::endl;
      return false;
    }
    for (auto& it: tags) {
      if (name == it.second.name) {
        std::cerr << "Error: Sequence name: " << name << " is present more than once" << std::endl;
        return false;
      }
    }

    std::vector<uint64_t> mismatches;
    generate_hamming_mismatches(seq, mismatch_dist, mismatches);
    for (int mismatch_seq_hash : mismatches) {
      if (tags.find(mismatch_seq_hash) != tags.end()) {
        std::cerr << "Error: Sequence: " << name << " collides with sequence: " << tags[mismatch_seq_hash].name << std::endl;
        return false;
      }
      tags.insert({mismatch_seq_hash,new_tag});
    }
    
    tags.insert({seq_hash,new_tag});
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
  
  std::unordered_map<uint64_t, SplitCodeTag> tags;
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
};


#endif // SPLITCODE_H

#ifndef SPLITCODE_H
#define SPLITCODE_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <sstream>

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
    std::string seq;
    int16_t file;
    int32_t pos_start;
    int32_t pos_end;
    bool discard_read_if_not_present;
    bool not_include_in_barcode;
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
    uint32_t new_tag_index = tags_vec.size();

    if (seq.length() > 0 && seq[0] == '*') {
      new_tag.initiator = true;
      seq.erase(0,1);
    }
    if (seq.length() > 0 && seq[seq.size()-1] == '*') {
      new_tag.terminator = true;
      seq.erase(seq.end()-1);
    }
    if (seq.length() == 0) {
      std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" is empty" << std::endl;
      return false;
    }
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    for (int i = 0; i < seq.size(); i++) {
      if (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C' && seq[i] != 'G') {
        std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" contains a non-ATCG character" << std::endl;
        return false;
      }
    }

    new_tag.seq = seq;
    new_tag.name = name;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.discard_read_if_not_present = discard_read_if_not_present;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    if (tags.find(seq) != tags.end()) { // If we've seen that sequence before
      const auto& v = tags[seq];
      uint32_t i;
      if (checkCollision(new_tag, v, i)) {
        std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << tags_vec[i].name << "\"" << std::endl;
        return false;
      }
    }

    std::vector<std::string> mismatches;
    generate_hamming_mismatches(seq, mismatch_dist, mismatches);
    for (std::string mismatch_seq : mismatches) {
      if (tags.find(mismatch_seq) != tags.end()) { // If we've seen that sequence (mismatch_seq) before
        auto& v = tags[mismatch_seq];
        uint32_t i;
        bool collision = checkCollision(new_tag, v, i);
        if (collision && tags_vec[i].seq == mismatch_seq) { // If there is a collision AND the sequence seen before is an original (i.e. non-mismatched-generated) sequence
          std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << tags_vec[i].name << "\"" << std::endl;
          return false;
        } else if (collision) { // Deal with collisions between mismatch-generated sequences from different original barcodes
          v.erase(std::remove(v.begin(), v.end(), i), v.end()); // Remove the previously found sequence from the vector stored in the map
          if (v.size() == 0) { // If the vector stored in the map is now empty, remove the entry from the map
            tags.erase(mismatch_seq);
          }
          continue; // Because of collision, don't want to update map
        }
      }
      addToMap(mismatch_seq, new_tag_index);
    }
    
    tags_vec.push_back(new_tag);
    addToMap(seq, new_tag_index);
    
    return true;
  }
  
  bool checkCollision(const SplitCodeTag& tag, const std::vector<uint32_t>& v, uint32_t& i) {
    for (auto x : v) {
      if (1==1) { // TODO: Replace with collision checks
        i = 0; // TODO: ^see above
        return true;
      }
    }
    return false;
  }
  
  void addToMap(const std::string& seq, uint32_t index) {
    if (tags.find(seq) != tags.end()) {
      auto& v = tags[seq];
      v.reserve(v.size()+1);
      v.push_back(index);
    } else {
      std::vector<uint32_t> v(0);
      v.reserve(1);
      v.push_back(index);
      tags.insert({seq,v});
    }
  }
  
  void addTags(/*file*/) {
    //for loop: addTag(...)
  }
  
  SplitCodeTag getTag(std::string& seq) { // TODO: Specify parameters (e.g. file number, location)
    return tags_vec[tags[seq][0]];
  }
  
  int getNumTags() {
    return tags_vec.size();
  }
  
  int getMapSize(bool unique = true) {
    if (unique) {
      return tags.size();
    } else {
      size_t map_size = 0;
      for (auto it : tags) {
        map_size += it.second.size();
      }
      return map_size;
    }
  }
  
  int getNumMapped() {
    int nummapped = 0;
    for (auto& n : idcount) {
      nummapped += n;
    }
    return nummapped;
  }
  
  static bool parseLocation(const std::string& location, int16_t& file, int32_t& pos_start, int32_t& pos_end, int nFiles = -1) {
    file = -1;
    pos_start = 0;
    pos_end = 0;
    if (location.empty()) {
      return true;
    }
    char delimeter = ':';
    if (location.find(',') < location.length()) {
      delimeter = ','; // If string contains commas, use commas as delimeter
    }
    std::stringstream ss_loc(location);
    std::string location_attribute;
    int i = 0;
    try {
      while (std::getline(ss_loc, location_attribute, delimeter)) {
        if (!location_attribute.empty()) {
          if (i == 0) {
            file = std::stoi(location_attribute);
          } else if (i == 1) {
            pos_start = std::stoi(location_attribute);
          } else if (i == 2) {
            pos_end = std::stoi(location_attribute);
          }
        }
        i++;
      }
      if (i > 3 || file < -1 || (file >= nFiles && nFiles != -1) || pos_start < 0 || pos_end < 0 || (pos_end <= pos_start && pos_end != 0)) {
        std::cerr << "Error: --locations is malformed; unable to parse \"" << location << "\"" << std::endl;
        return false;
      }
    } catch (std::invalid_argument &e) {
      std::cerr << "Error: Could not convert \"" << location_attribute << "\" to int in --locations" << std::endl;
      return false;
    }
    return true;
  }
  
  std::vector<SplitCodeTag> tags_vec;
  std::unordered_map<std::string, std::vector<uint32_t>> tags;
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
};


#endif // SPLITCODE_H

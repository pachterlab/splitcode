#ifndef SPLITCODE_H
#define SPLITCODE_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <sstream>

#include "hash.hpp"

struct SplitCode {
  SplitCode() {
    init = false;
  }
  
  void checkInit() { // Initialize if necessary (once initialized, can't add any more barcode tags)
    if (init) {
      return;
    }
    // Trim out some tags to create the final 'index'
    for (auto x : tags_to_remove) {
      auto& v = tags[x.first];
      v.erase(std::remove(v.begin(), v.end(), x.second), v.end()); // Remove the sequence from the vector stored in the map
      if (v.size() == 0) { // If the vector stored in the map is now empty, remove the entry from the map
        tags.erase(x.first);
      }
    }
    tags_to_remove.clear();
    init = true;
  }
  
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& IDs, const std::vector<std::vector<int>>& newIDs) {
    checkInit();
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
    uint16_t max_finds;
    uint16_t min_finds;
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
  
  bool addTag(std::string seq, std::string name, uint16_t mismatch_dist, 
              int16_t file, int32_t pos_start, int32_t pos_end,
              uint16_t max_finds, uint16_t min_finds, bool not_include_in_barcode) {
    if (init) {
      std::cerr << "Error: Already initialized" << std::endl;
      return false;
    }

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
    if (pos_end != 0 && pos_end - pos_start < seq.length()) {
      std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" is too long to fit in the supplied location" << std::endl;
      return false;
    }
    if (max_finds == 0) { // 0 = no restrictions
      max_finds = -1; // max unsigned int
    }
    if (max_finds < min_finds) {
      std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" -- max finds cannot be less than min finds" << std::endl;
      return false;
    }

    new_tag.seq = seq;
    new_tag.name = name;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.max_finds = max_finds;
    new_tag.min_finds = min_finds;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    if (tags.find(seq) != tags.end()) { // If we've seen that sequence before
      const auto& v = tags[seq];
      std::vector<uint32_t> vi;
      if (checkCollision(new_tag, v, vi)) {
        auto i = vi[0];
        std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << tags_vec[i].name << "\"" << std::endl;
        return false;
      }
    }

    std::vector<std::string> mismatches;
    generate_hamming_mismatches(seq, mismatch_dist, mismatches);
    for (std::string mismatch_seq : mismatches) {
      if (tags.find(mismatch_seq) != tags.end()) { // If we've seen that sequence (mismatch_seq) before
        auto& v = tags[mismatch_seq];
        std::vector<uint32_t> vi;
        bool collision = checkCollision(new_tag, v, vi);
        if (collision) {
          for (auto i : vi) {
            if (collision && tags_vec[i].seq == mismatch_seq) { // If there is a collision AND the sequence seen before is an original (i.e. non-mismatched-generated) sequence
              std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << tags_vec[i].name << "\"" << std::endl;
              return false;
            }
            tags_to_remove.insert(std::make_pair(mismatch_seq, i)); // Mark for removal
          }
          tags_to_remove.insert(std::make_pair(mismatch_seq, new_tag_index)); // Mark new_tag for removal (we remove later rather than now to account for future addTag(...) collision)
        }
      }
      addToMap(mismatch_seq, new_tag_index);
    }
    
    tags_vec.push_back(new_tag);
    addToMap(seq, new_tag_index);
    
    return true;
  }
  
  bool checkCollision(const SplitCodeTag& tag, const std::vector<uint32_t>& v, std::vector<uint32_t>& vi, bool fill_vi=true) {
    bool ret = false;
    for (auto x : v) {
      if (overlapRegion(tag, tags_vec[x])) {
        if (!fill_vi) { // If fill_vi is false, no need to check for every collision
          return true;
        }
        vi.push_back(x); // Add the indices of every tag in vector v that collided with the supplied tag to vector vi
        ret = true;
      }
    }
    return ret;
  }
  
  bool checkCollision(const SplitCodeTag& tag, const std::vector<uint32_t>& v) {
    std::vector<uint32_t> vi;
    return checkCollision(tag, v, vi, false);
  }
  
  bool overlapRegion(int16_t file_1, int32_t pos_start_1, int32_t pos_end_1, int16_t file_2, int32_t pos_start_2, int32_t pos_end_2) {
    if (file_1 != file_2 && !(file_1 == -1 || file_2 == -1)) {
      return false;
    }
    if (pos_start_1 < pos_start_2 && pos_end_1 <= pos_start_2 && pos_end_1 != 0) {
      return false;
    }
    if (pos_start_2 < pos_start_1 && pos_end_2 <= pos_start_1 && pos_end_2 != 0) {
      return false;
    }
    return true;
  }
  
  bool overlapRegion(const SplitCodeTag& tag1, const SplitCodeTag& tag2) {
    return overlapRegion(tag1.file, tag1.pos_start, tag1.pos_end, tag2.file, tag2.pos_start, tag2.pos_end);
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
    checkInit();
  }
  
  SplitCodeTag getTag(std::string& seq) { // TODO: Specify parameters (e.g. file number, location)
    checkInit();
    return tags_vec[tags[seq][0]];
  }
  
  int getNumTags() {
    return tags_vec.size();
  }
  
  int getMapSize(bool unique = true) {
    checkInit();
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
        std::cerr << "Error: Location string is malformed; unable to parse \"" << location << "\"" << std::endl;
        return false;
      }
    } catch (std::invalid_argument &e) {
      std::cerr << "Error: Could not convert \"" << location_attribute << "\" to int in location string" << std::endl;
      return false;
    }
    return true;
  }
  
  std::vector<SplitCodeTag> tags_vec;
  std::unordered_map<std::string, std::vector<uint32_t>> tags;
  std::set<std::pair<std::string, uint32_t>> tags_to_remove;
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
  
  bool init;
};


#endif // SPLITCODE_H

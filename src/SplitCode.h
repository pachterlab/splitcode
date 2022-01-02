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
#include <fstream>
#include <sys/stat.h>

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
  
  void generate_indels(std::string seq, int dist, std::unordered_map<std::string,int>& results, std::string original = "") {
    // Much less elegant than generate_hamming_mismatches(...) but oh well... 
    // Note: results contains the modified string and how many remaining modifications could be applied
    if (dist == 0) {
      return;
    }
    if (original.empty()) {
      original = seq;
    }
    size_t bc = seq.length();
    for (size_t i = 0; i <= bc; ++i) {
      char bases[] = {'A','T','C','G'};
      // Insertions: 
      for (int d = 0; d < 4; d++) {
        std::string y = seq;
        y.insert(i,1,bases[d]);
        if (y != original && (results.find(y) == results.end() || results[y] < dist-1)) {
          results[y] = dist-1;
          generate_indels(y, dist-1, results, original);
        }
      }
      // Deletions:
      if (i < bc) {
        std::string y = seq;
        y.erase(i,1);
        if (y != original && (results.find(y) == results.end() || results[y] < dist-1)) {
          results[y] = dist-1;
          generate_indels(y, dist-1, results, original);
        }
      }
    }
  }
  
  void generate_indels_hamming_mismatches(std::string seq, int mismatch_dist, int indel_dist, int total_dist, std::vector<std::string>& results) {
    mismatch_dist = std::min(mismatch_dist, total_dist);
    indel_dist = std::min(indel_dist, total_dist);
    if (indel_dist == 0) { // Handle hamming mismatches only
      generate_hamming_mismatches(seq, mismatch_dist, results);
      return;
    }
    std::unordered_map<std::string,int> indel_results;
    generate_indels(seq, indel_dist, indel_results);
    generate_hamming_mismatches(seq, mismatch_dist, results);
    for (auto r : indel_results) {
      results.push_back(r.first);
      int indels_dist_used = indel_dist - r.second;
      generate_hamming_mismatches(r.first, std::min(total_dist - indels_dist_used, mismatch_dist), results);
    }
    // Remove duplicates from vector:
    std::set<std::string> s;
    for(std::string r : results) {
      s.insert(r);
    }
    results.assign(s.begin(), s.end());
  }
  
  bool addTag(std::string seq, std::string name, int mismatch_dist, int indel_dist, int total_dist,
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
    generate_indels_hamming_mismatches(seq, mismatch_dist, indel_dist, total_dist, mismatches);
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
  
  bool addTags(std::string config_file) {
    if (init) {
      std::cerr << "Error: Already initialized" << std::endl;
      return false;
    }
    struct stat stFileInfo;
    auto intstat = stat(config_file.c_str(), &stFileInfo);
    if (intstat != 0) {
      std::cerr << "Error: file not found " << config_file << std::endl;
      return false;
    }
    std::ifstream cfile(config_file);
    std::string line;
    bool header_read = false;
    std::vector<std::string> h;
    while (std::getline(cfile,line)) {
      if (line.size() == 0) {
        continue;
      }
      if (line[0] == '#') {
        continue;
      }
      std::stringstream ss(line);
      std::string field;
      if (!header_read) {
        while (ss >> field) {
          std::transform(field.begin(), field.end(), field.begin(), ::toupper);
          h.push_back(field);
        }
        if (std::find(h.begin(), h.end(), "BARCODES") == h.end()) {
          std::cerr << "Error: The file \"" << config_file << "\" must contain a header with, minimally, a column header named barcodes" << std::endl;
          return false;
        }
        if (std::set<std::string>(h.begin(), h.end()).size() != h.size()) {
          std::cerr << "Error: The file \"" << config_file << "\" has a header with duplicate column names" << std::endl;
          return false;
        }
        header_read = true;
        continue;
      }
      std::string bc = "";
      std::string name = "";
      int mismatch, indel, total_dist;
      parseDistance("", mismatch, indel, total_dist); // Set up default values
      int16_t file;
      int32_t pos_start;
      int32_t pos_end;
      parseLocation("", file, pos_start, pos_end); // Set up default values
      uint16_t max_finds = 0;
      uint16_t min_finds = 0;
      bool exclude = false;
      bool ret = true;
      for (int i = 0; ss >> field; i++) {
        if (h[i] == "BARCODES") {
          bc = field;
        } else if (h[i] == "DISTANCES") {
          ret = ret && parseDistance(field, mismatch, indel, total_dist);
        } else if (h[i] == "LOCATIONS") {
          ret = ret && parseLocation(field, file, pos_start, pos_end);
        } else if (h[i] == "IDS") {
          name = field;
        } else if (h[i] == "MINFINDS") {
          std::stringstream(field) >> min_finds;
        } else if (h[i] == "MAXFINDS") {
          std::stringstream(field) >> max_finds;
        } else if (h[i] == "EXCLUDE") {
          std::stringstream(field) >> exclude;
        } else {
          std::cerr << "Error: The file \"" << config_file << "\" contains the invalid column header: " << h[i] << std::endl;
          return false;
        }
      }
      if (!ret || !addTag(bc, name.empty() ? bc : name, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude)) {
        std::cerr << "Error: The file \"" << config_file << "\" contains an error" << std::endl;
        return false;
      }
    }
    checkInit();
    return true;
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
  
  static bool parseDistance(const std::string& distance, int& mismatch, int& indel, int& total_dist) {
    mismatch = 0;
    indel = 0;
    total_dist = 0;
    if (distance.empty()) {
      return true;
    }
    char delimeter = ':';
    std::stringstream ss_dist(distance);
    std::string dist_attribute;
    int i = 0;
    try {
      while (std::getline(ss_dist, dist_attribute, delimeter)) {
        if (!dist_attribute.empty()) {
          if (i == 0) {
            mismatch = std::stoi(dist_attribute);
          } else if (i == 1) {
            indel = std::stoi(dist_attribute);
          } else if (i == 2) {
            total_dist = std::stoi(dist_attribute);
          }
        }
        i++;
      }
      if (i > 3 || mismatch < 0 || indel < 0 || total_dist < 0) {
        std::cerr << "Error: Distance string is malformed; unable to parse \"" << distance << "\"" << std::endl;
        return false;
      } else if (total_dist != 0 && (mismatch + indel < total_dist || mismatch > total_dist || indel > total_dist)) {
        std::cerr << "Error: Distance string is invalid: \"" << distance << "\"" << std::endl;
        return false;
      }
      if (total_dist == 0) {
        total_dist = mismatch + indel;
      }
    } catch (std::invalid_argument &e) {
      std::cerr << "Error: Could not convert \"" << dist_attribute << "\" to int in distance string" << std::endl;
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

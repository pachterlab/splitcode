#ifndef SPLITCODE_H
#define SPLITCODE_H

#define SPLITCODE_VERSION "0.10.0"

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <limits>
#include <random>

struct SplitCode {
  typedef std::pair<uint32_t,short> tval; // first element of pair is tag id, second is mismatch distance
  
  SplitCode() {
    init = false;
    discard_check = false;
    keep_check = false;
    setNFiles(0);
  }
  
  SplitCode(int nFiles) {
    init = false;
    discard_check = false;
    keep_check = false;
    setNFiles(nFiles);
  }
  
  void setNFiles(int nFiles) {
    this->nFiles = nFiles;
  }
  
  void checkInit() { // Initialize if necessary (once initialized, can't add any more barcode tags)
    if (init) {
      return;
    }
    if (nFiles <= 0) {
      std::cerr << "Error: nFiles must be set to a positive integer" << std::endl;
      exit(1);
    }
    // Trim out some tags to create the final 'index'
    for (auto x : tags_to_remove) {
      auto& v = tags[x.first];
      auto it = v.begin();
      while (it != v.end()) {
        if ((*it).first == x.second) {
          it = v.erase(it); // Remove the sequence from the vector stored in the map
        } else {
          ++it;
        }
      }
      if (v.size() == 0) { // If the vector stored in the map is now empty, remove the entry from the map
        tags.erase(x.first);
      }
    }
    tags_to_remove.clear();
    // Fill in k-mer sizes by location (e.g. search for k-mers of length n at positions a-b in file c):
    // (we have to be sure to merge overlapping intervals and having intervals in sorted order which is what most of what the code below does)
    int POS_MAX = std::numeric_limits<std::int32_t>::max();
    std::vector<std::map<int,std::vector<std::pair<int,int>>>> kmer_map_vec; // key = k-mer size, value = vector of position intervals; vector = one map for each file
    for (auto x : tags) {
      int kmer_size = x.first.length();
      for (auto y : x.second) {
        auto tag = tags_vec[y.first];
        if (tag.pos_end == 0) {
          tag.pos_end = POS_MAX;
        }
        std::vector<int> files(0);
        if (tag.file == -1) {
          kmer_map_vec.resize(nFiles);
          for (int i = 0; i < nFiles; i++) {
            files.push_back(i);
          }
          // TODO: What if tag.file = -1? resize vector to nfiles+1, then loop through nfiles; tag.file = 0...nfiles [what about -l 0 but -N 2]; e.g. -b ATCG,ATAA -l 0,-1 -d 1 completely gets rid of intersection ATCA (which should be allowed in files other than 0; same thing if small region overlapped by large region [should be allowed outside small region])
        } else {
          kmer_map_vec.resize(std::max((int)kmer_map_vec.size(),tag.file+1));
          files.push_back(tag.file);
        }
        for (int f : files) {
          auto &kmer_map = kmer_map_vec[f];
          if (kmer_map.find(kmer_size) == kmer_map.end()) {
            kmer_map[kmer_size] = std::vector<std::pair<int,int>>(0);
            kmer_map[kmer_size].push_back(std::make_pair(tag.pos_start, tag.pos_end));
          } else {
            // Take the union of the intervals:
            auto& curr_intervals = kmer_map[kmer_size];
            std::pair<int,int> new_interval = std::make_pair(tag.pos_start, tag.pos_end);
            bool modified = false;
            bool update_vector = true;
            for (auto &interval : curr_intervals) {
              if (new_interval.second >= interval.first && new_interval.first <= interval.second) {
                if (std::min(new_interval.first, interval.first) == interval.first && std::max(new_interval.second,interval.second) == interval.second) {
                  update_vector = false;
                } else {
                  modified = true;
                  interval = std::make_pair(std::min(new_interval.first, interval.first), std::max(new_interval.second, interval.second));
                }
              }
            }
            if (!modified) {
              if (update_vector) {
                curr_intervals.push_back(new_interval);
              }
            } else { // Existing intervals were modified so let's merge all overlapping intervals in the vector
              std::stack<std::pair<int,int>> s;
              std::sort(curr_intervals.begin(), curr_intervals.end(), [](std::pair<int,int> a, std::pair<int,int> b) {return a.first < b.first; });
              s.push(curr_intervals[0]);
              int n = curr_intervals.size();
              for (int i = 1; i < n; i++) {
                auto top = s.top();
                if (top.second < curr_intervals[i].first) {
                  s.push(curr_intervals[i]);
                } else if (top.second < curr_intervals[i].second) {
                  top.second = curr_intervals[i].second;
                  s.pop();
                  s.push(top);
                }
              }
              // Convert stack to vector
              curr_intervals.clear();
              while (!s.empty()) {
                curr_intervals.push_back(s.top());
                s.pop();
              }
              std::sort(curr_intervals.begin(), curr_intervals.end(), [](std::pair<int,int> a, std::pair<int,int> b) {return a.first < b.first; });
            }
          }
        }
      }
    }
    // Transfer kmer_map_vec into kmer_size_locations (which facilitates iteration while processing fastq reads in k-mers)
    kmer_size_locations.resize(kmer_map_vec.size());
    for (int i = 0; i < kmer_map_vec.size(); i++) {
      auto kmer_map = kmer_map_vec[i];
      for (auto x : kmer_map) {
        int kmer_size = x.first;
        for (auto v : x.second) {
          int start_pos = v.first;
          int end_pos = v.second;
          while (start_pos + kmer_size <= end_pos || end_pos == POS_MAX) {
            std::pair<int,int> kmer_location = std::make_pair(kmer_size, start_pos); // first = k-mer size; second = start position of interval
            kmer_size_locations[i].push_back(kmer_location);
            // DEBUG:
            // std::cout << "file=" << i << " k=" << kmer_location.first << " pos=" << kmer_location.second << std::endl;
            if (end_pos == POS_MAX) {
              kmer_location = std::make_pair(kmer_size, -1); // -1 = progress to end of read
              kmer_size_locations[i].push_back(kmer_location);
              break;
            }
            ++start_pos;
          }
        }
      }
    }
    for (int i = 0; i < tags_vec.size(); i++) { // Set up minFinds and maxFinds
      auto& tag = tags_vec[i];
      if (tag.min_finds != 0) {
        min_finds_map[i] = tag.min_finds;
      }
      if (tag.max_finds != 0) {
        max_finds_map[i] = tag.max_finds;
      }
    }
    init = true;
  }
  
  void update(const std::vector<int>& c, const std::vector<std::vector<int>>& IDs, const std::vector<std::vector<int>>& newIDs) {
    checkInit();
  }
  
  struct VectorHasher {
    size_t operator()(const std::vector<uint32_t>& v) const {
      uint64_t r = 0;
      int i=0;
      for (auto x : v) {
        r ^= (r<<2) + (r>>1) + (x<<2) + (x|i) + i;
        i = (i+1)%64;
      }
      return r;
    }
  };
  
  struct SplitCodeTag {
    bool initiator;
    bool terminator;
    uint32_t name_id;
    std::string seq;
    int16_t file;
    int32_t pos_start;
    int32_t pos_end;
    uint16_t max_finds;
    uint16_t min_finds;
    bool not_include_in_barcode;
  };

  struct Results {
    std::vector<uint32_t> name_ids;
    std::vector<std::pair<int,int>> pos;
    int id;
  };
  
  
  void generate_hamming_mismatches(std::string seq, int dist, std::unordered_map<std::string,int>& results, std::vector<size_t> pos = std::vector<size_t>()) {
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
            if (results.find(y) == results.end() || results[y] < dist-1) {
              results[y] = dist-1;
            }
            pos.push_back(i);
            generate_hamming_mismatches(y, dist-1, results, pos);
          }
        }
      }
    }
  }
  
  void generate_indels(std::string seq, int dist, std::unordered_map<std::string,int>& results, std::string original = "") {
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
        if (y != original && y != "" && (results.find(y) == results.end() || results[y] < dist-1)) {
          results[y] = dist-1;
          generate_indels(y, dist-1, results, original);
        }
      }
    }
  }
  
  void generate_indels_hamming_mismatches(std::string seq, int mismatch_dist, int indel_dist, int total_dist, std::unordered_map<std::string,int>& results) {
    mismatch_dist = std::min(mismatch_dist, total_dist);
    indel_dist = std::min(indel_dist, total_dist);
    if (indel_dist == 0) { // Handle hamming mismatches only
      generate_hamming_mismatches(seq, mismatch_dist, results);
      return;
    }
    std::unordered_map<std::string,int> indel_results; // Contains the modified string and how many remaining modifications could be applied
    generate_indels(seq, indel_dist, indel_results);
    results = indel_results;
    generate_hamming_mismatches(seq, mismatch_dist, results);
    for (auto r : indel_results) {
      int indels_dist_used = indel_dist - r.second;
      generate_hamming_mismatches(r.first, std::min(total_dist - indels_dist_used, mismatch_dist), results);
    }
    results.erase(seq); // Remove the original sequence in case it was generated
  }
  
  bool matchSequences(const SplitCodeTag& tag, const std::string& match_seq) {
    // Returns true if one of the sequences in tag is equal to seq
    // (Remember: tag.seq can have multiple sequences separated by '/')
    char delimeter = '/';
    std::stringstream ss(tag.seq);
    std::string seq;
    while (std::getline(ss, seq, delimeter)) {
      if (seq == match_seq) {
        return true;
      }
    }
    return false;
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
    if (max_finds < min_finds && max_finds != 0) {
      std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" -- max finds cannot be less than min finds" << std::endl;
      return false;
    }
    uint32_t name_id;
    const auto& itnames = find(names.begin(), names.end(), name);
    if (itnames == names.end()) {
      name_id = names.size();
      names.push_back(name);
    } else {
      name_id = itnames - names.begin();
    }
    
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    new_tag.seq = seq;
    new_tag.name_id = name_id;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.max_finds = max_finds;
    new_tag.min_finds = min_finds;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    
    // Now deal with adding the actual sequence:
    char delimeter = '/'; // Sequence can be delimited by '/' if the user gives multiple sequences for one tag record
    std::stringstream ss(new_tag.seq);
    int num_seqs = 0;
    tags_vec.push_back(new_tag);
    while (std::getline(ss, seq, delimeter)) {
      if (seq.empty()) {
        continue;
      }
      ++num_seqs;
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
      
      if (tags.find(seq) != tags.end()) { // If we've seen that sequence before
        const auto& v = tags[seq];
        std::vector<uint32_t> vi;
        if (checkCollision(new_tag, v, vi)) {
          for (auto i : vi) {
            if (i != new_tag_index) {
              std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << names[tags_vec[i].name_id] << "\"" << std::endl;
              return false;
            }
          }
        }
      }
      
      std::unordered_map<std::string,int> mismatches;
      generate_indels_hamming_mismatches(seq, mismatch_dist, indel_dist, total_dist, mismatches);
      for (auto mm : mismatches) {
        std::string mismatch_seq = mm.first;
        int error = total_dist - mm.second; // The number of substitutions, insertions, or deletions
        if (tags.find(mismatch_seq) != tags.end()) { // If we've seen that sequence (mismatch_seq) before
          auto& v = tags[mismatch_seq];
          std::vector<uint32_t> vi;
          bool collision = checkCollision(new_tag, v, vi);
          if (collision) {
            for (auto i : vi) {
              if (i == new_tag_index) {
                continue;
              }
              if (matchSequences(tags_vec[i], mismatch_seq)) { // If there is a collision AND the sequence seen before is an original (i.e. non-mismatched-generated) sequence
                std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" collides with sequence #" << i+1 << ": \"" << names[tags_vec[i].name_id] << "\"" << std::endl;
                return false;
              }
              tags_to_remove.insert(std::make_pair(mismatch_seq, i)); // Mark for removal
              tags_to_remove.insert(std::make_pair(mismatch_seq, new_tag_index)); // Mark new_tag for removal (we remove later rather than now to account for future addTag(...) collision)
            }
          }
        }
        addToMap(mismatch_seq, new_tag_index, error);
      }
      addToMap(seq, new_tag_index);
    }
    
    if (num_seqs == 0) {
      std::cerr << "Error: Sequence #" << new_tag_index+1 << ": \"" << name << "\" is empty" << std::endl;
      return false;
    }
    
    return true;
  }
  
  bool checkCollision(const SplitCodeTag& tag, const std::vector<tval>& v, std::vector<uint32_t>& vi, bool fill_vi=true) {
    bool ret = false;
    for (auto x : v) {
      if (overlapRegion(tag, tags_vec[x.first])) {
        if (!fill_vi) { // If fill_vi is false, no need to check for every collision
          return true;
        }
        vi.push_back(x.first); // Add the indices of every tag in vector v that collided with the supplied tag to vector vi
        ret = true;
      }
    }
    return ret;
  }
  
  bool checkCollision(const SplitCodeTag& tag, const std::vector<tval>& v) {
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
  
  void addToMap(const std::string& seq, uint32_t index, int dist = 0) {
    if (tags.find(seq) != tags.end()) {
      auto& v = tags[seq];
      for (auto i : v) {
        if (i.first == index) {
          return;
        }
      }
      v.reserve(v.size()+1);
      v.push_back(std::make_pair(index,dist));
    } else {
      std::vector<tval> v(0);
      v.reserve(1);
      v.push_back(std::make_pair(index,dist));
      tags.insert({seq,v});
    }
  }
  
  bool addTags(std::string config_file) {
    if (init) {
      std::cerr << "Error: Already initialized" << std::endl;
      return false;
    }
    if (nFiles <= 0) {
      std::cerr << "Error: nFiles must be set to a positive integer" << std::endl;
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
          ret = ret && parseLocation(field, file, pos_start, pos_end, nFiles);
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
  
  bool getTag(std::string& seq, uint32_t& tag_id) { // TODO: Specify parameters (e.g. file number, location)
    checkInit();
    const auto& it = tags.find(seq);
    if (it == tags.end()) {
      return false;
    }
    tag_id = (it->second)[0].first;
    return true;
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
  
  bool addFilterList(std::string keep_file, bool discard=false) {
    struct stat stFileInfo;
    auto intstat = stat(keep_file.c_str(), &stFileInfo);
    if (intstat != 0) {
      std::cerr << "Error: file not found " << keep_file << std::endl;
      return false;
    }
    std::ifstream kfile(keep_file);
    std::string line;
    while (std::getline(kfile,line)) {
      if (line.size() == 0) {
        continue;
      }
      std::stringstream ss(line);
      std::string name;
      char delimeter = ',';
      std::vector<uint32_t> u;
      while (std::getline(ss, name, delimeter)) {
        if (line.size() == 0) {
          continue;
        }
        const auto& itnames = find(names.begin(), names.end(), name);
        if (itnames == names.end()) {
          std::cerr << "Error: File " << keep_file << " contains the name \"" << name << "\" which does not exist" << std::endl;
          return false;
        }
        u.push_back(itnames - names.begin());
      }
      auto it = idmapinv.find(u);
      auto it1 = idmapinv_keep.find(u);
      auto it2 = idmapinv_discard.find(u);
      if (it1 != idmapinv_keep.end() || it2 != idmapinv_discard.end()) {
        std::cerr << "Error: In file " << keep_file << ", the following line is duplicated: " << line << std::endl;
        return false;
      } else if (discard && it != idmapinv.end()) {
        std::cerr << "Error: In file " << keep_file << ", the following line cannot be used: " << line << std::endl;
        return false;
      }
      if (discard) {
        idmapinv_discard.insert({u,0});
        discard_check = true;
      } else {
        idmapinv_keep.insert({u,0});
        keep_check = true;
      }
    }
    return true;
  }
  
  class Locations {
  public:
    Locations(std::vector<std::pair<int,int>>& kmers, int rlen) : kmers(kmers), size(kmers.size()), rlen(rlen) {
      invalid = false;
      jump_pos = 0;
      i = -1;
      operator++();
    };
    
    void operator++() {
      if (invalid) {
        return;
      }
      int kmer_size;
      int kmer_loc;
      if (i != -1) {
        kmer_size = kmers[i].first;
        kmer_loc = kmers[i].second;
        if (kmer_loc == -1) {
          if (pos < jump_pos) {
            pos = jump_pos;
          } else {
            pos++;
          }
          if (pos+kmer_size <= rlen) {
            return;
          }
        }
      }
      i++;
      while (i < size) {
        kmer_size = kmers[i].first;
        kmer_loc = kmers[i].second;
        if (kmer_loc == -1) {
          operator++();
          return;
        }
        pos = kmer_loc;
        if (pos+kmer_size <= rlen && pos >= jump_pos) {
          return;
        }
        i++;
      }
      invalid = true;
    }
    
    std::pair<int,int> get() {
      return std::make_pair(kmers[i].first, pos);
    }
    
    bool good() {
      return !invalid;
    }
    
    void setJump(int jump = 0) {
      jump = jump <= 0 ? kmers[i].first : jump;
      if (pos+jump > jump_pos) {
        jump_pos = pos+jump;
      }
      if (jump_pos >= rlen) {
        invalid = true;
      }
    }

  private:
    const std::vector<std::pair<int,int>>& kmers;
    const int size;
    const int rlen;
    int i;
    int pos;
    int jump_pos;
    bool invalid;
  };
  
  void processRead(std::vector<const char*>& s, std::vector<int>& l, int jmax, Results& results) {
    std::mt19937 gen;
    auto min_finds = min_finds_map; // copy
    auto max_finds = max_finds_map; // copy
    int n = std::min(jmax, (int)kmer_size_locations.size());
    for (int j = 0; j < n; j++) {
      int file = j;
      int readLength = l[file];
      std::string seq(s[file], readLength);
      bool found_weird_base = false; // non-ATCG bases
      for (auto& c: seq) {
        c &= 0xDF; // Convert a/t/c/g to upper case
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
          if (!found_weird_base) {
            int hash = hashSequence(seq);
            gen.seed(hash); // seed random number generator based on hash of read sequence
            found_weird_base = true;
          }
          c = "ATCG"[(gen()%4)]; // substitute non-ATCG base for pseudo-random base
        }
      }
      auto& kmers = kmer_size_locations[file];
      for (Locations locations(kmers, readLength); locations.good(); ++locations) {
        auto loc = locations.get();
        auto k = loc.first;
        auto pos = loc.second;
        std::string kmer = seq.substr(pos, k);
        // DEBUG:
        // std::cout << "file=" << file << " k=" << k << " pos=" << pos << " kmer=" << kmer << std::endl;
        uint32_t tag_id;
        if (getTag(kmer, tag_id)) {
          auto& tag = tags_vec[tag_id];
          if (tag.min_finds > 0) {
            min_finds[tag_id]--;
          }
          if (tag.max_finds > 0) {
            if (max_finds[tag_id]-- <= 0) {
              continue; // maxFinds exceeded; just continue
            }
          }
          if (!tag.not_include_in_barcode) {
            results.name_ids.push_back(tag.name_id);
          }
          locations.setJump();
        }
      }
    }
    for (auto& it : min_finds) {
      if (it.second > 0) {
        results.name_ids.clear(); // minFinds not met
        break;
      }
    }
  }
  
  void update(std::vector<Results>& rv) {
    // Should only be called under a lock (can't have multiple threads accessing a common container)
    for (auto& r : rv) {
      auto& u = r.name_ids;
      if (u.empty()) {
        continue;
      }
      auto it = idmapinv.find(u);
      int id;
      if (keep_check && idmapinv_keep.find(u) == idmapinv_keep.end()) {
        r.name_ids.clear();
        continue;
      }
      if (it != idmapinv.end()) {
        id = it->second;
        idcount[id]++;
      } else if (discard_check && idmapinv_discard.find(u) != idmapinv_discard.end()) {
        r.name_ids.clear();
        continue;
      } else {
        id = idmapinv.size();
        idmapinv.insert({u,id});
        idmap.push_back(u);
        idcount.push_back(1);
      }
      r.id = id;
    }
  }
  
  bool isAssigned(Results& r) {
    return !r.name_ids.empty();
  }
  
  std::string getNameString(Results& r) {
    std::string names_str = "";
    for (int n : r.name_ids) {
      names_str += "[" + names[n] + "]";
    }
    return names_str;
  }
  
  void writeBarcodeMapping(std::string fname) {
    std::ofstream of;
    of.open(fname);
    if (!of.is_open()) {
      std::cerr << "Error: Couldn't open file: " << fname << std::endl;
      exit(1);
    }
    for (int i = 0; i < idmap.size(); i++) {
      auto &u = idmap[i];
      auto it = idmapinv.find(u);
      int id = it->second;
      int n = idcount[i];
      std::string barcode_str = "";
      for (auto& tag_id : u) {
        barcode_str += names[tags_vec[tag_id].name_id] + ",";
      }
      if (!barcode_str.empty()) {
        barcode_str.resize(barcode_str.size()-1);
      }
      of << binaryToString(id, FAKE_BARCODE_LEN) << "\t" << barcode_str << "\t" << n << """\n";
    }
    of.close();
  }
  
  static std::string binaryToString(uint64_t x, size_t len) {
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
  
  static uint64_t hashKmer(const std::string& key) {
    uint64_t r = 0;
    const char* s = key.c_str();
    size_t k = key.length();
    if (k > MAX_K) {
      k = MAX_K;
    }
    for (size_t i = 0; i < k; ++i) {
      uint64_t x = ((*s) & 4) >> 1;
      r = r << 2;
      r |= (x + ((x ^ (*s & 2)) >>1));
      s++;
    }
    return r;
  }
  
  static uint64_t hashSequence(const std::string& s) {
    int len = s.length();
    uint64_t hash = 0;
    int n = len / MAX_K;
    int i = 0;
    while (i < n) {
      hash ^= hashKmer(s.substr(i*MAX_K,MAX_K));
      i++;
    }
    int remaining_size = len % MAX_K;
    hash ^= hashKmer(s.substr(i*MAX_K,remaining_size));
    return hash;
  }
  
  class KmerHasher {
  public:
    size_t operator() (const std::string& key) const {
      return hashKmer(key);
    }
  };
  
  std::vector<SplitCodeTag> tags_vec;
  std::unordered_map<std::string, std::vector<tval>, KmerHasher> tags;
  std::set<std::pair<std::string, uint32_t>> tags_to_remove;
  std::vector<std::string> names;
  
  std::vector<std::vector<uint32_t>> idmap;
  std::unordered_map<std::vector<uint32_t>, int, VectorHasher> idmapinv;
  std::vector<int> idcount;
  std::unordered_map<std::vector<uint32_t>, int, VectorHasher> idmapinv_keep;
  std::unordered_map<std::vector<uint32_t>, int, VectorHasher> idmapinv_discard;
  
  std::unordered_map<uint32_t,int> min_finds_map;
  std::unordered_map<uint32_t,int> max_finds_map;
  
  std::vector<std::vector<std::pair<int,int>>> kmer_size_locations;
  
  bool init;
  bool discard_check;
  bool keep_check;
  int nFiles;
  static const int MAX_K = 32;
  static const size_t FAKE_BARCODE_LEN = 16;
  static const char QUAL = 'K';
};


#endif // SPLITCODE_H

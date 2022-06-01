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
#include <stack>
#include "robin_hood.h"

struct SplitCode {
  typedef std::pair<uint32_t,short> tval; // first element of pair is tag id, second is mismatch distance
  enum dir {left, right, nodir};
  
  SplitCode() {
    init = false;
    discard_check = false;
    keep_check = false;
    discard_check_group = false;
    keep_check_group = false;
    always_assign = false;
    random_replacement = false;
    n_tag_entries = 0;
    curr_barcode_mapping_i = 0;
    curr_umi_id_i = 0;
    setNFiles(0);
  }
  
  SplitCode(int nFiles, bool trim_only = false, bool disable_n = true,
            std::string trim_5_str = "", std::string trim_3_str = "", std::string extract_str = "", std::string barcode_prefix = "") {
    init = false;
    discard_check = false;
    keep_check = false;
    discard_check_group = false;
    keep_check_group = false;
    n_tag_entries = 0;
    curr_barcode_mapping_i = 0;
    curr_umi_id_i = 0;
    this->trim_5_str = trim_5_str;
    this->trim_3_str = trim_3_str;
    this->extract_str = extract_str;
    this->barcode_prefix = barcode_prefix;
    setNFiles(nFiles);
    setTrimOnly(trim_only);
    setRandomReplacement(!disable_n);
  }
  
  void setNFiles(int nFiles) {
    if (init) {
      return;
    }
    this->nFiles = nFiles;
  }
  
  void setTrimOnly(bool trim_only) {
    if (init) {
      return;
    }
    this->always_assign = trim_only;
  }
  
  void setRandomReplacement(bool rand) {
    if (init) {
      return;
    }
    this->random_replacement = rand;
  }
  
  void checkInit() { // Initialize if necessary (once initialized, can't add any more barcode tags)
    if (init) {
      return;
    }
    if (nFiles <= 0) {
      std::cerr << "Error: nFiles must be set to a positive integer" << std::endl;
      exit(1);
    }
    // Process 5′/3′ end-trimming
    std::stringstream ss_trim_5(this->trim_5_str);
    std::stringstream ss_trim_3(this->trim_3_str);
    std::string trim_val;
    trim_5_3_vec.resize(nFiles, std::make_pair(0,0));
    try {
      int i = 0;
      while (std::getline(ss_trim_5, trim_val, ',')) {
        if (!trim_val.empty()) {
          trim_5_3_vec[i].first = std::max(0,std::stoi(trim_val));
        }
        i++;
      }
      if (i != nFiles && i != 0) {
        std::cerr << "Error: 5′/3′ trimming invalid; need to specify as many values as nFastqs" << std::endl;
        exit(1);
      }
      i = 0;
      while (std::getline(ss_trim_3, trim_val, ',')) {
        if (!trim_val.empty()) {
          trim_5_3_vec[i].second = std::max(0,std::stoi(trim_val));
        }
        i++;
      }
      if (i != nFiles && i != 0) {
        std::cerr << "Error: 5′/3′ trimming invalid; need to specify as many values as nFastqs" << std::endl;
        exit(1);
      }
    } catch (std::exception &e) {
      std::cerr << "Error: Could not convert \"" << trim_val << "\" to int in 5′/3′ trimming" << std::endl;
      exit(1);
    }
    // Process UMI string extraction
    if (!extract_str.empty()) {
      if (!parseExtractStr(extract_str)) {
        std::cerr << "Error: Could not parse extraction pattern: \"" << extract_str << "\"" << std::endl;
        exit(1);
      }
    }
    // Process barcode prefix
    if (!barcode_prefix.empty()) {
      std::transform(barcode_prefix.begin(), barcode_prefix.end(), barcode_prefix.begin(), ::toupper);
      if (barcode_prefix.length() > 8) {
        std::cerr << "Error: Barcode prefix is too long" << std::endl;
        exit(1);
      }
      for (int i = 0; i < barcode_prefix.size(); i++) {
        if (barcode_prefix[i] != 'A' && barcode_prefix[i] != 'T' && barcode_prefix[i] != 'C' && barcode_prefix[i] != 'G') {
          std::cerr << "Error: Barcode prefix contains a non-ATCG character" << std::endl;
          exit(1);
        }
      }
    }
    // Process before_after_vec
    for (auto& x: before_after_vec) {
      auto &tag = tags_vec[x.first];
      std::string& s = x.second.second;
      std::string name;
      uint32_t id;
      bool group = false;
      if (s[1] == '{') {
        name = s.substr(2,s.find_first_of('}')-2);
        group = true;
        const auto& itnames = std::find(group_names.begin(), group_names.end(), name);
        if (itnames == group_names.end()) {
          continue;
        } else {
          id = itnames - group_names.begin();
        }
      } else {
        name = s.substr(1,s.find_first_of('}')-1);
        const auto& itnames = std::find(names.begin(), names.end(), name);
        if (itnames == names.end()) {
          continue;
        } else {
          id = itnames - names.begin();
        }
      }
      uint16_t extra = 0;
      uint16_t extra2 = 0;
      if (s[s.length()-1] != '}') {
        std::string ssub = s.substr(s.find_last_of('}')+1);
        std::stringstream ss(ssub);
        std::string s2;
        int i = 0;
        while (std::getline(ss, s2, '-')) {
          if (i == 0) {
            extra = std::stoi(s2);
          } else if (i == 1) {
            extra2 = std::stoi(s2);
          }
          i++;
        }
      }
      if (x.second.first) {
        tag.extra_after = extra;
        tag.extra_after2 = extra2;
        tag.has_after = true;
        tag.has_after_group = group;
        tag.id_after = id;
      } else {
        tag.extra_before = extra;
        tag.extra_before2 = extra2;
        tag.has_before = true;
        tag.has_before_group = group;
        tag.id_before = id;
      }
    }
    before_after_vec.clear();
    before_after_vec.shrink_to_fit();
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
        } else {
          kmer_map_vec.resize(std::max((int)kmer_map_vec.size(),tag.file+1));
          files.push_back(tag.file);
        }
        for (int f : files) {
          auto &kmer_map = kmer_map_vec[f];
          if (kmer_map.find(kmer_size) == kmer_map.end()) {
            kmer_map[kmer_size] = std::vector<std::pair<int,int>>(0);
            kmer_map[kmer_size].push_back(std::make_pair(tag.pos_start < 0 ? 0 : tag.pos_start, tag.pos_end));
          } else {
            // Take the union of the intervals:
            auto& curr_intervals = kmer_map[kmer_size];
            std::pair<int,int> new_interval = std::make_pair(tag.pos_start < 0 ? 0 : tag.pos_start, tag.pos_end);
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
      // Sort kmer_size_locations[i] and ensure all elements are unique
      std::set<std::pair<int,int>> s(kmer_size_locations[i].begin(), kmer_size_locations[i].end());
      kmer_size_locations[i].assign(s.begin(), s.end()); // Make sure all elements are unique
      std::sort(kmer_size_locations[i].begin(), kmer_size_locations[i].end(), [](std::pair<int,int> a, std::pair<int,int> b) {
        return (a.first == b.first) ? (a.second == -1 || b.second == -1 ? a.second > b.second : a.second < b.second) : (a.first < b.first);
        });
    }
    initiator_files.resize(kmer_size_locations.size(), false);
    for (int i = 0; i < tags_vec.size(); i++) { // Set up minFinds and maxFinds and initiators
      auto& tag = tags_vec[i];
      if (tag.min_finds != 0) {
        min_finds_map[i] = tag.min_finds;
      }
      if (tag.max_finds != 0) {
        max_finds_map[i] = tag.max_finds;
      }
      if (tag.initiator && (tag.file < initiator_files.size() || tag.file == -1)) {
        if (tag.file == -1) {
          std::replace(initiator_files.begin(), initiator_files.end(), false, true); // All files have initiator sequences
        } else {
          initiator_files[tag.file] = true; // Identify which files have initiator sequences
        }
      }
    }
    // Set up expansions vectors (in a map with keys being k-mer sizes):
    struct Expansion {
      int kmer_size, kmer_size_2, start_pos, start_pos_2, file;
    };
    std::unordered_map<int,std::vector<Expansion>> emap;
    for (int i = 0; i < kmer_size_locations.size(); i++) { // Go through locations to figure out where expansions may 'potentially' occur (on the basis of location)
      int prev_kmer_size = -1;
      int prev_start_pos = -1;
      std::pair<int,int> max_loc = std::make_pair(-1,-1); // k-mer location with the rightmost non-(-1) location
      int smallest_kmer_unbound = -1; // the smallest k-mer size with a -1 location
      for (auto &loc : kmer_size_locations[i]) {
        if (smallest_kmer_unbound == -1 && loc.second == -1) {
          smallest_kmer_unbound = loc.first;
        }
        if (loc.second > max_loc.second) {
          max_loc.first = loc.first; // max_loc.first will be the smallest k if there are ties
          max_loc.second = loc.second;
        }
      }
      if (max_loc.second != -1) {
        // Extend all -1's to max_loc
        int prev_k = -1;
        int prev_pos = -1;
        for (int j = 0; j < kmer_size_locations[i].size(); j++) {
          if (kmer_size_locations[i][j].second == -1) {
            int k = kmer_size_locations[i][j].first;
            int pos = kmer_size_locations[i][j].second;
            if (pos == -1 && prev_k != -1 && k == prev_k && prev_pos != -1) {
              auto &v = kmer_size_locations[i];
              std::vector<std::pair<int,int>> new_v;
              for (int p = prev_pos+1; p <= max_loc.second; p++) {
                new_v.push_back(std::make_pair(k,p));
              }
              v.insert(v.begin()+j,new_v.begin(),new_v.end());
              j += new_v.size();
            }
          }
          prev_k = kmer_size_locations[i][j].first;
          prev_pos = kmer_size_locations[i][j].second;
        }
      }
      std::set<std::pair<int,int>> deletions;
      for (int j = 0; j < kmer_size_locations[i].size(); j++) {
        int kmer_size = kmer_size_locations[i][j].first;
        int start_pos = kmer_size_locations[i][j].second;
        int prev_kmer_size_2 = -1;
        int prev_start_pos_2 = -1;
        for (int k = j+1; k < kmer_size_locations[i].size(); k++) {
          int kmer_size_2 = kmer_size_locations[i][k].first;
          int start_pos_2 = kmer_size_locations[i][k].second;
          if (kmer_size_2 <= kmer_size) {
            continue;
          }
          int actual_start_pos_2 = start_pos_2 == -1 && prev_start_pos_2 != -1 && prev_kmer_size_2 != -1 && prev_kmer_size_2 == kmer_size_2 ? prev_start_pos_2 : start_pos_2;
          int actual_start_pos = start_pos == -1 && prev_start_pos != -1 && prev_kmer_size != -1 && prev_kmer_size == kmer_size ? prev_start_pos : start_pos;
          // If we can expand from kmer_size to kmer_size_2:
          if (((start_pos == -1 || start_pos_2 == -1) && ((actual_start_pos_2 == -1 || actual_start_pos == -1) || (actual_start_pos_2 == actual_start_pos))) 
                || (start_pos != -1 && start_pos_2 != -1 && start_pos_2 == start_pos)) {
            Expansion e;
            e.kmer_size = kmer_size;
            e.kmer_size_2 = kmer_size_2;
            e.start_pos = actual_start_pos;
            e.start_pos_2 = actual_start_pos_2;
            e.file = i;
            emap[kmer_size_2].push_back(e);
            if (!(smallest_kmer_unbound == e.kmer_size_2 && start_pos_2 == -1)) {
              deletions.insert(std::make_pair(e.kmer_size_2, start_pos_2)); // Mark location with larger k-mer for deletion (since we can expand to it)
            }
          }
          prev_kmer_size_2 = kmer_size_2;
          prev_start_pos_2 = start_pos_2;
        }
        prev_kmer_size = kmer_size;
        prev_start_pos = start_pos;
      }
      // Delete marked locations
      for (auto d : deletions) {
        kmer_size_locations[i].erase(std::remove(kmer_size_locations[i].begin(), kmer_size_locations[i].end(), d), kmer_size_locations[i].end());
      }
      // Find the -1 position (if it exists) and make sure it's preceded by max_loc.second+1
      auto it = std::find(kmer_size_locations[i].begin(), kmer_size_locations[i].end(), std::pair<int,int>(smallest_kmer_unbound, -1));
      if (it != kmer_size_locations[i].end()) {
        kmer_size_locations[i].insert(it, std::pair<int,int>(smallest_kmer_unbound, max_loc.second+1));
      }
      // Sort by location (aka the second element in the pair) rather than by k-mer size
      std::sort(kmer_size_locations[i].begin(), kmer_size_locations[i].end(), [](std::pair<int,int> a, std::pair<int,int> b) {
        return (a.second == -1 || b.second == -1 ? a.second > b.second : a.second < b.second);
      });
    }
    // For all sequences in map, decompose them into smaller substrings
    std::set<std::pair<std::string,int>> decomposed_kmers;
    for (auto& it: tags) {
      int kmer_size = it.first.length();
      if (emap.find(kmer_size) != emap.end()) { // emap[kmer_size]
        for (auto &e : emap[kmer_size]) {
          bool found = false;
          // Check if any tags associated with the current sequence are in the same location as the current expansion 
          for (auto x : it.second) {
            auto &tag = tags_vec[x.first];
            if (overlapRegion(e.file, e.start_pos_2 == -1 ? 0 : e.start_pos_2, e.start_pos_2 == -1 ? 0 : e.start_pos_2+e.kmer_size_2, tag.file, tag.pos_start, tag.pos_end)) {
              found = true;
              break;
            }
          }
          if (!found) {
            continue;
          }
          // Decompose kmer of kmer_size by substring'ing
          int decomposed_kmer_size = e.kmer_size;
          std::string s = it.first.s_;
          decomposed_kmers.insert(std::make_pair(s.substr(0, decomposed_kmer_size), kmer_size)); // to be added to tags map
        }
      }
    }
    for (auto& d : decomposed_kmers) { // Put decomposed k-mer strings into the tags map
      SeqString sstr(d.first);
      int k_expanded = d.second;
      if (tags.find(sstr) != tags.end()) {
        auto& tag_v = tags[sstr];
        if (tag_v.size() > 0 && tag_v[0].second == -1) {
          if (k_expanded < tag_v[0].first) {
            // If encounter duplicate expansions, use the one with the smaller k-mer size
            tag_v[0].first = k_expanded;
          }
        } else {
          tag_v.insert(tag_v.begin(), std::make_pair(k_expanded,-1)); // Put expansion at beginning of vector
        }
      } else { // String not previously seen in map (vector is empty)
        tags[sstr].push_back(std::make_pair(k_expanded,-1)); // Put expansion at beginning of vector
      }
    }
    // DEBUG: Print out final locations
    /*for (int i = 0; i < kmer_size_locations.size(); i++) {
      for (int j = 0; j < kmer_size_locations[i].size(); j++) {
        int kmer_size = kmer_size_locations[i][j].first;
        int start_pos = kmer_size_locations[i][j].second;
        std::cout << i << ":" << kmer_size << ":" << start_pos << std::endl;
      }
    }*/
    // DEBUG: Print out tags map
    /*for (auto& it: tags) {
      std::cout << it.first.s_ << "; k = " << it.first.length() << std::endl;
      auto &v = it.second;
      for (auto &x : v) {
        if (x.second != -1) {
          std::cout << "\t" << names[tags_vec[x.first].name_id] << " " << x.second << std::endl;
        } else {
          std::cout << "\t" << "Expansion: " << x.first << std::endl;
        }
      }
    }*/
    init = true;
  }
  
  struct VectorHasher {
    size_t operator()(const std::vector<uint32_t>& v) const {
      uint64_t r = v.size()-1;
      for (auto x : v) {
        r ^= x + 0x9e3779b9 + (r<<6) + (r>>2); // boost hash_combine method
      }
      return r;
    }
  };
  
  struct SplitCodeTag {
    bool initiator;
    bool terminator;
    uint32_t name_id;
    uint32_t group;
    std::string seq;
    int16_t file;
    int32_t pos_start;
    int32_t pos_end;
    uint16_t max_finds;
    uint16_t min_finds;
    bool not_include_in_barcode;
    dir trim;
    int trim_offset;
    bool has_before;
    bool has_after;
    bool has_before_group;
    bool has_after_group;
    uint32_t id_after;
    uint32_t id_before;
    uint16_t extra_before;
    uint16_t extra_after;
    uint16_t extra_before2;
    uint16_t extra_after2;
  };

  struct Results {
    std::vector<uint32_t> name_ids;
    std::vector<std::pair<int,std::pair<int,int>>> modtrim;
    int id;
    bool discard;
    std::string ofile;
  };
  
  struct SeqString {
    const char* p_;
    unsigned short l_;
    std::string s_;
    SeqString(const char* p, unsigned short l) : p_(p), l_(l) { }
    SeqString(const std::string& s) : p_(nullptr), l_(s.length()), s_(s) { }
    SeqString() : p_{nullptr}, s_(""), l_(0) { }
    bool operator==(const SeqString& ss) const {
      return p_ ? 
        (ss.p_ ? ss.l_ == l_ && std::strncmp(p_, ss.p_, l_) == 0 : ss.l_ == l_ && strncmp(p_, ss.s_.c_str(), l_) == 0) : 
        (ss.p_ ? ss.l_ == l_ && strncmp(ss.p_, s_.c_str(), l_) == 0 : s_ == ss.s_);
    }
    size_t length() const {
      return l_;
    }
  };
  
  struct UMI {
    uint32_t id1, id2;
    uint16_t length_range_start;
    uint16_t length_range_end;
    uint16_t padding_left;
    uint16_t padding_right;
    int16_t id;
    std::pair<int,int> location1;
    std::pair<int,int> location2;
    std::string name;
    bool group1, group2, id1_present, id2_present;
  };
  
  
  void generate_hamming_mismatches(std::string seq, int dist, std::unordered_map<std::string,int>& results, bool use_N = false, int initial_dist = 0, std::vector<size_t> pos = std::vector<size_t>()) {
    if (dist == 0) {
      return;
    }
    if (initial_dist == 0) {
      initial_dist = dist;
    }
    size_t bc = seq.length();
    for (size_t i = 0; i < bc; ++i) {
      if (std::find(pos.begin(), pos.end(),i)==pos.end()) {
        char bases[] = {'A','T','C','G','N'};
        int l = use_N ? 5 : 4;
        for (int d = 0; d < l; d++) {
          if (seq[i] != bases[d]) {
            std::string y = seq;
            y[i] = bases[d];
            if (results.find(y) == results.end() || results[y] > initial_dist-(dist-1)) {
              results[y] = initial_dist-(dist-1); // contains the hamming mismatch error
            }
            pos.push_back(i);
            generate_hamming_mismatches(y, dist-1, results, use_N, initial_dist, pos);
          }
        }
      }
    }
  }
  
  void generate_indels(std::string seq, int dist, std::unordered_map<std::string,int>& results, bool use_N = false, bool do_insertion = true, bool do_deletion = true, int initial_dist = 0) {
    if (dist == 0) {
      return;
    }
    if (initial_dist == 0) {
      initial_dist = dist;
    }
    size_t bc = seq.length();
    for (size_t i = 0; i <= bc; ++i) {
      char bases[] = {'A','T','C','G','N'};
      int l = use_N ? 5 : 4;
      // Insertions: 
      if (do_insertion) {
        for (int d = 0; d < l; d++) {
          std::string y = seq;
          y.insert(i,1,bases[d]);
          if (results.find(y) == results.end() || results[y] > initial_dist-(dist-1)) {
            results[y] = initial_dist-(dist-1);
            generate_indels(y, dist-1, results, use_N, true, false, initial_dist);
          }
        }
      }
      // Deletions:
      if (do_deletion) {
        if (i < bc) {
          std::string y = seq;
          y.erase(i,1);
          if (y != "" && (results.find(y) == results.end() || results[y] > initial_dist-(dist-1))) {
            results[y] = initial_dist-(dist-1);
            generate_indels(y, dist-1, results, use_N, false, true, initial_dist);
          }
        }
      }
    }
  }
  
  void generate_indels_hamming_mismatches(std::string seq, int mismatch_dist, int indel_dist, int total_dist, std::unordered_map<std::string,int>& results) {
    bool use_N = !random_replacement;
    mismatch_dist = std::min(mismatch_dist, total_dist);
    indel_dist = std::min(indel_dist, total_dist);
    if (indel_dist == 0) { // Handle hamming mismatches only
      generate_hamming_mismatches(seq, mismatch_dist, results, use_N);
      return;
    }
    std::unordered_map<std::string,int> indel_results; // Contains the modified string and how many remaining modifications could be applied
    generate_indels(seq, indel_dist, indel_results, use_N);
    results = indel_results;
    generate_hamming_mismatches(seq, mismatch_dist, results, use_N);
    for (auto r : indel_results) {
      int indels_dist_used = r.second;
      int dist_remaining = std::min(total_dist - indels_dist_used, mismatch_dist);
      generate_hamming_mismatches(r.first, dist_remaining, results, use_N, dist_remaining+indels_dist_used);
    }
    results.erase(seq); // Remove the original sequence in case it was generated
    // Remove all keys containing the original sequence as a substring:
    auto it = results.begin();
    while (it != results.end()) {
      if (it->first.length() > seq.length() && it->first.find(seq) != std::string::npos) {
        it = results.erase(it);
      } else {
        it++;
      }
    }
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
  
  bool addTag(std::string seq, std::string name, std::string group_name, int mismatch_dist, int indel_dist, int total_dist,
              int16_t file, int32_t pos_start, int32_t pos_end,
              uint16_t max_finds, uint16_t min_finds, bool not_include_in_barcode,
              dir trim, int trim_offset, std::string after_str, std::string before_str) {
    if (init) {
      std::cerr << "Error: Already initialized" << std::endl;
      return false;
    }

    SplitCodeTag new_tag;
    new_tag.initiator = false;
    new_tag.terminator = false;
    uint32_t new_tag_index = tags_vec.size();
    ++n_tag_entries;

    if (seq.length() > 0 && seq[0] == '*') {
      new_tag.initiator = true;
      seq.erase(0,1);
    }
    if (seq.length() > 0 && seq[seq.size()-1] == '*') {
      new_tag.terminator = true;
      seq.erase(seq.end()-1);
    }
    if (seq.length() == 0) {
      std::cerr << "Error: Sequence #" << n_tag_entries << ": \"" << name << "\" is empty" << std::endl;
      return false;
    }
    if (max_finds < min_finds && max_finds != 0) {
      std::cerr << "Error: Sequence #" << n_tag_entries << ": \"" << name << "\" -- max finds cannot be less than min finds" << std::endl;
      return false;
    }
    for (int i = 0; i < name.size(); i++) {
      if (name[i] == '#' || name[i] == '|' || name[i] == '(' || name[i] == ')' || name[i] == '[' || name[i] == ']' || name[i] == '{' || name[i] == '}' || name[i] == '@') {
        std::cerr << "Error: The name of sequence #" << n_tag_entries << ": \"" << name << "\" contains an invalid character" << std::endl;
        return false;
      }
    }
    uint32_t name_id;
    const auto& itnames = std::find(names.begin(), names.end(), name);
    if (itnames == names.end()) {
      name_id = names.size();
      names.push_back(name);
    } else {
      name_id = itnames - names.begin();
    }
    for (int i = 0; i < group_name.size(); i++) {
      if (group_name[i] == '#' || group_name[i] == '|' || group_name[i] == '(' || group_name[i] == ')' || group_name[i] == '[' || group_name[i] == ']' || group_name[i] == '{' || group_name[i] == '}' || group_name[i] == '@') {
        std::cerr << "Error: The group name: \"" << group_name << "\" contains an invalid character" << std::endl;
        return false;
      }
    }
    uint32_t group_name_id = -1;
    if (!group_name.empty()) {
      const auto& itgnames = std::find(group_names.begin(), group_names.end(), group_name);
      if (itgnames == group_names.end()) {
        group_name_id = group_names.size();
        group_names.push_back(group_name);
      } else {
        group_name_id = itgnames - group_names.begin();
      }
    }
    
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    new_tag.seq = seq;
    new_tag.name_id = name_id;
    new_tag.group = group_name_id;
    new_tag.file = file;
    new_tag.pos_start = pos_start;
    new_tag.pos_end = pos_end;
    new_tag.max_finds = max_finds;
    new_tag.min_finds = min_finds;
    new_tag.not_include_in_barcode = not_include_in_barcode;
    new_tag.trim = trim;
    new_tag.trim_offset = trim_offset;
    new_tag.has_after = false;
    new_tag.has_before = false;
    
    // Now deal with adding the actual sequence:
    if (nFiles <= 0 && new_tag.file == -1) { // Make sure we have nFiles set if tag can belong to any file
      std::cerr << "Error: nFiles must be set to a positive integer" << std::endl;
      return false;
    }
    // Only go through loop once; but if tag can belong to any file (new_tag.file == -1), iterate through nFiles (we'll make one new_tag per file)
    int start_file = new_tag.file == -1 ? 0 : new_tag.file;
    int end_file = new_tag.file == -1 ? nFiles : new_tag.file+1;
    for (file = start_file; file < end_file; file++,new_tag_index++) {
      new_tag.file = file;
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
            std::cerr << "Error: Sequence #" << n_tag_entries << ": \"" << name << "\" contains a non-ATCG character" << std::endl;
            return false;
          }
        }
        if (pos_end != 0 && pos_end - pos_start < seq.length()) {
          std::cerr << "Error: Sequence #" << n_tag_entries << ": \"" << name << "\" is too long to fit in the supplied location" << std::endl;
          return false;
        }
        
        std::unordered_map<std::string,int> mismatches;
        generate_indels_hamming_mismatches(seq, mismatch_dist, indel_dist, total_dist, mismatches);
        for (auto mm : mismatches) {
          std::string mismatch_seq = mm.first;
          int error = mm.second; // The number of substitutions, insertions, or deletions
          addToMap(mismatch_seq, new_tag_index, error);
          // DEBUG:
          // std::cout << seq << ": " << mismatch_seq << " " << error << " | " << total_dist << " " << mm.second << std::endl;
        }
        addToMap(seq, new_tag_index);
      }
      if (num_seqs == 0) {
        std::cerr << "Error: Sequence #" << n_tag_entries << ": \"" << name << "\" is empty" << std::endl;
        return false;
      }
      if (!after_str.empty()) {
        before_after_vec.push_back(std::make_pair(new_tag_index, std::make_pair(true, after_str)));
      }
      if (!before_str.empty()) {
        before_after_vec.push_back(std::make_pair(new_tag_index, std::make_pair(false, before_str)));
      }
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
  
  bool containsRegion(int16_t file_1, int32_t pos_start_1, int32_t pos_end_1, int16_t file_2, int32_t pos_start_2, int32_t pos_end_2, int l) {
    // Checks if 2 (the k-mer) is contained within 1 (the tag)
    if (file_1 != file_2 && !(file_1 == -1 || file_2 == -1)) {
      return false;
    }
    if (pos_start_1 < 0) { // If we have to calculate start position from the right-hand side of read
      pos_start_1 = l+pos_start_1;
      if (pos_start_1 < 0) {
        return false; // Length of read sequence too short to get start position from right-hand side of read
      }
    }
    if (pos_start_2 < pos_start_1) {
      return false;
    }
    if (pos_end_2 == 0 && pos_end_1 != 0) {
      return false;
    }
    if (pos_end_2 > pos_end_1 && pos_end_1 != 0) {
      return false;
    }
    return true;
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
    SeqString sstr(seq);
    if (tags.find(sstr) != tags.end()) {
      auto& v = tags[sstr];
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
      tags.insert({sstr,v});
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
      if (line[0] == '@') {
        std::stringstream ss(line);
        std::string field;
        std::string value;
        ss >> field >> value;
        if (value.empty()) {
          std::cerr << "Error: The file \"" << config_file << "\" contains an invalid line starting with @" << std::endl;
          return false;
        }
        if (field == "@trim-5") {
          if (!this->trim_5_str.empty()) {
            std::cerr << "Error: The file \"" << config_file << "\" specifies @trim-5 which was already previously set" << std::endl;
            return false;
          }
          this->trim_5_str = value;
        } else if (field == "@trim-3") {
          if (!this->trim_3_str.empty()) {
            std::cerr << "Error: The file \"" << config_file << "\" specifies @trim-3 which was already previously set" << std::endl;
            return false;
          }
          this->trim_3_str = value;
        } else if (field == "@prefix") {
          if (!this->barcode_prefix.empty()) {
            std::cerr << "Error: The file \"" << config_file << "\" specifies @prefix which was already previously set" << std::endl;
            return false;
          }
          this->barcode_prefix = value;
        } else if (field == "@extract") {
          if (!this->extract_str.empty()) {
            std::cerr << "Error: The file \"" << config_file << "\" specifies @extract which was already previously set" << std::endl;
            return false;
          }
          this->extract_str = value;
        }
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
      std::string group = "";
      std::string after_str = "";
      std::string before_str = "";
      int mismatch, indel, total_dist;
      parseDistance("", mismatch, indel, total_dist); // Set up default values
      int16_t file;
      int32_t pos_start;
      int32_t pos_end;
      parseLocation("", file, pos_start, pos_end); // Set up default values
      uint16_t max_finds = 0;
      uint16_t min_finds = 0;
      uint16_t max_finds_g = 0;
      uint16_t min_finds_g = 0;
      bool trim_left, trim_right;
      int trim_left_offset, trim_right_offset;
      parseTrimStr("", trim_left, trim_left_offset); // Set up default values
      parseTrimStr("", trim_right, trim_right_offset); // Set up default values
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
        } else if (h[i] == "GROUPS") {
          group = field;
        } else if (h[i] == "MINFINDS") {
          std::stringstream(field) >> min_finds;
        } else if (h[i] == "MAXFINDS") {
          std::stringstream(field) >> max_finds;
        } else if (h[i] == "MINFINDSG") {
          std::stringstream(field) >> min_finds_g;
        } else if (h[i] == "MAXFINDSG") {
          std::stringstream(field) >> max_finds_g;
        } else if (h[i] == "EXCLUDE") {
          std::stringstream(field) >> exclude;
        } else if (h[i] == "AFTER" || h[i] == "NEXT") {
          std::stringstream(field) >> after_str;
          ret = ret && validateBeforeAfterStr(after_str);
        } else if (h[i] == "BEFORE" || h[i] == "PREVIOUS") {
          std::stringstream(field) >> before_str;
          ret = ret && validateBeforeAfterStr(before_str);
        } else if (h[i] == "LEFT") {
          ret = ret && parseTrimStr(field, trim_left, trim_left_offset);
        } else if (h[i] == "RIGHT") {
          ret = ret && parseTrimStr(field, trim_right, trim_right_offset);
        } else {
          std::cerr << "Error: The file \"" << config_file << "\" contains the invalid column header: " << h[i] << std::endl;
          return false;
        }
      }
      if (trim_left && trim_right) {
        std::cerr << "Error: One of the barcodes has both left and right trimming specified" << std::endl;
        ret = false;
      }
      auto trim_dir = trim_left ? left : (trim_right ? right : nodir);
      auto trim_offset = trim_left ? trim_left_offset : (trim_right ? trim_right_offset : 0);
      if (!ret || !addTag(bc, name.empty() ? bc : name, group, mismatch, indel, total_dist, file, pos_start, pos_end, max_finds, min_finds, exclude, trim_dir, trim_offset, after_str, before_str)) {
        std::cerr << "Error: The file \"" << config_file << "\" contains an error" << std::endl;
        return false;
      }
      if (!group.empty() && (max_finds_g != 0 || min_finds_g != 0)) {
        if (!addGroupOptions(group, max_finds_g, min_finds_g)) {
          std::cerr << "Error: The file \"" << config_file << "\" contains an error" << std::endl;
          return false;
        }
      }
    }
    checkInit();
    return true;
  }
  
  bool getTag(std::string& seq, uint32_t& tag_id, int file, int pos, int& k, int l, bool look_for_initiator = false,
             bool search_tag_name_after = false, bool search_group_after = false, uint32_t search_id_after = -1,
             bool search_tag_before = false, uint32_t group_curr = -1, uint32_t name_id_curr = -1, int end_pos_curr = 0) {
    checkInit();
    int k_expanded = k;
    uint32_t updated_tag_id;
    uint32_t updated_name_id;
    int updated_k;
    int updated_error;
    bool found = false;
    while (k_expanded != -1) {
      bool found_curr = false;
      uint32_t tag_id_;
      int error_prev;
      uint32_t name_id_curr;
      uint32_t tag_id_curr;
      int curr_k = k_expanded;
      k_expanded = -1;
      const auto& it = tags.find(SeqString(seq.c_str()+pos, curr_k));
      if (it == tags.end()) {
        break;
      }
      for (auto &x : it->second) {
        if (x.second == -1) {
          k_expanded = x.first;
          continue;
        }
        tag_id_ = x.first;
        auto& tag = tags_vec[tag_id_];
        if (search_tag_name_after && tag.name_id != search_id_after) {
          continue;
        } else if (search_group_after && tag.group != search_id_after) {
          continue;
        }
        if (tag.has_before || tag.has_before_group) {
          if (!search_tag_before) {
            continue;
          }
          if (tag.has_before && tag.id_before != name_id_curr) {
            continue;
          } else if (tag.has_before_group && tag.id_before != group_curr) {
            continue;
          } else {
            if (pos-end_pos_curr < tag.extra_before) {
              continue;
            }
            if (tag.extra_before2 != 0 && pos-end_pos_curr >= tag.extra_before2) {
              continue;
            }
          }
        }
        if (containsRegion(tag.file, tag.pos_start, tag.pos_end, file, pos, pos+curr_k, l)) {
          if (!look_for_initiator || (look_for_initiator && tags_vec[tag_id_].initiator)) {
            if (found_curr && tag.name_id != name_id_curr) {
              found_curr = false; // seq of length curr_k maps to multiple tags of different names
              break;
            }
            if (!found_curr || (found_curr && error_prev > x.second)) {
              error_prev = x.second;
              tag_id_curr = tag_id_; // if tags have same name but different mismatch errors: choose the tag w/ smallest error
            }
            name_id_curr = tag.name_id;
            found_curr = true;
          }
        }
      }
      // Algorithm works as follows:
      // // for a given k, remove that k from consideration if there are multiple tag.name_id's for that k
      // // however, if there are multiple tags of the same name_id for that k, pick the tag with the smallest error
      // // afterwards, compare across all k's being considered: if multiple tag.name_id's across different k's, return false (-1), otherwise pick the tag associated with the largest k
      if (found_curr) {
        if (!found) { // First time identifying a tag
          found = true;
          updated_tag_id = tag_id_curr;
          updated_k = curr_k;
          updated_error = error_prev;
          updated_name_id = name_id_curr;
        } else { // Already previously identified a tag when looking at a smaller k
          if (updated_name_id != name_id_curr) {
            return false; // multiple tags of different names
          }
          if (updated_error >= error_prev) { // Choose smallest error first when deciding if to update to larger k
            updated_tag_id = tag_id_curr;
            updated_k = curr_k; // Update to larger k
            updated_error = error_prev;
            updated_name_id = name_id_curr;
          }
        }
      }
    }
    if (found) {
      tag_id = updated_tag_id;
      k = updated_k;
      return true;
    }
    return false;
  }
  
  int getNumTags() {
    return tags_vec.size();
  }
  
  int getNumTagsOriginallyAdded() {
    return n_tag_entries;
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
  
  static bool parseTrimStr(const std::string& s_trim, bool& trim, int& offset) {
    trim = false;
    offset = 0;
    char delimeter = ':';
    if (s_trim.find(',') < s_trim.length()) {
      delimeter = ','; // If string contains commas, use commas as delimeter
    }
    std::stringstream ss_trim(s_trim);
    std::string trim_attribute;
    int i = 0;
    try {
      while (std::getline(ss_trim, trim_attribute, delimeter)) {
        if (!trim_attribute.empty()) {
          if (i == 0) {
            trim = std::stoi(trim_attribute);
          } else if (i == 1) {
            offset = std::stoi(trim_attribute);
          }
        }
        i++;
      }
      if (i > 2 || (!trim && offset != 0)) {
        std::cerr << "Error: Trim string is malformed; unable to parse \"" << s_trim << "\"" << std::endl;
        return false;
      }
    } catch (std::exception &e) {
      std::cerr << "Error: Could not convert \"" << trim_attribute << "\" to int in trim string" << std::endl;
      return false;
    }
    return true;
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
      if (i > 3 || file < -1 || (file >= nFiles && nFiles != -1) || pos_end < 0 || (pos_end <= pos_start && pos_end != 0)) {
        std::cerr << "Error: Location string is malformed; unable to parse \"" << location << "\"" << std::endl;
        return false;
      }
    } catch (std::exception &e) {
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
    } catch (std::exception &e) {
      std::cerr << "Error: Could not convert \"" << dist_attribute << "\" to int in distance string" << std::endl;
      return false;
    }
    return true;
  }
  
  static bool validateBeforeAfterStr(const std::string& s) {
    bool ret = true;
    if (s.empty()) {
      return true;
    }
    if (s.length() < 3) {
      ret = false;
    } else if (s[0] != '{') {
      ret = false;
    } else if (s[1] == '{') {
      if (std::count(s.begin(), s.end(), '}') != 2 || s.find("}}") == std::string::npos) {
        ret = false;
      }
    } else if (std::count(s.begin(), s.end(), '}') != 1 || s.find("}") == std::string::npos) {
      ret = false;
    }
    if (ret && s[s.length()-1] != '}') {
      try {
        std::string ssub = s.substr(s.find_last_of('}')+1);
        std::stringstream ss(ssub);
        std::string s2;
        while (std::getline(ss, s2, '-')) {
          std::stoi(s2);
        }
      } catch (std::exception &e) {
        ret = false;
      }
    }
    if (!ret) {
      std::cerr << "Error: The following string is invalid: \"" << s << "\" for the next/previous fields" << std::endl;
    }
    return ret;
  }
  
  bool addExistingMapping(std::string mapping_file) {
    struct stat stFileInfo;
    auto intstat = stat(mapping_file.c_str(), &stFileInfo);
    if (intstat != 0) {
      std::cerr << "Error: file not found " << mapping_file << std::endl;
      return false;
    }
    std::ifstream mfile(mapping_file);
    std::string line;
    while (std::getline(mfile,line)) {
      if (line.size() == 0) {
        continue;
      }
      std::stringstream ss(line);
      std::string barcode;
      std::string names_list;
      uint32_t count;
      ss >> barcode >> names_list >> count;
      if (!barcode_prefix.empty()) {
        std::string p = barcode.substr(0,barcode_prefix.length());
        if (p != barcode_prefix) {
          std::cerr << "Error: File " << mapping_file << " contains barcode: " << barcode << ", which does not have prefix: " << barcode_prefix << std::endl;
          return false;
        }
      }
      std::string barcode_no_prefix = barcode.substr(barcode_prefix.length());
      if (barcode_no_prefix.length() != FAKE_BARCODE_LEN) {
        std::cerr << "Error: File " << mapping_file << " contains a barcode of invalid length: " << barcode << std::endl;
        return false;
      }
      barcode = barcode_no_prefix;
      if (hashKmer(barcode) != idmapinv.size()) { // Barcodes need to be ordered 0,1,2,3,... in terms of their binary values
        std::cerr << "Error: File " << mapping_file << " contains an invalid list of sequences in the first column" << std::endl;
        return false;
      }
      std::stringstream ss1(names_list);
      std::string name;
      char delimeter = ',';
      std::vector<uint32_t> u;
      while (std::getline(ss1, name, delimeter)) {
        if (name.size() == 0) {
          continue;
        }
        const auto& itnames = std::find(names.begin(), names.end(), name);
        if (itnames == names.end()) {
          std::cerr << "Error: File " << mapping_file << " contains the name \"" << name << "\" which does not exist" << std::endl;
          return false;
        }
        u.push_back(itnames - names.begin());
      }
      auto it = idmapinv.find(u);
      if (it != idmapinv.end()) {
        std::cerr << "Error: In file " << mapping_file << ", the following is duplicated: " << names_list << std::endl;
        return false;
      }
      idmapinv.insert({u,idmapinv.size()});
      idmap.push_back(std::move(u));
      idcount.push_back(0);
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
      std::string ofile = "";
      if (line.size() == 0) {
        continue;
      }
      std::stringstream s(line);
      if (!discard) {
        std::string a,b;
        s >> a >> b;
        line = a;
        if (!b.empty()) {
          ofile = b;
        }
      }
      std::stringstream ss(line);
      std::string name;
      char delimeter = ',';
      std::vector<uint32_t> u;
      while (std::getline(ss, name, delimeter)) {
        if (name.size() == 0) {
          continue;
        }
        const auto& itnames = std::find(names.begin(), names.end(), name);
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
        idmapinv_keep.insert({u,ofile});
        keep_check = true;
      }
    }
    return true;
  }
  
  bool addFilterListGroup(std::string keep_file, bool discard=false) {
    struct stat stFileInfo;
    auto intstat = stat(keep_file.c_str(), &stFileInfo);
    if (intstat != 0) {
      std::cerr << "Error: file not found " << keep_file << std::endl;
      return false;
    }
    std::ifstream kfile(keep_file);
    std::string line;
    while (std::getline(kfile,line)) {
      std::string ofile = "";
      if (line.size() == 0) {
        continue;
      }
      std::stringstream s(line);
      if (!discard) {
        std::string a,b;
        s >> a >> b;
        line = a;
        if (!b.empty()) {
          ofile = b;
        }
      }
      std::stringstream ss(line);
      std::string name;
      char delimeter = ',';
      std::vector<uint32_t> u;
      while (std::getline(ss, name, delimeter)) {
        if (name.size() == 0) {
          continue;
        }
        const auto& itnames = std::find(group_names.begin(), group_names.end(), name);
        if (itnames == group_names.end()) {
          std::cerr << "Error: File " << keep_file << " contains the group name \"" << name << "\" which does not exist" << std::endl;
          return false;
        }
        u.push_back(itnames - group_names.begin());
      }
      auto it1 = groupmapinv_keep.find(u);
      auto it2 = groupmapinv_discard.find(u);
      if (it1 != groupmapinv_keep.end() || it2 != groupmapinv_discard.end()) {
        std::cerr << "Error: In file " << keep_file << ", the following line is duplicated: " << line << std::endl;
        return false;
      }
      if (discard) {
        groupmapinv_discard.insert({u,0});
        discard_check_group = true;
      } else {
        groupmapinv_keep.insert({u,ofile});
        keep_check_group = true;
      }
    }
    return true;
  }
  
  bool addGroupOptions(std::string group_name, uint16_t max_finds, uint16_t min_finds) {
    const auto& itnames = std::find(group_names.begin(), group_names.end(), group_name);
    if (itnames == group_names.end()) {
      std::cerr << "Error: Group name \"" << group_name << "\" does not exist" << std::endl;
      return false;
    }
    auto i = itnames - group_names.begin();
    if (max_finds > 0) {
      if (max_finds_group_map.find(i) != max_finds_group_map.end() && max_finds_group_map[i] != max_finds) {
        std::cerr << "Error: Group name \"" << group_name << "\" had max finds specified multiple times" << std::endl;
        return false;
      }
      max_finds_group_map[i] = max_finds;
    }
    if (min_finds > 0) {
      if (min_finds_group_map.find(i) != min_finds_group_map.end() && min_finds_group_map[i] != min_finds) {
        std::cerr << "Error: Group name \"" << group_name << "\" had min finds specified multiple times" << std::endl;
        return false;
      }
      min_finds_group_map[i] = min_finds;
    }
    return true;
  }
  
  bool parseExtractStr(std::string extract_str) {
    std::string &s = extract_str;
    s.erase(remove(s.begin(),s.end(),' '),s.end()); // remove spaces from string
    if (s.find(',') != std::string::npos) { // If we are supplied multiple extraction patterns (comma-separated)
      bool ret = false;
      std::string s_ = "";
      for (int i = 0; i < s.length(); i++) {
        if (s[i] == ',') {
          if (!parseExtractStr(s_)) {
            return false;
          }
          s_ = "";
        }
        else {
          s_.push_back(s[i]);
        }
      }
      return parseExtractStr(s_);
    }
    UMI umi;
    umi.id = curr_umi_id_i++;
    auto& umi_name = umi.name;
    auto& group1 = umi.group1;
    auto& group2 = umi.group2;
    auto& name1_present = umi.id1_present;
    auto& name2_present = umi.id2_present;
    auto& length_range_start = umi.length_range_start;
    auto& length_range_end = umi.length_range_end;
    auto& padding_left = umi.padding_left;
    auto& padding_right = umi.padding_right;
    length_range_start = 0;
    length_range_end = 0;
    padding_left = 0;
    padding_right = 0;
    int file1 = -1, file2 = -1, pos1 = -1, pos2 = -1;
    try {
      // Find the UMI name and position:
      auto umi_open = s.find_first_of('<');
      auto umi_close = s.find_first_of('>');
      if (umi_open == std::string::npos || umi_close == std::string::npos || umi_open > umi_close ||
          s.find_last_of('<') != umi_open || s.find_last_of('>') != umi_close) {
        return false; // malformed
      }
      umi_name = s.substr(umi_open+1,umi_close-umi_open-1);
      // Find the UMI length range (e.g. <umi[number]> or <umi[start-end]>):
      auto length_range_open = umi_name.find_first_of('[');
      auto length_range_close = umi_name.find_first_of(']');
      if ((length_range_open == std::string::npos || length_range_close == std::string::npos || length_range_open > length_range_close ||
          umi_name.find_last_of('[') != length_range_open || umi_name.find_last_of(']') != length_range_close || length_range_close != umi_name.length()-1) && 
          !(length_range_open == std::string::npos && length_range_close == std::string::npos)) {
          return false; // malformed
      }
      if (length_range_open != std::string::npos) {
        std::string length_range_str = umi_name.substr(length_range_open+1,length_range_close-length_range_open-1); // extract length range
        umi_name = umi_name.substr(0, length_range_open); // remove length range from umi_name
        auto dash_pos = length_range_str.find_first_of('-');
        if (dash_pos == std::string::npos) { // no range supplied, just a number
          auto n = std::stoi(length_range_str);
          if (n <= 0) {
            return false;
          }
          length_range_start = n;
          length_range_end = length_range_start;
        } else {
          auto n1 = std::stoi(length_range_str.substr(0, dash_pos));
          auto n2 = std::stoi(length_range_str.substr(dash_pos+1, length_range_str.length()-dash_pos));
          if (n1 <= 0 || n2 <= 0 || n1 > n2) {
            return false;
          }
          length_range_start = n1;
          length_range_end = n2;
        }
      }
      if (std::count_if(umi_name.begin(),umi_name.end(),[](char c) { return !(std::isalnum(c) || c == '_'); }) > 1) {
        return false; // malformed; non-alphanumeric (and non-underscore) character found in umi name
      }
      // Function to parse location strings: "file:pos"
      auto parse_location = [](std::string str) { 
        std::stringstream ss(str);
        std::string location_attribute;
        int file = -1;
        int pos = -1;
        int i = 0;
        while (std::getline(ss, location_attribute, ':')) {
          if (!location_attribute.empty()) {
            if (i == 0) {
              file = std::stoi(location_attribute);
            } else if (i == 1) {
              pos = std::stoi(location_attribute);
            }
          }
          i++;
        }
        if (i > 2 || file < 0 || pos < -1) {
          return std::pair<int,int>(-1,-1);
        }
        return std::pair<int,int>(file,pos);
      };
      // Find the barcode or location string before the UMI:
      std::string name1 = "";
      std::string s1 = s.substr(0,umi_open);
      group1 = false; // Whether it's a barcode name {bc} or a barcode group {{group}}
      auto pos1_open = s1.find_first_of('{');
      auto pos1_close = s1.find_first_of('}');
      name1_present = !(pos1_open == std::string::npos && pos1_close == std::string::npos);
      if (!name1_present) {
        if (!s1.empty()) {
          // The string is formatted as [location]<umi> or [location]<umi>{bc} or [location]<umi>[padding]{bc}
          auto loc = parse_location(s1);
          file1 = loc.first;
          pos1 = loc.second;
          if (pos1 < 0 || file1 < 0 || file1 >= nFiles) {
            return false;
          }
        }
      } else if (pos1_open == std::string::npos || pos1_close == std::string::npos || pos1_open > pos1_close || pos1_open != 0) {
        return false; // malformed
      } else {
        // Extract barcode or group name:
        if (s1[pos1_open+1] == '{') {
          if (s1[pos1_close+1] != '}') {
            return false; // malformed
          }
          group1 = true;
          name1 = s1.substr(pos1_open+2,pos1_close-pos1_open-2);
          pos1_close++;
        } else {
          name1 = s1.substr(pos1_open+1,pos1_close-pos1_open-1);
        }
        // Extract padding: {bc}[padding]<umi>...
        auto ss = s1.substr(pos1_close+1);
        if (!ss.empty()) {
          padding_left = std::stoi(ss);
        }
      }
      // Find the barcode or location string after the UMI:
      std::string s2 = s.substr(umi_close+1);
      std::string name2 = "";
      group2 = false; // Whether it's a barcode name {bc} or a barcode group {{group}}
      auto pos2_open = s2.find_first_of('{');
      auto pos2_close = s2.find_first_of('}');
      name2_present = !(pos2_open == std::string::npos && pos2_close == std::string::npos);
      if (!name2_present) {
        if (!s2.empty()) {
          // The string is formatted as {bc}<umi>[location] or {bc}[padding]<umi>[location] or <umi>[location] or [location]<umi>[location]
          auto loc = parse_location(s2);
          file2 = loc.first;
          pos2 = loc.second;
          if (pos2 < -1 || file2 < 0 || file2 >= nFiles) { // pos=-1 means we anchor at the end of the read
            return false;
          }
        }
      } else if (pos2_close == std::string::npos || pos2_close == std::string::npos || pos2_open > pos2_close || (pos2_close != s2.length()-1 && !(pos2_close+1 == s2.length()-1 && s2[pos2_close+1] == '}'))) {
        return false; // malformed
      } else {
        // Extract barcode or group name:
        if (s2[pos2_open+1] == '{') {
          if (s2[pos2_close+1] != '}') {
            return false; // malformed
          }
          group2 = true;
          name2 = s2.substr(pos2_open+2,pos2_close-pos2_open-2);
          pos2_close++;
        } else {
          name2 = s2.substr(pos2_open+1,pos2_close-pos2_open-1);
        }
        // Extract padding: ...<umi>[padding]{bc}
        auto ss = s2.substr(0, pos2_open);
        if (!ss.empty()) {
          padding_right = std::stoi(ss);
        }
      }
      // Process the barcodes, locations, paddings, and UMIs:
      auto process_umi = [this, &umi](std::string name, bool first) {
        auto name_present = first ? umi.id1_present : umi.id2_present;
        auto group = first ? umi.group1 : umi.group2;
        auto file = first ? umi.location1.first : umi.location2.first;
        auto pos = first ? umi.location1.second : umi.location2.second;
        auto& id = first ? umi.id1 : umi.id2;
        if (name_present) {
          if (group) {
            const auto& itnames = std::find(this->group_names.begin(), this->group_names.end(), name);
            if (itnames == this->group_names.end()) {
              return false;
            }
            id = itnames - this->group_names.begin();
            this->umi_group_map[id].push_back(umi);
          } else {
            const auto& itnames = std::find(this->names.begin(), this->names.end(), name);
            if (itnames == this->names.end()) {
              return false;
            }
            id = itnames - this->names.begin();
            this->umi_name_map[id].push_back(umi);
          }
        } else {
          this->umi_loc_map[std::make_pair(file,pos)].push_back(umi);
        }
        return true;
      };
      bool length_range_supplied = !(length_range_start == 0 && length_range_end == 0); // If a length range was supplied (e.g. <umi[number]> or <umi[start-end]>)
      if ((name1_present || file1 >= 0) && (name2_present || file2 >= 0)) { // <umi> is sandwiched between {bc}'s and/or locations
        // Length range is optional (if length range is provided and extracted UMI doesn't fit in length range, it's not identified)
        // i.e. <umi[number]> or <umi[start-end]> or <umi> are all acceptable
        if (!name1_present && !name2_present && file1 != file2) {
          return false; // malformed; different file numbers supplied
        } else if (!name1_present && !name2_present && pos1 >= pos2 && pos2 != -1) {
          return false; // malformed; <umi> is sandwiched between two locations but the locations are invalid
        }
      } else if ((name1_present || file1 >= 0) || (name2_present || file2 >= 0)) { // <umi> has {bc} or location on either left side or right side of it
        // Length range is necessary but only whole numbers, not ranges, allowed (number [length_range_start] dictates how many characters to extract on the right or left of the UMI)
        // i.e. <umi[number]> is necessary; <umi[start-end]> or <umi> are not acceptable
        if (!length_range_supplied || length_range_start == 0 || length_range_end != length_range_start) {
          return false;
        }
      } else {
        return false; // malformed; <umi> must have at least one {bc} or location next to it
      }
      umi.location1 = std::make_pair(file1, pos1);
      umi.location2 = std::make_pair(file2, pos2);
      bool ret = process_umi(name1, true) && process_umi(name2, false);
      if (!ret) {
        return false;
      }
      // DEBUG:
      /*std::cout << "UMI " << umi.id << ": " << umi.name 
                << "; length: (" << length_range_start << "," << length_range_end 
                << "); padding: (" << padding_left << "," << padding_right << ")";
      if (umi.id1_present) {
        std::cout << " --- " << (umi.group1 ? "group1: " : "barcode1: ") << name1 << " " << umi.id1;
      } else {
        std::cout << " --- " << "location1: (" << umi.location1.first << "," << umi.location1.second << ")";
      }
      if (umi.id2_present) {
        std::cout << ", " << (umi.group2 ? "group2: " : "barcode2: ") << name2 << " " << umi.id2;
      } else {
        std::cout << ", " << "location2: (" << umi.location2.first << "," << umi.location2.second << ")";
      }
      std::cout << std::endl;*/
      // Add UMI name to vector:
      if (std::find(umi_names.begin(), umi_names.end(), umi.name) == umi_names.end()) {
        umi_names.push_back(umi.name);
      }
      return true;
    } catch (std::exception &e) {
      return false;
    }
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
    // Note: s and l may end up being trimmed/modified (even if the read ends up becoming unassigned)
    results.id = -1;
    results.discard = false;
    auto min_finds = min_finds_map; // copy
    auto max_finds = max_finds_map; // copy
    auto min_finds_group = min_finds_group_map; // copy
    auto max_finds_group = max_finds_group_map; // copy
    bool check_group = keep_check_group || discard_check_group;
    std::vector<uint32_t> group_v(0);
    if (check_group) {
      group_v.reserve(16);
    }
    int n = std::min(jmax, (int)kmer_size_locations.size());
    for (int j = 0; j < jmax; j++) {
      int file = j;
      int readLength = l[file];
      // First, let's do end-trimming
      int trim_5 = std::min(trim_5_3_vec[file].first, l[file]);
      s[file] += trim_5;
      l[file] -= trim_5;
      int trim_3 = std::min(trim_5_3_vec[file].second, l[file]);
      l[file] = l[file] - trim_3 <= 0 ? 0 : l[file] - trim_3;
      if (l[file] != readLength) {
        results.modtrim.push_back(std::make_pair(file, std::make_pair(trim_5, l[file])));
      }
      if (j >= n) {
        continue;
      }
      readLength = l[file];
      bool look_for_initiator = initiator_files[file];
      std::string seq(s[file], readLength);
      bool found_weird_base = false; // non-ATCG bases
      uint32_t rando;
      for (auto& c: seq) {
        c &= 0xDF; // Convert a/t/c/g to upper case
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
          if (random_replacement) {
            if (!found_weird_base) {
              rando = hashSequence(seq);
              found_weird_base = true;
            }
            c = "ATCG"[rando%4]; // substitute non-ATCG base for pseudo-random base
            rando = ((rando ^ (rando >> 3)) ^ (rando << 20)) ^ (rando >> 9);
          } else {
            c = 'N';
          }
        }
      }
      int left_trim = 0;
      int right_trim = 0;
      bool right_trim_found = false;
      auto& kmers = kmer_size_locations[file];
      bool search_tag_before = false;
      uint32_t group_curr;
      uint32_t name_id_curr;
      bool search_tag_name_after = false;
      bool search_group_after = false;
      uint32_t search_id_after;
      uint16_t search_extra_after;
      uint16_t search_extra_after2;
      int search_after_start;
      for (Locations locations(kmers, readLength); locations.good(); ++locations) {
        auto loc = locations.get();
        auto k = loc.first;
        auto pos = loc.second;
        // DEBUG:
        // std::cout << "file=" << file << " k=" << k << " pos=" << pos << std::endl;
        if (search_tag_name_after || search_group_after) {
          if (pos-search_after_start < search_extra_after) {
            continue;
          }
          if (search_extra_after2 != 0 && pos-search_after_start >= search_extra_after2) {
            continue; // Could break here?
          }
        }
        uint32_t tag_id;
        if (getTag(seq, tag_id, file, pos, k, readLength, look_for_initiator, 
                   search_tag_name_after, search_group_after, search_id_after,
                   search_tag_before, group_curr, name_id_curr, search_after_start)) {
          look_for_initiator = false;
          auto& tag = tags_vec[tag_id];
          if (tag.min_finds > 0) {
            min_finds[tag_id]--;
          }
          if (min_finds_group.find(tag.group) != min_finds_group.end()) {
            min_finds_group[tag.group]--;
          }
          if (tag.max_finds > 0) {
            if (max_finds[tag_id]-- <= 0) {
              continue; // maxFinds exceeded; just continue
            }
          }
          if (max_finds_group.find(tag.group) != max_finds_group.end()) {
            if (max_finds_group[tag.group]-- <= 0) {
              continue; // maxFindsG exceeded; just continue
            }
          }
          search_group_after = false; // reset
          search_tag_name_after = false; // reset
          search_tag_before = true;
          name_id_curr = tag.name_id;
          group_curr = tag.group;
          search_after_start = pos+k; // aka end_pos_curr
          if (tag.has_after) {
            if (tag.has_after_group) {
              search_group_after = true;
            } else {
              search_tag_name_after = true;
            }
            search_id_after = tag.id_after;
            search_extra_after = tag.extra_after;
            search_extra_after2 = tag.extra_after2;
          }
          if (!tag.not_include_in_barcode) {
            results.name_ids.push_back(tag.name_id);
            if (check_group) {
              group_v.push_back(tag.group);
            }
          }
          if (tag.trim == left) {
            left_trim = pos+k+tag.trim_offset;
          } else if (tag.trim == right && !right_trim_found) {
            right_trim = (readLength-pos)+tag.trim_offset;
            right_trim_found = true;
          }
          
          if (tag.terminator) {
            break; // End the search for the current (j'th) read file's sequence
          }
          locations.setJump(k);
        }
      }
      // Modify (trim) the reads
      if (left_trim != 0) {
        left_trim = std::min(left_trim, readLength);
        s[file] += left_trim;
        l[file] -= left_trim;
      }
      if (right_trim != 0) {
        right_trim = std::min(right_trim, readLength);
        if (l[file] - right_trim <= 0) {
          l[file] = 0;
        } else {
          l[file] -= right_trim;
        }
      }
      if (l[file] != readLength) {
        results.modtrim.push_back(std::make_pair(file, std::make_pair(left_trim, l[file])));
      }
    }
    for (auto& it : min_finds) {
      if (it.second > 0) {
        results.name_ids.clear(); // minFinds not met
        break;
      }
    }
    for (auto& it : min_finds_group) {
      if (it.second > 0) {
        results.name_ids.clear(); // minFindsG not met
        break;
      }
    }
    auto &u = results.name_ids;
    if (u.empty()) {
      return;
    }
    if (keep_check) {
      auto it = idmapinv_keep.find(u);
      if (it == idmapinv_keep.end()) {
        results.discard = true;
        return;
      }
      results.ofile = it->second;
    }
    if (discard_check && idmapinv_discard.find(u) != idmapinv_discard.end()) {
      results.discard = true;
      return;
    }
    if (keep_check_group) {
      auto it = groupmapinv_keep.find(group_v);
      if (it == groupmapinv_keep.end()) {
        results.discard = true;
        return;
      }
      if (results.ofile.empty()) { // We prioritize name IDs over groups for user-specified output files
        results.ofile = it->second;
      }
    }
    if (discard_check_group && groupmapinv_discard.find(group_v) != groupmapinv_discard.end()) {
      results.discard = true;
      return;
    }
  }
  
  static void modifyRead(std::vector<std::pair<const char*, int>>& seqs, std::vector<std::pair<const char*, int>>& quals, int i, const Results& results) {
    // Modify (trim) the reads in the actual read buffer itself
    size_t q_size = quals.size();
    for (auto& mt : results.modtrim) {
      int j = mt.first;
      int index = i+j;
      int leftOffset = mt.second.first;
      int readLength = mt.second.second;
      seqs[index].first += leftOffset;
      seqs[index].second = readLength;
      if (q_size > index) { // Just in case we decided not to store quality scores
        quals[index].first += leftOffset;
        quals[index].second = readLength;
      }
    }
  }
  
  void update(std::vector<Results>& rv) {
    // Should only be called under a lock (can't have multiple threads accessing a common container)
    for (auto& r : rv) {
      auto& u = r.name_ids;
      if (u.empty() || !isAssigned(r)) {
        continue;
      }
      auto it = idmapinv.find(u);
      int id;
      if (it != idmapinv.end()) {
        id = it->second;
        idcount[id]++;
      } else {
        id = idmapinv.size();
        // DEBUG hash function:
        //VectorHasher h;
        //std::cout << h(u) << std::endl;
        idmapinv.insert({u,id});
        idmap.push_back(u);
        idcount.push_back(1);
      }
      r.id = id;
    }
  }
  
  bool isAssigned(Results& r) {
    return (!r.discard && !r.name_ids.empty()) || always_assign;
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
    std::string o;
    while ((o = fetchNextBarcodeMapping()) != "") {
      of << o;
    }
    of.close();
  }
  
  std::string fetchNextBarcodeMapping() {
    int i = curr_barcode_mapping_i;
    if (i >= idmap.size()) {
      curr_barcode_mapping_i = 0;
      return "";
    }
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
    std::string o = binaryToString(getID(id), getBarcodeLength()) + "\t" + barcode_str + "\t" + std::to_string(n) + "\n";
    ++curr_barcode_mapping_i;
    return o;
  }
  
  uint64_t getID(uint64_t id) { // Get the "real ID" (aka results ID merged with prefix)
    if (barcode_prefix.empty()) {
      return id;
    }
    return ((hashKmer(barcode_prefix) << (2*FAKE_BARCODE_LEN)) | id);
  }
  
  int getBarcodeLength() {
    return FAKE_BARCODE_LEN+barcode_prefix.length();
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
  
  static uint64_t hashKmer(const char* s, size_t k) {
    uint64_t r = 0;
    if (k > MAX_K) {
      k = MAX_K;
    }
    size_t k1 = k/4;
    size_t k2 = k%4;
    for (size_t i = 0; i < k1; s += 4, ++i) {
      uint32_t x = (*(const uint32_t*)(s));
      x = (x ^ (x >> 1)) & 101058054;
      r = r << 8;
      r |= (((x<<5) | (x>>5) | (x>>15) | (x>>25)) & 255);
    }
    for (size_t i = 0; i < k2; ++i, ++s) {
      r = r << 2;
      r |= ((((*s) ^ ((*s) >> 1))) >> 1) & 3;
    }
    return r;
  }
  
  static uint64_t hashKmer(const std::string& key) {
    const char* s = key.c_str();
    size_t k = key.length();
    return hashKmer(s, k);
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
  
  class SeqStringHasher {
  public:
    size_t operator()(const SeqString& key) const {
      return key.p_ ? hashKmer(key.p_, key.l_) : hashKmer(key.s_.c_str(), key.l_);
    }
  };
  
  std::vector<SplitCodeTag> tags_vec;
  robin_hood::unordered_flat_map<SeqString, std::vector<tval>, SeqStringHasher> tags;
  std::vector<std::string> names;
  std::vector<std::string> group_names;
  
  std::vector<std::pair<uint32_t,std::pair<bool,std::string>>> before_after_vec;
  
  std::vector<std::vector<uint32_t>> idmap;
  robin_hood::unordered_flat_map<std::vector<uint32_t>, int, VectorHasher> idmapinv;
  std::vector<int> idcount;
  std::unordered_map<std::vector<uint32_t>, std::string, VectorHasher> idmapinv_keep;
  std::unordered_map<std::vector<uint32_t>, int, VectorHasher> idmapinv_discard;
  std::unordered_map<std::vector<uint32_t>, std::string, VectorHasher> groupmapinv_keep;
  std::unordered_map<std::vector<uint32_t>, int, VectorHasher> groupmapinv_discard;
  
  std::unordered_map<uint32_t,int> min_finds_map;
  std::unordered_map<uint32_t,int> max_finds_map;
  std::unordered_map<uint32_t,int> min_finds_group_map;
  std::unordered_map<uint32_t,int> max_finds_group_map;
  std::vector<bool> initiator_files;
  
  std::unordered_map<uint32_t,std::vector<UMI>> umi_name_map;
  std::unordered_map<uint32_t,std::vector<UMI>> umi_group_map;
  std::map<std::pair<int,int>,std::vector<UMI>> umi_loc_map; // hash for pair not defined for unordered_map
  std::vector<std::string> umi_names;
  
  std::vector<std::vector<std::pair<int,int>>> kmer_size_locations;
  
  std::string barcode_prefix;
  std::string trim_5_str, trim_3_str;
  std::string extract_str;
  std::vector<std::pair<int,int>> trim_5_3_vec;
  
  bool init;
  bool discard_check;
  bool keep_check;
  bool discard_check_group;
  bool keep_check_group;
  bool always_assign;
  bool random_replacement;
  int nFiles;
  int n_tag_entries;
  int curr_barcode_mapping_i;
  int curr_umi_id_i;
  static const int MAX_K = 32;
  static const size_t FAKE_BARCODE_LEN = 16;
  static const char QUAL = 'K';
};


#endif // SPLITCODE_H

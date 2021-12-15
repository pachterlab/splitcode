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
  
  std::vector<std::vector<int>> idmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> idmapinv;
  std::vector<int> idcount;
};


#endif // SPLITCODE_H

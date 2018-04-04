#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <utility>
#include "zz_util_const.h"
#include "phg.h"

#pragma once

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

class RS_set{
public:
	RS_set() : map() { _size = 0; } 
	void insert(std::pair<ZOLTAN_GNO_TYPE, size_t> a) { 
		map[a.first] = a.second; 
		_size += a.second;
	}
	void emplace(ZOLTAN_GNO_TYPE a, size_t s){ 
		map[a] = s; 
		_size += s;
	}
	void print_all(){
		for(auto &el : map){
			std::cout << "(" << el.first << ":" << el.second << ") ";
		}
		std::cout << std::endl;
	}
	void erase(ZOLTAN_GNO_TYPE a){
		_size -= map[a]; 
		map.erase(a); 
	}
	void erase(ZOLTAN_GNO_TYPE a, size_t s){ 
		map.erase(a); 
		_size -= s;
	}
	void erase(std::pair<ZOLTAN_GNO_TYPE, size_t> a){ 
		map.erase(a.first); 
		_size -= a.second;
	}
	size_t size() { return _size; }
	std::unordered_map<ZOLTAN_GNO_TYPE, size_t>::iterator begin(){ return map.begin(); }
	std::unordered_map<ZOLTAN_GNO_TYPE, size_t>::iterator end(){ return map.end(); }

private:
	size_t _size;
	std::unordered_map<ZOLTAN_GNO_TYPE, size_t> map;
};


std::unordered_map<ZOLTAN_GNO_TYPE, RS_set> RS_compute_stable_matching(HGraph* hg, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_c, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_f, const std::unordered_map<ZOLTAN_GNO_TYPE, unsigned int> &capacities_c);

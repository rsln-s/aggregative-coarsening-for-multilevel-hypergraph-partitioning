#include <iostream>
#include "phg.h"
#include "zz_util_const.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#pragma once

int RS_get_vertex_degree(HGraph* hg, int vtx);
int RS_amg_matching(ZZ *zz, HGraph *hg, PHGPartParams *hgp, ZOLTAN_GNO_TYPE *match, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R);

/* From https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes */

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

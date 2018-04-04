#include "rs_amg_match.hpp"
#include "rs_stable_assignment.hpp"
#include <map>
#include <numeric>

float rs_total_vwgt = -1.0;

//#define RS_INTEGER_FUTURE_VOLUMES
//#define RS_STABLE_MATCHING

int RS_get_vertex_degree(HGraph* hg, int vtx){ return hg->vindex[vtx+1] - hg->vindex[vtx]; }

void RS_print_neighbors(HGraph* hg, int vtx, const std::unordered_set<ZOLTAN_GNO_TYPE>& coarse, const std::unordered_set<ZOLTAN_GNO_TYPE>& fine){
	std::unordered_set<ZOLTAN_GNO_TYPE> neigh_c, neigh_f;
	std::cout << "Printing neighbors for vertex " << vtx << std::endl;
	 /* for every hyperedge containing the vertex vtx */
	for (ZOLTAN_GNO_TYPE i = hg->vindex[vtx]; i < hg->vindex[vtx+1]; i++) {
		ZOLTAN_GNO_TYPE e = hg->vedge[i];
		/* for every other vertex in the hyperedge */
		for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
			ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];				
			if(vtx == v2) { continue; } 
			if(coarse.find(v2) != coarse.end()){
				neigh_c.insert(v2);
				continue;
			}
			if(fine.find(v2) != fine.end()){
				neigh_f.insert(v2);
				continue;
			}
			std::cout << "RS_print_neighbors :: OOOOOOPS" << std::endl;
		}	
	}	
	for (auto &v: neigh_c){
		std::cout << v << ":c "; 
	}
	for (auto &v: neigh_f){
		std::cout << v << ":f "; 
	}
	
	std::cout << std::endl << "Total neighbors for "<< vtx << " : " << neigh_c.size() + neigh_f.size() << std::endl;
}

int RS_get_conditional_degree(HGraph* hg, int vtx, const std::unordered_set<ZOLTAN_GNO_TYPE>& coarse, const std::unordered_set<ZOLTAN_GNO_TYPE>& fine, bool in_coarse){
	std::unordered_set<ZOLTAN_GNO_TYPE> neigh;
	 /* for every hyperedge containing the vertex vtx */
	for (ZOLTAN_GNO_TYPE i = hg->vindex[vtx]; i < hg->vindex[vtx+1]; i++) {
		ZOLTAN_GNO_TYPE e = hg->vedge[i];
		/* for every other vertex in the hyperedge */
		for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
			ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];				
			if(in_coarse){
				// vtx in coarse
				if(vtx == v2 || coarse.find(v2) != coarse.end()){ continue; }
				neigh.insert(v2);
			}else{
				// vtx in fine
				if(vtx == v2 || fine.find(v2) != fine.end()){ continue; }
				neigh.insert(v2);
			}
		}	
	}	
	return neigh.size();
}

unsigned int RS_get_max_vertex_weight(HGraph* hg){
	ZOLTAN_GNO_TYPE v;	
	float max = -1.0;
	for (v = 0; v < hg->nVtx; v++){
		if (hg->vwgt[v] > max) { max = hg->vwgt[v]; }
	}
	return (unsigned int) max;
}

// sees if there's enough neighbors from a given set
bool RS_strongly_coupled_to_set(HGraph *hg, ZOLTAN_GNO_TYPE v1, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, const std::unordered_set<ZOLTAN_GNO_TYPE>& coarse){
	ZOLTAN_GNO_TYPE v2, e;
	RS_ALG_COORD_TYPE sum_to_c = 0.0; // sum of alg_dist weights to coarse vertices
	RS_ALG_COORD_TYPE sum_all = 0.0; // sum of alg_dist weights to coarse vertices
	int i, j;
	
	/* for every hyperedge containing the vertex v1 */
	for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
		e = hg->vedge[i];
		/* for every other vertex in the hyperedge */
		for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
			v2 = hg->hvertex[j];
			if (v2 == v1) { continue; }
			sum_all += new_hedge_weights[e];
			if (coarse.find(v2) != coarse.end()){ 
				sum_to_c += new_hedge_weights[e];
			}
		}	
	}
	if ((sum_to_c / sum_all) > 0.5){
		return true;
	}else{
		return false;
	}
}

// assumes vertex_alg_coords is an R by hg->nVtx sized 2d array
int RS_compute_future_volumes(HGraph *hg, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, std::vector<RS_ALG_COORD_TYPE>& future_volumes, const std::unordered_set<ZOLTAN_GNO_TYPE>& coarse){
	char* yo = "RS_compute_future_volumes";
#ifdef RS_INTEGER_FUTURE_VOLUMES 
	std::cout << yo << "-> Using integer future volumes" << std::endl;
	ZOLTAN_GNO_TYPE v1, v2, e;
	RS_ALG_COORD_TYPE d_sum;

	/* Now each vertex just gives itself to the most connected neighbor */

	for (v1 = 0; v1 < hg->nVtx; v1++){
		if (coarse.find(v1) != coarse.end()){
				continue;
		}
		int best = -1;
		float best_w = 0;
		std::unordered_map<ZOLTAN_GNO_TYPE, float> ips;
		/* for every hyperedge containing the vertex v1 */
		for (ZOLTAN_GNO_TYPE i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			ZOLTAN_GNO_TYPE e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];
				if (v2 == v1 || coarse.find(v2) != coarse.end()) { continue; }
				ips[v2] += (hg->ewgt ? hg->ewgt[e] : 1.0);
			}	
		}
		for (auto &i : ips){
			if (i.second > best_w){
				best = i.first;
				best_w = i.second;
			}
		}
		if (best == -1 || best_w == 0){
			//std::cout << "Couldn't find best partner for " << v1 << " coarse.size() " << coarse.size() << " fine.size() " << fine.size() << std::endl;
			//return ZOLTAN_FATAL;
			continue;
		}
		future_volumes[best] += (hg->vwgt ? hg->vwgt[v1] : 1.0);
	}
	
	return ZOLTAN_OK;	
#else
	std::cout << yo << "-> Using fractional future volumes" << std::endl;
	int i, j;
	ZOLTAN_GNO_TYPE v1, v2, e;
	RS_ALG_COORD_TYPE d_sum;

	/* first, precompute \Sigma_i alg_dist_ij for each vertex j */

	std::vector<RS_ALG_COORD_TYPE> precomp_sums;
	precomp_sums.reserve(hg->nVtx);
	for (v1 = 0; v1 < hg->nVtx; v1++){
		if (coarse.find(v1) != coarse.end()){
			precomp_sums.push_back(-1.0);
			continue;
		}
		d_sum = 0.0;
		/* for every hyperedge containing the vertex v1 */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1 || coarse.find(v2) != coarse.end()) { continue; }
				d_sum += (hg->ewgt ? hg->ewgt[e] : 1.0);
			}	
		}
		precomp_sums.push_back(d_sum);
	}

	for (v1 = 0; v1 < hg->nVtx; v1++){
		if (coarse.find(v1) != coarse.end()){
			future_volumes.push_back(-1.0);
			continue;
		}
		RS_ALG_COORD_TYPE fv = (RS_ALG_COORD_TYPE)(hg->vwgt ? hg->vwgt[v1] : 1.0);
		/* for every hyperedge containing the vertex v */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1) { continue; }
				/* TODO put more code here */
				fv += (hg->vwgt ? hg->vwgt[v2] : 1.0)*((hg->ewgt ? hg->ewgt[e] : 1.0) / precomp_sums[v2]); // if made [1] weighted, add weights here too
			}
		}
		future_volumes.push_back(fv);
	}
	
	return ZOLTAN_OK;	
#endif
}

int RS_compute_preferences(ZZ *zz, HGraph *hg, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_f, std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_c, const std::unordered_set<ZOLTAN_GNO_TYPE> &coarse, std::unordered_set<ZOLTAN_GNO_TYPE> &fine){
	ZOLTAN_GNO_TYPE i,j,v2,e;
	for (auto &v1 : coarse){
		std::unordered_map<ZOLTAN_GNO_TYPE, RS_ALG_COORD_TYPE> rm;
		std::multimap<RS_ALG_COORD_TYPE, ZOLTAN_GNO_TYPE> m;
		/* for every hyperedge containing the vertex v1 */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1 || coarse.find(v2) != coarse.end()) { continue; }
				rm[v2] -= new_hedge_weights[e];
			}	
		}
		for (auto &i : rm){
			m.emplace(i.second, i.first);
		}
		std::vector<ZOLTAN_GNO_TYPE> v;
		v.reserve(m.size());
		for (auto &i : m){
			v.push_back(i.second);
		}
		pref_dict_c[v1] = v;
	}
	
	for (auto &v1 : fine){
		std::unordered_map<ZOLTAN_GNO_TYPE, RS_ALG_COORD_TYPE> rm;
		std::multimap<RS_ALG_COORD_TYPE, ZOLTAN_GNO_TYPE> m;
		/* for every hyperedge containing the vertex v1 */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1 || fine.find(v2) != fine.end()) { continue; }
				rm[v2] -= new_hedge_weights[e];
			}	
		}
		for (auto &i : rm){
			m.emplace(i.second, i.first);
		}
		std::vector<ZOLTAN_GNO_TYPE> v;
		v.reserve(m.size());
		for (auto &i : m){
			v.push_back(i.second);
		}
		pref_dict_f[v1] = v;
	}	
	return ZOLTAN_OK;	
}

int RS_amg_matching(ZZ *zz, HGraph *hg, PHGPartParams *hgp, ZOLTAN_GNO_TYPE *match, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R){
	int err;
	char* yo = "RS_amg_matching";

	if (rs_total_vwgt == -1.0){
			std::cout << "Computing total vertex weight...";
			if (hg->vwgt){
					std::cout << " the slow way" << std::endl;
					rs_total_vwgt = 0.0;
					for (int i = 0; i < hg->nVtx; i++){
							rs_total_vwgt += hg->vwgt[i];
					}
			}else{
					std::cout << std::endl;
					rs_total_vwgt = hg->nVtx;
			}
			std::cout << "Total vertex weight = " << rs_total_vwgt << std::endl;
	}

	/* ----------- stage one: compute seeds ------------- */

	std::unordered_set<ZOLTAN_GNO_TYPE> coarse, fine;

#ifdef RS_INTEGER_FUTURE_VOLUMES 
	std::vector<RS_ALG_COORD_TYPE> future_volumes(hg->nVtx, 0.0);
#else
	std::vector<RS_ALG_COORD_TYPE> future_volumes;
	future_volumes.reserve(hg->nVtx);
#endif	
	err = RS_compute_future_volumes(hg, vertex_alg_coords, new_hedge_weights, R, future_volumes, coarse);
	
	if (err != ZOLTAN_OK) { return ZOLTAN_FATAL; }

	/* Put top two standard deviations over mean in seeds */

	RS_ALG_COORD_TYPE fv_sum = std::accumulate(future_volumes.begin(), future_volumes.end(), 0.0);
	RS_ALG_COORD_TYPE fv_mean = fv_sum / future_volumes.size();

	std::vector<RS_ALG_COORD_TYPE> fv_diff(future_volumes.size());
	std::transform(future_volumes.begin(), future_volumes.end(), fv_diff.begin(), [fv_mean](double x) { return x - fv_mean; });
	RS_ALG_COORD_TYPE fv_sq_sum = std::inner_product(fv_diff.begin(), fv_diff.end(), fv_diff.begin(), 0.0);
	RS_ALG_COORD_TYPE fv_stdev = std::sqrt(fv_sq_sum / future_volumes.size());

	//std::cout << "Future volumes mean: " << fv_mean << " standard deviation: " << fv_stdev << std::endl;

	RS_ALG_COORD_TYPE fv_threshold = fv_mean + 2 * fv_stdev;

	for (int v = 0; v < future_volumes.size(); v++){
			if (future_volumes[v] > fv_threshold){
				coarse.insert(v);
			}else{
				fine.insert(v);
			}
	}
	//std::cout << "After two standard deviations over mean: " << future_volumes.size() << " coarse: " << coarse.size()  << " fine: " << fine.size() << std::endl;	

	/* Recompute future volumes */

#ifdef RS_INTEGER_FUTURE_VOLUMES 
	std::fill(future_volumes.begin(), future_volumes.end(), 0.0);
#else
	future_volumes.clear();
#endif
	
	err = RS_compute_future_volumes(hg, vertex_alg_coords, new_hedge_weights, R, future_volumes, coarse);
	
	if (err != ZOLTAN_OK) { return ZOLTAN_FATAL; }

	/* Now, traverse fine vertices and add them to coarse if needed */

	for (auto v : sort_indexes(future_volumes)){ // order can be different
		if(coarse.find(v) != coarse.end()) { continue; }
		if(!RS_strongly_coupled_to_set(hg, v, vertex_alg_coords, new_hedge_weights, R, coarse)){
			coarse.insert(v);
			fine.erase(fine.find(v));	
		}	
	}

	/* ------------- stage two ---------------------- */

	//std::cout << "After strongly connected: " << future_volumes.size() << " coarse: " << coarse.size()  << " fine: " << fine.size() << std::endl;	

#ifdef RS_STABLE_MATCHING
	/* first, compute preferences and capacities */
	std::cout << yo << "-> Using stable matching" << std::endl;

	std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> pref_dict_f, pref_dict_c;
	
	RS_compute_preferences(zz, hg, vertex_alg_coords, new_hedge_weights, R, pref_dict_f, pref_dict_c, coarse, fine);
	
	std::unordered_map<ZOLTAN_GNO_TYPE, unsigned int> capacities_c; 
	
	unsigned int capacity = RS_get_max_vertex_weight(hg)*((fine.size() / coarse.size()) + 2) + 10;
	
	//std::cout << "Assigned capacity " << capacity << " with max vertex_weight " << RS_get_max_vertex_weight(hg) << std::endl;
	
	for (auto &i : coarse){
		capacities_c[i] = capacity;
	}


	/* second, compute stable assignment */

	std::unordered_map<ZOLTAN_GNO_TYPE, RS_set> assignment = RS_compute_stable_matching(hg, pref_dict_c, pref_dict_f, capacities_c);
	
	unsigned int matched_num = 0;
	std::unordered_map<ZOLTAN_GNO_TYPE, float> clwgts; // cluster weights

	for (auto &i : assignment){
		for(auto &v : i.second){
			float clwgt = hg->vwgt[v.first] + hg->vwgt[i.first]; // prospective cluster weight
			if (clwgts.find(i.first) != clwgts.end()){
				clwgt += clwgts[i.first];
			}
			if (clwgt < rs_total_vwgt / zz->LB.Num_Global_Parts){
				match[v.first] = i.first; // have to unpack pair (vertex, weight)
				matched_num++;
				clwgts[i.first] += hg->vwgt[v.first];
			}else{
				match[v.first] = v.first;
				std::cout << "Rejected match for " << v.first << " due to imbalance constraint" <<  std::endl;
			}
		}
	}

	for(auto &v: coarse){
        if (assignment.find(v) == assignment.end()){
            match[v] = v;
        }
    }

#else
	/* compute agglomerative matching */
	std::cout << yo << "-> Using agglomerative matching" << std::endl;
	
	unsigned int matched_num = 0;

	for (int i = 0; i < hg->nVtx; i++){
		match[i] = i;
	}

	std::unordered_map<ZOLTAN_GNO_TYPE, float> clwgts; // cluster weights

	for (auto &v1 : fine){

		int best = -1;
		float best_w = 0;
		std::unordered_map<ZOLTAN_GNO_TYPE, float> ips;
		/* for every hyperedge containing the vertex v1 */
		for (ZOLTAN_GNO_TYPE i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			ZOLTAN_GNO_TYPE e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];
				if (v2 == v1) { continue; }
				if (fine.find(v2) != fine.end()){
				//	if (match[v2] == v2) { continue; }
				//	ips[match[v2]] += new_hedge_weights[e];
					continue;
				}
				ips[v2] += new_hedge_weights[e];
			}	
		}
		for (auto &i : ips){
			float clwgt = hg->vwgt[v1] + hg->vwgt[i.first]; // prospective cluster weight
			if (clwgts.find(i.first) != clwgts.end()){
					clwgt += clwgts[i.first];
			}
			//i.second /= clwgt;
			if (i.second > best_w && (clwgt < rs_total_vwgt / zz->LB.Num_Global_Parts)){
				best = i.first;
				best_w = i.second;
			}
		}
		if (best == -1 || best_w == 0){
			std::cout << "Couldn't find a match for " << v1 << " with weight " << hg->vwgt[v1] << " out of " << ips.size() << " candidates" << std::endl;
			match[v1] = v1;
		}else{
			// all is well
			match[v1] = best;
			clwgts[best] += hg->vwgt[v1];
		}
	}

	//std::cout << "clwgts.size() = " << clwgts.size() << " coarse.size() = " << coarse.size() << std::endl;

	//std::cout << "Number of vertices matched: " << matched_num << std::endl;

#endif

	return ZOLTAN_OK;
}

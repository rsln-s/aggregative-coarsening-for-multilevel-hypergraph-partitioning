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

int RS_get_conditional_degree(HGraph* hg, int vtx, const std::vector<bool>& is_coarse, bool in_coarse){
	std::unordered_set<ZOLTAN_GNO_TYPE> neigh;
	 /* for every hyperedge containing the vertex vtx */
	for (ZOLTAN_GNO_TYPE i = hg->vindex[vtx]; i < hg->vindex[vtx+1]; i++) {
		ZOLTAN_GNO_TYPE e = hg->vedge[i];
		/* for every other vertex in the hyperedge */
		for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
			ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];				
			if(in_coarse){
				// vtx in coarse
				if(vtx == v2 || is_coarse[v2]){ continue; }
				neigh.insert(v2);
			}else{
				// vtx in fine
				if(vtx == v2 || !is_coarse[v2]){ continue; }
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
bool RS_strongly_coupled_to_set(HGraph *hg, ZOLTAN_GNO_TYPE v1, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, const std::vector<bool>& is_coarse){
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
			if (is_coarse[v2]){ 
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
int RS_compute_future_volumes(HGraph *hg, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, std::vector<RS_ALG_COORD_TYPE>& future_volumes, const std::vector<bool>& is_coarse){
	char* yo = "RS_compute_future_volumes";
#ifdef RS_INTEGER_FUTURE_VOLUMES 
	std::cout << yo << "-> Using integer future volumes" << std::endl;
	ZOLTAN_GNO_TYPE v1, v2, e;
	RS_ALG_COORD_TYPE d_sum;

	/* Now each vertex just gives itself to the most connected neighbor */

	for (v1 = 0; v1 < hg->nVtx; v1++){
		if (is_coarse[v1]){
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
				if (v2 == v1 || is_coarse[v2]) { continue; }
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
		if (is_coarse[v1]){
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
				if (v2 == v1 || is_coarse[v2]) { continue; }
				d_sum += (hg->ewgt ? hg->ewgt[e] : 1.0);
			}	
		}
		precomp_sums.push_back(d_sum);
	}
	
	for (e = 0; e < hg->nEdge; e++){
		float curr_edgew = (hg->ewgt ? hg->ewgt[e] : 1.0);
		/* for every vertex in the hyperedge */
		for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
			v1 = hg->hvertex[j];
			if (is_coarse[v1]){ continue; }
			/* for every other vertex in the hyperedge */
			for (i = hg->hindex[e]; i < hg->hindex[e+1]; i++) {
				v2 = hg->hvertex[i];
				if (v2 == v1) { continue; }
				future_volumes[v2] += (hg->vwgt ? hg->vwgt[v1] : 1.0) * (curr_edgew / precomp_sums[v1]);
			}	
		}	
	}

	return ZOLTAN_OK;	
#endif
}

int RS_compute_preferences(ZZ *zz, HGraph *hg, RS_ALG_COORD_TYPE** vertex_alg_coords, float* new_hedge_weights, int R, std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_f, std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_c, const std::vector<bool>& is_coarse){
	ZOLTAN_GNO_TYPE i,j,v2,e;
	for (ZOLTAN_GNO_TYPE v1 = 0; v1 < hg->nVtx; v1++){
		if (!is_coarse[v1]){ continue; } 
		std::unordered_map<ZOLTAN_GNO_TYPE, RS_ALG_COORD_TYPE> rm;
		std::multimap<RS_ALG_COORD_TYPE, ZOLTAN_GNO_TYPE> m;
		/* for every hyperedge containing the vertex v1 */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1 || is_coarse[v2]) { continue; }
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
	
	for (ZOLTAN_GNO_TYPE v1 = 0; v1 < hg->nVtx; v1++){
		if (is_coarse[v1]){ continue; }
		std::unordered_map<ZOLTAN_GNO_TYPE, RS_ALG_COORD_TYPE> rm;
		std::multimap<RS_ALG_COORD_TYPE, ZOLTAN_GNO_TYPE> m;
		/* for every hyperedge containing the vertex v1 */
		for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				v2 = hg->hvertex[j];
				if (v2 == v1 || !is_coarse[v2]) { continue; }
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
	int num_coarse = -1;
	int num_fine = -1;
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

	std::vector<bool> is_coarse(hg->nVtx);

	std::vector<RS_ALG_COORD_TYPE> future_volumes(hg->nVtx, 0.0);
	err = RS_compute_future_volumes(hg, vertex_alg_coords, new_hedge_weights, R, future_volumes, is_coarse);
	
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
				is_coarse[v] = true;
			}else{
				is_coarse[v] = false;
			}
	}
//	num_coarse = std::accumulate(is_coarse.begin(), is_coarse.end(), 0);	
//	num_fine = hg->nVtx - num_coarse;
//	std::cout << "After two standard deviations over mean: " << future_volumes.size() << " coarse: " << num_coarse  << " fine: " << num_fine << std::endl;	

	/* Recompute future volumes */

	std::fill(future_volumes.begin(), future_volumes.end(), 0.0);
	
	err = RS_compute_future_volumes(hg, vertex_alg_coords, new_hedge_weights, R, future_volumes, is_coarse);
	
	if (err != ZOLTAN_OK) { return ZOLTAN_FATAL; }

	/* Now, traverse fine vertices and add them to coarse if needed */

	for (auto v : sort_indexes(future_volumes)){ // order can be different
		if(is_coarse[v]) { continue; }
		if(!RS_strongly_coupled_to_set(hg, v, vertex_alg_coords, new_hedge_weights, R, is_coarse)){
			is_coarse[v] = true;
		}	
	}

	/* ------------- stage two ---------------------- */

//	num_coarse = std::accumulate(is_coarse.begin(), is_coarse.end(), 0);	
//	num_fine = hg->nVtx - num_coarse;
//	std::cout << "After strongly connected: " << future_volumes.size() << " coarse: " << num_coarse  << " fine: " << num_fine << std::endl;	

#ifdef RS_STABLE_MATCHING
	/* first, compute preferences and capacities */
	std::cout << yo << "-> Using stable matching" << std::endl;

	std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> pref_dict_f, pref_dict_c;
	
	RS_compute_preferences(zz, hg, vertex_alg_coords, new_hedge_weights, R, pref_dict_f, pref_dict_c, is_coarse);
	
	std::unordered_map<ZOLTAN_GNO_TYPE, unsigned int> capacities_c; 

	int num_coarse = std::accumulate(is_coarse.begin(), is_coarse.end(), 0);	
	int num_fine = hg->nVtx - num_coarse;
	unsigned int capacity = RS_get_max_vertex_weight(hg)*(( num_fine / num_coarse) + 2) + 10;
	
	//std::cout << "Assigned capacity " << capacity << " with max vertex_weight " << RS_get_max_vertex_weight(hg) << std::endl;
	
	for (ZOLTAN_GNO_TYPE i = 0; i < hg->nVtx; i++){
		if (!is_coarse[i]) { continue; }
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

	for(ZOLTAN_GNO_TYPE v = 0; v < hg->nVtx; v++){
		if (!is_coarse[v]){ continue; }
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

	std::vector<float> clwgts(hg->nVtx); // cluster weights
	std::vector<float> ips(hg->nVtx); // inner products
	std::vector<float> adj(hg->nVtx);
	ZOLTAN_GNO_TYPE n; // number of neighbors

	for (ZOLTAN_GNO_TYPE v1 = 0; v1 < hg->nVtx; v1++){
		if (is_coarse[v1]) { continue; }
		int best = -1;
		float best_w = 0;
		n = 0;

		/* for every hyperedge containing the vertex v1 */
		for (ZOLTAN_GNO_TYPE i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
			ZOLTAN_GNO_TYPE e = hg->vedge[i];
			/* for every other vertex in the hyperedge */
			for (ZOLTAN_GNO_TYPE j = hg->hindex[e]; j < hg->hindex[e+1]; j++) {
				ZOLTAN_GNO_TYPE v2 = hg->hvertex[j];
				if (v2 == v1) { continue; }
				if (!is_coarse[v2]){
				//	if (match[v2] == v2) { continue; }
				//	ips[match[v2]] += new_hedge_weights[e];
					continue;
				}
				if (ips[v2] == 0.0)
						adj[n++] = v2;
				ips[v2] += new_hedge_weights[e];
			}	
		}

		for (int i = 0; i < n; i++){
			ZOLTAN_GNO_TYPE v2 = adj[i];
			float clwgt = hg->vwgt[v1] + hg->vwgt[v2]; // prospective cluster weight
			if (clwgts[v2] != 0){
					clwgt += clwgts[v2];
			}
			if (ips[v2] > best_w && (clwgt < rs_total_vwgt / zz->LB.Num_Global_Parts)){
				best = v2;
				best_w = ips[v2];
			}
			ips[v2] = 0.0;
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

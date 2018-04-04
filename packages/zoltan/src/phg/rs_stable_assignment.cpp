#include "rs_stable_assignment.hpp"

int propose(HGraph* hg, const ZOLTAN_GNO_TYPE key, std::unordered_map<ZOLTAN_GNO_TYPE, RS_set> &propos_c, std::unordered_map<ZOLTAN_GNO_TYPE, ZOLTAN_GNO_TYPE> &propos_f, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_c, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_f, const std::unordered_map<ZOLTAN_GNO_TYPE, unsigned int> &capacities_c){
	/* propose to everyone in the preference list pref_dict[key] in the order of preference */
	for (auto &pref : pref_dict_c.at(key)){
		if((propos_c[key].size() + hg->vwgt[pref]) > capacities_c.at(key)){
			return 0;
		}
		/* current proposal held by the pref, to whom we are proposing */
		auto curr_propos = propos_f.find(pref);
		if (curr_propos == propos_f.end()){
			/* if no proposal help by pref, key is now held by pref */
			propos_f[pref] = key;
			propos_c[key].emplace(pref, hg->vwgt[pref]);
			continue;
		}
		/* if pref had a proposal, check if it's better than what we're offering */
		if (std::distance(pref_dict_f.at(pref).begin(), std::find(pref_dict_f.at(pref).begin(), pref_dict_f.at(pref).end(), key)) < std::distance(pref_dict_f.at(pref).begin(), std::find(pref_dict_f.at(pref).begin(), pref_dict_f.at(pref).end(), propos_f[pref]))){
			ZOLTAN_GNO_TYPE rejected = propos_f[pref];
			propos_c[rejected].erase(pref);
			propos_f[pref] = key;
			propos_c[key].emplace(pref, hg->vwgt[pref]);
			propose(hg, rejected, propos_c, propos_f, pref_dict_c, pref_dict_f, capacities_c);
		}
	}	
	//std::cout << key << " got rejected by everybody!" << std::endl;
	return -1;
}

std::unordered_map<ZOLTAN_GNO_TYPE, RS_set> RS_compute_stable_matching(HGraph* hg, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_c, const std::unordered_map<ZOLTAN_GNO_TYPE, std::vector<ZOLTAN_GNO_TYPE>> &pref_dict_f, const std::unordered_map<ZOLTAN_GNO_TYPE, unsigned int> &capacities_c){

	std::unordered_map<ZOLTAN_GNO_TYPE, RS_set> propos_c; // start with an empty dictionary and only include els if somebody proposed to them
	std::unordered_map<ZOLTAN_GNO_TYPE, ZOLTAN_GNO_TYPE> propos_f;
	
	/* fine vertices propose to coarse ones */
	for (const auto &i : pref_dict_c){
		propose(hg, i.first, propos_c, propos_f, pref_dict_c, pref_dict_f, capacities_c);
	}

	return propos_c;
}


#include "Hash_indeces.h"
#include <NTL/tools.h>

//std::vector<vec_ZZ> get_hash_indeces(const vec_ZZ& x) {
//	std::vector < vec_ZZ > indeces;
//	indeces.push_back(vec_ZZ());
//	std::vector < vec_ZZ > indeces2;
//	for (auto& x_i : x) {
//		if (x_i > ZZ(0)) {
//			for (auto& ind : indeces) {
//				ind.append(ZZ(1));
//			}
//			continue;
//		}
//		if (x_i < ZZ(0)) {
//			for (auto& ind : indeces) {
//				ind.append(ZZ(0));
//			}
//			continue;
//		}
//		indeces2 = indeces;
//		for (auto& ind : indeces)
//			ind.append(ZZ(1));
//		for (auto& ind : indeces2)
//			ind.append(ZZ(0));
//		//indeces = indeces + indeces2 -> concatenate
//		for (auto& ind : indeces2) {
//			indeces.push_back(ind);
//		}
//	}
//
//	indeces2.clear();
//	indeces2.reserve(indeces.size());
//	for (auto& ind : indeces) {
//		vec_ZZ vec(INIT_SIZE, ind.length());
//		for (int i = 0; i < ind.length(); ++i) {
//			vec[i] = ZZ(1) - ind[i];
//		}
//		indeces2.push_back(std::move(vec));
//	}
//
//	for (auto& ind : indeces2) {
//		indeces.push_back(ind);
//	}
//	return indeces;
//}

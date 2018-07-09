#pragma once
#include <vector>
#include <NTL/vec_RR.h>
#include <bitset>
#include "Utils.h"
using namespace NTL;

template<typename T, typename S>
std::vector<std::bitset<DIM>> get_hash_indeces(const T& x) {
	std::vector<T> indeces;
	indeces.push_back(T());
	std::vector<T> indeces2;
	for (auto& x_i : x) {
		if (x_i > S(0)) {
			for (auto& ind : indeces) {
				ind.append(S(1));
			}
			continue;
		}
		if (x_i < S(0)) {
			for (auto& ind : indeces) {
				ind.append(S(0));
			}
			continue;
		}
		indeces2 = indeces;
		for (auto& ind : indeces)
			ind.append(S(1));
		for (auto& ind : indeces2)
			ind.append(S(0));
		//indeces = indeces + indeces2 -> concatenate
		for (auto& ind : indeces2) {
			indeces.push_back(ind);
		}
	}

	indeces2.clear();
	indeces2.reserve(indeces.size());
	for (auto& ind : indeces) {
		T vec(INIT_SIZE, ind.length());
		for (int i = 0; i < ind.length(); ++i) {
			vec[i] = S(1) - ind[i];
		}
		indeces2.push_back(std::move(vec));
	}

	for (auto& ind : indeces2) {
		indeces.push_back(ind);
	}

//convert to bitset
	std::vector < std::bitset < DIM >> res;
	for (auto&it : indeces) {
		std::bitset < DIM > vec;
		for (size_t i = 0; i < DIM; ++i) {
			if (it[i] == S(1))
				vec[i] = 1;
		}
		res.push_back(vec);
	}
	return res;
}

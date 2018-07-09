#pragma once
#include <NTL/vec_ZZ.h>
#include <vector>
using namespace NTL;

std::vector<vec_ZZ> get_index_sets(int r, int c, int min_val = 0, const vec_ZZ& index_set = vec_ZZ());
std::vector<vec_ZZ> get_v_vals(int r, int c);
vec_ZZ pad_with_zeros(const vec_ZZ& v);
std::vector<vec_ZZ> get_v_vals_padded(int r, int c);
std::vector<vec_ZZ> get_index_sets_new(int c, const std::vector<vec_ZZ>& indeces, int min_val = 0,
		const std::vector<vec_ZZ>& index_set = std::vector<vec_ZZ>());
std::vector<vec_ZZ> get_valid_v_1(const vec_ZZ& v);

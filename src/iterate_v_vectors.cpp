#include "iterate_v_vectors.h"
#include <utility>

std::vector<vec_ZZ> get_index_sets(int r, int c, int min_val,
		const vec_ZZ& index_set) {
	std::vector < vec_ZZ > res;
	if (c == index_set.length()) {
		res.push_back(index_set);
	} else {
		for (int i = min_val; i < r - c + index_set.length() + 1; ++i) {
			std::vector<vec_ZZ> res_recursive = get_index_sets(r, c, i + 1);
			std::move(res_recursive.begin(), res_recursive.end(), std::back_inserter(res));
			res_recursive.erase(res_recursive.begin(), res_recursive.end());
		}
	}
	return res;
}

std::vector<vec_ZZ> get_v_vals(int r, int c){
	std::vector<vec_ZZ> index_sets = get_index_sets(r,c);
	std::vector<vec_ZZ> v_vals;
	for (auto& ind : index_sets){
		vec_ZZ v_val(INIT_SIZE, r);
		for (auto& i : ind){
			int j = conv<int>(i);
			v_val[j] = ZZ(1);
		}
		v_vals.push_back(v_val);
	}
	return v_vals;
}

vec_ZZ pad_with_zeros(const vec_ZZ& v){
	vec_ZZ vals(INIT_SIZE, 2*v.length());
	int pos = 1;
	for (auto& i : v){
		vals[pos] = i;
		pos += 2;
	}
	return vals;
}

std::vector<vec_ZZ> get_v_vals_padded(int r, int c){
	std::vector<vec_ZZ> v_vals = get_v_vals(r/2, c);
	std::vector<vec_ZZ> res;
	res.reserve(v_vals.size());
	for (auto& v: v_vals){
		res.push_back(pad_with_zeros(v));
	}
	return res;
}

std::vector<vec_ZZ> get_index_sets_new(int c, const std::vector<vec_ZZ>& indeces, int min_val,
		const std::vector<vec_ZZ>& index_set){
	int r = indeces.size();
	if (c == index_set.size()){
		return index_set;
	}
	std::vector<vec_ZZ> results;
	int len;
	// TODO: attention!
	//len(indeces[min_val:r-c+len(index_set)+1])
	if (r-c+index_set.size()+1<= indeces.size()){
		len = r-c+index_set.size()+1 - min_val + 1;
	} else {
		len = indeces.size() - min_val +1;
	}

	for (int i = 0; i < len; ++i){
		vec_ZZ val(INIT_SIZE, len);
		for (int j = 0; j < len; ++j){
			val[j] = indeces[j][i];
		}
		//TODO :optimize init of index_set_
		std::vector<vec_ZZ> index_set_ = index_set;
		index_set_.push_back(val);
		auto tmp = get_index_sets_new(c, indeces, min_val + i + 1, index_set_);
		results.insert(results.end(), tmp.begin(), tmp.end());
	}
	return results;
}

std::vector<vec_ZZ> get_valid_v_1(const vec_ZZ& v){
	vec_ZZ indeces_;
	for (int i = 0; i < v.length(); i++){
		if (v[i] == ZZ(1)){
			indeces_.append(ZZ(i));
		}
	}
	int r = indeces_.length();
	int c = r/2;
	std::vector<vec_ZZ> indeces = get_index_sets_new(c, std::vector<vec_ZZ>(1,indeces_));
	std::vector<vec_ZZ> v_vals;
	v_vals.reserve(indeces.size());
	for (auto& ind: indeces){
		vec_ZZ v_val(INIT_SIZE, v.length());
		for (auto& i : ind){
			int j = conv<int>(i);
			v_val[j] = ZZ(1);
		}
		v_vals.push_back(v_val);
	}
	return v_vals;
}

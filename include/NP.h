#pragma once
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <limits.h>
#include "Utils.h"
using namespace NTL;

extern int n;
extern int m;
extern int r;
extern double NP_call_duration;

static mat_RR zero_mat(INIT_SIZE, m - r, m - r);
/**
 * B: B_1_reduced_transposed in attack
 */
//template<typename T, typename S>
//S NP(const T& B, const S& t, mat_RR B_gs = zero_mat, int k = INT_MAX) {
//	if (B_gs == zero_mat) {
//		B_gs = GSO(B);
//	}
//	if (k == INT_MAX) {
//		k = B.NumRows() - 1;
//	}
//	if (k == -1) {
//		return S(INIT_SIZE, B.NumRows());
//	}
//	RR B_gs_t_innerProd;
//	InnerProduct(B_gs_t_innerProd, B_gs[k], t);
//	RR B_gs_self_innerProd;
//	InnerProduct(B_gs_self_innerProd, B_gs[k], B_gs[k]);
//	RR c = round(B_gs_t_innerProd / B_gs_self_innerProd);
//	return (c * B[k]) + NP<T, S>(B, t - (c * B[k]), B_gs, k - 1);
//}


// iterative version
//template<typename T, typename S>
//S NP(const T& B, const S& t, mat_RR B_gs = zero_mat, int k = INT_MAX) {
//	if (B_gs == zero_mat) {
//		B_gs = GSO(B);
//	}
//	if (k == INT_MAX) {
//		k = B.NumRows() - 1;
//	}
//
//	S v = t;
//	for (int j = k; j >= 0; j--) {
//		RR B_gs_t_innerProd;
//		InnerProduct(B_gs_t_innerProd, B_gs[j], v);
//		RR B_gs_self_innerProd;
//		InnerProduct(B_gs_self_innerProd, B_gs[j], B_gs[j]);
//		RR c = round(B_gs_t_innerProd / B_gs_self_innerProd);
//
//		v = v - c*B[j];
//	}
//	return t-v;
//}

vec_RR NP(const mat_RR& B, const vec_RR& t, const mat_RR& B_gs, int k = INT_MAX);

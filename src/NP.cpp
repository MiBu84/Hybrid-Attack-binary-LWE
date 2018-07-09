#include "NP.h"
#include <omp.h>

vec_RR NP(const mat_RR& B, const vec_RR& t, const mat_RR& B_gs, int k) {
//	if (B_gs == zero_mat) {
//		B_gs = GSO(B);
//	}
	double startTimer = omp_get_wtime();

	if (k == INT_MAX) {
		k = B.NumRows() - 1;
	}

	vec_RR v = t;

	for (int j = k; j >= 0; j--) {
		RR B_gs_t_innerProd;
		InnerProduct(B_gs_t_innerProd, B_gs[j], v);

		RR B_gs_self_innerProd;
		InnerProduct(B_gs_self_innerProd, B_gs[j], B_gs[j]);

		RR c = round(B_gs_t_innerProd / B_gs_self_innerProd);

		v = v - c * B[j];
	}

#pragma omp atomic
	NP_call_duration+=omp_get_wtime()-startTimer;

	return t-v;
}


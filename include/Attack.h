#pragma once
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_GF2.h>
#include <vector>
using namespace NTL;

extern int m;
extern int n;
extern int q;
extern int r;
extern int c;
extern long beta;
extern mat_ZZ A;
extern vec_ZZ b;
extern vec_ZZ error;
extern int testcase;
extern long prune;
extern int numThread;
extern std::string input;
extern std::string output;
extern std::vector<int> permutation;

#ifdef USING_MPI
#include "mpi.h"
extern int world_rank;
extern int world_size;
extern int precomputing_reduction_iterCount;
const int PRECOMPUTING_FINISH_TAG = 9990;
const int ATTACK_FINISH_TAG = 9999;
#include <future>
#endif

/**
 * returns a binary vector of length r with exactly c bits 1.
 */
vec_RR sample_bin_vector_RR(int r, int c);

/**
 * output A and b.
 */
void get_instance(mat_ZZ& A, vec_ZZ& b, int n, int m, int q);

/**
 * attack is divided into 2 subroutines: precomputing and hybridAttack
 */
void attack(mat_ZZ& A, vec_ZZ& b, int q, int r, int c);

/**
 * @return
 * B
 * B_not_updated_upperleftCorner
 * B_reduced_gs is GS-basis obtained from B_1_red (B_1_red is B' in paper)
 */
//void precomputing(mat_ZZ& A, vec_ZZ& b, int q, int r, long beta,
//		mat_RR& B, mat_RR& B_1_reduced_transposed, mat_RR& B_reduced_gs);
//
//vec_RR hybridAttack(const mat_RR& B, const mat_RR& B_1_reduced_transposed,
//		const mat_RR& B_reduced_gs, int r, int c);




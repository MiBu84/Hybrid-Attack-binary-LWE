#pragma once
#define DIM 140 // = m -r

#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <bitset>
#include <random>
using namespace NTL;

// for setRow
#include <climits>

// for permutation
#include <algorithm>

// for cpu freq
#include <thread>

extern int m;
extern int n;
extern int q;
extern int r;
extern int c;
extern long beta;
extern mat_ZZ A;
extern vec_ZZ b;
extern vec_ZZ error;
extern std::vector<int> permutation;
extern std::string output;
extern std::string input;

// return GS-orthogonal basis from the input basis
mat_RR GSO(const mat_RR& in);
mat_RR GSO(const mat_ZZ& in);

void read(std::string file, mat_ZZ& A, vec_ZZ& b); // read input matrix/vector from a text
void read(mat_RR& mat, std::string file);
void read_v1(std::vector<vec_RR>& vec, std::string file); // read v1 from file without sampling

void setRow(vec_ZZ& a, std::string line, int len);
void setRow(vec_RR& a, std::string line, int len); // read output of sage/cpp NTL
void setRow(std::vector<int>& permutation, std::string line, int len);

//TODO: optimize it
template<typename T>
T submatrix(const T& input, int beginRow, int beginCol, int numRows,
		int numCols) {
	T res(INIT_SIZE, numRows, numCols);
	for (int i = 0; i < numRows; ++i)
		for (int j = 0; j < numCols; ++j) {
			res[i][j] = input[i + beginRow][j + beginCol];
		}
	return res;
}

#ifndef USING_TBB
template<typename T>
using HashMap = std::unordered_map<std::bitset<DIM 96-4-4>, std::vector<T>>;
#else
#include "tbb/concurrent_unordered_map.h"
template<typename T>
using HashMap = tbb::concurrent_unordered_map<size_t, std::vector<T>>;
#endif

//auxiliary functions
mat_ZZ augment(const mat_ZZ& A, const vec_ZZ& b);
mat_ZZ stack(const mat_ZZ& A, const vec_ZZ& b);

template<typename T, typename S>
bool is_binary(const T& x) {
	for (auto& val : x)
		if (val != S(0) && val != S(1))
			return false;
	return true;
}

/*
 * A in top rows, B below
 * */
template<typename T, typename S>
T stack(const T& A, const vec_ZZ& b) {
	T res;
	res.SetDims(A.NumRows() + 1, A.NumCols());
	for (int i = 0; i < A.NumRows(); ++i) {
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
	}

	for (int j = 0; j < A.NumCols(); ++j) {
		res[res.NumRows() - 1][j] = conv < S > (b[j]);
	}
	return res;
}

/*
 * A in top rows, B below
 * */
template<typename T>
T stack(const T& A, const T& B) {
	T res;
	res.SetDims(A.NumRows() + B.NumRows(), A.NumCols());
	int i = 0;
	for (; i < A.NumRows(); ++i)
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
	for (; i < res.NumRows(); ++i)
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = B[i - A.NumRows()][j];
		}
	return res;
}

/**
 *  A in ���rst columns, matrix B to the right
 */
template<typename T>
T augment(const T& A, const T& B) {
	T res;
	res.SetDims(A.NumRows(), A.NumCols() + B.NumCols());
	for (int i = 0; i < res.NumRows(); i++) {
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
		for (int j = 0; j < B.NumCols(); ++j) {
			res[i][j + A.NumCols()] = B[i][j];
		}
	}
	return res;
}

template<typename T, typename S>
T augment(const T& A, const mat_ZZ& B) {
	T res;
	res.SetDims(A.NumRows(), A.NumCols() + B.NumCols());
	for (int i = 0; i < res.NumRows(); i++) {
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
		for (int j = 0; j < B.NumCols(); ++j) {
			res[i][j + A.NumCols()] = conv < S > (B[i][j]);
		}
	}
	return res;
}

RR MyInnerProduct(const vec_RR& a, const vec_RR& b);

// https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
#include <algorithm>
#include <functional>

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
	std::vector<T> result;
	result.reserve(a.size());

	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result),
			std::plus<T>());
	return result;
}

//https://nathanbrixius.wordpress.com/2013/09/13/generating-random-unimodular-matrices-with-a-column-of-ones/
mat_RR random_unimodular_matrix(int d);

// only permutate the first m entries (rows)
/*
 * return true: if this permutation has 2c binary ones in the last r entries, ie attack will succeed
 * */
void permutate(mat_ZZ& A, vec_ZZ& b, vec_ZZ& error,
		std::vector<int>& permutation, std::mt19937& g);

void repermutate(vec_RR& result, const std::vector<int>& permutation);

// for evaluating
RR euclidean_norm(const vec_RR& vec);

// for monitoring cpu freq
void cpu_freq();

void output_PrecomputedBasis_ToFile(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs,
		const vec_RR& euclid_norm);

bool test_necessity_condition_for_attack(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs,
		const vec_RR& norm);

// output reduced matrices to file, one matrix per file
void output_individual_matrix(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs);

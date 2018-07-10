#include "Attack.h"
#include "Utils.h"
#include <cstring>
#include <NTL/LLL.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/tools.h>
#include <NTL/vec_GF2.h>

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

// for random
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "NP.h"
#include "Hash_indeces.h"
#include <omp.h>

// for monitoring cpu
#include "thread_guard.h"

// for output to text
#include <iostream>
#include <fstream>
#include <sstream>

// for permutation
#include <set>

#include <functional>

/*
 * read output file given by -output
 * */
void readBasis_fromFile(mat_RR& B, mat_RR& B_1_reduced_transposed,
		mat_RR& B_reduced_gs) {

	// read ZZ matrix then convert to mat_RR
	mat_ZZ B_(INIT_SIZE, m + 1, m + 1);
	mat_ZZ B_1_reduced_transposed_(INIT_SIZE, m - r, m - r);
	B_reduced_gs.SetDims(m - r, m - r);

	std::ifstream fileHandler;
	std::string output;
	std::string testName =
			input.substr(input.find_last_of("/") + 1,
					input.find_last_of(".") - input.find_last_of("/") - 1); // cut testx out of input string
#ifdef USING_MPI
							output = "reduced_basis/"+testName+"_r_" + std::to_string(r)+"_beta_"+std::to_string(beta)+".txt";
							std::cout << "MPI proc. " << world_rank << " reading basis file " << output
							<< std::endl;
#else
	output = "reduced_basis/" + testName + "_process0.txt";
	output = "reduced_basis/" + testName + "_beta_" + std::to_string(beta)
			+ ".txt";
	std::cout << "reading basis file " << output << std::endl;
#endif

	fileHandler.open(output);
	if (fileHandler.fail()) {
		std::cout << "No basis file\n";
	} else {

		std::string line;
		char toRead = '\0';

		int row = 0;
		std::map<std::string, char> tokens = { { "B", 'a' }, {
				"B_1_reduced_transposed", 'b' }, { "B_reduced_gs", 'c' }, {
				"permutation", 'd' } }; // convention: read case 'a', 'b', 'c'
		while (std::getline(fileHandler, line)) {
			if (toRead != '\0') {
				switch (toRead) {
				case 'a': { // read B
					// read a new row
					setRow(B_[row], line, m + 1); // row length = m+1
					++row;
					if (row == m + 1) {
						toRead = '\0';
						row = 0; // to read next matrix
					}
					break;
				}
				case 'b': { // read B_1_reduced_transposed
					// read a new row
					setRow(B_1_reduced_transposed_[row], line, m - r); // row length = m-r
					++row;
					if (row == m - r) {
						toRead = '\0';
						row = 0; // to read next matrix
					}
					break;
				}
				case 'c': { // read B_reduced_gs
					// read a new row
					setRow(B_reduced_gs[row], line, m - r); // row length = m-r
					++row;
					if (row == m - r) {
						toRead = '\0';
					}
					break;
				}
				case 'd': {
					// read a new row
					setRow(permutation, line, m); // row length = m
					toRead = '\0';
					break;
				}
				}
			} else {
				auto found = tokens.find(line);
				if (found != tokens.end()) {
					toRead = found->second;
				} else {
					toRead = '\0';
				}
			}
		}
	}

	B = conv < mat_RR > (B_);
	B_1_reduced_transposed = conv < mat_RR > (B_1_reduced_transposed_);
}

/**
 * returns a binary vector of length r with exactly c bits 1.
 */
vec_RR sample_bin_vector_RR(int r, int c) {
	vec_RR res(INIT_SIZE, r);
	int i = 0;
	while (i < c) {
		int pos = rand() % r;
		if (res[pos] == RR(0)) {
			res[pos] = RR(1);
			++i;
		}
	}
	return res;
}

std::vector<float> sample_bin_vector_float(int r, int c, std::mt19937& g,
		std::uniform_int_distribution<>& dis) {
	std::vector<float> res(r, 0);
	int i = 0;
	while (i < c) {
		int pos = dis(g);
		if (res[pos] == 0) {
			res[pos] = 1;
			++i;
		}
	}
	return res;
}

vec_ZZ sample_bin_vector_ZZ(int r, int c) {
	vec_ZZ res(INIT_SIZE, r);
	int i = 0;
	while (i < c) {
		int pos = rand() % r;
		if (res[pos] == ZZ(0)) {
			res[pos] = ZZ(1);
			++i;
		}
	}
	return res;
}

/**
 * returns a binary vector of length r with bit 1 at the odd positions.
 */
vec_ZZ bin_vector_with_odd_pos_1(int r) {
	vec_ZZ res(INIT_SIZE, r, ZZ(0));
	for (int i = 0; i < r; ++i) {
		if (i % 2) {
			res[i] = ZZ(1);
		}
	}
	return res;
}

double euclidian_norm(const vec_ZZ& vec) {
	int res;
	for (int i = 0; i < vec.length(); ++i) {
		//TODO: precision?
		res += conv<int>(vec[i] * vec[i]);
	}
	return sqrt(res);
}

void precomputing(mat_ZZ& A, vec_ZZ& b, int q, int r, long beta, mat_RR& B,
		mat_RR& B_1_reduced_transposed, mat_RR& B_reduced_gs, std::mt19937& g) {
	double total_precomputing_time = 0;
	vec_RR euclid_norm(INIT_SIZE, m - r);

#ifdef USING_MPI
	int buffer[world_size];
	int recvBuffer;

	MPI_Request recvRequest;
	MPI_Request sendReq[world_size];
	bool setRecv = false;
	int flag;
#endif

	do {
		double startTimer = omp_get_wtime();
#ifdef MODE_2
		permutate(A, b, error, permutation, g);

		std::cout << "Permutation: \n";
		for (const auto& it : permutation) {
			std::cout << it << " ";
		}
		std::cout << std::endl;
#endif
		//	A = A.stack(matrix(1,n))  						# NEW!!!!
		vec_ZZ zero_vec_length_n(INIT_SIZE, A.NumCols());
		mat_ZZ local_A = stack(A, zero_vec_length_n);

		//	b=matrix(b).transpose()   						# NEW!!!!
		vec_ZZ local_b = b;

		//		b = b.stack(matrix([1]))  						# NEW!!!!
		local_b.append(ZZ(1));
		// m_plus1 = m+1
		int m_plus1 = m + 1;

		// C=A.augment(b)
		mat_ZZ C = augment(local_A, local_b);

		// C=C.change_ring(Zmod(q))
		ZZ_p::init(ZZ(q));
		mat_ZZ_p C_q(INIT_SIZE, C.NumRows(), C.NumCols());

		for (int i = 0; i < C.NumRows(); i++)
			for (int j = 0; j < C.NumCols(); ++j)
				C_q[i][j] = C[i][j] % q;

		// A_2=C.submatrix(m_plus1-n-1,0) -> A_2 takes the last (n+1) rows of C
		mat_ZZ_p A_2 = submatrix<mat_ZZ_p>(C_q, m_plus1 - n - 1, 0, n + 1,
				n + 1);

		mat_ZZ AA;
		try {
			// mat_ZZ_p AA_q = C_q * inv(A_2); // in ZZ_q
			AA = conv < mat_ZZ > (C_q * inv(A_2));
		} catch (const std::exception& e) {
			std::cout << e.what() << std::endl;
			continue;
		}
		// AA.change_ring(ZZ)
		B = q * ident_mat_RR(m_plus1 - n - 1);

		// B=B.stack(matrix(ZZ,n+1,m-n-1))
		mat_RR zero1(INIT_SIZE, n + 1, m_plus1 - n - 1);
		B = stack<mat_RR>(B, zero1);
		B = augment<mat_RR, RR>(B, AA);

// B_1
		mat_RR B_1(INIT_SIZE, m - r, m - r);
		B_1 = submatrix<mat_RR>(B, 0, 0, m - r, m - r); // column vector

#ifdef MODE_2
// B_1_reduced
// replace B_1 with B_1,i unimodular
		B_1 = B_1 * random_unimodular_matrix(m - r);
		std::cout << "Unimodular matrix generated!\n";
#endif
		mat_ZZ B_1_t = transpose(conv < mat_ZZ > (B_1));	// row vector
		double delta = 0.99;
		std::cout << "r = " << r << " c = " << c << " beta = " << beta
				<< std::endl;
		std::cout << "prune = " << prune << std::endl;

#ifdef MODE_2
		if (!setRecv) {
			// waiting for complete signal from others
			MPI_Irecv(&recvBuffer, 1, MPI_INT, MPI_ANY_SOURCE, PRECOMPUTING_FINISH_TAG,
					MPI_COMM_WORLD, &recvRequest);// only one receive pro mpi process
			setRecv = true;
		} else {
			MPI_Test(&recvRequest, &flag, MPI_STATUS_IGNORE);
			if (flag) {
				std::cout << "MPI proc " << world_rank << " terminated by proc " << recvBuffer << std::endl;
				break;
			}
		}
#endif
		double startTimerBKZ = omp_get_wtime();
		BKZ_FP(B_1_t, delta, beta, prune);

#ifdef USING_MPI
		++precomputing_reduction_iterCount;
		std::cout << " proc " << world_rank << " time BKZ: "
		<< omp_get_wtime() - startTimerBKZ << " s\n";
#else
		std::cout << "time BKZ: " << omp_get_wtime() - startTimerBKZ << " s\n";
#endif
		mat_ZZ B_1_reduced = transpose(B_1_t);	// column vector

// replace the upper left part of B by B_1_reduced
		for (int i = 0; i < m - r; ++i) {
			for (int j = 0; j < m - r; ++j) {
				B[i][j] = conv < RR > (B_1_reduced[i][j]);
			}
		}

// 2. output matrix B_reduced_gs
// B_1_reduced_transposed
		B_1_reduced_transposed = conv < mat_RR > (transpose(B_1_reduced));// row vector

// B_reduced_gs
		double startTimerGS = omp_get_wtime();
		B_reduced_gs = GSO(B_1_reduced_transposed);	// row vector
		std::cout << "time GS: " << omp_get_wtime() - startTimerGS << " s\n";

		// for evaluation
		// compute euclidean norms of row vectors in B_reduced_gs

		for (int i = 0; i < m - r; ++i) {
			euclid_norm[i] = euclidean_norm(B_reduced_gs[i]);
		}

		total_precomputing_time += omp_get_wtime() - startTimer;
		bool cont = test_necessity_condition_for_attack(B,
				B_1_reduced_transposed, B_reduced_gs, euclid_norm);

#ifdef MODE_2
		if (cont) {
			// signaling others
			// signaling all processes
			for (int i = 0; i < world_size; ++i) {
				if (i!=world_rank) {
					buffer[i] = world_rank;
					MPI_Isend(&buffer[i], 1, MPI_INT, i, PRECOMPUTING_FINISH_TAG,
							MPI_COMM_WORLD, &sendReq[i]);
				}
			}
			for (int i = 0; i < world_size; ++i) {
				if (i!=world_rank) {
					MPI_Wait(&sendReq[i], MPI_STATUS_IGNORE);
				}
			}
		}
//		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if (cont) {
			break;
		}
	} while (true);

	//	output_PrecomputedBasis_ToFile(B, B_1_reduced_transposed, B_reduced_gs,euclid_norm);
//	return;

	std::cout << "Total precomputing time: " << total_precomputing_time
			<< std::endl;

}

vec_RR hybridAttack(const mat_RR& B, const mat_RR& B_1_reduced_transposed,
		const mat_RR& B_reduced_gs, int r, int c, std::mt19937& g) {

// hash_map
	double startTimer = omp_get_wtime();
	HashMap<std::vector<float>> hash_map;

	std::cout << "Running " << numThread << " attack processes\n";

// each thread has its own result vector
	vec_RR res_[numThread];
	for (int i = 0; i < numThread; ++i) {
		res_[i].SetLength(m);
	}
	bool breakFlag = false;
	int resIdx = 0;

#ifdef USING_MPI
	int buffer[world_size];
	int recvBuffer;

	MPI_Request recvRequest;
	MPI_Request sendReq[world_size];
	int flag;
	bool setRecv = false;
	size_t iteration_count = 0;
	const int check_period = 20;
#endif

// for uniform distribution in [0,r-1]
	std::uniform_int_distribution<> dis(0, r - 1);

// for evaluation
	size_t guess_vector_counter = 0;

#pragma omp parallel reduction(+:guess_vector_counter) shared(hash_map, breakFlag, res_, resIdx) num_threads(numThread)
	{
		vec_RR w(INIT_SIZE, m - r);
		vec_RR x(INIT_SIZE, m - r);

		std::hash < std::bitset < DIM >> hash_fn;

		while (!breakFlag) {
#ifdef USING_MPI
#pragma omp master
			{
				++iteration_count;
				if (!setRecv) {
					// waiting for complete signal from others
					MPI_Irecv(&recvBuffer, 1, MPI_INT, MPI_ANY_SOURCE, ATTACK_FINISH_TAG,
							MPI_COMM_WORLD, &recvRequest);// only one receive pro mpi process
					setRecv = true;
				} else {
					if (iteration_count % check_period == 0) {
						MPI_Test(&recvRequest, &flag, MPI_STATUS_IGNORE);
						if (flag) {
							std::cout << "MPI proc " << world_rank << " terminated by proc " << recvBuffer << std::endl;
#pragma omp atomic write
							breakFlag = true;
						}
					}
				}
			} // end of master-construct
#endif

//			vec_RR v_1 = sample_bin_vector_RR(r, c);
			std::vector<float> v_1 = sample_bin_vector_float(r, c, g, dis);
			++guess_vector_counter;

			//v_1_padded
			vec_RR v_1_padded(INIT_SIZE, m - r);
//			v_1_padded.append(v_1);
			for (auto& i : v_1) {
				v_1_padded.append(RR(i));
			}
			v_1_padded.append(RR(0.5));

			// w
			vec_RR w_ = B * v_1_padded - v_1_padded;
			w_.SetLength(m - r);
			w = w_;

// x
			x = w - NP(B_1_reduced_transposed, w, B_reduced_gs);

// indeces
			std::vector < std::bitset<DIM> > indeces = get_hash_indeces<vec_RR,
					RR>(x);

			for (auto& ind : indeces) {
				if (breakFlag) {
					break;
				}

				auto it = hash_map.find(hash_fn(ind));
				if (it == hash_map.end()) {
					// hash_map[ind].append(v_1)

					hash_map.insert(
							HashMap<std::vector<float>>::value_type(
									hash_fn(ind),
									std::vector<std::vector<float>>(1, v_1)));
				} else {
					for (auto& v_2 : it->second) {
						if (breakFlag) {
							break;
						}
//						if (v_2 != v_1) {
//							std::cout << "\tCollision: " << v_2 << std::endl;
//						}
						std::vector<float> sum_v12 = v_1 + v_2;
						if (!is_binary<std::vector<float>, float>(sum_v12)) {
//							std::cout << ".";
							continue;
						}
						// v_padded
						vec_RR v_padded(INIT_SIZE, m - r);
//						v_padded.append(sum_v12);
						for (auto& i : sum_v12) {
							v_padded.append(RR(i));
						}
						v_padded.append(RR(1));

						// w
						w_ = B * v_padded - v_padded;
						w_.SetLength(m - r);
						w = w_;

						// x
						x = w - NP(B_1_reduced_transposed, w, B_reduced_gs);

						if (is_binary<vec_RR, RR>(x)) {
							vec_RR res = x;
//							res.append(sum_v12);
							for (auto& i : sum_v12) {
								res.append(RR(i));
							}

#pragma omp atomic write
							breakFlag = true;
#pragma omp critical
							{
#ifdef USING_MPI
								std::cout << "\n MPI-proc " << world_rank
								<< " time attack: "
								<< omp_get_wtime() - startTimer
								<< " s\n";

								std::cout << " MPI-proc " << world_rank
								<< " thread " << omp_get_thread_num()
								<< " finished." << std::endl;
#else
								std::cout << "\ntime attack: "
										<< omp_get_wtime() - startTimer
										<< " s\n";

								std::cout << "\nthread " << omp_get_thread_num()
										<< " finished." << std::endl;
#endif
								resIdx = omp_get_thread_num();
								res_[resIdx] = res;
#ifdef USING_MPI
								// signaling all processes
								for (int i = 0; i < world_size; ++i) {
									buffer[i] = world_rank;
									MPI_Isend(&buffer[i], 1, MPI_INT, i, ATTACK_FINISH_TAG,
											MPI_COMM_WORLD, &sendReq[i]);
								}
								for (int i = 0; i < world_size; ++i) {
									MPI_Wait(&sendReq[i], MPI_STATUS_IGNORE);
								}
#endif
							}
						}
					}
					// hash_map[ind].append(v_1)
					it->second.push_back(v_1);
				} // end else
			}	// end for-loop

		}	// end while-loop
	} // end omp-section
	std::cout << "number of guessed vectors: " << guess_vector_counter
			<< std::endl;
	return res_[resIdx];
}

void attack(mat_ZZ& A, vec_ZZ& b, int q, int r, int c) {
	mat_RR B;
	mat_RR B_1_reduced_transposed;
	mat_RR B_reduced_gs;

// random generator for permutation
	std::random_device rd;
	std::mt19937 g(rd());

#ifdef USING_MPI
#ifdef MODE_1
	precomputing(A, b, q, r, beta, B, B_1_reduced_transposed, B_reduced_gs, g);
#endif
#ifdef MODE_2
	precomputing(A, b, q, r, beta, B, B_1_reduced_transposed, B_reduced_gs, g);
	return;
#endif
#ifdef MODE_3
	readBasis_fromFile(B, B_1_reduced_transposed, B_reduced_gs);
#endif
#ifdef MODE_4
	if (world_rank % 2 == 0) {
		precomputing(A, b, q, r, beta, B, B_1_reduced_transposed, B_reduced_gs, g);
	} else {
		readBasis_fromFile(B, B_1_reduced_transposed, B_reduced_gs);
	}
#endif
#else
	precomputing(A, b, q, r, beta, B, B_1_reduced_transposed, B_reduced_gs, g);
#endif

	// init the random seed
	vec_RR res = hybridAttack(B, B_1_reduced_transposed, B_reduced_gs, r, c, g);

#ifdef MODE_2
	repermutate(res, permutation);
#endif

#ifdef USING_MPI
	std::cout << "\nMPI-proc " << world_rank << " returns " << res << std::endl;
#else
	std::cout << "\nfound " << res << std::endl;
#endif
}


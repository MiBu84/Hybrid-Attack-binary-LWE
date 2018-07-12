#include <string>
#include <iostream>
#include <fstream>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>
#include "Attack.h"
#include "Utils.h"
#include <omp.h>

using namespace NTL;

int m;
int n;
int q;
int r;
int c;
long beta; // blocksize
long prune = 10;
mat_ZZ A;
vec_ZZ b;
vec_ZZ error;
int testcase;
std::string input; // A,b,error
std::string output; // output reduced basis
int timeout; // timeout for attack, given in seconds

// for evaluation
double NP_call_duration = 0;

std::vector<int> permutation;
int numThread;

#include "mpi.h"
int world_rank;
int world_size;
int precomputing_reduction_iterCount = 0;

int main(int argc, char** argv) {
	if (argc >= 2) {
		for (int i = 1; i < argc; i++) {
			std::string s_par = argv[i];
			// m
			if (s_par == "-m") {
				m = atoi(argv[i + 1]);
				i++;
			}

			// r
			if (s_par == "-r") {
				r = atoi(argv[i + 1]);
				i++;
			}

			// c
			if (s_par == "-c") {
				c = atoi(argv[i + 1]);
				i++;
			}

			// input
			if (s_par == "-input") {
				input = argv[i + 1];
				i++;
			}

			// output
			if (s_par == "-output") {
				output = argv[i + 1];
				i++;
			}

			// numthread
			if (s_par == "-numthread") {
				numThread = atoi(argv[i + 1]);
				i++;
			}

			// blocksize
			if (s_par == "-beta") {
				beta = atol(argv[i + 1]);
				i++;
			}

			// timeout
			if (s_par == "-timeout") {
				timeout = atol(argv[i + 1]);
				i++;
			}

		}
	}
	if (MPI_Init(NULL, NULL) != MPI_SUCCESS) {
		std::cout << "MPI init failed\n";
		return -1;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// redirect output of MPI proc to text file
	std::string output_redirect = "/output_proc_" + std::to_string(world_rank)
			+ ".txt";
	std::ofstream out("MPI-output/" + output + output_redirect);
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	read(input, A, b);
	// id permutation
	for (int i = 0; i < m; ++i) {
		permutation.push_back(i);
	}

	double start = omp_get_wtime();
	attack(A, b, q, r, c);
	std::cout << "\nMPI-proc " << world_rank << " total runtime: "
			<< omp_get_wtime() - start << " s.\n";

	std::cout.rdbuf(coutbuf); //reset to standard output again

	MPI_Finalize();
	return 0;
}

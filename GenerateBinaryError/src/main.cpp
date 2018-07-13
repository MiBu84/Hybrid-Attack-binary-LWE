//============================================================================
// Name        : GenerateInputHybridAttack.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

// for output to text
#include <iostream>
#include <fstream>
#include <sstream>

using namespace NTL;

// for random
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

int m;
int n;
int q;
int r = 0;
int c = 0;
long beta = 0; // blocksize
long prune = 10;
mat_ZZ A;
vec_ZZ b;
vec_ZZ error;
std::string output;

void get_instance(int n, int m, int q) {
	error = vec_ZZ(INIT_SIZE, m);
	//# for test files --> e = (0,1,0,1,0,1, ...)
	for (int i = 0; i < m; ++i) {
		error[i] = i % 2;
	}

	A.SetDims(m, n);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			A[i][j] = rand() % q;

	vec_ZZ s(INIT_SIZE, n);
	for (int i = 0; i < n; ++i) {
		s[i] = rand() % q;
	}

	b = A * s + error;
	for (int i = 0; i < m; ++i) {
		b[i] = b[i] % q;
		if (b[i] < 0) {
			b[i] += q;
		}
	}

	// print to file
	std::ofstream myfile;
	std::cout << "output: " << output << std::endl;

	myfile.open(output);
	//m
	myfile << "m\n";
	myfile << m << std::endl;

	//n
	myfile << "n\n";
	myfile << n << std::endl;

	//q
	myfile << "q\n";
	myfile << q << std::endl;

	//r
	if (r != 0) {
		myfile << "r\n";
		myfile << r << std::endl;
	}
	//c
	if (c != 0) {
		myfile << "c\n";
		myfile << c << std::endl;
	}
	//beta
	if (beta != 0) {
		myfile << "beta\n";
		myfile << beta << std::endl;
	}
	// dim
	if (r != 0) {
		myfile << "dim\n";
		myfile << m - r << std::endl;
	}
	// A
	myfile << "A\n";
	myfile << A << std::endl;
	// b
	myfile << "b\n";
	myfile << b << std::endl;

	// error
	myfile << "error\n";
	myfile << error << std::endl;

	myfile.close();
	std::cout << "output finished...\n";
}

int main(int argc, char** argv) {
//if (argc != 8) {
	for (int i = 1; i < argc; i++) {
		std::string s_par = argv[i];

		// m
		if (s_par == "-m") {
			m = atoi(argv[i + 1]);
			i++;
		}

		// n
		if (s_par == "-n") {
			n = atoi(argv[i + 1]);
			i++;
		}

		// q
		if (s_par == "-q") {
			q = atoi(argv[i + 1]);
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

		// output
		if (s_par == "-output") {
			output = argv[i + 1];
			i++;
		}

		// blocksize
		if (s_par == "-beta") {
			beta = atol(argv[i + 1]);
			//std::cout << "beta = " << beta << std::endl;
			i++;
		}
	}
	srand(time(NULL));
	get_instance(n, m, q);
	return EXIT_SUCCESS;
//} else {
//	std::cout << "Arguments missing\n";
//	return -1;
//}
}

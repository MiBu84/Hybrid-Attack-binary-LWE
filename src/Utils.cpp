#include "Utils.h"
#include "Attack.h"
#include <NTL/LLL.h>

// handling input file
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <algorithm>
#include <map>
#include "NP.h"
// for randomize unimodular matrix
#include <NTL/ZZ.h>

mat_RR GSO(const mat_ZZ& input) {
	mat_RR mu;
	vec_RR c;
	ComputeGS(input, mu, c);
	// set the diagonal
	for (int i = 0; i < input.NumRows(); ++i)
		mu[i][i] = RR(1);

	mat_RR res = inv(mu) * conv < mat_RR > (input);
	return res;
}

mat_RR GSO(const mat_RR& in) {
	mat_RR mu;
	vec_RR c;

	mat_ZZ input(INIT_SIZE, in.NumRows(), in.NumCols());
	for (int i = 0; i < input.NumRows(); ++i) {
		for (int j = 0; j < input.NumCols(); ++j) {
			input[i][j] = conv < ZZ > (in[i][j]);
		}
	}

	ComputeGS(input, mu, c);
	// set the diagonal
	for (int i = 0; i < input.NumRows(); ++i)
		mu[i][i] = RR(1);

	mat_RR res = inv(mu) * in;
	return res;
}

void setRow(vec_ZZ& a, std::string line, int len) {
	// remove [ ]
	size_t pos;
	do {
		pos = line.find('[');
		if (pos != std::string::npos) {
			line.erase(pos, 1);
		} else {
			break;
		}
	} while (true);

	do {
		pos = line.find(']');
		if (pos != std::string::npos) {
			line.erase(pos, 1);
		} else {
			break;
		}
	} while (true);

	std::vector < std::string > strs;
	if (line.find(",") != std::string::npos) {
		boost::split(strs, line, boost::is_any_of(","));
	} else {
		boost::split(strs, line, boost::is_any_of(" "));
	}

	if (strs.size() != len) { // remove EOL
		strs.pop_back();
	}

	// remove blank spaces
	for (auto& str : strs) {
		str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
	}
	for (int i = 0; i < len; ++i) {
		a[i] = ZZ(std::stoi(strs[i]));
	}
}

void setRow(std::vector<int>& permutation, std::string line, int len) {
	std::vector < std::string > strs;
	boost::split(strs, line, boost::is_any_of(" "));
	for (int i = 0; i < m; ++i) {
		permutation.push_back(std::stoi(strs[i]));
	}
}

void setRow(vec_RR& a, std::string line, int len) {
	// remove [ ]
	size_t pos;
	do {
		pos = line.find('[');
		if (pos != std::string::npos) {
			line.erase(pos, 1);
		} else {
			break;
		}
	} while (true);

	do {
		pos = line.find(']');
		if (pos != std::string::npos) {
			line.erase(pos, 1);
		} else {
			break;
		}
	} while (true);

	std::vector < std::string > strs;
	if (line.find(",") != std::string::npos) {
		boost::split(strs, line, boost::is_any_of(","));
	} else {
		boost::split(strs, line, boost::is_any_of(" "));
	}

	if (strs.size() != len) {
		strs.pop_back();
	}

	// remove blank spaces
	for (auto& str : strs) {
		str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
	}

	for (int i = 0; i < len; ++i) {
		size_t found = strs[i].find('/'); // sage output

		if (found != std::string::npos) { // reading fraction
			// nominator
			RR nom = RR(std::stod(strs[i].substr(0, found)));
			// denominator
			RR denom = RR(
					std::stod(
							strs[i].substr(found + 1,
									strs[i].size() - 1 - found)));
			a[i] = nom / denom;
		} else {
			a[i] = RR(std::stod(strs[i]));
		}
	}
}

void setDim() {
	A.SetDims(m, n);
	b.SetLength(m);
	error.SetLength(m);
}

void read(std::string file, mat_ZZ& A, vec_ZZ& b) {
	std::ifstream input(file);
	std::string line;
	char toRead = '\0';

	int d; // dimension in BKZ

	int row = 0;
	std::map<std::string, char> tokens = { { "m", 'm' }, { "n", 'n' }, { "q",
			'q' }, { "r", 'r' }, { "c", 'c' }, { "beta", 's' }, { "A", 'A' }, {
			"b", 'b' }, { "error", 'e' }, { "dim", 'd' } };
	while (std::getline(input, line)) {
		if (toRead != '\0') {
			switch (toRead) {
			case 'e':
				setRow(error, line, m);
				toRead = '\0';
				break;
			case 'm':
				m = std::stoi(line);
				toRead = '\0';
				break;
			case 'n':
				n = std::stoi(line);
				toRead = '\0';
				setDim();
				break;
			case 'q':
				q = std::stoi(line);
				toRead = '\0';
				break;
			case 'r':
				r = std::stoi(line);
				toRead = '\0';
				break;
			case 'c':
				c = std::stoi(line);
				toRead = '\0';
				break;
//			case 's':
//				beta = std::stol(line);
//				toRead = '\0';
//				break;
			case 'A': {
				// read a new row
				setRow(A[row], line, n);
				++row;
				if (row == m) {
					toRead = '\0';
					//	std::cout << "A = \n" << A << std::endl;
				}
				break;
			}
			case 'b': {
				setRow(b, line, m);
				toRead = '\0';
				//	std::cout << "b = \n" << b << std::endl;
				break;
			}
			case 'd':
				d = std::stoi(line);
				toRead = '\0';
				std::cout << "dimension in BKZ: " << d << std::endl;
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
	std::cout << "read input finished\n";
}

void read(mat_RR& mat, std::string file) {
	int i = 0;
	std::ifstream input(file);
	std::string line;
	vec_RR row(INIT_SIZE, mat.NumCols());
	while (std::getline(input, line)) {
		setRow(row, line, mat.NumCols());
		mat[i] = row;
		i++;
	}
}

void read_v1(std::vector<vec_RR>& vec, std::string file) {
	std::ifstream input(file);
	std::string line;
	vec_RR row(INIT_SIZE, r); // v1 has length r
	while (std::getline(input, line)) {
		setRow(row, line, r);
		vec.push_back(row);
	}
}

//TODO: optimize it
mat_ZZ augment(const mat_ZZ& A, const vec_ZZ& b) {
	mat_ZZ res;
	res.SetDims(A.NumRows(), A.NumCols() + 1);
	for (int i = 0; i < res.NumRows(); i++) {
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
		res[i][res.NumCols() - 1] = b[i];
	}
	return res;
}

//TODO: optimize it
mat_ZZ stack(const mat_ZZ& A, const vec_ZZ& b) {
	mat_ZZ res;
	res.SetDims(A.NumRows() + 1, A.NumCols());
	for (int i = 0; i < A.NumRows(); ++i) {
		for (int j = 0; j < A.NumCols(); ++j) {
			res[i][j] = A[i][j];
		}
	}

	for (int j = 0; j < A.NumCols(); ++j) {
		res[res.NumRows() - 1][j] = b[j];
	}
	return res;
}

RR MyInnerProduct(const vec_RR& a, const vec_RR& b) {
	RR res; // = 0
	for (size_t i = 0; i < a.length(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}

mat_RR random_unimodular_matrix(int d) {
	SetSeed (ZZ(time (NULL)));ZZ bound(q);

	mat_ZZ upper(INIT_SIZE, d, d);
	for (int i = 0; i < d; ++i) {
		upper[i][i] = ZZ(1);
		for (int j = i+1; j < d; ++j) upper[i][j] = RandomBnd(bound);
	}

	mat_ZZ lower(INIT_SIZE, d, d);
	for (int i = 0; i < d; ++i) {
		lower[i][i] = ZZ(1);
		for (int j = 0; j < i; ++j) lower[i][j] = RandomBnd(bound);
	}
	return conv<mat_RR>(upper * lower); // convert to mat_RR
}

void permutate(mat_ZZ& A, vec_ZZ& b, vec_ZZ& error,
		std::vector<int>& permutation, std::mt19937& g) {
	std::shuffle(permutation.begin(), permutation.end(), g);

	// do permutation on A,b
	for (int i = 0; i < m; ++i) {
		swap(error[i], error[permutation[i]]);
		swap(A[i], A[permutation[i]]);
		swap(b[i], b[permutation[i]]);
	}
}

void repermutate(vec_RR& result, const std::vector<int>& permutation) {
	for (int i = 0; i < result.length(); ++i) {
		swap(result[i], result[permutation[i]]);
	}
}

RR euclidean_norm(const vec_RR& vec) {
	RR sum = RR(0);
	for (auto& it : vec) {
		sum += it * it;
	}
	return sqrt(sum);
}

void cpu_freq() {
	do {
		system("~/Projects/HybridAttack/cpuinfo.sh");
		std::this_thread::sleep_for(std::chrono::seconds(15));
	} while (true);
}

void output_PrecomputedBasis_ToFile(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs,
		const vec_RR& euclid_norm) {
	std::ofstream myfile;
	std::string testName =
			input.substr(input.find_last_of("/") + 1,
					input.find_last_of(".") - input.find_last_of("/") - 1); // cut testx out of input string
	output = "reduced_basis/"+testName+"_r_" + std::to_string(r)+"_beta_"+std::to_string(beta)+".txt";

	myfile.open(output);
	// B
	myfile << "B\n";
	myfile << B;
	// B_1_reduced_t
	myfile << "\nB_1_reduced_transposed\n";
	myfile << B_1_reduced_transposed;
	// B_reduced_gs
	myfile << "\nB_reduced_gs\n";
	myfile << B_reduced_gs;
	myfile.close();
	std::cout << "output finished...\n";
}

bool test_necessity_condition_for_attack(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs,
		const vec_RR& norm) {
	std::cout << "MPI-Proc " << world_rank << " checks necessity condition: ";
	// count number of ones
	int count = 0;
	for (int i = m - r; i < m; ++i) {
		if (error[i] == ZZ(1)) {
			++count;
		}
	}

	if (count != 2 * c) {
		std::cout << "Number of ones not correct\n";
		std::cout << "---------------------------------------------\n";
		return false;
	}

	vec_RR v_g(INIT_SIZE, 0);
	v_g.append(conv < vec_RR > (error));
	v_g.append(RR(1));

	for (size_t i = 0; i < m - r; ++i) {
		v_g[i] = RR(0);
	}

// w_test
	vec_RR w_test = B * v_g - v_g;
	w_test.SetLength(m - r);

// x_test
	vec_RR x_test = w_test - NP(B_1_reduced_transposed, w_test, B_reduced_gs);

	if (is_binary<vec_RR, RR>(x_test)) {
		std::cout << "brute force guessing would work\n";
		std::cout << "---------------------------------------------\n";
//		output_PrecomputedBasis_ToFile(B, B_1_reduced_transposed, B_reduced_gs,
//				norm);
		return true;
	} else {
		std::cout << "attack can not work\n";
		std::cout << "---------------------------------------------\n";
		return false;
	}
}

void output_individual_matrix(const mat_RR& B,
		const mat_RR& B_1_reduced_transposed, const mat_RR& B_reduced_gs) {
	std::string inputNumber = input.substr(input.find('/') + 1,
			input.find('.') - input.find('/') - 1);
	// A

	// b

	// error

	std::ofstream myfile;
	std::string fileName;
	//mat_RR B;
	fileName = "B_" + inputNumber + "_r_" + std::to_string(r) + "_beta_" + std::to_string(beta) + ".txt";
	myfile.open(fileName);
	myfile << B;
	myfile.close();
	//mat_RR B_1_reduced_transposed;
	fileName = "B_1_reduced_transposed_" + inputNumber + "_r_" + std::to_string(r) + "_beta_"
			+ std::to_string(beta) + ".txt";
	myfile.open(fileName);
	myfile << B_1_reduced_transposed;
	myfile.close();
	//mat_RR B_reduced_gs;
	fileName = "B_reduced_gs_" + inputNumber + "_r_" + std::to_string(r) + "_beta_" + std::to_string(beta)
			+ ".txt";
	myfile.open(fileName);
	myfile << B_reduced_gs;
	myfile.close();
}

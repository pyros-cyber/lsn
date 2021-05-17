#include "myStatFunc.h"
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  // Declaring and initializing "global" variables
  Random rnd("Primes", "seed.in");
  unsigned int M = 1000000, N = 100;
  unsigned int L = M / N;
  vector<double> distribution(M);
  for (auto &elem : distribution) {
    elem = rnd.Rannyu();
  }

  /******************************************/
  /* Exercise 01.1
  estimate <r> and its uncertainty
  */
  vector<double> ave(N);
  // Calculating the average of the N blocks
  auto distr_begin = distribution.begin();
  for (int i{}; i < N; i++) {
    ave[i] = accumulate(distr_begin + i * L, distr_begin + (i + 1) * L, 0.) /
             static_cast<double>(L);
  }
  // Calculating the average and the error as a function of the number of blocks
  vector<double> error = blocking_error(ave);
  // Printing the results
  ofstream out("results01.1.dat");
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      out << (i + 1) << " " << ave[i] << " " << error[i] << endl;
    }
  } else {
    cerr << "ERROR: error while opening output file." << endl;
    return 2;
  }
  out.close();

  return 0;
}

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
  int L = M / N;
  vector<double> distribution(M);
  for (auto &elem : distribution) {
    elem = rnd.Rannyu();
  }

  /******************************************/
  /* Exercise 01.1.1
   * estimate <r> and its uncertainty
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
  ofstream out("results01.1.1.dat");
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      out << (i + 1) << " " << ave[i] << " " << error[i] << endl;
    }
  } else {
    cerr << "ERROR: error while opening output file." << endl;
    return 2;
  }
  out.close();

  /******************************************/
  /* Exercise 01.1.2
   * estimate <sigma^2> and its uncertainty
   */
  fill(ave.begin(), ave.end(), 0);
  fill(error.begin(), error.end(), 0);
  // Same as above, but this time we accumelate sigma^2 = (<r>-0.5)^2
  for (int i{}; i < N; i++) {
    ave[i] =
        accumulate(distr_begin + i * L, distr_begin + (i + 1) * L, 0.,
                   [](auto lhs, auto rhs) { return lhs + pow(rhs - 0.5, 2); }) /
        static_cast<double>(L);
  }
  error = blocking_error(ave);
  // Printing the results
  out.open("results01.1.2.dat");
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      out << (i + 1) << " " << ave[i] << " " << error[i] << endl;
    }
  } else {
    cerr << "ERROR: error while opening output file." << endl;
    return 2;
  }
  out.close();

  /******************************************/
  /* Exercise 01.1.3
   * estimate the chiquadro^2
   */

  // Assign new values to M,N and thus L
  // and re-initializing the random distribution
  N = 100;
  M = 10000;
  L = M / N;

  distribution.resize(M * N);
  for (auto &elem : distribution) {
    elem = rnd.Rannyu();
  }

  // Opening the output stream
  out.open("results01.1.3.dat");

  vector<int> events(N);
  double chi_sq = 0.;

  // Computing the chi^2 100 times
  for (int i{}; i < N; i++) {
    fill(events.begin(), events.end(), 0);
    for (int j{}; j < M; j++) {
      events[static_cast<int>(distribution[i * M + j] * N)]++;
    }
    for (int j{}; j < N; j++) {
      chi_sq += static_cast<double>(pow(events[j] - L, 2) / L);
    }
    out << chi_sq << endl;
    chi_sq = 0.;
  }
  return 0;
}

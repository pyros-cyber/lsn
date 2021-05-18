#include "myStatFunc.h"
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Integrand function to evaluate
inline double integrand(double x) { return M_PI_2 * cos(M_PI_2 * x); }
inline double reject_distribution(double x) { return 1 + M_PI_2 * (0.5 - x); }
int main(int argc, char *argv[]) {
  // Variables
  Random rnd("Primes", "seed.in");
  unsigned int N = 100, M = 100000;
  int throws_per_block = M / N;

  /*************************/
  /* Exercise 02.1.1
  estimanate the integral with uniform sampling
  */
  vector<double> integral(N);
  vector<double> samples(M);

  for (auto &elem : samples) {
    elem = integrand(rnd.Rannyu());
  }

  for (int i{}; i < N; i++) {
    integral[i] = accumulate(samples.begin() + i * throws_per_block,
                             samples.begin() + (i + 1) * throws_per_block, 0.) /
                  static_cast<double>(throws_per_block);
  }
  // Error estimation
  vector<double> error = blocking_error(integral);

  // Out the results to file
  ofstream out("results02.1.uniform.dat");
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      out << (i + 1) << " " << integral[i] << " " << error[i] << endl;
    }
  } else {
    cerr << "ERROR: can't open output file" << endl;
    return 2;
  }

  out.close();

  /*************************/
  /* Exercise 02.1.2
  estimate the integral with non uniform sampling
  with good approx function, using rejection method
  */
  fill(integral.begin(), integral.end(), 0.);
  fill(samples.begin(), samples.end(), 0.);

  // upper bound of reject_distribution
  double max = 1. + M_PI * 0.25;
  for (auto &elem : samples) {
    double x = rnd.Accept_Reject(reject_distribution, max);
    elem = integrand(x) / reject_distribution(x);
  }

  for (int i{}; i < N; i++) {
    integral[i] = accumulate(samples.begin() + i * throws_per_block,
                             samples.begin() + (i + 1) * throws_per_block, 0.) /
                  static_cast<double>(throws_per_block);
  }
  // Error estimation
  fill(error.begin(), error.end(), 0.);
  error = blocking_error(integral);

  out.open("results02.1.reject.dat");
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      out << (i + 1) << " " << integral[i] << " " << error[i] << endl;
    }
  } else {
    cerr << "ERROR: can't open output file" << endl;
    return 2;
  }

  out.close();

  return 0;
}
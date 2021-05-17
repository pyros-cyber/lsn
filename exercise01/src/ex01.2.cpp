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

  // Declaring variables
  Random rnd("Primes", "seed.in");
  unsigned int M = 10000;
  double lambda = 1., mu = 1., gamma = 1.;
  vector<int> N = {1, 2, 10, 100};
  ofstream out_gauss("results01.2.gauss.dat"), out_exp("results01.2.exp.dat"),
      out_lorentz("results01.2.lorentz.dat");

  if (out_gauss.is_open() && out_exp.is_open() && out_lorentz.is_open()) {
    // Generating the results
    for (int i{}; i < N.size(); i++) {
      for (int j{}; j < M; j++) {
        double gauss = 0., exp = 0., lorentz = 0.;
        for (int k{}; k < N[i]; k++) {
          gauss += rnd.Rannyu();
          exp += rnd.Exp(lambda);
          lorentz += rnd.Lorentz(gamma, mu);
        }
        gauss /= static_cast<double>(N[i]);
        exp /= static_cast<double>(N[i]);
        lorentz /= static_cast<double>(N[i]);
        out_gauss << gauss << endl;
        out_exp << exp << endl;
        out_lorentz << lorentz << endl;
      }
    }
  } else {
    cerr << "ERROR: can't open output files." << endl;
    return 2;
  }

  out_gauss.close();
  out_exp.close();
  out_lorentz.close();
  rnd.SaveSeed();
  return 0;
}
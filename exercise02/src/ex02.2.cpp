#include "myStatFunc.h"
#include "random.h"
#include "randomWalk.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  // Variables
  Random rnd("Primes", "seed.in");
  unsigned int N = 100, M = 100000, steps = 100;
  int throws_per_block = M / N;

  /*************/
  /*Exercise 02.2.1
  simulate a 3d random walk on a discrete lattice
  */
  ofstream out("results02.2.discrete.dat");
  ofstream out1("results02.2.discrete.error.dat");

  if (out.is_open() && out1.is_open()) {
    for (int i{}; i < steps; i++) {
      vector<double> r_sqDiscrete(N);
      for (int j{}; j < N; j++) {
        for (int k{}; k < throws_per_block; k++) {
          shared_ptr<randomWalk> w = make_shared<DiscreteWalk>();
          for (int l{}; l < i; l++) {
            w->walk(&rnd);
          }
          r_sqDiscrete[j] += pow(w->distance_from_origin, 2);
        }
        r_sqDiscrete[j] /= static_cast<double>(throws_per_block);
      }
      vector<double> error = blocking_error(r_sqDiscrete);
      out << i << " " << sqrt(r_sqDiscrete[N - 1]) << " " << error[N - 1]
          << endl;
      if (i == steps - 1) {
        for (int j{}; j < N; j++) {
          out1 << j << " " << r_sqDiscrete[j] << " " << error[j] << endl;
        }
      }
    }
  } else {
    cerr << "ERROR: can't open output files" << endl;
  }

  out.close();
  out1.close();

  /*************/
  /*Exercise 02.2.2
  simulate a 3d random walk on a continuos lattice
  */

  out.open("results02.2.continuous.dat");
  out1.open("results02.2.continuous.error.dat");

  if (out.is_open() && out1.is_open()) {
    for (int i{}; i < steps; i++) {
      vector<double> r_sqContinuous(N);
      for (int j{}; j < N; j++) {
        for (int k{}; k < throws_per_block; k++) {
          shared_ptr<randomWalk> w = make_shared<ContinuousWalk>();
          for (int l{}; l < i; l++) {
            w->walk(&rnd);
          }
          r_sqContinuous[j] += pow(w->distance_from_origin, 2);
        }
        r_sqContinuous[j] /= static_cast<double>(throws_per_block);
      }
      vector<double> error = blocking_error(r_sqContinuous);
      out << i << " " << sqrt(r_sqContinuous[N - 1]) << " " << error[N - 1]
          << endl;
      if (i == steps - 1) {
        for (int j{}; j < N; j++) {
          out1 << j << " " << r_sqContinuous[j] << " " << error[j] << endl;
        }
      }
    }
  } else {
    cerr << "ERROR: can't open output files" << endl;
  }

  out.close();
  out1.close();

  return 0;
}
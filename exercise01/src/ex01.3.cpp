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
  // L: length of the needle, d: space between lines
  double L = 0.5, d = 1.;
  double av = 0., av2 = 0., sum = 0., sum2 = 0., err = 0.;
  int hit = 0;
  // M: number of throws, N: number of blocks
  unsigned int M = 100000, N = 100;
  int throws_per_block = M / N;
  ofstream out("results01.3.dat");

  // Estimating PI
  if (out.is_open()) {
    for (int i{}; i < N; i++) {
      hit = 0;
      for (int j{}; j < throws_per_block; j++) {
        double x = rnd.Rannyu(0., d);
        double angle = rnd.Buffon_Angle();
        if ((x - L / 2. * cos(angle)) <= 0. || (x + L / 2. * cos(angle)) >= d)
          hit++;
      }
      av = static_cast<double>(2 * L * throws_per_block / (hit * d));
      av2 = pow(av, 2);
      sum = static_cast<double>((sum * i + av) / (i + 1));
      sum2 = static_cast<double>((sum2 * i + av2) / (i + 1));
      if (i == 0) {
        err = 0;
      } else {
        err = error(sum, sum2, i);
      }
      out << i + 1 << " " << sum << " " << err << endl;
    }
  } else {
    cerr << "ERROR: can't open output file." << endl;
    return 2;
  }

  out.close();
  rnd.SaveSeed();

  return 0;
}
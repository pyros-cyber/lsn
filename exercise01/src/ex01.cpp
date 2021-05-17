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
  Random rnd("Primes", "seed.in");
  unsigned int M = 1000000, N = 100;
  unsigned int L = M / N;
  vector<double> distribution(M);
  for (auto &elem : distribution)
    elem = rnd.Rannyu();

  return 0;
}

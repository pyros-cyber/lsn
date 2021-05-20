#include "myStatFunc.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  Random rnd("../Primes", "../seed.in");
  /* # of estimates, # of blocks, # of discrete time steps */
  int M = 100000, N = 100, steps = 100;
  /* # of estimates per block */
  int throws_per_block = M / N;
  /* initial price of asset, strike price, delivery time, risk-free interest
   * rate, volatility */
  double S0 = 100., K = 100., T = 1., r = 0.1, sigma = 0.25;
  vector<double> callprice(N);
  vector<double> putprice(N);

  /****************/
  /* Exercise 03.1.1
   * Calculate the call and pull price by direct sampling of S(t=T)
   */

  ofstream out("results03.1.1.call.direct_sampling.dat");
  ofstream out1("results03.1.1.put.direct_sampling.dat");
  if (out.is_open() && out1.is_open()) {
    for (int i{}; i < N; i++) {
      double call = 0., put = 0.;
      for (int j{}; j < throws_per_block; j++) {
        double WT = rnd.Gauss(0., T);
        double ST = S0 * exp((r - sigma * sigma * 0.5) * T + sigma * WT);
        if (ST > K) {
          call += (ST - K) * exp(-r * T);
        } else {
          put += (K - ST) * exp(-r * T);
        }
      }
      callprice[i] = static_cast<double>(call / throws_per_block);
      putprice[i] = static_cast<double>(put / throws_per_block);
    }

    vector<double> call_error = blocking_error(callprice);
    vector<double> put_error = blocking_error(putprice);

    for (int i{}; i < N; i++) {
      out << i << " " << callprice[i] << " " << call_error[i] << endl;
      out1 << i << " " << putprice[i] << " " << put_error[i] << endl;
    }
  } else {
    cerr << "ERROR: can't open output file" << endl;
    return 2;
  }

  out.close();
  out1.close();

  /****************/
  /* Exercise 03.1.2
   * Calculate the call and pull price by discretized sampling of S(t=T)
   */

  out.open("results03.1.1.call.discretized_sampling.dat");
  out1.open("results03.1.1.put.discretized_sampling.dat");
  if (out.is_open() && out1.is_open()) {
    for (int i{}; i < N; i++) {
      double call = 0., put = 0.;
      for (int j{}; j < throws_per_block; j++) {
        double delta = static_cast<double>(T / steps);
        double ST = S0;
        for (int k{}; k < steps; k++) {
          double Z = rnd.Gauss(0., 1);
          ST = ST *
               exp((r - sigma * sigma * 0.5) * delta + sigma * Z * sqrt(delta));
        }
        if (ST > K) {
          call += (ST - K) * exp(-r * T);
        } else {
          put += (K - ST) * exp(-r * T);
        }
      }
      callprice[i] = static_cast<double>(call / throws_per_block);
      putprice[i] = static_cast<double>(put / throws_per_block);
    }

    vector<double> call_error = blocking_error(callprice);
    vector<double> put_error = blocking_error(putprice);

    for (int i{}; i < N; i++) {
      out << i << " " << callprice[i] << " " << call_error[i] << endl;
      out1 << i << " " << putprice[i] << " " << put_error[i] << endl;
    }
  } else {
    cerr << "ERROR: can't open output file" << endl;
    return 2;
  }

  out.close();
  out1.close();

  rnd.SaveSeed();
  return 0;
}

#ifndef ISING_H
#define ISING_H

#include "random.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

class Ising1D {
private:
  Random rnd;
  // keys for the map corresponding to properties:
  // - energy
  // - capacity
  // - magnetization
  // - susceptibility
  vector<string> props;
  int seed[4]; // to initialize the random number generator
  map<string, double> walker;
  // for blocking method averages:
  map<string, double> block_average;
  map<string, double> glob_average;
  map<string, double> glob_average2;

  double blk_norm, accepted, attempted;
  double stima_u, stima_c, stima_m, stima_x;
  double err_u, err_c, err_m, err_x;

  // configuration
  vector<double> spin_conf;
  int n_spin;
  // thermodynamical state
  double beta, temp, J, h;
  // simulation
  int nstep, nblk;
  // output files
  //  std::ofstream Ene, Heat, Mag, Chi;

  void Input();
  void Reset(int);
  void Accumulate();
  void BlockAverages(int);
  void ConfFinal();
  void Measure(int step = 0);

  // inline methods
  inline double EnergyGap(int i) {
    return 2. * J * spin_conf[i] *
               (spin_conf[Pbc(i - 1)] + spin_conf[Pbc(i + 1)]) +
           2. * h * spin_conf[i];
  }
  // Algorithm for periodic boundary conditions
  inline int Pbc(int i) const {
    if (i >= n_spin) {
      i -= n_spin;
    } else if (i < 0) {
      i += n_spin;
    }
    return i;
  }
  inline double Error(double sum, double sum2, int iblk) const {
    return sqrt((sum2 / iblk - pow(sum / iblk, 2)) / iblk);
  }

public:
  // The constructor of the model can start from a previous configuration;
  // the default option is that no configuration file is provided: in this
  // case the values of the spins are selected randomly.
  // the constructor makes use of a file "input.dat" where the Simulation
  // parameters are specified.
  Ising1D(string old_configuration = "");

  void MetropolisMove();
  void GibbsMove();
  // run the simulation (default is NOT to print instant values)
  void Run(function<void()>, bool instant = false);
};

#endif
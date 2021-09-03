#ifndef MCMD_H
#define MCMD_H

#include "random.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class mcmd {
private:
  /* VARIABLES */
  Random rnd;
  // Coordinates of the molecules
  vector<double> x, y, z;
  // r range in each histogram bin
  vector<double> r_range;
  // Keys for the maps of the observables
  vector<string> props;
  // Maps for blocking method averages
  map<string, double> walker;
  map<string, double> block_average;
  map<string, double> glob_average;
  map<string, double> glob_average_sq;
  // g(r) histogram and related variables for blocking method averages
  vector<double> g_histo;
  vector<double> g_histo_block_ave;
  vector<double> g_histo_glob_ave;
  vector<double> g_histo_glob_ave_sq;

  // simulation parameters
  double accepted, attempted, box, temp, rho, rcut, delta, beta, vol, vtail,
      ptail;
  unsigned int npart, nblk, nstep;
  const unsigned int nbins;
  bool instantenous;

  // output files
  ofstream Gave, Gerr, Epot, Pres, Binn, ist_pot, ist_pres;

  /* METHODS */
  // to reset the block average
  void Reset(unsigned int);
  // get block average
  void Accumulate();
  // blocking averages
  void Averages(unsigned int);
  // one MC step
  void Move();
  void ConfFinal() const;
  void ConfXYZ(unsigned int);
  void Measure(unsigned int istep = 0);
  double Boltzmann(double, double, double, unsigned int);
  void getRadiusRange(vector<double> &);

  /* INLINE METHODS */
  inline double Pbc(double r) const { return r - box * rint(r / box); }
  inline double Error(double sum, double sum_sq, unsigned int iblk) const {
    return sqrt((sum_sq / iblk - pow(sum / iblk, 2)) / iblk);
  }

public:
  mcmd(string, string, bool _instantenous = false);
  ~mcmd();
  void Run();
};

#endif

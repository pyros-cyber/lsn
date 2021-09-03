#ifndef MOLDYN_NVE_H
#define MOLDYN_NVE_H

#include "myStatFunc.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

struct MolDyn_NVE {
private:
  // thermodynamical state
  // # of particles in the Simulation
  unsigned int npart;
  // physical parameters
  double energy, temp, vol, rho, box, rcut;

  // simulation parameters:
  unsigned int nstep, iprint;
  double delta;

  Random rand;

  // time interval that separates each measurement:
  // compute block size (for the blocking averages) and setting block index
  // (L = M/N of the previous exercises, M is nstep/measure_time_interval, N is
  // nblocks, L is block_size)
  unsigned int measure_time_interval, n_blocks, block_size, imeasure, iblock;

  // output streams
  ofstream Epot, Ekin, Etot, Temp, Press, Gave, Gerr, Binn;

  // configuration:
  // positions, old positions, velocities, forces acting on each particle
  vector<double> x, y, z, xold, yold, zold, vx, vy, vz, fx, fy, fz;
  // observable properties:
  double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
  // vectors to get the blocking average:
  vector<double> est_pot, est_kin, est_etot, est_temp, est_press;
  // histogram for g(r) estimation
  vector<vector<double>> g_histo;
  // # of bins for estimation of g(r)
  const int nbins;
  // r range in each bin
  vector<double> r_range;

  // functions: I declare them as private because I will only call them
  // from methods inside the class itself
  void Input(string);
  // computes the forces on the particles (from the potential):
  void Force();
  // computes Periodic Boundary Conditions for box of length L=box
  inline double Pbc(double r) const { return r - box * rint(r / box); }
  // Moves the particles with Verlet algorithm
  void Move();
  // writes the final configuration
  void ConfFinal(string) const;
  // writes the configuration in .xyz format
  void ConfXYZ(int) const;
  // print values of properties to file
  void Measure();
  // computes and prints averages and statistical uncertaintes with blocking
  // method
  void BlockingResults();
  void getRadiusRange(vector<double> &);

public:
  /**
   * @brief Construct a new MolDyn_NVE object without old position
   * - first string is the filename with input parameters
   * - second string is the filename with initial molecular configuration
   */
  MolDyn_NVE(string, string);
  /**
   * @brief Construct a new MolDyn_NVE object with old position
   * - first and second string are as above
   * - the third string is the filename with old molecular configuration, used
   * to estrapolate the velocities
   */
  MolDyn_NVE(string, string, string);
  /**
   * @brief Runs a simulation using the Verlet algorithm
   *
   */
  ~MolDyn_NVE();
  void RunSimulation();
};

#endif

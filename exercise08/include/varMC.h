#ifndef VARMC_H
#define VARMC_H

#include "random.h"

using namespace std;

class Particle {
private:
  /* variational parameters of the wavefunction */
  double mu, mu_opt, sigma, sigma_opt;
  /* Metropolis parameters */
  double x, x_trial, x_jump;
  int accepted, attempted;
  /* Random variable of our class */
  Random rnd;

public:
  /**
   * @brief Default constructor of a Particle object
   */
  Particle(double _mu, double _sigma, double _x_start = 0., double _x_jump = 0.)
      : rnd{"Primes", "seed.in"} {
    mu = _mu;
    sigma = _sigma;
    mu_opt = NAN;
    sigma_opt = NAN;
    x = _x_start;
    x_jump = _x_jump;
    accepted = 0;
    attempted = 0;
  }

  /**
   * @brief Reset the values of the simulation
   */
  inline void Reset(double _mu, double _sigma, double _x_start = 0.,
                    double _x_jump = 0.) {
    mu = _mu;
    sigma = _sigma;
    x = _x_start;
    x_jump = _x_jump;
    accepted = 0;
    attempted = 0;
  }
  /**
   * @brief Density distribution function of the particle
   */
  inline double Distribution(double _x) {
    return pow(exp(-(_x - mu) * (_x - mu) / (sigma * sigma * 2.)) +
                   exp(-(_x + mu) * (_x + mu) / (sigma * sigma * 2.)),
               2);
  }
  /**
   * @brief Compute the local energy of the particle
   */
  inline double LocEnergy(double _x) {
    double s = sigma * sigma;
    return -0.5 / s *
               ((_x * _x + mu * mu) / s - 1. -
                (2. * _x * mu / s) * tanh(_x * mu / s)) +
           (_x * _x - 2.5) * _x * _x;
  }
  /**
   * @brief Return the Metropolis algorithm acceptance
   */
  inline double Acceptance() {
    return static_cast<double>(accepted) / attempted;
  }
  /**
   * @brief Getter function for the position
   */
  inline double getX() { return x; }
  /**
   * @brief Set the results of the simulation
   */
  inline void Results(double &_mu, double &_sigma) {
    _mu = mu_opt;
    _sigma = sigma_opt;
  }
  /**
   * @brief Make a step of the Metropolis algorithm, which samples the density
   * distribution function of the particle, and use a uniform transition
   * probability
   */
  void MetropolisMove();
};
#endif

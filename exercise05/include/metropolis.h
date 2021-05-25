#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "random.h"
#include <algorithm>
#include <functional>
#include <string>

using namespace std;
/**
 * @brief My implementation of the Metropolis Algorithm
 *
 */
struct Metropolis {
  Random rnd;
  double x_start, y_start, z_start, step;
  int accepted, total;

  /**
   * @brief Construct a new Metropolis object
   *
   * @param _step : step of the algorithm
   * @param _x_start : starting point
   * @param _y_start : starting point
   * @param _z_start : starting point
   */
  Metropolis(double _step, double _x_start = 0., double _y_start = 0.,
             double _z_start = 0.)
      : rnd{"../Primes", "../seed.in"} {
    x_start = _x_start;
    y_start = _y_start;
    z_start = _z_start;
    step = _step;
    accepted = 0;
    total = 0;
  }
  /**
   * @brief Makes a uniformly distributed step
   *
   * @param initial_coordinate
   * @return double
   */
  inline double uniform_step(double initial_coordinate) {
    return initial_coordinate + rnd.Rannyu(-step, step);
  }
  /**
   * @brief Makes a gaussian distributed step
   *
   * @param initial_coordinate
   * @return double
   */
  inline double gaussian_step(double initial_coordinate) {
    return rnd.Gauss(initial_coordinate, step / 2.);
  }
  /**
   * @brief One step of the algorithm. Accept to functions, representing the
   * transition probability and the distribution function used to decide if the
   * point is accepted or no
   *
   */
  void Run(function<double(double)>, function<double(double, double, double)>);

  /**
   * @brief Return the algorithm acceptance
   *
   * @return double
   */
  inline double acceptance() const {
    return static_cast<double>(accepted / total);
  }
};

#endif
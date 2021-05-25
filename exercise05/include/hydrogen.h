#ifndef HYDROGEN_H
#define HYDROGEN_H

#include <cmath>
#include <string>
#include <vector>

using namespace std;
/**
 * @brief Header only library representing a hydrogen atom
 *
 */
struct HydrogenAtom {
  double x, y, z;
  string state;
  /**
   * @brief Construct a new Hydrogen Atom object
   *
   * @param _x : coordinate in the space
   * @param _y : coordinate in the space
   * @param _z : coordinate in the space
   * @param _state : state of the atom, can either be `ground state` or `first
   * excited`
   */
  HydrogenAtom(double _x, double _y, double _z, string _state) {
    x = _x;
    y = _y;
    z = _z;
    state = _state;
  }

  /**
   * @brief Get the position of the atom
   *
   * @return vector<double>
   */
  inline vector<double> get_position() {
    vector<double> position{x, y, z};
    return position;
  }

  /**
   * @brief Get the radius of the atom
   *
   * @return double
   */
  inline double get_radius() { return sqrt(x * x + y * y + z * z); }
};

/**
 * @brief Probability distribution of the hydrogen atom ground state.
 *
 * @param u
 * @param v
 * @param w
 * @return double
 */
inline double ground_state(double u, double v, double w) {
  return exp(-2. * sqrt(u * u + v * v + w * w)) / M_PI;
}

/**
 * @brief Probability distribution of the hydrogen atom first excited state.
 *
 * @param u
 * @param v
 * @param w
 * @return double
 */
inline double first_excited_state(double u, double v, double w) {
  double r = u * u + v * v + w * w;
  return exp(-sqrt(r)) * r * r / (32. * M_PI);
}
#endif
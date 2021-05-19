#include "randomWalk.h"

using namespace std;

void DiscreteWalk::set_distance() {
  double d_sq = 0.;
  for (auto &elem : position) {
    d_sq += elem * elem;
  }
  distance_from_origin = sqrt(d_sq);
}

void DiscreteWalk::walk(Random *rnd) {
  int chosen_nn = static_cast<int>(rnd->Rannyu(1, 7));
  if (chosen_nn % 2 == 0) {
    position[chosen_nn / 2 - 1] += 1;
  } else {
    position[(chosen_nn + 1) / 2 - 1] += -1;
  }
  set_distance();
  n_step++;
}

void ContinuousWalk::set_distance() {
  double d_sq = 0.;
  for (auto &elem : position) {
    d_sq += elem * elem;
  }
  distance_from_origin = sqrt(d_sq);
}

void ContinuousWalk::walk(Random *rnd) {
  double theta = 0., phi = 0.;
  rnd->Solid_Angle(theta, phi);
  position[0] += sin(theta) * cos(phi);
  position[1] += sin(theta) * sin(phi);
  position[2] += cos(theta);
  set_distance();
  n_step++;
}
#ifndef RANDOMWALK_H
#define RANDOMWALK_H

#include "random.h"
#include <memory>

struct randomWalk {
  int n_step;
  double distance_from_origin;

  virtual void set_distance() = 0;
  virtual void walk(Random *) = 0;
};

struct DiscreteWalk : public randomWalk {
  int position[3];

  DiscreteWalk() {
    for (auto &elem : position)
      elem = 0;
  }
  void set_distance();
  void walk(Random *);
};

struct ContinuousWalk : public randomWalk {
  double position[3];

  ContinuousWalk() {
    for (auto &elem : position)
      elem = 0.;
  }
  void set_distance();
  void walk(Random *);
};

#endif
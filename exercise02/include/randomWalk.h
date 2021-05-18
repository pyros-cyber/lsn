#ifndef RANDOMWALK_H
#define RANDOMWALK_H

#include "random.h"
#include <memory>

struct randomWalk {
  int n_step;
  double distance_from_origin;

  virtual void set_distance() = 0;
  virtual void walk(std::shared_ptr<Random>) = 0;
};

struct DiscreteWalk : public randomWalk {
  int position[3];

  DiscreteWalk() {
    for (auto &elem : position)
      elem = 0;
  }
  void set_distance();
  void walk(std::shared_ptr<Random>);
};

struct ContinuosWalk : public randomWalk {
  double position[3];

  ContinuosWalk() {
    for (auto &elem : position)
      elem = 0.;
  }
  void set_distance();
  void walk(std::shared_ptr<Random>);
};

#endif
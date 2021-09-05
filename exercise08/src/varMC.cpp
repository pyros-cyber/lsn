#include "varMC.h"

void Particle::MetropolisMove() {
  x_trial = (rnd.Rannyu() - 0.5) * 2. * x_jump + x;
  double p = min(1., Distribution(x_trial) / Distribution(x));

  if (p > rnd.Rannyu()) {
    x = x_trial;
    accepted++;
  }
  attempted++;
}

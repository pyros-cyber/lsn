#include "metropolis.h"

void Metropolis::Run(function<double(double)> trans_prob,
                     function<double(double, double, double)> distr_prob) {
  double x_trial = trans_prob(x_start);
  double y_trial = trans_prob(y_start);
  double z_trial = trans_prob(z_start);
  double alpha = min(1., distr_prob(x_trial, y_trial, z_trial) /
                             distr_prob(x_start, y_start, z_start));
  if (rnd.Rannyu() <= alpha) {
    x_start = x_trial;
    y_start = y_trial;
    z_start = z_trial;
    accepted++;
  }
  total++;
}
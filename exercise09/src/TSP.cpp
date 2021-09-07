#include "TSP.h"

Salesman::Salesman(shared_ptr<Random> _rnd, string _shape, int _Ncities,
                   double _pair_prob = 0.5, double _multi_prob = 0.5,
                   double _inv_prob = 0.5, double _shift_prob = 0.5)
    : rnd{_rnd}, Cities{_Ncities} {

  Ncities = _Ncities;
  pair_prob = _pair_prob;
  multi_prob = _multi_prob;
  inv_prob = _inv_prob;
  shift_prob = _shift_prob;

  for (int i{1}, i <= Ncities; ++i) {
    Path.push_back(i);
  }

  random_shuffle(Path.begin() + 1, Path.end());

  if (shape == "circumference") {
    double theta;
    for (int i{}; i < Ncities; ++i) {
      theta = rnd->Rannyu(0, 2 * M_PI);
      Cities[i].first = cos(theta);
      Cities[i].second = sin(theta);
    }
  } else if (shape == "square") {
    for (int i{}; i < Ncities; ++i) {
      Cities[i].first = 2 * rnd->Rannyu() - 1.;
      Cities[i].second = 2 * rnd->Rannyu() - 1.;
    }
  } else {
    cerr << "ERROR: unknown shape. Exiting." << endl;
    exit(1);
  }
}

bool Salesman::Check() {
  bool check = true;
  int count;
  if (Path[0] != 1) {
    check = false;
  } else {
    for (int i{1}; i <= Ncities; ++i) {
      count = 0;
      for (int j{}; j < Ncities; ++j) {
        if (Path[j] == i) {
          count++;
        }
      }
      if (count > 1) {
        check = false;
        break;
      }
    }
  }
  return check;
}

double Salesman::AbsoluteLoss() {
  double loss = 0;
  for (int i{}; i < Ncities; ++i) {
    loss += sqrt(
        pow(Cities[Path[i] - 1].first - Cities[Path[PbcPath(i + 1)] - 1].first,
            2) +
        pow(Cities[Path[i] - 1].second -
                Cities[Path[PbcPath(i + 1)] - 1].second,
            2));
  }

  return loss;
}

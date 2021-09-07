#include "TSP.h"

Salesman::Salesman(shared_ptr<Random> _rnd, string _shape, int _Ncities,
                   double _pair_prob = 0.5, double _multi_prob = 0.5,
                   double _inv_prob = 0.5, double _shift_prob = 0.5)
    : rnd{_rnd} {

  Ncities = _Ncities;
  Cities.resize(Ncities);
  pair_prob = _pair_prob;
  multi_prob = _multi_prob;
  inv_prob = _inv_prob;
  shift_prob = _shift_prob;

  for (int i{1}; i <= Ncities; ++i) {
    Path.push_back(i);
  }

  random_shuffle(Path.begin() + 1, Path.end());

  if (_shape == "circumference") {
    double theta;
    for (int i{}; i < Ncities; ++i) {
      theta = rnd->Rannyu(0, 2 * M_PI);
      Cities[i].first = cos(theta);
      Cities[i].second = sin(theta);
    }
  } else if (_shape == "square") {
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

void Salesman::PairPermutation() {
  int i1, i2;
  if (rnd->Rannyu() < pair_prob) {
    i1 = rnd->RandInt(1, Ncities - 1);
    do {
      i2 = rnd->RandInt(1, Ncities - 1);
    } while (i2 == i1);
    swap(Path[i1], Path[i2]);
  }
}

void Salesman::Shift() {
  int i1, i2, n;
  if (rnd->Rannyu() < shift_prob) {
    i1 = rnd->RandInt(1, Ncities - 1);
    i2 = rnd->RandInt(1, Ncities - 1);
    n = rnd->RandInt(1, i2);
    for (int i{}; i < n; ++i) {
      for (int j{}; j < i2 - 1; ++j) {
        swap(Path[PbcMutation(i1 + j)], Path[PbcMutation(i1 + i2 - 1)]);
      }
    }
  }
}

void Salesman::MultiPermutation() {
  int m, i1, i2;
  if (rnd->Rannyu() < multi_prob) {
    m = rnd->RandInt(1, ((Ncities + 1) / 2) - 1); // 1 <= m < Ncities/2
    i1 = rnd->RandInt(1, Ncities - 1);
    i2 = PbcMutation(rnd->RandInt(i1 + m, Ncities + i1 - m - 1));
    for (int i{}; i < m; ++i) {
      swap(Path[PbcMutation(i1 + i)], Path[PbcMutation(i2 + i)]);
    }
  }
}

void Salesman::Inversion() {
  int i1, m;
  if (rnd->Rannyu() < inv_prob) {
    i1 = rnd->RandInt(1, Ncities - 1);
    m = rnd->RandInt(1, Ncities);
    for (int i{}; i < abs(m / 2); ++i) {
      swap(Path[PbcMutation(i1 + i)], Path[PbcMutation(i1 + (m - 1) - i)]);
    }
  }
}

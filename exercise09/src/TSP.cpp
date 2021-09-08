#include "TSP.h"

Salesman::Salesman(shared_ptr<Random> _rnd, string _shape, int _Ncities,
                   double _pair_prob, double _multi_prob, double _inv_prob,
                   double _shift_prob)
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
  // cout << loss << endl;

  return loss;
}

vector<pair<double, double>> Salesman::GetCities() {
  vector<pair<double, double>> ordered_cities(Ncities);
  for (int i{}; i < Ncities; ++i) {
    ordered_cities[i] = Cities[Path[i] - 1];
  }
  return ordered_cities;
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

Population::Population(shared_ptr<Random> _rnd, string _shape, int _Npop,
                       int _Ncities, double _cross_prob)
    : rnd{_rnd} {
  Npop = _Npop;
  Ncities = _Ncities;
  cross_prob = _cross_prob;

  Losses.resize(Npop);
  fill(Losses.begin(), Losses.end(), 0.);

  for (int i{}; i < Npop; ++i) {
    Pop.push_back(Salesman(_rnd, _shape, Ncities));
  }

  vector<pair<double, double>> cities(Ncities);

  if (_shape == "circumference") {
    double theta;
    for (int i{}; i < Ncities; ++i) {
      theta = rnd->Rannyu(0, 2 * M_PI);
      cities[i].first = cos(theta);
      cities[i].second = sin(theta);
    }

    for (int i{}; i < Npop; ++i) {
      Pop[i].Cities = cities;
    }
  } else if (_shape == "square") {

    for (int i{}; i < Ncities; ++i) {
      cities[i].first = 2 * rnd->Rannyu() - 1.;
      cities[i].second = 2 * rnd->Rannyu() - 1.;
    }

    for (int i{}; i < Npop; ++i) {
      Pop[i].Cities = cities;
    }
  } else {
    cerr << "ERROR: unknown shape. Exiting." << endl;
    exit(1);
  }
}

double Population::LossesAverage() {
  double av = 0.;
  for (int i{Npop - 1}; i > Npop / 2; --i) {
    av += Losses[i];
  }
  av /= (Npop / 2 - 1);
  return av;
}

void Population::OrderPop() {

  vector<Mileage> mileage(Npop);
  for (int i{}; i < Npop; ++i) {
    mileage[i].Chromo = Pop[i];
    mileage[i].Miles = Pop[i].AbsoluteLoss();
  }

  /* We order the population on a fitness basis:
   * the higher AbsoluteLoss, the lower the rank
   */
  sort(mileage.begin(), mileage.end(),
       [](Mileage const &a, Mileage const &b) { return a.Miles > b.Miles; });

  for (int i{}; i < Npop; ++i) {
    Pop[i] = mileage[i].Chromo;
    Losses[i] = mileage[i].Miles;
  }
}

vector<vector<int>> Population::Crossover() {
  vector<vector<int>> offspring(2);
  vector<int> offspring1, offspring2;
  int i1, i2, cut, start1 = 1, start2 = 1, elem;

  i1 = Selection();
  do {
    i2 = Selection();
  } while (i2 == i1);
  if (rnd->Rannyu() < cross_prob) {
    cut = rnd->RandInt(1, Ncities - 1);

    offspring1 = Pop[i1].Path;
    offspring2 = Pop[i2].Path;

    for (int i{cut}; i < Ncities; ++i) {
      for (int j{start1}; j < Ncities; ++j) {
        elem = Pop[i2].Path[j];
        if (find(offspring1.begin(), offspring1.begin() + cut, elem) ==
            offspring1.begin() + cut) {
          offspring1[i] = elem;
          start1++;
          break;
        } else {
          start1++;
        }
      }
      for (int j{start2}; j < Ncities; ++j) {
        elem = Pop[i1].Path[j];
        if (find(offspring2.begin(), offspring2.begin() + cut, elem) ==
            offspring2.begin() + cut) {
          offspring2[i] = elem;
          start2++;
          break;
        } else {
          start2++;
        }
      }
    }
    offspring[0] = offspring1;
    offspring[1] = offspring2;
  } else {
    offspring[0] = Pop[i1].Path;
    offspring[1] = Pop[i2].Path;
  }
  return offspring;
}

void Population::Mutations() {
  for (int i{}; i < Npop; ++i) {
    Pop[i].PairPermutation();
    Pop[i].Shift();
    Pop[i].MultiPermutation();
    Pop[i].Inversion();
  }
}

void Population::Evolve(int gen, string shape) {
  vector<vector<int>> offsprings(2);
  vector<vector<int>> pop(Npop);
  ofstream loss("results/loss_" + shape + ".dat");
  ofstream lossave("results/loss_ave_" + shape + ".dat");
  if (loss.is_open() && lossave.is_open()) {
    for (int j{}; j < gen; ++j) {
      for (int i{}; i < Npop / 2; ++i) {
        offsprings = this->Crossover();
        pop[2 * i] = offsprings[0];
        pop[2 * i + 1] = offsprings[1];
      }

      for (int i{}; i < Npop; ++i) {
        this->Pop[i].Path = pop[i];
      }
      this->Mutations();
      this->OrderPop();
      loss << j + 1 << " " << Losses[Npop - 1] << endl;
      lossave << j + 1 << " " << LossesAverage() << endl;
    }
  } else {
    cerr << "ERROR: can't open output files! Exiting." << endl;
    exit(1);
  }
  loss.close();
  lossave.close();
}

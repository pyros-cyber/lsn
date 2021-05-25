/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"

using namespace std;

Random ::Random() {}
Random ::Random(string primes, string seed) {
  int s[4];
  int p1 = 0, p2 = 0;
  ifstream Primes(primes);

  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else {
    cerr << "ERROR: unable to open " << primes << " file" << endl;
  }
  Primes.close();

  ifstream input(seed);
  string property;

  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> s[0] >> s[1] >> s[2] >> s[3];
        SetRandom(s, p1, p2);
      }
    }
    input.close();
  } else {
    std::cerr << "ERROR: unable to open " << seed << " file" << endl;
  }
}

Random ::~Random() {}

void Random ::SetRandom(int *s, int p1, int p2) {
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0] % 4096;
  l2 = s[1] % 4096;
  l3 = s[2] % 4096;
  l4 = s[3] % 4096;
  l4 = 2 * (l4 / 2) + 1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;
}

void Random ::SaveSeed() {
  ofstream WriteSeed;
  WriteSeed.open("seed.out");
  if (WriteSeed.is_open()) {
    WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
  } else
    cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
}

double Random ::Rannyu(void) {
  const double twom12 = 0.000244140625;
  int i1, i2, i3, i4;
  double r;

  i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
  i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
  i3 = l3 * m4 + l4 * m3 + n3;
  i4 = l4 * m4 + n4;
  l4 = i4 % 4096;
  i3 = i3 + i4 / 4096;
  l3 = i3 % 4096;
  i2 = i2 + i3 / 4096;
  l2 = i2 % 4096;
  l1 = (i1 + i2 / 4096) % 4096;
  r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

  return r;
}

double Random ::Rannyu(double min, double max) {
  return min + (max - min) * Rannyu();
}

double Random ::Gauss(double mean, double sigma) {
  double s = Rannyu();
  double t = Rannyu();
  double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
  return mean + x * sigma;
}

double Random::Exp(double lambda) { return -log(1. - Rannyu()) / lambda; }

double Random::Lorentz(double gamma, double mu) {
  return mu + gamma * tan(M_PI * (Rannyu() - 0.5));
}

double Random::Buffon_Angle() {
  double x, y, r;
  do {
    x = this->Rannyu();
    y = this->Rannyu();
    r = sqrt(x * x + y * y);
  } while (r >= 1.);
  return acos(x / r);
}

void Random::Solid_Angle(double &theta, double &phi) {
  double x, y, z, r_sq;
  do {
    x = this->Rannyu(-1., 1.);
    y = this->Rannyu(-1., 1.);
    z = this->Rannyu(-1., 1.);
    r_sq = x * x + y * y + z * z;
  } while (r_sq >= 1.);
  theta = acos(z / sqrt(r_sq));
  if (y < 0.) {
    phi = -acos(x / (sqrt(r_sq) * sin(theta)));
  } else {
    phi = acos(x / (sqrt(r_sq) * sin(theta)));
  }
}
double Random::Accept_Reject(std::function<double(double)> f, double max) {
  double x = 0., r = 0.;
  do {
    x = this->Rannyu();
    r = this->Rannyu();
  } while (r >= f(x) / max);
  return x;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

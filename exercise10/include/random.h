/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

class Random {

private:
  int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;

protected:
public:
  // constructors
  Random();
  Random(int);
  Random(std::string, std::string);
  // destructor
  ~Random();
  // methods
  void SetRandom(int *, int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double, double);
  int RandInt(int, int);
  double Gauss(double, double);
  double Exp(double);
  double Lorentz(double, double);
  double Buffon_Angle();
  void Solid_Angle(double &, double &);
  double Accept_Reject(std::function<double(double)>, double);
};

#endif // RANDOM_H

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

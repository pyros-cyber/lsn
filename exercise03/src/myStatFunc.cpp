#include "myStatFunc.h"

using namespace std;

vector<double> blocking_error(vector<double> &average) {
  vector<double> err(average.size(), 0.);
  double sum = 0., sum2 = 0.;
  sum = average[0];
  sum2 = pow(sum, 2);

  for (int i{1}; i < average.size(); i++) {
    sum += average[i];
    sum2 += pow(average[i], 2);
    average[i] = sum / static_cast<double>(i + 1);
    err[i] = sqrt((sum2 / static_cast<double>(i + 1) - pow(average[i], 2)) /
                  static_cast<double>(i));
  }
  return err;
}

double error(double av, double av2, int n) {
  if (n == 0)
    return 0;
  else
    return sqrt((av2 - pow(av, 2)) / static_cast<double>(n));
}
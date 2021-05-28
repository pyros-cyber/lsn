#ifndef MYSTATFUNC_H
#define MYSTATFUNC_H
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

/**
 * @brief Takes a vector as input, where each element is the average of a block
 * and returns which elements are the progressive blocking error of the input
 * vector, and modifies the input vector to obtain the progressive average of
 * the vector
 *
 * @param average
 * @return std::vector<double>
 */
std::vector<double> blocking_error(std::vector<double> &average);

/**
 * @brief uncertainty calculated as Standard Deviation
 *
 * @param av
 * @param av2
 * @param n
 * @return double
 */
double error(double av, double av2, int n);

#endif
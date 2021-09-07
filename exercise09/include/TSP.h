#ifndef TSP_H
#define TSP_H

#include "random.h"
#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

using namespace std;

/**
 * @brief Class representing a Travelling Salesman
 */
struct Salesman {
private:
  int Ncities;
  double pair_prob, multi_prob, inv_prob, shift_prob;
  shared_ptr<Random> rnd;

public:
  vector<int> Path;
  vector<pair<double, double>> Cities;

  /**
   * @brief Construct a new Salesman object
   * @param _rnd: shared_ptr to Random object
   * @param _shape: shape of the path ("circular" or "square")
   * @param _Ncities: number of cities to visit
   * @param _pair_prob: probability of pair permutation
   * @param _multi_prob: probability of permutation
   * @param _inv_prob: probability of inversion
   * @param _shift_prob: probability of shift
   */
  Salesman(shared_ptr<Random>, string, int, double, double, double, double);

  /**
   * @brief Check function to check if Salesman fulfills the bonds
   */
  bool Check();
  /**
   * @brief Loss function defined as the sum of the absolute distance
   * between the cities visited by the Salesman
   */
  double AbsoluteLoss();
  /**
   * @brief Pbc for the Salesman's path
   */
  inline int PbcPath(int i) { return i % Ncities; }
  /**
   * @brief Pbc used when applying mutations
   */
  inline int PbcMutation(int i) {
    if (abs(i) < Ncities)
      return i;
    else
      return PbcPath(i) + 1;
  }
  /**
   * @brief Permutes a pair of elements
   */
  void PairPermutation();
  /**
   * @brief Shifts elements
   */
  void Shift();
  /**
   * @brief Permute elements
   */
  void MultiPermutation();
  /**
   * @brief Invert elements
   */
  void Inversion();
};
#endif

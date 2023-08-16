#pragma once

#include "BasisGenerator.hpp"
#include "container.hpp"

#include <Eigen/Dense>
#include <vector>

class OperatorGenerator {
  BasisGenerator* basis_generator;

public:
  OperatorGenerator() = default;

  explicit OperatorGenerator(BasisGenerator* generator)
      : basis_generator(generator) {}

  /**
   * @brief Calculates the F matrix based on input values.
   *
   * The F matrix is a 2x2 matrix that depends on four integer parameters (a1,
   * a2, a3, and outcome). The values of these parameters determine the values
   * in the F matrix as per the provided conditions.
   *
   * @param a1 Integer parameter a1.
   * @param a2 Integer parameter a2.
   * @param a3 Integer parameter a3.
   * @param outcome Integer parameter outcome.
   * @return Eigen::Matrix<double, 2, 2> The F matrix.
   */
  Eigen::Matrix<double, 2, 2> F(int a1, int a2, int a3, int outcome);

  /**
   * @brief Calculates the R matrix based on input values.
   *
   * The R matrix is a 2x2 matrix that depends on two integer parameters (a1 and
   * a2). The values of these parameters determine the values in the R matrix as
   * per the provided conditions.
   *
   * @param a1 Integer parameter a1.
   * @param a2 Integer parameter a2.
   * @return Eigen::Matrix<std::complex<double>, 2, 2> The R matrix.
   */
  Eigen::Matrix<std::complex<double>, 2, 2> R(int a1, int a2);

  /**
   * @brief Calculates the B matrix based on input values.
   * Braid matrix
   *
   *  @param a0 Matrix representing a0
   *  @param a1 Matrix representing a1
   *  @param a2 Matrix representing a2
   *  @param outcome: Matrix representing outcome
   *  @return b_matrix The braid matrix computed using Eigen.
   */
  Eigen::Matrix<std::complex<double>, 2, 2> B(int a0, int a1, int a2,
                                              int outcome);

  /**
   * Amplitude of getting state_f by applying the braiding operator
   * sigma_{index} on state_i.
   *
   * @param index The index of the braiding operator.
   * @param state_f The final state represented as a vector of integers.
   * @param state_i The initial state represented as a vector of integers.
   *
   * @return The component (state_f, state_i) of the sigma_{index} matrix.
   *
   * @throws std::invalid_argument If the index value is not valid (less than or
   * equal to 0, or greater than the size of state_i).
   */
  std::complex<double> sigma(int index, const std::vector<int>& state_f,
                             const std::vector<int>& state_i);

  /**
   * L matrix component used in the calculation of braiding between two anyons
   * separated in two qudits.
   *
   * @param k int: k
   * @param h int: i_{m(q-1)}
   * @param i_ int: i'_{mq}
   * @param i int: i_{mq}
   * @param jj_ std::vector<int>: [i'_{(m+1)1},....i'_{(m+1)q}]
   * @param jj std::vector<int>: [i_{(m+1)1},....i_{(m+1)q}]
   *
   * @return The L matrix component.
   */
  std::complex<double> L(int k, int h, int i_, int i,
                         const std::vector<int>& jj_,
                         const std::vector<int>& jj);

  /**
   * S matrix or sewing matrix used in the calculation of the braiding operator
   * between two anyons separated between two qudits not fused immediately.
   *
   * @param jm int: j_m
   * @param jmo int: j_{m-1}
   * @param jmoo int: j_{m-2}
   * @param jmo_ int: j'_{m-1}
   * @param h int: i_{m(q-1)}
   * @param i_ int: i'_{mq}
   * @param i int: i_{mq}
   * @param jj_ std::vector<int>: [i'_{(m+1)1},....i'_{(m+1)q}]
   * @param jj std::vector<int>: [i_{(m+1)1},....i_{(m+1)q}]
   *
   * @return The S matrix component.
   */
  std::complex<double> S(int jm, int jmo, int jmoo, int jmo_, int h, int i_,
                         int i, const std::vector<int>& jj_,
                         const std::vector<int>& jj);

  /**
   * Generate the Sigma matrix component used in the calculation of the braiding
   * operator between two anyons separated between two qudits not fused
   * immediately.
   *
   * @param index int: The index of the Sigma matrix component.
   * @param state_i std::unordered_map<std::string,
   * std::vector<std::vector<int>>>: The initial state represented as a
   * dictionary containing "qudits" and "roots".
   * @param state_f std::unordered_map<std::string,
   * std::vector<std::vector<int>>>: The final state represented as a dictionary
   * containing "qudits" and "roots".
   *
   * @return The Sigma matrix component.
   */
  std::complex<double> gen_sigma(int index, const State& state_i,
                                 const State& state_f);

  /**
   * Generate the braiding operator of index 'index' for a system of
   * a given number of qudits and anyons per qudit. This operator braids
   * anyons at positions 'index' and 'index'+1.
   *
   * @param index The operator's index.
   * @param nb_qudits Number of qudits in the circuit.
   * @param nb_anyons_per_qudit Number of anyons in each qudit.
   * @return Matrix representation of the braiding operator as a list of lists
   * (sigmas).
   */
  Eigen::MatrixXcd generate_braiding_operator(int index, int nb_qudits,
                                              int nb_anyons_per_qudit);
};

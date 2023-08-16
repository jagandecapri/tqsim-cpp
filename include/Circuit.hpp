#pragma once

#include "BasisGenerator.hpp"
#include "CircuitInterface.hpp"
#include "OperatorGenerator.hpp"
#include "utils.hpp"

#include <iostream>
#include <map>
#include <pcg/pcg_extras.hpp>
#include <pcg/pcg_random.hpp>
#include <random>

class Circuit : public CircuitInterface {
private:
  /* data */
  OperatorGeneratorInterface* operator_generator;
  BasisGeneratorInterface* basis_generator;
  int nb_qudits;
  int nb_anyons_per_qudit;
  int nb_anyons;
  int nb_braids;
  bool measured;
  Basis basis;
  size_t dim;
  std::vector<std::tuple<int, int>> braids_history;
  Eigen::VectorXcd initial_state;
  std::vector<Eigen::MatrixXcd> sigmas;
  Eigen::MatrixXcd unitary;

public:
  Circuit(int nb_qudits, int nb_anyons_per_qudit);

  int get_nb_qudits() override { return nb_qudits; }

  int get_nb_anyons_per_qudit() override { return nb_anyons_per_qudit; }

  int get_dim() override { return dim; }

  Basis get_basis() override { return basis; }

  std::vector<std::tuple<int, int>> get_braids_history() override {
    return braids_history;
  }

  std::vector<Eigen::MatrixXcd> get_braiding_operators() override {
    return sigmas;
  }

  Eigen::MatrixXcd get_unitary() override { return unitary; }

  void generate_basis();

  void get_sigmas();

  /**
   * Initializes the circuit in the state input_state.
   *
   * @param input_state A normalized quantum state with the same dimensions as
   * the fusion space. For example, for a 1-qudit circuit with 3 anyons,
   *                    input_state must be a 3 dimensional vector with norm 1.
   *
   * @throws std::exception Will be thrown if initialization is attempted after
   * performing braiding operations.
   * @throws std::invalid_argument Will be thrown if the input state has the
   * wrong dimension or is not normalized.
   *
   * @return A reference to the same circuit.
   */
  CircuitInterface& initialize(const Eigen::VectorXcd& input_state) override;

  CircuitInterface& braid(int n, int m) override;

  /**
   * Takes a sequence of [sigma operator, power], and applies the successive
   * operators to the 'power'. The first operator in the sequence is the first
   * to be applied.
   *
   * @param braid A sequence of pairs of integers representing a braiding
   * operator and an exponent.
   *
   * @throws std::invalid_argument Is thrown if one of the operators' indices is
   * not an integer greater or equal to 1, or is an incorrect index.
   *
   * @return A reference to the same circuit.
   */
  Circuit& braid_sequence(Sequence& braid) override;

  Circuit& measure() override;

  // Constructor and other methods

  /**
   * Computes and returns the current state vector of the circuit.
   *
   * @return Eigen::VectorXcd The state vector of the circuit.
   */
  Eigen::VectorXcd statevector() override { return unitary * initial_state; }

  /**
   * Simulates the quantum circuit for a specified number of shots and returns
   * the measurement results.
   *
   * @param shots Number of times the circuit is simulated.
   * @return std::map<int, int> Contains the number of measurements for each
   * measured state.
   * @throws std::runtime_error Is thrown if the circuit is run without a
   * measurement.
   */
  Result run(int shots) override;
};
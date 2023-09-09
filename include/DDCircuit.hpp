#pragma once

#include "BasisGenerator.hpp"
#include "Circuit.hpp"
#include "OperatorGenerator.hpp"
#include "container.hpp"
#include "dd/Package.hpp"
#include "utils.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <pcg/pcg_extras.hpp>
#include <pcg/pcg_random.hpp>
#include <random>

class DDCircuit {
private:
  /* data */
  std::unique_ptr<OperatorGenerator> operatorGenerator;
  std::unique_ptr<BasisGenerator> basisGenerator;
  int nbQudits{};
  int nbAnyonsPerQudit{};
  int nbAnyons{};
  int nbBraids{};
  bool measured{};
  Basis basis;
  size_t dim{};
  std::vector<std::tuple<int, int>> braidsHistory;
  dd::vEdge currentState;
  std::vector<Eigen::MatrixXcd> sigmas;
  std::unique_ptr<dd::Package<>> circuitDD;
  BraidingOperators braidingOperators;
  Eigen::MatrixXcd unitary;

public:
  DDCircuit() = default;
  DDCircuit(int nbQudits, int nbAnyonsPerQudit);

  [[nodiscard]] int getNbQudits() const { return nbQudits; }

  [[nodiscard]] int getNbAnyonsPerQudit() const { return nbAnyonsPerQudit; }

  [[nodiscard]] size_t getDim() const { return dim; }

  Basis getBasis() { return basis; }

  std::vector<std::tuple<int, int>> getBraidsHistory() { return braidsHistory; }

  std::vector<Eigen::MatrixXcd> getBraidingOperators() { return sigmas; }

  Eigen::MatrixXcd getUnitary() { return unitary; }

  void generateBasis();

  void getSigmas();

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
  void initialize(const std::vector<std::complex<double>>& inputState);

  void braid(int n, int m);
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
  void braidSequence(Sequence& braid);

  void measure();

  // Constructor and other methods

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
  ResultDD run(int shots, std::size_t seed);
};
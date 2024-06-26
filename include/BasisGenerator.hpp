#pragma once

#include "container.hpp"

#include <vector>

class BasisGenerator {
public:
  /**
   * Checks whether the Fibonacci fusion rules are obeyed for given anyon
   * charges and outcome charge.
   *
   * @param anyon1 Anyon charge of the 1st anyon.
   * @param anyon2 Anyon charge of the 2nd anyon.
   * @param outcome Anyon charge of the fusion result.
   * @return True if the Fibonacci fusion rules are obeyed, False otherwise.
   */
  static bool checkRule(int anyon1, int anyon2, int outcome);

  /**
   * Checks whether the list of outcomes obeys a specific rule.
   *
   * @param outcomes List of outcomes.
   * @return True if all outcomes obey the rule, False otherwise.
   */
  static bool checkOutcomes(const std::vector<int>& outcomes);

  /**
   * Checks whether the state obeys specific rules.
   *
   * @param state Input state represented by a Basis struct.
   * @return True if the state obeys the rules, False otherwise.
   */
  static bool checkState(const State& state);

  static State genState(const std::vector<int>& comb, int nbQudits, int quditLen);

  /**
   * Generate all the basis states for a system of a given number of qudits
   * and a given number of anyons per qudit.
   *
   * @param nbQudits Number of qudits in the circuit.
   * @param nbAnyonsPerQudit Number of anyons in each qudit.
   * @return A list of basis states represented by a vector of State objects.
   */
  static Basis generateBasis(int nbQudits, int nbAnyonsPerQudit);
};

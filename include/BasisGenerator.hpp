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
  bool check_rule(int anyon1, int anyon2, int outcome);

  /**
   * Checks whether the list of outcomes obeys a specific rule.
   *
   * @param outcomes List of outcomes.
   * @return True if all outcomes obey the rule, False otherwise.
   */
  bool check_outcomes(const std::vector<int>& outcomes);

  /**
   * Checks whether the state obeys specific rules.
   *
   * @param state Input state represented by a Basis struct.
   * @return True if the state obeys the rules, False otherwise.
   */
  bool check_state(const State& state);

  State gen_state(const std::vector<int>& comb, int nb_qudits, int qudit_len);

  /**
   * Generate all the basis states for a system of a given number of qudits
   * and a given number of anyons per qudit.
   *
   * @param nb_qudits Number of qudits in the circuit.
   * @param nb_anyons_per_qudit Number of anyons in each qudit.
   * @return A list of basis states represented by a vector of State objects.
   */
  Basis generate_basis(int nb_qudits, int nb_anyons_per_qudit);
};

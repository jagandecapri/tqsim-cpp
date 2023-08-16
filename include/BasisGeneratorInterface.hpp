#pragma once

#include "container.hpp"

#include <vector>

class BasisGeneratorInterface {
private:
  virtual bool check_rule(int anyon1, int anyon2, int outcome) = 0;

  virtual bool check_outcomes(std::vector<int> outcomes) = 0;

  virtual bool check_state(const State& state) = 0;

  virtual State gen_state(const std::vector<int>& comb, int nb_qubits,
                          int qudit_len) = 0;

public:
  virtual ~BasisGeneratorInterface() = default;

  virtual Basis generate_basis(int nb_qudits, int nb_anyons_per_qudit) = 0;
};

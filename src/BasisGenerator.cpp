#include "BasisGenerator.hpp"

bool BasisGenerator::checkRule(int anyon1, int anyon2, int outcome) {
  if (anyon1 != 0 && anyon2 != 0) {
    return true;
  }

  if ((anyon1 != 0 || anyon2 != 0) && outcome == 1) {
    return true;
  }

  if (anyon1 == 0 && anyon2 == 0 && outcome == 0) {
    return true;
  }
  return false;
}

bool BasisGenerator::checkOutcomes(const std::vector<int>& outcomes) {
  int previousOutcome = 1;
  for (const auto outcome : outcomes) {
    if (checkRule(previousOutcome, 1, outcome)) {
      previousOutcome = outcome;
    } else {
      return false;
    }
  }
  return true;
}

bool BasisGenerator::checkState(const State& state) {
  std::size_t nb_qudits = state.qudits.size();
  std::size_t qudit_len = state.qudits[0].size();

  for (const std::vector<int>& qudit : state.qudits) {
    if (qudit.size() == qudit_len) {
      if (!checkOutcomes(qudit)) {
        return false;
      }
    } else {
      return false;
    }
  }

  if (nb_qudits != static_cast<int>(state.roots.size()) + 1) {
    return false;
  }

  int previous_outcome = state.qudits[0].back();

  for (size_t i = 0; i < state.roots.size(); ++i) {
    if (!checkRule(previous_outcome, state.qudits[i + 1].back(),
                   state.roots[i])) {
      return false;
    }
    previous_outcome = state.roots[i];
  }

  return true;
}

State BasisGenerator::genState(const std::vector<int>& comb, int nbQudits,
                                int quditLen) {
  State state;

  for (size_t i = 0; i < comb.size(); ++i) {
    int label = comb[i];
    if (i < (nbQudits * quditLen)) {
      if (i % quditLen) {
        state.qudits.back().push_back(label);
      } else {
        state.qudits.push_back({label});
      }
    } else {
      state.roots.push_back(label);
    }
  }

  return state;
}

Basis BasisGenerator::generateBasis(int nbQudits, int nbAnyonsPerQudit) {
  int nb_roots = nbQudits - 1;
  int qudit_len = nbAnyonsPerQudit - 1;
  int nb_labels = nbQudits * qudit_len + nb_roots;

  std::vector<State> basis;
  std::vector<int> curr_comb(nb_labels, 0);
  std::vector<int> final_comb(nb_labels, 1);

  State curr_state = genState(curr_comb, nbQudits, qudit_len);

  if (checkState(curr_state)) {
    basis.push_back(curr_state);
  }

  while (curr_comb != final_comb) {
    for (int i = 0; i < nb_labels; ++i) {
      if (curr_comb[i] == 0) {
        curr_comb[i] = 1;
        break;
      }
      curr_comb[i] = 0;
    }

    curr_state = genState(curr_comb, nbQudits, qudit_len);

    if (checkState(curr_state)) {
      basis.push_back(curr_state);
    }
  }

  return basis;
}

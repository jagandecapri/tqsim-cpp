#include "BasisGenerator.hpp"

bool BasisGenerator::check_rule(int anyon1, int anyon2, int outcome) {
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

bool BasisGenerator::check_outcomes(const std::vector<int>& outcomes) {
  int previousOutcome = 1;
  for (const auto outcome : outcomes) {
    if (check_rule(previousOutcome, 1, outcome)) {
      previousOutcome = outcome;
    } else {
      return false;
    }
  }
  return true;
}

bool BasisGenerator::check_state(const State& state) {
  std::size_t nb_qudits = state.qudits.size();
  std::size_t qudit_len = state.qudits[0].size();

  for (const std::vector<int>& qudit : state.qudits) {
    if (qudit.size() == qudit_len) {
      if (!check_outcomes(qudit)) {
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
    if (!check_rule(previous_outcome, state.qudits[i + 1].back(),
                    state.roots[i])) {
      return false;
    }
    previous_outcome = state.roots[i];
  }

  return true;
}

State BasisGenerator::gen_state(const std::vector<int>& comb, int nb_qudits,
                                int qudit_len) {
  State state;

  for (size_t i = 0; i < comb.size(); ++i) {
    int label = comb[i];
    if (i < (nb_qudits * qudit_len)) {
      if (i % qudit_len) {
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

Basis BasisGenerator::generate_basis(int nb_qudits, int nb_anyons_per_qudit) {
  int nb_roots = nb_qudits - 1;
  int qudit_len = nb_anyons_per_qudit - 1;
  int nb_labels = nb_qudits * qudit_len + nb_roots;

  std::vector<State> basis;
  std::vector<int> curr_comb(nb_labels, 0);
  std::vector<int> final_comb(nb_labels, 1);

  State curr_state = gen_state(curr_comb, nb_qudits, qudit_len);

  if (check_state(curr_state)) {
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

    curr_state = gen_state(curr_comb, nb_qudits, qudit_len);

    if (check_state(curr_state)) {
      basis.push_back(curr_state);
    }
  }

  return basis;
}

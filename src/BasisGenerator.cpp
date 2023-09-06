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
  std::size_t const nbQudits = state.qudits.size();
  std::size_t const quditLen = state.qudits[0].size();

  for (const std::vector<int>& qudit : state.qudits) {
    if (qudit.size() == quditLen) {
      if (!checkOutcomes(qudit)) {
        return false;
      }
    } else {
      return false;
    }
  }

  if (nbQudits != static_cast<int>(state.roots.size()) + 1) {
    return false;
  }

  int previousOutcome = state.qudits[0].back();

  for (size_t i = 0; i < state.roots.size(); ++i) {
    if (!checkRule(previousOutcome, state.qudits[i + 1].back(),
                   state.roots[i])) {
      return false;
    }
    previousOutcome = state.roots[i];
  }

  return true;
}

State BasisGenerator::genState(const std::vector<int>& comb, int nbQudits,
                                int quditLen) {
  State state;

  for (size_t i = 0; i < comb.size(); ++i) {
    int const label = comb[i];
    if (i < (static_cast<size_t>(nbQudits) * quditLen)) {
      if ((i % quditLen) != 0U) {
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
  int const nbRoots = nbQudits - 1;
  int const quditLen = nbAnyonsPerQudit - 1;
  int const nbLabels = nbQudits * quditLen + nbRoots;

  std::vector<State> basis;
  std::vector<int> currComb(nbLabels, 0);
  std::vector<int> const finalComb(nbLabels, 1);

  State currState = genState(currComb, nbQudits, quditLen);

  if (checkState(currState)) {
    basis.push_back(currState);
  }

  while (currComb != finalComb) {
    for (int i = 0; i < nbLabels; ++i) {
      if (currComb[i] == 0) {
        currComb[i] = 1;
        break;
      }
      currComb[i] = 0;
    }

    currState = genState(currComb, nbQudits, quditLen);

    if (checkState(currState)) {
      basis.push_back(currState);
    }
  }

  return basis;
}

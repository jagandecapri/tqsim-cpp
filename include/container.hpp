#pragma once

#include <complex>
#include <map>
#include <vector>

struct State {
  std::vector<std::vector<int>> qudits;
  std::vector<int> roots;

  // Define the comparison operator (operator==) for State
  bool operator==(const State& other) const {
    return qudits == other.qudits && roots == other.roots;
  }

  // Define the comparison operator (operator!=) for State
  bool operator!=(const State& other) const { return !(*this == other); }
};

using Basis = std::vector<State>;

using Sequence = std::vector<std::vector<int>>;

struct Result {
  std::map<int, int> counts_dict;
  std::vector<int> memory;
};

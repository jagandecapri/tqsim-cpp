#pragma once

#include "dd/Package.hpp"

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
  std::map<int, int> countsDict;
  std::vector<int> memory;
};

using ResultDD = std::map<std::string, std::size_t>;

struct BraidingOperator {
  std::unique_ptr<dd::Package<>> dd;
  dd::mEdge ddMatrix;
  dd::mEdge ddAdjointMatrix;

  // Constructor to initialize the member using the std::unique_ptr<T> parameter
  BraidingOperator(std::unique_ptr<dd::Package<>> uniqueDd, dd::mEdge a,
                   dd::mEdge b)
      : dd(std::move(uniqueDd)), ddMatrix(a), ddAdjointMatrix(b) {}
  BraidingOperator(dd::mEdge a, dd::mEdge b)
      : ddMatrix(a), ddAdjointMatrix(b) {}
};

using BraidingOperators = std::vector<BraidingOperator>;

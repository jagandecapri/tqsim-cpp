//
// Created by Jagan on 06/08/2023.
//
#include <map>
#include <complex>
#include <vector>

#ifndef TQSIM_CPP_BASIS_H
#define TQSIM_CPP_BASIS_H

struct State {
    std::vector<std::vector<int>> qudits;
    std::vector<int> roots;

    // Define the comparison operator (operator==) for State
    bool operator==(const State &other) const {
        return qudits == other.qudits && roots == other.roots;
    }

    // Define the comparison operator (operator!=) for State
    bool operator!=(const State &other) const {
        return !(*this == other);
    }
};

typedef std::vector<State> Basis;

typedef std::vector<std::vector<int>> Sequence;

struct Result {
    std::map<int, int> counts_dict;
    std::vector<int> memory;
};
#endif //TQSIM_CPP_BASIS_H

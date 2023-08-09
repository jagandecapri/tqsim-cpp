//
// Created by Jagan on 06/08/2023.
//
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
};

typedef std::vector<State> Basis;

typedef std::vector<std::vector<int>> Sequence;
#endif //TQSIM_CPP_BASIS_H

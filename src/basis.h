//
// Created by Jagan on 06/08/2023.
//
#include <vector>

#ifndef TQSIM_CPP_BASIS_H
#define TQSIM_CPP_BASIS_H
struct State {
    std::vector<std::vector<int>> qudits;
    std::vector<int> roots;
};
#endif //TQSIM_CPP_BASIS_H

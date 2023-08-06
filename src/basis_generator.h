//
// Created by Jagan on 06/08/2023.
//
#include <vector>
#include "container.h"

#ifndef TQSIM_CPP_BASIS_GENERATOR_H
#define TQSIM_CPP_BASIS_GENERATOR_H

class BasisGeneratorInterface {
    virtual bool check_rule(int, int, int) = 0;

    virtual bool check_outcomes(std::vector<int>) = 0;

    virtual bool check_state(const State &) = 0;

    virtual State gen_state(const std::vector<int> &, int, int) = 0;

    virtual Basis generate_basis(int, int) = 0;
};

#endif //TQSIM_CPP_BASIS_GENERATOR_H

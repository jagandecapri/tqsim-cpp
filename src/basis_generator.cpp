//
// Created by Jagan on 06/08/2023.
//
#include "basis_generator.h"

class BasisGenerator : public BasisGeneratorInterface {
    bool check_rule(int, int, int) override {
        return true;
    };

    bool check_outcomes(std::vector<int>) override {
        return true;
    };

    bool check_state(std::vector<std::vector<int>>) override {
        return true;
    };

    Basis gen_state(std::vector<int>, int, int) override {
        Basis basis;
        return basis;
    };

    std::vector<Basis> generate_basis(int, int) override {
        std::vector<Basis> basis;
        return basis;
    };
};
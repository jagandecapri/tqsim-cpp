//
// Created by Jagan on 06/08/2023.
//
#include "basis_generator.h"

class BasisGenerator : public BasisGeneratorInterface {
public:
    /**
     * Checks whether the Fibonacci fusion rules are obeyed for given anyon charges and outcome charge.
     *
     * @param anyon1 Anyon charge of the 1st anyon.
     * @param anyon2 Anyon charge of the 2nd anyon.
     * @param outcome Anyon charge of the fusion result.
     * @return True if the Fibonacci fusion rules are obeyed, False otherwise.
     */
    bool check_rule(int anyon1, int anyon2, int outcome) override {
        if (anyon1 && anyon2) {
            return true;
        } else if ((anyon1 || anyon2) && outcome == 1) {
            return true;
        } else if (!(anyon1 || anyon2) && outcome == 0) {
            return true;
        } else {
            return false;
        }
    }

    /**
 * Checks whether the list of outcomes obeys a specific rule.
 *
 * @param outcomes List of outcomes.
 * @return True if all outcomes obey the rule, False otherwise.
 */
    bool check_outcomes(std::vector<int> outcomes) override {
        int previousOutcome = 1;

        for (int outcome: outcomes) {
            if (this->check_rule(previousOutcome, 1, outcome)) {
                previousOutcome = outcome;
            } else {
                return false;
            }
        }

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
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

    bool check_state(const Basis &state) override {
        std::size_t nb_qudits = state.qudits.size();
        std::size_t qudit_len = state.qudits[0].size();

        for (const std::vector<int> &qudit: state.qudits) {
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
            if (!check_rule(previous_outcome, state.qudits[i + 1].back(), state.roots[i])) {
                return false;
            }
            previous_outcome = state.roots[i];
        }

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
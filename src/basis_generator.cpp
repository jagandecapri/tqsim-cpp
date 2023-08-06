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

    /**
 * Checks whether the state obeys specific rules.
 *
 * @param state Input state represented by a Basis struct.
 * @return True if the state obeys the rules, False otherwise.
 */
    bool check_state(const State &state) override {
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

    State gen_state(const std::vector<int> &comb, int nb_qudits, int qudit_len) override {
        State state;

        for (size_t i = 0; i < comb.size(); ++i) {
            int label = comb[i];
            if (i < (nb_qudits * qudit_len)) {
                if (i % qudit_len) {
                    state.qudits.back().push_back(label);
                } else {
                    state.qudits.push_back({label});
                }
            } else {
                state.roots.push_back(label);
            }
        }

        return state;

    };

    /**
 * Generate all the basis states for a system of a given number of qudits
 * and a given number of anyons per qudit.
 *
 * @param nb_qudits Number of qudits in the circuit.
 * @param nb_anyons_per_qudit Number of anyons in each qudit.
 * @return A list of basis states represented by a vector of State objects.
 */
    Basis generate_basis(int nb_qudits, int nb_anyons_per_qudit) override {
        int nb_roots = nb_qudits - 1;
        int qudit_len = nb_anyons_per_qudit - 1;
        int nb_labels = nb_qudits * qudit_len + nb_roots;

        std::vector<State> basis;
        std::vector<int> curr_comb(nb_labels, 0);
        std::vector<int> final_comb(nb_labels, 1);

        State curr_state = this->gen_state(curr_comb, nb_qudits, qudit_len);

        if (this->check_state(curr_state)) {
            basis.push_back(curr_state);
        }

        while (curr_comb != final_comb) {
            for (int i = 0; i < nb_labels; ++i) {
                if (curr_comb[i] == 0) {
                    curr_comb[i] = 1;
                    break;
                } else {
                    curr_comb[i] = 0;
                }
            }

            curr_state = this->gen_state(curr_comb, nb_qudits, qudit_len);

            if (this->check_state(curr_state)) {
                basis.push_back(curr_state);
            }
        }

        return basis;
    };
};
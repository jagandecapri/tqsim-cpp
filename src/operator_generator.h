//
// Created by Jagan on 05/08/2023.
//
#include <cmath>
#include <Eigen/Dense>

#ifndef TQSIM_CPP_OPERATOR_GENERATOR_H
#define TQSIM_CPP_OPERATOR_GENERATOR_H

/**
 * @brief Calculates the F matrix based on input values.
 *
 * The F matrix is a 2x2 matrix that depends on four integer parameters (a1, a2, a3, and outcome).
 * The values of these parameters determine the values in the F matrix as per the provided conditions.
 *
 * @param a1 Integer parameter a1.
 * @param a2 Integer parameter a2.
 * @param a3 Integer parameter a3.
 * @param outcome Integer parameter outcome.
 * @return Eigen::Matrix<double, 2, 2> The F matrix.
 */
Eigen::Matrix<double, 2, 2> F(int a1, int a2, int a3, int outcome) {
    const double inv_phi = (std::sqrt(5) - 1) / 2;  // inverse of golden number
    Eigen::Matrix<double, 2, 2> f_matrix;

    // a1 + a2 + a3 + outcome = 4
    if (a1 + a2 + a3 + outcome == 4) {
        f_matrix << inv_phi, std::sqrt(inv_phi),
                std::sqrt(inv_phi), -inv_phi;
    }
        // a1 + a2 + a3 + outcome = 3
    else if (a1 + a2 + a3 + outcome == 3) {
        f_matrix << 0, 0,
                0, 1;
    }
        // a1 + a2 + a3 + outcome = 2
    else if (a1 + a2 + a3 + outcome == 2) {
        if (a1 + a2 == 2) {
            f_matrix << 0, 1,
                    0, 0;
        } else if (a2 + a3 == 2) {
            f_matrix << 0, 0,
                    1, 0;
        } else if (a1 + a3 == 2) {
            f_matrix << 0, 0,
                    0, 1;
        } else if (a3 + outcome == 2) {
            f_matrix << 0, 1,
                    0, 0;
        } else if (a1 + outcome == 2) {
            f_matrix << 0, 0,
                    1, 0;
        } else if (a2 + outcome == 2) {
            f_matrix << 0, 0,
                    0, 1;
        }
    }
        // a1 + a2 + a3 + outcome = 1 or 0
    else if (a1 + a2 + a3 + outcome <= 1) {
        f_matrix << 1, 0,
                0, 0;
    }

    return f_matrix;
}
#endif //TQSIM_CPP_OPERATOR_GENERATOR_H

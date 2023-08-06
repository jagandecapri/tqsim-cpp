//
// Created by Jagan on 05/08/2023.
//
#include <complex>
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

/**
 * @brief Calculates the R matrix based on input values.
 *
 * The R matrix is a 2x2 matrix that depends on two integer parameters (a1 and a2).
 * The values of these parameters determine the values in the R matrix as per the provided conditions.
 *
 * @param a1 Integer parameter a1.
 * @param a2 Integer parameter a2.
 * @return Eigen::Matrix<std::complex<double>, 2, 2> The R matrix.
 */
Eigen::Matrix<std::complex<double>, 2, 2> R(int a1, int a2) {
    using complex_t = std::complex<double>;
    Eigen::Matrix<complex_t, 2, 2> r_matrix;

    if (a1 + a2 == 2) {
        const double pi = std::acos(-1);
        complex_t exp_val1 = std::exp(complex_t(0.0, -4 * pi / 5.0));
        complex_t exp_val2 = std::exp(complex_t(0.0, 3 * pi / 5.0));
        r_matrix << exp_val1, complex_t(0, 0),
                complex_t(0, 0), exp_val2;
    } else {
        r_matrix << complex_t(1, 0), complex_t(0, 0),
                complex_t(0, 0), complex_t(1, 0);
    }

    return r_matrix;
}

// Define the B function (equivalent to Python implementation)
Eigen::Matrix<std::complex<double>, 2, 2> B(int a0, int a1,
                                            int a2, int outcome) {
    /*
    Braid matrix

    Parameters:
        a0: Matrix representing a0
        a1: Matrix representing a1
        a2: Matrix representing a2
        outcome: Matrix representing outcome

    Returns:
        b_matrix: The braid matrix computed using Eigen.
    */

    // Compute the intermediate matrices
    Eigen::Matrix<double, 2, 2> f_result = F(a0, a1, a2, outcome);
    Eigen::Matrix<std::complex<double>, 2, 2> r_result = R(a1, a2);
    Eigen::Matrix<double, 2, 2> f_transpose_conjugate_result = F(a0, a2, a1, outcome).transpose().conjugate();

    // Compute the braid matrix using Eigen matrix operations
    Eigen::Matrix<std::complex<double>, 2, 2> b_matrix = f_result * r_result * f_transpose_conjugate_result;

    return b_matrix;
}

#endif //TQSIM_CPP_OPERATOR_GENERATOR_H

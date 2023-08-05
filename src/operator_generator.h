//
// Created by Jagan on 05/08/2023.
//
#include <cmath>
#include <Eigen/Dense>

#ifndef TQSIM_CPP_OPERATOR_GENERATOR_H
#define TQSIM_CPP_OPERATOR_GENERATOR_H

Eigen::Matrix<double, 2, 2> F(int a1, int a2, int a3, int outcome) {
    /**
    F matrix
    **/
    //inverse of golden number
    float inv_phi = (sqrt(5) - 1) / 2;
    Eigen::Matrix<double, 2, 2> f_matrix;

    // Initialize all elements to zero
    f_matrix.setZero();


    // a1 + a2 + a3 + outcome = 4
    if (a1 + a2 + a3 + outcome == 4) {
        f_matrix << inv_phi, sqrt(inv_phi),
                    sqrt(inv_phi), -inv_phi;
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
        }
        else if (a2 + a3 == 2) {
            f_matrix << 0, 0,
                        1, 0;
        }
        else if (a1 + a3 == 2){
            f_matrix << 0, 0,
                        0, 1;
        }
        else if (a3 + outcome == 2) {
            f_matrix << 0, 1,
                        0, 0;
        }
        else if (a1 + outcome == 2) {
            f_matrix << 0, 0,
                        1, 0;
        }
        else if (a2 + outcome == 2) {
            f_matrix << 0, 0,
                        0, 1;
        }
    // a1 + a2 + a3 + outcome = 1
    // a1 + a2 + a3 + outcome = 0
    }
    else if (a1 + a2 + a3 + outcome == 0) {
        f_matrix << 1, 0,
                    0, 0;
    }

    return f_matrix;
}
#endif //TQSIM_CPP_OPERATOR_GENERATOR_H

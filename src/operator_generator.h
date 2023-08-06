//
// Created by Jagan on 05/08/2023.
//
#include <Eigen/Dense>

#ifndef TQSIM_CPP_OPERATOR_GENERATOR_H
#define TQSIM_CPP_OPERATOR_GENERATOR_H

class OperatorGeneratorInterface {
public:
    virtual Eigen::Matrix<double, 2, 2> F(int, int, int, int) = 0;

    virtual Eigen::Matrix<std::complex<double>, 2, 2> R(int, int) = 0;

    virtual Eigen::Matrix<std::complex<double>, 2, 2> B(int, int, int, int) = 0;

    virtual std::complex<double> sigma(int, const std::vector<int> &, const std::vector<int> &) = 0;

    virtual std::complex<double> L(int, int, int, int, const std::vector<int> &, const std::vector<int> &) = 0;

    virtual std::complex<double> S(int, int, int, int, int, int, int, const std::vector<int> &,
                                   const std::vector<int> &) = 0;

    virtual std::complex<double> gen_sigma(int,
                                           const State &,
                                           const State &) = 0;
};

#endif //TQSIM_CPP_OPERATOR_GENERATOR_H

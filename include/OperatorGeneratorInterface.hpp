#pragma once

#include <Eigen/Dense>
#include <complex>
#include <vector>

struct State;

class OperatorGeneratorInterface {
private:
  virtual Eigen::Matrix<double, 2, 2> F(int, int, int, int) = 0;

  virtual Eigen::Matrix<std::complex<double>, 2, 2> R(int, int) = 0;

  virtual Eigen::Matrix<std::complex<double>, 2, 2> B(int, int, int, int) = 0;

  virtual std::complex<double> sigma(int, const std::vector<int>&,
                                     const std::vector<int>&) = 0;

  virtual std::complex<double> L(int, int, int, int, const std::vector<int>&,
                                 const std::vector<int>&) = 0;

  virtual std::complex<double> S(int, int, int, int, int, int, int,
                                 const std::vector<int>&,
                                 const std::vector<int>&) = 0;

  virtual std::complex<double> gen_sigma(int, const State&, const State&) = 0;

public:
  virtual ~OperatorGeneratorInterface() = default;

  virtual Eigen::MatrixXcd generate_braiding_operator(int, int, int) = 0;
};

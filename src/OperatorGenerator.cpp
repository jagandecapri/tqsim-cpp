#include "OperatorGenerator.hpp"

#include "container.hpp"

Eigen::Matrix<double, 2, 2> OperatorGenerator::F(int a1, int a2, int a3,
                                                 int outcome) {
  const double inv_phi = (std::sqrt(5) - 1) / 2; // inverse of golden number
  Eigen::Matrix<double, 2, 2> f_matrix;

  // a1 + a2 + a3 + outcome = 4
  if (a1 + a2 + a3 + outcome == 4) {
    f_matrix << inv_phi, std::sqrt(inv_phi), std::sqrt(inv_phi), -inv_phi;
  }
  // a1 + a2 + a3 + outcome = 3
  else if (a1 + a2 + a3 + outcome == 3) {
    f_matrix << 0, 0, 0, 1;
  }
  // a1 + a2 + a3 + outcome = 2
  else if (a1 + a2 + a3 + outcome == 2) {
    if (a1 + a2 == 2) {
      f_matrix << 0, 1, 0, 0;
    } else if (a2 + a3 == 2) {
      f_matrix << 0, 0, 1, 0;
    } else if (a1 + a3 == 2) {
      f_matrix << 0, 0, 0, 1;
    } else if (a3 + outcome == 2) {
      f_matrix << 0, 1, 0, 0;
    } else if (a1 + outcome == 2) {
      f_matrix << 0, 0, 1, 0;
    } else if (a2 + outcome == 2) {
      f_matrix << 0, 0, 0, 1;
    }
  }
  // a1 + a2 + a3 + outcome = 1 or 0
  else if (a1 + a2 + a3 + outcome <= 1) {
    f_matrix << 1, 0, 0, 0;
  }

  return f_matrix;
}

Eigen::Matrix<std::complex<double>, 2, 2> OperatorGenerator::R(int a1, int a2) {
  using complex_t = std::complex<double>;
  Eigen::Matrix<complex_t, 2, 2> r_matrix;

  if (a1 + a2 == 2) {
    const double pi = std::acos(-1);
    complex_t exp_val1 = std::exp(complex_t(0.0, -4 * pi / 5.0));
    complex_t exp_val2 = std::exp(complex_t(0.0, 3 * pi / 5.0));
    r_matrix << exp_val1, complex_t(0, 0), complex_t(0, 0), exp_val2;
  } else {
    r_matrix << complex_t(1, 0), complex_t(0, 0), complex_t(0, 0),
        complex_t(1, 0);
  }

  return r_matrix;
}

Eigen::Matrix<std::complex<double>, 2, 2>
OperatorGenerator::B(int a0, int a1, int a2, int outcome) {
  // Compute the intermediate matrices
  Eigen::Matrix<double, 2, 2> f_result = this->F(a0, a1, a2, outcome);
  Eigen::Matrix<std::complex<double>, 2, 2> r_result = this->R(a1, a2);
  Eigen::Matrix<double, 2, 2> f_transpose_conjugate_result =
      this->F(a0, a2, a1, outcome).transpose().conjugate();

  // Compute the braid matrix using Eigen matrix operations
  Eigen::Matrix<std::complex<double>, 2, 2> b_matrix =
      f_result * r_result * f_transpose_conjugate_result;

  return b_matrix;
}

std::complex<double> OperatorGenerator::sigma(int index,
                                              const std::vector<int>& state_f,
                                              const std::vector<int>& state_i) {
  if (index <= 0 || index > static_cast<int>(state_i.size())) {
    throw std::invalid_argument("index value is not valid!");
  }

  std::vector<int> stt_f = {1};
  stt_f.insert(stt_f.end(), state_f.begin(), state_f.end());

  std::vector<int> stt_i = {1};
  stt_i.insert(stt_i.end(), state_i.begin(), state_i.end());

  int a0;
  if (index - 2 < 0) {
    a0 = 0;
  } else if (index - 2 == 0) {
    a0 = 1;
  } else {
    a0 = state_i[index - 3];
  }

  int outcome = state_i[index - 1];
  int a = stt_i[index - 1];
  int b = stt_f[index - 1];
  Eigen::Matrix<std::complex<double>, 2, 2> b_matrix =
      this->B(a0, 1, 1, outcome);
  std::complex<double> amplitude = b_matrix(a, b);

  std::vector<int> ket = stt_i;
  ket[index - 1] = b;
  std::vector<int> bra = stt_f;
  double braket = (ket == bra) ? 1 : 0;

  return amplitude * braket;
}

std::complex<double> OperatorGenerator::L(int k, int h, int i_, int i,
                                          const std::vector<int>& jj_,
                                          const std::vector<int>& jj) {
  std::complex<double> component(0.0, 0.0);

  int qudit_len = static_cast<int>(jj.size());
  std::vector<int> jjj_ = jj_;
  std::vector<int> jjj = jj;
  jjj_.insert(jjj_.begin(), 1);
  jjj.insert(jjj.begin(), 1);

  std::vector<int> init_p(qudit_len, 0);
  std::vector<int> final_p(qudit_len, 1);
  std::vector<int> new_p = init_p;

  while (new_p != final_p) {
    std::vector<int> pp = new_p;
    pp.push_back(k);

    std::complex<double> product(1.0, 0.0);
    for (int ii = 0; ii < qudit_len; ++ii) {
      product *= F(i, jjj[ii], 1, pp[ii + 1])
                     .transpose()
                     .conjugate()(jjj[ii + 1], pp[ii]) *
                 F(i_, jjj_[ii], 1, pp[ii + 1])(pp[ii], jjj_[ii + 1]);
    }

    product *= B(h, 1, 1, pp[0])(i, i_);
    component += product;

    // iterate
    for (int ii = 0; ii < qudit_len; ++ii) {
      if (new_p[ii] == 0) {
        new_p[ii] = 1;
        break;
      } else {
        new_p[ii] = 0;
      }
    }
  }

  // final iteration
  std::vector<int> pp = new_p;
  pp.push_back(k);
  std::complex<double> product(1.0, 0.0);
  for (int ii = 0; ii < qudit_len; ++ii) {
    product *= F(i, jjj[ii], 1, pp[ii + 1])
                   .transpose()
                   .conjugate()(jjj[ii + 1], pp[ii]) *
               F(i_, jjj_[ii], 1, pp[ii + 1])(pp[ii], jjj_[ii + 1]);
  }

  product *= B(h, 1, 1, pp[0])(i, i_);
  component += product;

  return component;
}

std::complex<double> OperatorGenerator::S(int jm, int jmo, int jmoo, int jmo_,
                                          int h, int i_, int i,
                                          const std::vector<int>& jj_,
                                          const std::vector<int>& jj) {
  std::complex<double> component(0.0, 0.0);

  for (int kk : {0, 1}) {
    component += F(jmoo, i, jj.back(), jm)(jmo, kk) * L(kk, h, i_, i, jj_, jj) *
                 F(jmoo, i_, jj_.back(), jm).transpose().conjugate()(kk, jmo_);
  }

  return component;
}

std::complex<double> OperatorGenerator::gen_sigma(int index,
                                                  const State& state_i,
                                                  const State& state_f) {
  int qudit_len = static_cast<int>(state_i.qudits[0].size());
  int nb_anyons_per_qudit = qudit_len + 1;

  std::complex<double> amplitude(0.0, 0.0);
  std::complex<double> braket(1.0, 0.0);

  if (index % nb_anyons_per_qudit > 0) {
    int m = index / nb_anyons_per_qudit;
    int idx = index % nb_anyons_per_qudit;
    amplitude = sigma(idx, state_f.qudits.at(m), state_i.qudits.at(m));

    for (size_t i = 0; i < state_i.qudits.size(); ++i) {
      if (i == static_cast<size_t>(m)) {
        continue;
      } else if (state_i.qudits[i] != state_f.qudits[i]) {
        braket = {0.0, 0.0};
      }
    }

    for (size_t i = 0; i < state_i.roots.size(); ++i) {
      if (state_i.roots[i] != state_f.roots[i]) {
        braket = {0.0, 0.0};
      }
    }

  } else {
    int m = (index / nb_anyons_per_qudit) - 1;

    State new_state_i = state_i;
    new_state_i.qudits.at(m).back() = state_f.qudits.at(m).back();
    new_state_i.qudits.at(m + 1) = state_f.qudits.at(m + 1);

    int jm = 0, jmo = 0, jmoo = 0, jmo_ = 0, h = 0, i_ = 0, i = 0;
    std::vector<int> jj_, jj;

    if (m + 1 > 2) {
      new_state_i.roots[m - 1] = state_f.roots[m - 1];

      jj_ = new_state_i.qudits[m + 1];
      jj = state_i.qudits[m + 1];
      h = state_i.qudits.at(m).at(qudit_len - 2);
      i = state_i.qudits.at(m).at(qudit_len - 1);
      i_ = new_state_i.qudits.at(m).at(qudit_len - 1);

      jmo_ = new_state_i.roots[m - 1];
      jmoo = state_i.roots[m - 2];
      jmo = state_i.roots[m - 1];
      jm = state_i.roots[m];

    } else if (m + 1 == 2) {
      new_state_i.roots[m - 1] = state_f.roots[m - 1];

      jj_ = new_state_i.qudits[m + 1];
      jj = state_i.qudits[m + 1];
      h = state_i.qudits.at(m).at(qudit_len - 2);
      i = state_i.qudits.at(m).at(qudit_len - 1);
      i_ = new_state_i.qudits.at(m).at(qudit_len - 1);

      jmo_ = new_state_i.roots[m - 1];
      jmoo = state_i.qudits.at(0).at(qudit_len - 1);
      jmo = state_i.roots[m - 1];
      jm = state_i.roots[m];

    } else if (m + 1 == 1) {

      jj_ = new_state_i.qudits[m + 1];
      jj = state_i.qudits[m + 1];
      h = state_i.qudits.at(m).at(qudit_len - 2);
      i = state_i.qudits.at(m).at(qudit_len - 1);
      i_ = new_state_i.qudits.at(m).at(qudit_len - 1);

      jmo_ = new_state_i.qudits.at(0).at(qudit_len - 1);
      jmoo = 0;
      jmo = state_i.qudits.at(0).at(qudit_len - 1);
      jm = state_i.roots[m];
    }

    amplitude += S(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj);

    if (new_state_i != state_f) {
      braket = {0.0, 0.0};
    }
  }

  return braket * amplitude;
}

Eigen::MatrixXcd
OperatorGenerator::generate_braiding_operator(int index, int nb_qudits,
                                              int nb_anyons_per_qudit) {
  // Generate the basis states
  Basis basis =
      this->basis_generator->generate_basis(nb_qudits, nb_anyons_per_qudit);

  size_t basisSize = basis.size();
  Eigen::MatrixXcd sigma(basisSize, basisSize);

  for (size_t f = 0; f < basisSize; ++f) {
    auto f_index = static_cast<Eigen::Index>(f);
    for (size_t base_i = 0; base_i < basisSize; ++base_i) {
      auto base_i_index = static_cast<Eigen::Index>(base_i);
      sigma(f_index, base_i_index) = gen_sigma(index, basis[base_i], basis[f]);
    }
  }

  return sigma;
}

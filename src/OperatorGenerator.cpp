#include "OperatorGenerator.hpp"

#include "container.hpp"

Eigen::Matrix<double, 2, 2> OperatorGenerator::F(int a1, int a2, int a3,
                                                 int outcome) {
  const double invPhi = (std::sqrt(5) - 1) / 2; // inverse of golden number
  Eigen::Matrix<double, 2, 2> fMatrix;

  // a1 + a2 + a3 + outcome = 4
  if (a1 + a2 + a3 + outcome == 4) {
    fMatrix << invPhi, std::sqrt(invPhi), std::sqrt(invPhi), -invPhi;
  }
  // a1 + a2 + a3 + outcome = 3
  else if (a1 + a2 + a3 + outcome == 3) {
    fMatrix << 0, 0, 0, 1;
  }
  // a1 + a2 + a3 + outcome = 2
  else if (a1 + a2 + a3 + outcome == 2) {
    if (a1 + a2 == 2) {
      fMatrix << 0, 1, 0, 0;
    } else if (a2 + a3 == 2) {
      fMatrix << 0, 0, 1, 0;
    } else if (a1 + a3 == 2) {
      fMatrix << 0, 0, 0, 1;
    } else if (a3 + outcome == 2) {
      fMatrix << 0, 1, 0, 0;
    } else if (a1 + outcome == 2) {
      fMatrix << 0, 0, 1, 0;
    } else if (a2 + outcome == 2) {
      fMatrix << 0, 0, 0, 1;
    }
  }
  // a1 + a2 + a3 + outcome = 1 or 0
  else if (a1 + a2 + a3 + outcome <= 1) {
    fMatrix << 1, 0, 0, 0;
  }

  return fMatrix;
}

Eigen::Matrix<std::complex<double>, 2, 2> OperatorGenerator::R(int a1, int a2) {
  using complex_t = std::complex<double>;
  Eigen::Matrix<complex_t, 2, 2> rMatrix;

  if (a1 + a2 == 2) {
    const double pi = std::acos(-1);
    complex_t const expVal1 = std::exp(complex_t(0.0, -4 * pi / 5.0));
    complex_t const expVal2 = std::exp(complex_t(0.0, 3 * pi / 5.0));
    rMatrix << expVal1, complex_t(0, 0), complex_t(0, 0), expVal2;
  } else {
    rMatrix << complex_t(1, 0), complex_t(0, 0), complex_t(0, 0),
        complex_t(1, 0);
  }

  return rMatrix;
}

Eigen::Matrix<std::complex<double>, 2, 2>
OperatorGenerator::B(int a0, int a1, int a2, int outcome) {
  // Compute the intermediate matrices
  Eigen::Matrix<double, 2, 2> const fResult = F(a0, a1, a2, outcome);
  Eigen::Matrix<std::complex<double>, 2, 2> const rResult = R(a1, a2);
  // TODO: Can reuse previous computation of Result?
  Eigen::Matrix<double, 2, 2> const fTransposeConjugateResult =
      F(a0, a2, a1, outcome).transpose().conjugate();

  // Compute the braid matrix using Eigen matrix operations
  Eigen::Matrix<std::complex<double>, 2, 2> bMatrix =
      fResult * rResult * fTransposeConjugateResult;

  return bMatrix;
}

std::complex<double> OperatorGenerator::sigma(int index,
                                              const std::vector<int>& stateF,
                                              const std::vector<int>& stateI) {
  if (index <= 0 || index > static_cast<int>(stateI.size())) {
    throw std::invalid_argument("index value is not valid!");
  }

  std::vector<int> sttF = {1};
  sttF.insert(sttF.end(), stateF.begin(), stateF.end());

  std::vector<int> sttI = {1};
  sttI.insert(sttI.end(), stateI.begin(), stateI.end());

  int a0 = 0;
  if (index - 2 < 0) {
    a0 = 0;
  } else if (index - 2 == 0) {
    a0 = 1;
  } else {
    a0 = stateI[index - 3];
  }

  int const outcome = stateI[index - 1];
  int const a = sttI[index - 1];
  int const b = sttF[index - 1];
  Eigen::Matrix<std::complex<double>, 2, 2> bMatrix =
      B(a0, 1, 1, outcome);
  std::complex<double> const amplitude = bMatrix(a, b);

  std::vector<int> ket = sttI;
  ket[index - 1] = b;
  std::vector<int> const bra = sttF;
  double const braket = (ket == bra) ? 1 : 0;

  return amplitude * braket;
}

std::complex<double> OperatorGenerator::L(int k, int h, int i_, int i,
                                          const std::vector<int>& jj_,
                                          const std::vector<int>& jj) {
  std::complex<double> component(0.0, 0.0);

  int const quditLen = static_cast<int>(jj.size());
  std::vector<int> jjj_ = jj_;
  std::vector<int> jjj = jj;
  jjj_.insert(jjj_.begin(), 1);
  jjj.insert(jjj.begin(), 1);

  std::vector<int> const initP(quditLen, 0);
  std::vector<int> const finalP(quditLen, 1);
  std::vector<int> newP = initP;

  while (newP != finalP) {
    std::vector<int> pp = newP;
    pp.push_back(k);

    std::complex<double> product(1.0, 0.0);
    for (int ii = 0; ii < quditLen; ++ii) {
      product *= F(i, jjj[ii], 1, pp[ii + 1])
                     .transpose()
                     .conjugate()(jjj[ii + 1], pp[ii]) *
                 F(i_, jjj_[ii], 1, pp[ii + 1])(pp[ii], jjj_[ii + 1]);
    }

    product *= B(h, 1, 1, pp[0])(i, i_);
    component += product;

    // iterate
    for (int ii = 0; ii < quditLen; ++ii) {
      if (newP[ii] == 0) {
        newP[ii] = 1;
        break;
      }
      newP[ii] = 0;
    }
  }

  // final iteration
  std::vector<int> pp = newP;
  pp.push_back(k);
  std::complex<double> product(1.0, 0.0);
  for (int ii = 0; ii < quditLen; ++ii) {
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

  for (int const kk : {0, 1}) {
    component += F(jmoo, i, jj.back(), jm)(jmo, kk) * L(kk, h, i_, i, jj_, jj) *
                 F(jmoo, i_, jj_.back(), jm).transpose().conjugate()(kk, jmo_);
  }

  return component;
}

std::complex<double> OperatorGenerator::genSigma(int index,
                                                  const State& stateI,
                                                  const State& stateF) {
  int const quditLen = static_cast<int>(stateI.qudits[0].size());
  int const nbAnyonsPerQudit = quditLen + 1;

  std::complex<double> amplitude(0.0, 0.0);
  std::complex<double> braket(1.0, 0.0);

  if (index % nbAnyonsPerQudit > 0) {
    int const m = index / nbAnyonsPerQudit;
    int const idx = index % nbAnyonsPerQudit;
    amplitude = sigma(idx, stateF.qudits.at(m), stateI.qudits.at(m));

    for (size_t i = 0; i < stateI.qudits.size(); ++i) {
      if (i == static_cast<size_t>(m)) {
        continue;
      }
      if (stateI.qudits[i] != stateF.qudits[i]) {
        braket = {0.0, 0.0};
      }
    }

    for (size_t i = 0; i < stateI.roots.size(); ++i) {
      if (stateI.roots[i] != stateF.roots[i]) {
        braket = {0.0, 0.0};
      }
    }

  } else {
    int const m = (index / nbAnyonsPerQudit) - 1;

    State newStateI = stateI;
    newStateI.qudits.at(m).back() = stateF.qudits.at(m).back();
    newStateI.qudits.at(m + 1) = stateF.qudits.at(m + 1);

    int jm = 0;
    int jmo = 0;
    int jmoo = 0;
    int jmo_ = 0;
    int h = 0;
    int i_ = 0;
    int i = 0;
    std::vector<int> jj_;
    std::vector<int> jj;

    if (m + 1 > 2) {
      newStateI.roots[m - 1] = stateF.roots[m - 1];

      jj_ = newStateI.qudits[m + 1];
      jj = stateI.qudits[m + 1];
      h = stateI.qudits.at(m).at(quditLen - 2);
      i = stateI.qudits.at(m).at(quditLen - 1);
      i_ = newStateI.qudits.at(m).at(quditLen - 1);

      jmo_ = newStateI.roots[m - 1];
      jmoo = stateI.roots[m - 2];
      jmo = stateI.roots[m - 1];
      jm = stateI.roots[m];

    } else if (m + 1 == 2) {
      newStateI.roots[m - 1] = stateF.roots[m - 1];

      jj_ = newStateI.qudits[m + 1];
      jj = stateI.qudits[m + 1];
      h = stateI.qudits.at(m).at(quditLen - 2);
      i = stateI.qudits.at(m).at(quditLen - 1);
      i_ = newStateI.qudits.at(m).at(quditLen - 1);

      jmo_ = newStateI.roots[m - 1];
      jmoo = stateI.qudits.at(0).at(quditLen - 1);
      jmo = stateI.roots[m - 1];
      jm = stateI.roots[m];

    } else if (m + 1 == 1) {

      jj_ = newStateI.qudits[m + 1];
      jj = stateI.qudits[m + 1];
      h = stateI.qudits.at(m).at(quditLen - 2);
      i = stateI.qudits.at(m).at(quditLen - 1);
      i_ = newStateI.qudits.at(m).at(quditLen - 1);

      jmo_ = newStateI.qudits.at(0).at(quditLen - 1);
      jmoo = 0;
      jmo = stateI.qudits.at(0).at(quditLen - 1);
      jm = stateI.roots[m];
    }

    amplitude += S(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj);

    if (newStateI != stateF) {
      braket = {0.0, 0.0};
    }
  }

  return braket * amplitude;
}

Eigen::MatrixXcd
OperatorGenerator::generateBraidingOperator(int index, int nbQudits,
                                              int nbAnyonsPerQudit) {
  // Generate the basis states
  Basis basis = basisGenerator->generateBasis(nbQudits, nbAnyonsPerQudit);

  size_t const basisSize = basis.size();
  Eigen::MatrixXcd sigma(basisSize, basisSize);

  for (size_t f = 0; f < basisSize; ++f) {
    auto fIndex = static_cast<Eigen::Index>(f);
    for (size_t baseI = 0; baseI < basisSize; ++baseI) {
      auto baseIIndex = static_cast<Eigen::Index>(baseI);
      sigma(fIndex, baseIIndex) = genSigma(index, basis[baseI], basis[f]);
    }
  }

  return sigma;
}

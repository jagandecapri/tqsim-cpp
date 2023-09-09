#include "DDCircuit.hpp"

#include "dd/DDDefinitions.hpp"

#include <cmath>
#include <memory>
#include <vector>

DDCircuit::DDCircuit(int nbQudits, int nbAnyonsPerQudit)
    : basisGenerator(std::make_unique<BasisGenerator>()), nbQudits(nbQudits),
      nbAnyonsPerQudit(nbAnyonsPerQudit), nbAnyons(nbQudits * nbAnyonsPerQudit),
      currentState() {
  circuitDD = std::make_unique<dd::Package<>>(nbQudits);
  operatorGenerator = std::make_unique<OperatorGenerator>(basisGenerator.get());

  generateBasis();
  dim = static_cast<size_t>(pow(2, nbQudits)); // NOLINT

  getSigmas();
}

void DDCircuit::generateBasis() {
  basis = basisGenerator->generateBasis(nbQudits, nbAnyonsPerQudit);
}

void DDCircuit::getSigmas() {
  sigmas.reserve(static_cast<size_t>(nbAnyons - 1));
  braidingOperators.reserve(static_cast<size_t>(nbAnyons - 1));
  circuitDD = std::make_unique<dd::Package<>>(nbQudits);

  size_t const startRow = 1; // First row
  auto const endRow = static_cast<size_t>(pow(2, nbQudits));
  size_t const startCol = 1; // First column
  auto const endCol = static_cast<size_t>(pow(2, nbQudits));
  //  std::cout << "Operator matrix endRow: " << endRow << '\n';
  //  std::cout << "Operator matrix endCol: " << endCol << '\n';

  for (int index = 1; index < nbAnyons; ++index) {
    const Eigen::MatrixXcd sigma = operatorGenerator->generateBraidingOperator(
        index, nbQudits, nbAnyonsPerQudit);
    const Eigen::MatrixXcd complexMatrix =
        sigma
            .block(startRow, startCol, endRow - startRow + 1,
                   endCol - startCol + 1)
            .cast<std::complex<double>>();

    dd::CMat sigmaMatrix(static_cast<size_t>(complexMatrix.rows()),
                         dd::CVec(static_cast<size_t>(complexMatrix.cols())));

    for (int i = 0; i < complexMatrix.rows(); ++i) {
      for (int j = 0; j < complexMatrix.cols(); ++j) {
        auto row = static_cast<size_t>(i);
        auto col = static_cast<size_t>(j);
        sigmaMatrix[row][col] = complexMatrix(i, j);
      }
    }

    auto matDD = circuitDD->makeDDFromMatrix(sigmaMatrix);
    auto matAdjointDD = circuitDD->conjugateTranspose(matDD);
    BraidingOperator braidingOperator{matDD, matAdjointDD};
    braidingOperators.push_back(std::move(braidingOperator));

    // TODO: Not necessary to store sigma since DD are used.
    sigmas.push_back(sigma);
  }

  assert(braidingOperators[2].ddMatrix.p->e[0].w.r != nullptr);
  for (int index = 1; index < nbAnyons; ++index) {
    assert(braidingOperators[static_cast<size_t>(index - 1)]
               .ddMatrix.p->e[0]
               .w.r != nullptr);
  }
}

void DDCircuit::initialize(
    const std::vector<std::complex<dd::fp>>& inputState) {
  if (nbBraids > 0) {
    throw InitializationException("Initialization should happen before any "
                                  "braiding operation is performed!");
  }

  if (inputState.size() != dim) {
    throw std::invalid_argument("The state has wrong dimension. Should be " +
                                std::to_string(dim));
  }

  // TODO: Check whether statevector is normalized
  //  double const norm = inputState.dot(inputState.conjugate()).real();
  //  if (Eigen::NumTraits<double>::epsilon() <= std::abs(norm - 1)) {
  //    throw std::invalid_argument("The input state is not normalized
  //    correctly!");
  //  }
  auto dd = std::make_unique<dd::Package<>>(nbQudits);
  currentState = dd->makeStateFromVector(inputState);
}

void DDCircuit::braid(int n, int m) {
  if (m < 1 || n < 1) {
    throw std::invalid_argument("n, m must be higher than 0!");
  }

  if (m > nbAnyons || n > nbAnyons) {
    throw std::invalid_argument("The system has only " +
                                std::to_string(nbAnyons) +
                                " anyons! n, m are erroneous!");
  }

  if (std::abs(n - m) != 1) {
    throw std::runtime_error("You can only braid adjacent anyons!");
  }

  circuitDD->incRef(currentState);

  int index = 0;
  dd::Edge<dd::vNode> e{};
  if (n < m) {
    // auto dd = sigmas[n - 1] * unitary;
    index = n - 1;
    assert(braidingOperators[static_cast<size_t>(index)].ddMatrix.p->e[0].w.r !=
           nullptr);
    e = circuitDD->multiply(
        braidingOperators[static_cast<size_t>(index)].ddMatrix, currentState);
  } else {
    // unitary = sigmas[m - 1].adjoint() * unitary;
    index = m - 1;
    // TODO: Need to adjoint the matrix for the second case
    assert(braidingOperators[static_cast<size_t>(index)]
               .ddAdjointMatrix.p->e[0]
               .w.r != nullptr);
    e = circuitDD->multiply(
        braidingOperators[static_cast<size_t>(index)].ddAdjointMatrix,
        currentState);
  }
  circuitDD->decRef(currentState);
  currentState = e;
  circuitDD->garbageCollect();

  braidsHistory.emplace_back(n, m);

  nbBraids++;
}

void DDCircuit::braidSequence(Sequence& braid) {
  for (const auto& step : braid) {
    if (step.size() != 2 ||
        !std::all_of(step.begin(), step.end(),
                     [](int /*val*/) { return std::is_integral_v<int>; })) {
      throw std::invalid_argument("Indices and powers must be integers!");
    }

    int const ind = step[0];
    int const power = step[1];

    if (ind < 1 || ind >= nbAnyons) {
      throw std::invalid_argument("Invalid braiding operator index!");
    }

    // Computing m and n
    int m = 0;
    int n = 0;
    if (power > 0) {
      n = ind;
      m = ind + 1;
    } else if (power < 0) {
      m = ind;
      n = ind + 1;
    } else { // if power=0, do nothing (identity)
      continue;
    }

    for (int _ = 0; _ < std::abs(power); ++_) {
      this->braid(n, m); // Assuming the braid function exists
    }
  }
}

void DDCircuit::measure() {
  if (measured) {
    throw std::runtime_error("Cannot carry the measurements twice!");
  }
  measured = true;
  // Assuming drawer.measure function exists
  // drawer.measure();
}

ResultDD DDCircuit::run(int shots, std::size_t seed) {
  std::mt19937_64 mt{};
  if (seed != 0U) {
    mt.seed(seed);
  } else {
    // create and properly seed rng
    std::array<std::mt19937_64::result_type, std::mt19937_64::state_size>
        randomData{};
    std::random_device rd;
    std::generate(std::begin(randomData), std::end(randomData),
                  [&rd]() { return rd(); });
    std::seed_seq seeds(std::begin(randomData), std::end(randomData));
    mt.seed(seeds);
  }

  // measure all qubits
  ResultDD counts{};
  auto shotsSt = static_cast<size_t>(shots);
  for (std::size_t i = 0U; i < shotsSt; ++i) {
    // measure all returns a string of the form "q(n-1) ... q(0)"
    auto measurement = circuitDD->measureAll(currentState, false, mt);
    counts.operator[](measurement) += 1U;
  }
  // reduce reference count of measured state
  // circuitDD->decRef(currentState);

  return counts;
}
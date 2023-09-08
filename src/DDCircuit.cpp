#include <cmath>
#include <memory>
#include <vector>
#include "DDCircuit.hpp"
#include "dd/DDDefinitions.hpp"

DDCircuit::DDCircuit(int nbQudits, int nbAnyonsPerQudit)
    : basisGenerator(std::make_unique<BasisGenerator>()), nbQudits(nbQudits),
      nbAnyonsPerQudit(nbAnyonsPerQudit), nbAnyons(nbQudits * nbAnyonsPerQudit){
  circuitDD = std::make_unique<dd::Package<>>(nbQudits);
  operatorGenerator = std::make_unique<OperatorGenerator>(basisGenerator.get());

  generateBasis();
  dim = static_cast<size_t>(pow(2, nbQudits)); // NOLINT

  getSigmas();
//  //        for (int i = 0; i < sigmas.size(); ++i) {
//  //            std::cout << "Sigmas: " << i << " " << sigmas[i] << std::endl;
//  //        }
//  unitary = Eigen::MatrixXcd::Identity(dim, dim);
}

void DDCircuit::generateBasis() {
  basis = basisGenerator->generateBasis(nbQudits, nbAnyonsPerQudit);
}

void DDCircuit::getSigmas() {
  sigmas.reserve(
      static_cast<uint64_t>(nbAnyons - 1)); // Reserve space for sigmas

  size_t const startRow = 1;  // First row
  auto const endRow =
      static_cast<size_t>(pow(2, nbQudits));    // Fifth row (0-based indexing)
  size_t const startCol = 1;  // First column
  auto const endCol = static_cast<size_t>(
      pow(2, nbQudits));    // Fifth column (0-based indexing)
  std::cout << "endRow: " << endRow << '\n';
  std::cout << "endCol: " << endCol << '\n';

  for (int index = 1; index < nbAnyons; ++index) {
    const Eigen::MatrixXcd sigma = operatorGenerator->generateBraidingOperator(
        index, nbQudits, nbAnyonsPerQudit);
    const Eigen::MatrixXcd complexMatrix = sigma.block(startRow, startCol, endRow - startRow + 1, endCol - startCol + 1).cast<std::complex<double>>();
    // Initialize a std::vector<std::vector<std::complex<double>>> with the same dimensions as the Eigen matrix
    dd::CMat sigmaMatrix(
        static_cast<size_t>(complexMatrix.rows()), dd::CVec(static_cast<size_t>(complexMatrix.cols())));
    // Copy elements from the Eigen matrix to the std::vector<std::vector<std::complex<double>>>
    for (int i = 0; i < complexMatrix.rows(); ++i) {
      for (int j = 0; j < complexMatrix.cols(); ++j) {
        auto row = static_cast<size_t>(i);
        auto col = static_cast<size_t>(j);
        sigmaMatrix[row][col] = complexMatrix(i, j);
        //std::cout << sigmaMatrix[row][col] << " ";
      }
      //std::cout << '\n';
    }
    //std::cout << '\n';

    const auto dd = std::make_unique<dd::Package<>>(nbQudits);
    const auto matDD = dd->makeDDFromMatrix(sigmaMatrix);
    braidingOperatorsDD.push_back(matDD);
    assert(braidingOperatorsDD[static_cast<size_t>(index - 1)].p->e[0].w.r != nullptr);
    //TODO: Not necessary to store sigma since DD are used.
    //For debugging purpose
    sigmas.push_back(sigma);
  }

  for(auto& braidingOperator: braidingOperatorsDD){
    assert(braidingOperator.p->e[0].w.r != nullptr);
  }
}

void DDCircuit::initialize(const std::vector<std::complex<dd::fp>>& inputState) {
  if (nbBraids > 0) {
    throw InitializationException("Initialization should happen before any "
                                  "braiding operation is performed!");
  }

  if (inputState.size() != dim) {
    throw std::invalid_argument("The state has wrong dimension. Should be " +
                                std::to_string(dim));
  }

  //TODO: Check whether statevector is normalized
//  double const norm = inputState.dot(inputState.conjugate()).real();
//  if (Eigen::NumTraits<double>::epsilon() <= std::abs(norm - 1)) {
//    throw std::invalid_argument("The input state is not normalized correctly!");
//  }
  // std::cout << "Input state: " << input_state << std::endl;
  auto dd = std::make_unique<dd::Package<>>(nbQudits);
  currentState = dd->makeStateFromVector(inputState);
}

void DDCircuit::braid(int n, int m) {
  if (measured) {
    throw std::runtime_error(
        "System already measured! Cannot perform further braiding!");
  }

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

  //        std::cout << "Unitary before modificaition: " << unitary <<
  //        std::endl;

  if (n < m) {
    unitary = sigmas[n - 1] * unitary;
  } else {
    unitary = sigmas[m - 1].adjoint() * unitary;
  }

  braidsHistory.emplace_back(n, m);

  nbBraids++;

  //        std::cout << "n = " << n << " m = " << m << std::endl;
  //        std::cout << "nb_braids: " << nb_braids << std::endl;
  //        if (n < m) {
  //            std::cout << "n<m sigmas: " << sigmas[n - 1] << std::endl;
  //        } else {
  //            std::cout << "n >= m sigmas: " << sigmas[m - 1].adjoint() <<
  //            std::endl;
  //        }
  //        std::cout << "unitary: " << unitary << std::endl;

  // Assuming drawer.braid function exists
  // drawer.braid(m, n);
}

void DDCircuit::braidDD(int n, int m){
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
  dd::Edge<dd::vNode> e;
  if (n < m) {
    //auto dd = sigmas[n - 1] * unitary;
    index = n - 1;
    assert(braidingOperatorsDD[static_cast<size_t>(index)].p->e[0].w.r != nullptr);
    e = circuitDD->multiply(braidingOperatorsDD[static_cast<size_t>(index)], currentState);
  } else {
    //unitary = sigmas[m - 1].adjoint() * unitary;
    index = m - 1;
    //TODO: Need to adjoint the matrix for the second case
    assert(braidingOperatorsDD[static_cast<size_t>(index)].p->e[0].w.r != nullptr);
    e = circuitDD->multiply(braidingOperatorsDD[static_cast<size_t>(index)], currentState);
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
      this->braidDD(n, m); // Assuming the braid function exists
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

Result DDCircuit::run(int shots, std::size_t seed) {
//  if (!measured) {
//    throw std::runtime_error("The system was not measured!");
//  }
//
//  const Eigen::VectorXcd sv = statevector();
//  Eigen::VectorXd probs = (sv.cwiseProduct(sv.conjugate())).real();
//  //  std::discrete_distribution<int> distribution(probs.data(),
//  //                                               probs.data() + probs.size());
//  std::discrete_distribution<int> distribution(std::begin(probs),
//                                               std::end(probs));
//
//  // Create a random number generator using PCG
//  // Seed with a real random value, if available
//  pcg_extras::seed_seq_from<std::random_device> seedSource;
//  pcg64 generator(seedSource);
//
//  std::vector<int> memory(shots);
//  for (int i = 0; i < shots; ++i) {
//    memory[i] = distribution(generator);
//  }
//
//  std::map<int, int> countsDict;
//  for (int const value : memory) {
//    countsDict[value]++;
//  }
//
//  Result result = {countsDict, memory};
//  return result;

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
  std::map<std::string, std::size_t> counts{};
  for (std::size_t i = 0U; i < shots; ++i) {
    // measure all returns a string of the form "q(n-1) ... q(0)"
    auto measurement = circuitDD->measureAll(currentState, false, mt);
    counts.operator[](measurement) += 1U;
  }
  // reduce reference count of measured state
  circuitDD->decRef(currentState);

//  std::map<std::string, std::size_t> actualCounts{};
//  for (const auto& [bitstring, count] : counts) {
//    std::string measurement(qc->getNcbits(), '0');
//    if (hasMeasurements) {
//      // if the circuit contains measurements, we only want to return the
//      // measured bits
//      for (const auto& [qubit, bit] : measurementMap) {
//        // measurement map specifies that the circuit `qubit` is measured into
//        // a certain `bit`
//        measurement[qc->getNcbits() - 1U - bit] =
//            bitstring[bitstring.size() - 1U - qubit];
//      }
//    } else {
//      // otherwise, we consider the output permutation for determining where
//      // to measure the qubits to
//      for (const auto& [qubit, bit] : qc->outputPermutation) {
//        measurement[qc->getNcbits() - 1 - bit] =
//            bitstring[bitstring.size() - 1U - qubit];
//      }
//    }
//    actualCounts[measurement] += count;
//  }
  return Result{};
}
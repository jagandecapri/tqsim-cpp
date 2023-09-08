#include <memory>
#include "DDCircuit.hpp"
#include "dd/DDDefinitions.hpp"
#include "dd/Package.hpp"

DDCircuit::DDCircuit(int nbQudits, int nbAnyonsPerQudit)
    : basisGenerator(std::make_unique<BasisGenerator>()), nbQudits(nbQudits),
      nbAnyonsPerQudit(nbAnyonsPerQudit), nbAnyons(nbQudits * nbAnyonsPerQudit){
  operatorGenerator = std::make_unique<OperatorGenerator>(basisGenerator.get());

  generateBasis();
  dim = static_cast<Eigen::Index>(basis.size()); // NOLINT
  initialState = Eigen::VectorXcd::Zero(dim);
  initialState(0) = 1.0;

  getSigmas();
  //        for (int i = 0; i < sigmas.size(); ++i) {
  //            std::cout << "Sigmas: " << i << " " << sigmas[i] << std::endl;
  //        }
  unitary = Eigen::MatrixXcd::Identity(dim, dim);
}

void DDCircuit::generateBasis() {
  basis = basisGenerator->generateBasis(nbQudits, nbAnyonsPerQudit);
}

void DDCircuit::getSigmas() {
  sigmas.reserve(nbAnyons - 1); // Reserve space for sigmas

  int const startRow = 1;  // First row
  int const endRow = pow(2, nbQudits);    // Fifth row (0-based indexing)
  int const startCol = 1;  // First column
  int const endCol = pow(2, nbQudits);    // Fifth column (0-based indexing)
  for (int index = 1; index < nbAnyons; ++index) {
    const Eigen::MatrixXcd sigma = operatorGenerator->generateBraidingOperator(
        index, nbQudits, nbAnyonsPerQudit).block(startRow, startCol, endRow - startRow + 1, endCol - startCol + 1).cast<std::complex<double>>();
    // Initialize a std::vector<std::vector<std::complex<double>>> with the same dimensions as the Eigen matrix
    dd::CMat sigmaMatrix(
        sigma.rows(), std::vector<std::complex<double>>(sigma.cols()));
    const auto dd = std::make_unique<dd::Package<>>(nbQudits);
    const auto matDD = dd->makeDDFromMatrix(sigmaMatrix);
    decisionDiagrams.push_back(matDD);

    //TODO: Not necessary to store sigma since DD are used.
    //For debugging purpose
    sigmas.push_back(sigma);
  }
}

void DDCircuit::initialize(const Eigen::VectorXcd& inputState) {
  if (nbBraids > 0) {
    throw InitializationException("Initialization should happen before any "
                                  "braiding operation is performed!");
  }

  if (inputState.size() != dim) {
    throw std::invalid_argument("The state has wrong dimension. Should be " +
                                std::to_string(dim));
  }

  double const norm = inputState.dot(inputState.conjugate()).real();
  if (Eigen::NumTraits<double>::epsilon() <= std::abs(norm - 1)) {
    throw std::invalid_argument("The input state is not normalized correctly!");
  }
  // std::cout << "Input state: " << input_state << std::endl;
  initialState = inputState;
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

Result DDCircuit::run(int shots) {
  if (!measured) {
    throw std::runtime_error("The system was not measured!");
  }

  const Eigen::VectorXcd sv = statevector();
  Eigen::VectorXd probs = (sv.cwiseProduct(sv.conjugate())).real();
  //  std::discrete_distribution<int> distribution(probs.data(),
  //                                               probs.data() + probs.size());
  std::discrete_distribution<int> distribution(std::begin(probs),
                                               std::end(probs));

  // Create a random number generator using PCG
  // Seed with a real random value, if available
  pcg_extras::seed_seq_from<std::random_device> seedSource;
  pcg64 generator(seedSource);

  std::vector<int> memory(shots);
  for (int i = 0; i < shots; ++i) {
    memory[i] = distribution(generator);
  }

  std::map<int, int> countsDict;
  for (int const value : memory) {
    countsDict[value]++;
  }

  Result result = {countsDict, memory};
  return result;
}
#include "Circuit.hpp"

Circuit::Circuit(int nb_qudits, int nb_anyons_per_qudit)
    : nb_qudits(nb_qudits), nb_anyons_per_qudit(nb_anyons_per_qudit),
      nb_anyons(nb_qudits * nb_anyons_per_qudit), nb_braids(0),
      measured(false) {
  basis_generator = new BasisGenerator();
  operator_generator = new OperatorGenerator(basis_generator);

  generate_basis();
  dim = basis.size();
  initial_state = Eigen::VectorXcd::Zero(dim);
  initial_state(0) = 1.0;

  get_sigmas();
  //        for (int i = 0; i < sigmas.size(); ++i) {
  //            std::cout << "Sigmas: " << i << " " << sigmas[i] << std::endl;
  //        }
  unitary = Eigen::MatrixXcd::Identity(dim, dim);
}

void Circuit::generate_basis() {
  basis = basis_generator->generate_basis(nb_qudits, nb_anyons_per_qudit);
}

void Circuit::get_sigmas() {
  sigmas.reserve(nb_anyons - 1); // Reserve space for sigmas

  for (int index = 1; index < nb_anyons; ++index) {
    Eigen::MatrixXcd sigma = operator_generator->generate_braiding_operator(
        index, nb_qudits, nb_anyons_per_qudit);
    sigmas.push_back(sigma);
  }
}

void Circuit::initialize(const Eigen::VectorXcd& input_state) {
  if (nb_braids > 0) {
    throw InitializationException("Initialization should happen before any "
                                  "braiding operation is performed!");
  }

  if (input_state.size() != dim) {
    throw std::invalid_argument("The state has wrong dimension. Should be " +
                                std::to_string(dim));
  }

  double norm = input_state.dot(input_state.conjugate()).real();
  if (Eigen::NumTraits<double>::epsilon() <= std::abs(norm - 1)) {
    throw std::invalid_argument("The input state is not normalized correctly!");
  }
  // std::cout << "Input state: " << input_state << std::endl;
  initial_state = input_state;
}

void Circuit::braid(int n, int m) {
  if (measured) {
    throw std::runtime_error(
        "System already measured! Cannot perform further braiding!");
  }

  if (m < 1 || n < 1) {
    throw std::invalid_argument("n, m must be higher than 0!");
  }

  if (m > nb_anyons || n > nb_anyons) {
    throw std::invalid_argument("The system has only " +
                                std::to_string(nb_anyons) +
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

  braids_history.emplace_back(n, m);

  nb_braids++;

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

void Circuit::braid_sequence(Sequence& braid_sequence) {
  for (const auto& step : braid_sequence) {
    if (step.size() != 2 || !std::all_of(step.begin(), step.end(), [](int val) {
          return std::is_integral<int>::value;
        })) {
      throw std::invalid_argument("Indices and powers must be integers!");
    }

    int ind = step[0];
    int power = step[1];

    if (ind < 1 || ind >= nb_anyons) {
      throw std::invalid_argument("Invalid braiding operator index!");
    }

    // Computing m and n
    int m = 0, n = 0;
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
      braid(n, m); // Assuming the braid function exists
    }
  }
}

void Circuit::measure() {
  if (measured) {
    throw std::runtime_error("Cannot carry the measurements twice!");
  }
  measured = true;
  // Assuming drawer.measure function exists
  // drawer.measure();
}

Result Circuit::run(int shots) {
  if (!measured) {
    throw std::runtime_error("The system was not measured!");
  }

  Eigen::VectorXcd sv = statevector();
  Eigen::VectorXd probs = (sv.cwiseProduct(sv.conjugate())).real();
  std::discrete_distribution<int> distribution(probs.data(),
                                               probs.data() + probs.size());

  // Create a random number generator using PCG
  // Seed with a real random value, if available
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  pcg64 generator(seed_source);

  std::vector<int> memory(shots);
  for (int i = 0; i < shots; ++i) {
    memory[i] = distribution(generator);
  }

  std::map<int, int> counts_dict;
  for (int value : memory) {
    counts_dict[value]++;
  }

  Result result = {counts_dict, memory};
  return result;
}

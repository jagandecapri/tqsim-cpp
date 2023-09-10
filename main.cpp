#include "dd/DDDefinitions.hpp"
#include "dd/Export.hpp"
#include "dd/Package.hpp"
#include "src/Circuit.cpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int /*unused*/, char** /*unused*/) {
  cout << "Hello, from tqsim-cpp!\n";
  // H_2
  Sequence hadSequence2 = {{4, 2}, {5, 2},  {4, -2}, {5, -2}, {4, 2},
                           {5, 4}, {4, -2}, {5, 2},  {4, 2},  {5, -2},
                           {4, 2}, {5, -2}, {4, 4}};

  // CNOT_2->1
  Sequence cnotSequence = {
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, -1}, {3, -1}, {3, -1}, {4, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, -1}, {2, -1}, {2, -1}, {1, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, -1}, {2, -1}, {2, -1}, {1, -1}, {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {1, -1}, {2, -1}, {2, -1}, {1, -1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, -1}, {2, -1}, {2, -1}, {1, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1}};

  // Start the timer
  auto start = std::chrono::high_resolution_clock::now();
  const auto nrQubits = 2U;

  Circuit circuit = Circuit(nrQubits, 3);
  Eigen::VectorXcd initSequence(13);
  initSequence << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  circuit.initialize(initSequence);
  circuit.braidSequence(hadSequence2);
  circuit.braidSequence(cnotSequence);

  circuit.measure();
  int const shots = 10000;
  Result const res = circuit.run(shots);
  // Stop the timer
  auto end = std::chrono::high_resolution_clock::now();

  // Calculate the elapsed time
  std::chrono::duration<double> const elapsed = end - start;

  // Output the elapsed time in seconds
  std::cout << "Number of shots: " << shots << '\n';
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << '\n';
  for (const auto& count : res.countsDict) {
    std::cout << count.first << ": " << count.second << '\n';
  }

  std::vector<Eigen::MatrixXcd> const braidingOperators =
      circuit.getBraidingOperators();

  // Define the range for rows and columns you want to select
  size_t const startRow = 1; // First row
  auto const endRow =
      static_cast<size_t>(pow(2, nrQubits)); // Fifth row (0-based indexing)
  int const startCol = 1;                    // First column
  auto const endCol =
      static_cast<size_t>(pow(2, nrQubits)); // Fifth column (0-based indexing)

  std::cout << "endRow: " << endRow << '\n';
  std::cout << "endCol: " << endCol << '\n';

  for (int index = 0; index < static_cast<int>(braidingOperators.size());
       ++index) {

    // Create a std::ostringstream object
    std::ostringstream oss;

    // Format the string
    oss << std::to_string(nrQubits) << "_qubit_"
        << "sigma_" << index + 1 << "_operator.dot";

    // Get the formatted string
    std::string const formattedString = oss.str();

    // Print the formatted string
    std::cout << formattedString << '\n';

    // Select the sub-matrix
    Eigen::MatrixXcd complexMatrix =
        braidingOperators[static_cast<size_t>(index)]
            .block(startRow, startCol, endRow - startRow + 1,
                   endCol - startCol + 1)
            .cast<std::complex<double>>();

    std::cout << "Braiding operator: \n"
              << braidingOperators[static_cast<uint64_t>(index)] << '\n';

    // Initialize a std::vector<std::vector<std::complex<double>>> with the same
    // dimensions as the Eigen matrix
    dd::CMat sigmaMatrix(static_cast<size_t>(complexMatrix.rows()),
                         dd::CVec(static_cast<size_t>(complexMatrix.cols())));

    // Copy elements from the Eigen matrix to the
    // std::vector<std::vector<std::complex<double>>>
    for (int i = 0; i < complexMatrix.rows(); ++i) {
      for (int j = 0; j < complexMatrix.cols(); ++j) {
        auto row = static_cast<size_t>(i);
        auto col = static_cast<size_t>(j);
        sigmaMatrix[row][col] = complexMatrix(i, j);
        std::cout << sigmaMatrix[row][col] << " ";
      }
      std::cout << '\n';
    }
    std::cout << '\n';

    const auto dd = std::make_unique<dd::Package<>>(nrQubits);
    const auto matDD = dd->makeDDFromMatrix(sigmaMatrix);
    export2Dot(matDD, formattedString, true, true, false, false, true);
  }

  return 0;
}

#include "dd/DDDefinitions.hpp"
#include "dd/Export.hpp"
#include "dd/Package.hpp"
#include "src/DDCircuit.cpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int /*unused*/, char** /*unused*/) {
  cout << "Hello, from tqsim-cpp DD!\n";

  // X
  Sequence x_seq = {{1, -2}, {2, -4}, {1, 4},  {2, -2}, {1, 2}, {2, 2},
                    {1, -2}, {2, 4},  {1, -2}, {2, 4},  {1, 2}, {2, -4},
                    {1, 2},  {2, -2}, {1, 2},  {2, -2}, {1, -2}};

  // Hadamard on second qubit
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

  // Sequence hadSequence2 = {{2, 1}};
  //  Start the timer
  auto start = std::chrono::high_resolution_clock::now();
  const auto nrQubits = 2;

  DDCircuit circuit = DDCircuit(nrQubits, 3);
  std::vector<std::complex<dd::fp>> const initSequence{1, 0, 0, 0};
  circuit.initialize(initSequence);

  circuit.braidSequence(x_seq);
  circuit.braidSequence(hadSequence2);

  // TODO: CNOT not working properly as Warning below is thrown:
  // WARNING in MAll: numerical instability occurred during simulation:
  // |alpha|^2 + |beta|^2 = xyz, but should be 1!
  circuit.braidSequence(cnotSequence);

  circuit.measure();
  int const shots = 50;
  size_t const seeds = 123;
  ResultDD const res = circuit.run(shots, seeds);
  // Stop the timer
  auto end = std::chrono::high_resolution_clock::now();

  // Calculate the elapsed time
  std::chrono::duration<double> const elapsed = end - start;

  // Output the elapsed time in seconds
  std::cout << "Number of shots: " << shots << '\n';
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << '\n';
  for (const auto& [bitstring, count] : res) {
    std::cout << bitstring << ": " << count << '\n';
  }
  return 0;
}

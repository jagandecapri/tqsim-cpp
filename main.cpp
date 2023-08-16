#include "src/Circuit.cpp"

#include <chrono>
#include <iostream>

using namespace std;

int main(int, char**) {
  cout << "Hello, from tqsim-cpp!\n";
  // H_2
  Sequence had_sequence_2 = {{4, 2}, {5, 2},  {4, -2}, {5, -2}, {4, 2},
                             {5, 4}, {4, -2}, {5, 2},  {4, 2},  {5, -2},
                             {4, 2}, {5, -2}, {4, 4}};

  // CNOT_2->1
  Sequence cnot_sequence = {
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

  Circuit circuit = Circuit(2, 3);
  Eigen::VectorXcd init_sequence(13);
  init_sequence << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  circuit.initialize(init_sequence);
  circuit.braidSequence(had_sequence_2);
  circuit.braidSequence(cnot_sequence);
  circuit.measure();
  int shots = 100000000;
  Result res = circuit.run(shots);
  // Stop the timer
  auto end = std::chrono::high_resolution_clock::now();

  // Calculate the elapsed time
  std::chrono::duration<double> elapsed = end - start;

  // Output the elapsed time in seconds
  std::cout << "Number of shots: " << shots << std::endl;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
  for (auto& count : res.counts_dict) {
    std::cout << count.first << ": " << count.second << std::endl;
  }
  return 0;
}

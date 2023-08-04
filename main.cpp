#include <iostream>
#include "src/circuit.h"

using namespace std;

int main(int, char**){
    Circuit circuit = Circuit(1);
    circuit.print();
    cout << "Hello, from tqsim-cpp!\n";
    return 0;
}

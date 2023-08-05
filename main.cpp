#include <iostream>
#include "src/circuit.h"

using namespace std;

int main(int, char**){
    Circuit circuit = Circuit(2,3);
    circuit.print();
    cout << "Hello, from tqsim-cpp!\n";
    return 0;
}

# tqsim-cpp

A minimal C++ implementation of a Topological Quantum Computing (TQC) simulator adapted
from [arXiv:2307.01892](https://arxiv.org/abs/2307.01892) ([Python impl.](https://github.com/Constantine-Quantum-Tech/tqsim))

Rudimentary execution time test for `100 million` shots gives a ~6 times execution time improvement over the Python
implementation (~1.1 seconds vs ~6.5751 seconds).

No drawing functionality is implemented as it is more friendly to be implemented using Python.

Engineering improvements can be done to parallelize random number generation which takes the most time when program is
profiled and potentially other improvement.

Algorithm wise, potentially better improvement can be done by using better data structures as quoted from the paper:

> it is important to note that this method does not guarantee
> a reduction in computational complexity since the size of the braid generator grows exponentially with the number
> of anyons

Entry point of the program is `main.cpp` which has a simple example.

Tested with GNU GCC/G++-13.1.0 with C++17 standard.

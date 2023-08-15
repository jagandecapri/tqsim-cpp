# tqsim-cpp

A minimal C++ implementation of Topological Quantum Computing (TQC) simulator adapted
from [arXiv:2307.01892](https://arxiv.org/abs/2307.01892) ([Python impl.](https://github.com/Constantine-Quantum-Tech/tqsim))

Rudimentary execution time test for `100 million` shots gives a ~5 times execution time improvement over the Python
implementation (1.43034 seconds vs 6.5751 seconds).

No drawing functionality is implemented as it is more friendly to be implemented using Python.

Engineering improvements can be done to parallelize random number generation which takes the most time when program is
profiled.

Potentially bigger improvement can be done by using better data structures.

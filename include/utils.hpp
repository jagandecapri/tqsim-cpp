#pragma once

#include <exception>

class InitializationException : public std::exception {
public:
  explicit InitializationException(const char* message) : msg(message) {}

  [[nodiscard]] const char* what() const noexcept override { return msg; }

private:
  const char* msg;
};

// https://noobtuts.com/cpp/compare-float-values
bool cmpf(double A, double B, double epsilon = 0.005f) {
  return (fabs(A - B) < epsilon);
}

#pragma once

#include <exception>

class InitializationException : public std::exception {
public:
  explicit InitializationException(const char* message) : msg(message) {}

  [[nodiscard]] const char* what() const noexcept override { return msg; }

private:
  const char* msg;
};

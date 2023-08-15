//
// Created by Jagan on 09/08/2023.
//

#ifndef TQSIM_CPP_UTILS_H
#define TQSIM_CPP_UTILS_H

#include <exception>

class InitializationException : public std::exception {
public:
    explicit InitializationException(const char *message) : msg(message) {}

    const char *what() const noexcept override {
        return msg;
    }

private:
    const char *msg;
};

#endif //TQSIM_CPP_UTILS_H

//
// Created by luca on 13/07/22.
//
#include <exception>
#include <string>
#ifndef QECC_QECCEXCEPTION_HPP
#define QECC_QECCEXCEPTION_HPP

struct QeccException : public std::exception {
private:
    std::string message;

public:
    explicit QeccException(const char* msg) : message(msg) {}
    std::string getMessage() {
        return message;
    }
    [[nodiscard]] const char* what() const noexcept override {
        return message.c_str();
    }
};
#endif // QECC_QECCEXCEPTION_HPP

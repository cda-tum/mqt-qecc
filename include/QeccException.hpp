#pragma once

#include <exception>
#include <string>

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

//
// Created by luca on 13/07/22.
//
#include <exception>
#include <string>
#ifndef QECC_QECCEXCEPTION_HPP
#define QECC_QECCEXCEPTION_HPP

struct QeccException: public std::exception {
private:
    const char* message;
public:
    explicit QeccException(const char * msg) : message(msg) {}
    const char * getMessage () {
        return message;
    }
};
#endif //QECC_QECCEXCEPTION_HPP

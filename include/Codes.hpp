//
// Created by lucas on 17/05/2022.
//

#ifndef QUNIONFIND_CODES_HPP
#define QUNIONFIND_CODES_HPP
#include "Code.hpp"

#include <utility>

class SteaneXCode : public Code {
public:
    SteaneXCode() : Code("./resources/codes/testCode.txt") {
        K = 3;
        N = 7;
    }
};

class SteaneCode : public Code {
public:
    SteaneCode() : Code("./resources/codes/testCode.txt", "./resources/codes/testCode.txt") {
        K = 3;
        N = 7;
    }
};

class HGPcode : public Code {
public:
    HGPcode(const std::string& inFile, const std::size_t K) : Code(inFile) { this->K = K; };
    HGPcode() : Code("./resources/codes/hgp_(4,7)-[[900,36,10]]_hz.txt") {
        K = 36;
        N = 900;
    }
    HGPcode(const std::string& hxIn, const std::string& hzIn, const std::size_t K) : Code(hxIn, hzIn) { this->K = K; };
};

class ToricCode_8 : public Code {
public:
    ToricCode_8() : Code("./resources/codes/toric_(nan,nan)-[[8,2,2]]_hx.txt") {
        K = 2;
        N = 8;
    }
};

class ToricCode_18 : public Code {
public:
    ToricCode_18() : Code("./resources/codes/toric_(nan,nan)-[[18,2,3]]_hx.txt") {
        K = 2;
        N = 18;
    }
};

class ToricCode_32 : public Code {
public:
    ToricCode_32() : Code("./resources/codes/toric_(nan,nan)-[[32,2,4]]_hx.txt") {
        K = 2;
        N = 32;
    }
};
#endif // QUNIONFIND_CODES_HPP

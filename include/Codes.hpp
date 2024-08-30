#pragma once

#include "Code.hpp"

#include <cstddef>
#include <string>

class SteaneXCode : public Code {
public:
    SteaneXCode() : Code("./resources/codes/testCode.txt") {
        k = 3;
        n = 7;
    }
};

class SteaneCode : public Code {
public:
    SteaneCode() : Code("./resources/codes/testCode.txt", "./resources/codes/testCode.txt") {
        k = 3;
        n = 7;
    }
};

class HGPcode : public Code {
public:
    HGPcode(const std::string& inFile, const std::size_t size) : Code(inFile) { k = size; };
    HGPcode() : Code("./resources/codes/hgp_(4,7)-[[900,36,10]]_hz.txt") {
        k = 36;
        n = 900;
    }
    HGPcode(const std::string& hxIn, const std::string& hzIn, const std::size_t size) : Code(hxIn, hzIn) { k = size; };
};

class ToricCode8 : public Code {
public:
    ToricCode8() : Code("./resources/codes/toric_(nan,nan)-[[8,2,2]]_hx.txt") {
        k = 2;
        n = 8;
    }
};

class ToricCode18 : public Code {
public:
    ToricCode18() : Code("./resources/codes/toric_(nan,nan)-[[18,2,3]]_hx.txt") {
        k = 2;
        n = 18;
    }
};

class ToricCode32 : public Code {
public:
    ToricCode32() : Code("./resources/codes/toric_(nan,nan)-[[32,2,4]]_hx.txt") {
        k = 2;
        n = 32;
    }
};

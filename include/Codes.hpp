//
// Created by luca on 17/05/2022.
//

#ifndef QUNIONFIND_CODES_HPP
#define QUNIONFIND_CODES_HPP
#include "Code.hpp"

#include <utility>

class SteaneXCode: public Code {
public:
    SteaneXCode():
        Code(ParityCheckMatrix({{1, 0, 0, 1, 0, 1, 1},
                                {0, 1, 0, 1, 1, 0, 1},
                                {0, 0, 1, 0, 1, 1, 1}})) {K = 3;}
};
class LargeCode: public Code {
public: LargeCode() :
        Code(ParityCheckMatrix({{1, 0},
                                {0, 1}})) {K = 100;}
};
#endif //QUNIONFIND_CODES_HPP

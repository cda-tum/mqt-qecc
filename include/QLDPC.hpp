//
// Created by luca on 09/06/22.
//

#ifndef QUNIONFIND_QLDPC_HPP
#define QUNIONFIND_QLDPC_HPP
#include "Code.hpp"
#include "Decoder.hpp"
/*
 * API
 */
class QLDPC {
public:
    /*
     * Code Factory
     */
    Code getCode(const std::string name);

    /*
     * Decoder Factory
     */
    Decoder getDecoder();

    /**
     * Computes estimated X error given a Z syndrome measurement
     * @param error
     * @return
     */
    DecodingResult decode(std::vector<bool> syndrome);

};
#endif //QUNIONFIND_QLDPC_HPP

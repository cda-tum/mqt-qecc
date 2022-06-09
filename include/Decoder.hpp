//
// Created by lucas on 21/04/2022.
//

#include "Code.hpp"
#include "Codes.hpp"
#include "TreeNode.hpp"

#include <chrono>
#include <nlohmann/json.hpp>
#include <utility>
#include <vector>

using json = nlohmann::json;
#ifndef QUNIONFIND_DECODER_HPP
#define QUNIONFIND_DECODER_HPP

enum DecodingResultStatus {
    SUCCESS,
    FLAGGED_ERROR,
    FAILURE //only to be set from outside when logical operator is introduced
};

struct DecodingResult {
    DecodingResultStatus     status{};
    std::size_t              decodingTime = 0U;
    std::vector<std::size_t> estimNodeIdxVector;
    std::vector<bool>        estimBoolVector;
    [[nodiscard]] json       to_json() const {
              return json{
                {"staus", status},
                {"decodingTime(ms)", decodingTime}};
    }
};

class Decoder {
public:
    DecodingResult result;
    explicit Decoder(Code& code):
        code(code) {}
    virtual void decode(std::set<std::shared_ptr<TreeNode>>& syndrome){}; // todo set or vector more efficient?
protected:
    Code code;
};
#endif //QUNIONFIND_DECODER_HPP

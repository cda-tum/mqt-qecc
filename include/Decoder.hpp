//
// Created by lucas on 21/04/2022.
//
#ifndef QUNIONFIND_DECODER_HPP
#define QUNIONFIND_DECODER_HPP
#include "Code.hpp"
#include "Codes.hpp"
#include "TreeNode.hpp"

#include <chrono>
#include <nlohmann/json.hpp>
#include <utility>
#include <vector>

using json = nlohmann::json;
struct DecodingResult {
    std::size_t              decodingTime       = 0U; // in ms
    std::vector<std::size_t> estimNodeIdxVector = {};
    gf2Vec                   estimBoolVector    = {};
    [[nodiscard]] json       to_json() const {
              return json{
                {"decodingTime(ms)", decodingTime},
                {"estimate", Utils::getStringFrom(estimBoolVector)}};
    }
    void from_json(const json& j) {
        j.at("decodingTime(ms)").get_to(decodingTime);
        j.at("estimate").get_to(estimBoolVector);
        j.at("estimatedNodes").get_to(estimNodeIdxVector);
    }
};
class Decoder {
public:
    DecodingResult result;

    explicit Decoder(Code code):
        code(std::move(code)) {}
    virtual void decode(std::vector<bool>&){};

    Decoder(const Decoder& other):
        Decoder(Code(other.code)) {
        this->result = DecodingResult();
    }
    [[nodiscard]] Code getCode() const {
        return this->code;
    }

protected:
    Code code;
};
#endif //QUNIONFIND_DECODER_HPP

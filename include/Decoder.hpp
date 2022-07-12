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

enum GrowthVariant {
    ALL_COMPONENTS, // standard growth
    INVALID_COMPONENTS,
    SINGLE_SMALLEST,
    SINGLE_RANDOM
};
NLOHMANN_JSON_SERIALIZE_ENUM(GrowthVariant, {{ALL_COMPONENTS, "all components"},
                                             {INVALID_COMPONENTS, "invalid components"},
                                             {SINGLE_SMALLEST, "smallest component only"},
                                             {SINGLE_RANDOM, "single random component"}})
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
    GrowthVariant  growth = ALL_COMPONENTS; // standard
    explicit Decoder(Code* code): code(code){}
    virtual void decode(std::vector<bool>&){};

    Decoder(const Decoder& other): code(other.code) {
        this->result = DecodingResult();
    }
    [[nodiscard]] Code* getCode() const {
        return this->code;
    }
    [[nodiscard]] GrowthVariant getGrowth() const {
        return growth;
    }
    void setGrowth(GrowthVariant g) {
        Decoder::growth = g;
    }
protected:
    Code* code;
};
#endif //QUNIONFIND_DECODER_HPP

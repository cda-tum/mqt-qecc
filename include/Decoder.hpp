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

enum class GrowthVariant {
    ALL_COMPONENTS, // standard growth
    INVALID_COMPONENTS,
    SINGLE_SMALLEST,
    SINGLE_RANDOM,
    SINGLE_QUBIT_RANDOM
};

[[maybe_unused]] static GrowthVariant growthVariantFromString(const std::string& architecture) {
    if (architecture == "ALL_COMPONENTS" || architecture == "0") {
        return GrowthVariant::ALL_COMPONENTS;
    } else if (architecture == "INVALID_COMPONENTS" || architecture == "1") {
        return GrowthVariant::INVALID_COMPONENTS;
    } else if (architecture == "SINGLE_SMALLEST" || architecture == "2") {
        return GrowthVariant::SINGLE_SMALLEST;
    } else if (architecture == "SINGLE_RANDOM" || architecture == "3") {
        return GrowthVariant::SINGLE_RANDOM;
    } else if (architecture == "SINGLE_QUBIT_RANDOM" || architecture == "4") {
        return GrowthVariant::SINGLE_QUBIT_RANDOM;
    } else {
        throw std::invalid_argument("Invalid growth variant: " + architecture);
    }
}

NLOHMANN_JSON_SERIALIZE_ENUM(GrowthVariant, {{GrowthVariant::ALL_COMPONENTS, "all components"},
                                             {GrowthVariant::INVALID_COMPONENTS, "invalid components"},
                                             {GrowthVariant::SINGLE_SMALLEST, "smallest component only"},
                                             {GrowthVariant::SINGLE_QUBIT_RANDOM, "single random qubit only"},
                                             {GrowthVariant::SINGLE_RANDOM, "single random component"}})
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
    [[nodiscard]] std::string toString() const{
        return this->to_json().dump(2U);
    }
};
class Decoder {
private:
    std::unique_ptr<Code> code;

public:
    DecodingResult result{};
    GrowthVariant  growth = GrowthVariant::ALL_COMPONENTS; // standard

    Decoder() = default;
    virtual void decode(const std::vector<bool>&){};
    virtual ~Decoder() = default;

    [[nodiscard]] const std::unique_ptr<Code>& getCode() const {
        return code;
    }
    [[nodiscard]] GrowthVariant getGrowth() const {
        return growth;
    }
    void setGrowth(GrowthVariant g) {
        Decoder::growth = g;
    }
    void setCode(Code& c) {
        this->code = std::make_unique<Code>(*c.getHz()->pcm);
    }
    virtual void reset(){};
};
#endif //QUNIONFIND_DECODER_HPP

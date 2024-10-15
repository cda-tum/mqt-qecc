#pragma once

#include "Code.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using json = nlohmann::basic_json<>;

enum class GrowthVariant : std::uint8_t {
    AllComponents, // standard growth
    InvalidComponents,
    SingleSmallest,
    SingleRandom,
    SingleQubitRandom
};

[[maybe_unused]] static GrowthVariant growthVariantFromString(const std::string& architecture) {
    if (architecture == "ALL_COMPONENTS" || architecture == "0") {
        return GrowthVariant::AllComponents;
    }
    if (architecture == "INVALID_COMPONENTS" || architecture == "1") {
        return GrowthVariant::InvalidComponents;
    }
    if (architecture == "SINGLE_SMALLEST" || architecture == "2") {
        return GrowthVariant::SingleSmallest;
    }
    if (architecture == "SINGLE_RANDOM" || architecture == "3") {
        return GrowthVariant::SingleRandom;
    }
    if (architecture == "SINGLE_QUBIT_RANDOM" || architecture == "4") {
        return GrowthVariant::SingleQubitRandom;
    }
    throw std::invalid_argument("Invalid growth variant: " + architecture);
}

NLOHMANN_JSON_SERIALIZE_ENUM(GrowthVariant, {{GrowthVariant::AllComponents, "all components"}, // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,misc-include-cleaner)
                                             {GrowthVariant::InvalidComponents, "invalid components"},
                                             {GrowthVariant::SingleSmallest, "smallest component only"},
                                             {GrowthVariant::SingleQubitRandom, "single random qubit only"},
                                             {GrowthVariant::SingleRandom, "single random component"}})
struct DecodingResult {
    std::size_t              decodingTime = 0U; // in ms
    std::vector<std::size_t> estimNodeIdxVector;
    gf2Vec                   estimBoolVector;

    [[nodiscard]] json to_json() const { // NOLINT(readability-identifier-naming)
        return json{{"decodingTime(ms)", decodingTime},
                    {"estimate", Utils::getStringFrom(estimBoolVector)}};
    }
    void from_json(const json& j) { // NOLINT(readability-identifier-naming)
        j.at("decodingTime(ms)").get_to(decodingTime);
        j.at("estimate").get_to(estimBoolVector);
        j.at("estimatedNodes").get_to(estimNodeIdxVector);
    }
    [[nodiscard]] std::string toString() const {
        return this->to_json().dump(2U);
    }
};
class Decoder {
private:
    std::unique_ptr<Code> code;

public:
    DecodingResult result{};
    GrowthVariant  growth = GrowthVariant::AllComponents; // standard

    Decoder() = default;
    virtual void decode(const std::vector<bool>&) {}; // NOLINT(readability-named-parameter)
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
        if (c.gethX() == nullptr) {
            this->code = std::make_unique<Code>(*c.gethZ()->pcm);
        } else {
            this->code = std::make_unique<Code>(*c.gethZ()->pcm, *c.gethX()->pcm);
        }
    }
    virtual void reset() {};
};

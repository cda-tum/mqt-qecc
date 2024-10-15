#pragma once

#include "Decoder.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <utility>

using json = nlohmann::basic_json<>;

enum DecodingResultStatus : std::uint8_t {
    SUCCESS, // estimated correctly up to stabilizer
    FAILURE  // logical operator introduced
    // FLAGGED_ERROR // not implemented
};

[[maybe_unused]] static DecodingResultStatus decodingResultStatusFromString(const std::string& status) {
    if (status == "SUCCESS" || status == "0") {
        return DecodingResultStatus::SUCCESS;
    }
    if (status == "FAILURE" || status == "1") {
        return DecodingResultStatus::FAILURE;
    }
    throw std::invalid_argument("Invalid decoding result status: " + status);
}

NLOHMANN_JSON_SERIALIZE_ENUM(DecodingResultStatus, {{SUCCESS, "success"}, // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,misc-include-cleaner)
                                                    {FAILURE, "failure"}})
/**
 * Contains information about a single run of a decoder
 * status determines whether the estimate is valid or not
 * result contains information obtained from the decoder
 */
struct DecodingRunInformation {
    DecodingRunInformation(double               physicalErrorRate,
                           std::size_t          size,
                           gf2Vec               err,
                           gf2Vec               syn,
                           DecodingResultStatus resultStatus,
                           DecodingResult       res) : physicalErrR(physicalErrorRate),
                                                 codeSize(size), error(std::move(err)), syndrome(std::move(syn)), status(resultStatus), result(std::move(res)) {}
    DecodingRunInformation(double         physicalErrorRate,
                           std::size_t    size,
                           gf2Vec         err,
                           gf2Vec         syn,
                           DecodingResult res) : physicalErrR(physicalErrorRate),
                                                 codeSize(size), error(std::move(err)), syndrome(std::move(syn)), result(std::move(res)) {}
    DecodingRunInformation() = default;

    double               physicalErrR = 0.0;
    std::size_t          codeSize     = 0U;
    gf2Vec               error;
    gf2Vec               syndrome;
    DecodingResultStatus status{};
    DecodingResult       result{};

    [[nodiscard]] json to_json() const { // NOLINT(readability-identifier-naming)
        return json{{"physicalErrRate", physicalErrR},
                    {"codeSize", codeSize},
                    {"error", Utils::getStringFrom(error)},
                    {"syndrome", Utils::getStringFrom(syndrome)},
                    {"decodingResult", result.to_json()},
                    {"decodingStatus", status}};
    }

    void from_json(const json& j) { // NOLINT(readability-identifier-naming)
        j.at("physicalErrRate").get_to(physicalErrR);
        j.at("codeSize").get_to(codeSize);
        j.at("error").get_to(error);
        j.at("syndrome").get_to(syndrome);
        j.at("decodingStatus").get_to(status);
    }

    [[nodiscard]] std::string toString() const {
        return to_json().dump(2U);
    }
    void print() const {
        std::cout << toString() << '\n';
    }
};

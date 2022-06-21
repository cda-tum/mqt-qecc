//
// Created by lucas on 21/06/22.
//
#ifndef QUNIONFIND_DECODINGRUNINFORMATION_HPP
#define QUNIONFIND_DECODINGRUNINFORMATION_HPP
#include "Decoder.hpp"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

enum DecodingResultStatus {
    SUCCESS, // estimated correctly up to stabilizer
    FAILURE  // logical operator introduced
    //FLAGGED_ERROR // not implemented
};
NLOHMANN_JSON_SERIALIZE_ENUM(DecodingResultStatus, {{SUCCESS, "success"},
                                                    {FAILURE, "failure"}})
/**
 * Contains information about a single run of a decoder
 * status determines whether the estimate is valid or not
 * result contains information obtained from the decoder
 */
struct DecodingRunInformation {
    DecodingResultStatus status{};
    DecodingResult       result{};
    [[nodiscard]] json   to_json() const {
          return json{
                {"decodingResult", result.to_json(),
                   "decodingStatus", status}};
    }
};
#endif //QUNIONFIND_DECODINGRUNINFORMATION_HPP

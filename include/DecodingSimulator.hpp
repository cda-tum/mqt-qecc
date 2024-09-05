#pragma once

#include "Code.hpp"

#include <cstddef>
#include <cstdint>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

using json = nlohmann::basic_json<>;

enum DecoderType : std::uint8_t {
    UfHeuristic,
    UfDecoder
};
[[maybe_unused]] static DecoderType decoderTypeFromString(const std::string& status) {
    if (status == "UF_HEURISTIC" || status == "0") {
        return DecoderType::UfHeuristic;
    }
    if (status == "ORIGINAL_UF" || status == "1") {
        return DecoderType::UfDecoder;
    }
    throw std::invalid_argument("Invalid decodinger type: " + status);
}
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,misc-include-cleaner)
NLOHMANN_JSON_SERIALIZE_ENUM(DecoderType, {{UfHeuristic, "UF_HEURISTIC"},
                                           {UfDecoder, "UF_DECODER"}})

class DecodingSimulator {
public:
    /**
     * Runs a simulation for the given decoder starting with the minimum physical error rate.
     * For each error rate there is a number of runs made to calculate the WER for the error rate.
     * The error rate between runs is increased by stepsize until the maximum physical error rate is reached.
     * It is assumed that the minimum physical error rate is smaller than the maximum and that
     * the field code is set in the given @decoder.
     * Results are written to the two output files specified.
     * @param rawDataOutputFilepath path to file to write raw output data in the form [physicalErrorRate:WER] to
     * @param statsOutputFilepath path to file to write statistics of decoding run to (@DecodingRunInformation)
     * @param minPhysicalErrRate starting physical error rate
     * @param maxPhysicalErrRate maximum physical error rate
     * @param physErrRateStepSize stepsize between error rates
     * @param nrRunsPerRate number of runs to average WER over
     * @param decoder
     */
    static void simulateWER(const std::string& rawDataOutputFilepath,
                            const std::string& statsOutputFilepath,
                            double             minPhysicalErrRate,
                            double             maxPhysicalErrRate,
                            std::size_t        nrRunsPerRate,
                            Code&              code,
                            double             perStepSize,
                            const DecoderType& decoderType); // code is field of decoder

    /**
     * Runs the specified number of decoding runs for each physical error rate on each code and
     * computes the average runtime (in ms) per code. For this to reflect runtime scaling reasonably
     * the codes should be from the same family (e.g. toric codes) with increasing size.
     * @param rawDataOutputFilepath
     * @param decodingInfoOutfilePath
     * @param physicalErrRates
     * @param nrRuns
     * @param codes
     */
    static void simulateAverageRuntime(const std::string& rawDataOutputFilepath,
                                       const std::string& decodingInfoOutfilePath,
                                       const double&      physicalErrRate,
                                       std::size_t        nrRuns,
                                       const std::string& codesPath,
                                       std::size_t        nrSamples,
                                       const DecoderType& decoderType);
};

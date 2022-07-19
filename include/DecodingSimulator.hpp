//
// Created by lucas on 21/06/22.
//

#ifndef QUNIONFIND_DECODINGSIMULATOR_HPP
#define QUNIONFIND_DECODINGSIMULATOR_HPP
#include "Decoder.hpp"
#include "ImprovedUFD.hpp"

#include <string>
#include <utility>

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
     * @param nrOfRunsPerErrRate number of runs to average WER over
     * @param decoder
     */
    static void simulateWER(const std::string& rawDataOutputFilepath,
                            const std::string& statsOutputFilepath,
                            double             minPhysicalErrRate,
                            double             maxPhysicalErrRate,
                            double             physErrRateStepSize,
                            std::size_t        nrOfRunsPerErrRate,
                            Decoder&           decoder); // code is field of decoder

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
                                       std::size_t        nrSamples);

private:
    static std::string generateOutFileName(const std::string& filepath);
};

#endif //QUNIONFIND_DECODINGSIMULATOR_HPP

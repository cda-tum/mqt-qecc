/**
 * Simulate average runtime for codes with growing nr of n for several physical err rates (err rates should only increase slope of curve)
 */
#include "Code.hpp"
#include "DecodingRunInformation.hpp"
#include "UFDecoder.hpp"
#include "Utils.hpp"

#include <cmath> /* pow */
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void runtime(const std::string& codename) {
    /**
     * ***************** Comment out accordingly *****************
     */
    //**** server:
    const std::string inPath = "in/toricCodes2/";
    //**** local:
    // const std::string outPath = "/home/luca/Documents/uf-simulations/runtime/original/";
    // const std::string inPath  = "/home/luca/Documents/codeRepos/qecc/examples/toricCodes2/";
    // ***************** config end *****************
    const std::size_t nrRuns    = 1000;
    const std::size_t nrSamples = 50;
    const double      per       = 0.01;
    auto              code      = Code(inPath + codename);
    const auto        coden     = code.getN();

    std::size_t samplesSum = 0;

    for (std::size_t i = 0; i < nrSamples; i++) {
        auto runsSum = 0U;
        for (std::size_t j = 0; j < nrRuns; j++) {
            auto decoder = UFDecoder();
            decoder.setCode(code);
            std::vector<bool> error;
            while (error.empty() || std::none_of(error.begin(), error.end(), [](bool c) { return c; })) {
                error = Utils::sampleErrorIidPauliNoise(coden, per);
            }
            auto syndrome = code.getXSyndrome(error);
            decoder.decode(syndrome);
            runsSum += decoder.result.decodingTime;
            decoder.reset();
        }
        std::cout << runsSum << std::endl;
        samplesSum += runsSum;
    }
    std::cout << codename << ":" << samplesSum / nrSamples << std::endl;
}

void decodingPerformance(const double per) {
    const std::string rootPath  = "qecc/examples/lp_(4,8)-[[1024,18,nan]]_hx.txt";
    const std::string rootPath2 = "qecc/examples/lp_(4,8)-[[1024,18,nan]]_hz.txt";
    const std::size_t codeK     = 18;
    const std::size_t nrRuns    = 1000;

    auto       code = HGPcode(rootPath, rootPath2, codeK);
    const auto n    = code.getN();

    std::size_t nrOfFailedRuns   = 0U;
    std::size_t nrSuccessfulRuns = 0U;
    std::size_t j                = 0;
    while (nrOfFailedRuns < 1 || j < nrRuns) {
        auto decoder = std::make_unique<UFDecoder>();
        decoder->setCode(code);
        auto error    = Utils::sampleErrorIidPauliNoise(n, per);
        auto syndrome = code.getXSyndrome(error);
        decoder->decode(syndrome);
        auto const&       decodingResult = decoder->result;
        std::vector<bool> residualErr    = decodingResult.estimBoolVector;
        Utils::computeResidualErr(error, residualErr);
        auto success = code.isXStabilizer(residualErr);
        if (!success) {
            nrOfFailedRuns++;
        } else {
            nrSuccessfulRuns++;
        }
        decoder->reset();
        j++;
    }
    // compute word error rate WER
    const auto logicalErrRate = 1 - static_cast<double>(nrSuccessfulRuns) / static_cast<double>(nrRuns);
    const auto wordErrRate    = 1.0 - std::pow(1 - logicalErrRate, (1.0 / codeK)); // rate of codewords for decoder does not give correct answer (fails or introduces logical operator)
    std::cout << "per:wer = " << per << ":" << wordErrRate << std::endl;
    std::cout.flush();
}

int main(int argc, char* argv[]) {         // NOLINT(clang-diagnostic-unused-parameter, bugprone-exception-escape,misc-unused-parameters)
    const double per = std::stod(argv[1]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    decodingPerformance(per);
}

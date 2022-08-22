/**
* Simulate average runtime for codes with growing nr of N for several physical err rates (err rates should only increase slope of curve)
*/
#include "Code.hpp"
#include "DecodingRunInformation.hpp"
#include "UFDecoder.hpp"
#include "UFHeuristic.hpp"
#include "Utils.hpp"

#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void runtime(const std::string& codeName) {
    /**
     * ***************** Comment out accordingly *****************
     */
    //**** server:
    //const std::string rootPath = "/home/berent/ufpaper/simulations/montecarlo/final/";
    //const std::string inPath   = rootPath + "in/toricCodes2/";
    //**** local:
    const std::string outPath = "/home/luca/Documents/uf-simulations/runtime/original/";
    const std::string inPath  = "/home/luca/Documents/codeRepos/qecc/examples/toricCodes2/";
    // ***************** config end *****************
    const std::size_t nrRuns = 100'000;
    const std::size_t nrSamples = 3;
    const double      per    = 0.02;
    auto       code    = Code(inPath + codeName);
    const auto codeN   = code.getN();

    std::size_t runsSum = 0;
    std::size_t samplesSum = 0;

    for (std::size_t i = 0; i < nrSamples; i++) {
        runsSum = 0;
        for (std::size_t j = 0; j < nrRuns; j++) {
            auto       decoder = UFDecoder();
            decoder.setCode(code);
            auto error    = Utils::sampleErrorIidPauliNoise(codeN, per);
            auto syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            runsSum += decoder.result.decodingTime;
            decoder.reset();
        }
        std::cout << runsSum << std::endl;
        samplesSum += runsSum;
    }
    std::cout << codeName << ":" << samplesSum/nrSamples  << std::endl;
    flint_cleanup();
}

void decodingPerformance(const double per) {
    const std::string rootPath   = "/home/luca/Documents/codeRepos/qecc/examples/lp_(4,8)-[[1024,18,nan]]_hz.txt"; // TODO adapt
    const std::size_t code_K     = 18;

    const std::size_t nrOfRunsPerRate    = 10;

    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;
    //    decodingResOutput << "{ \"runs\" : [ ";
    //    rawIntermediateOut << "{ ";

    auto        code = HGPcode(rootPath, code_K);
    const auto  K    = code.getK();
    const auto  N    = code.getN();

    std::size_t nrOfFailedRuns = 0U;
    for (std::size_t j = 0; j < nrOfRunsPerRate; j++) {
        std::cout << "run " << j << std::endl;
        auto* decoder = new UFHeuristic();
        decoder->setCode(code);
        decoder->setGrowth(GrowthVariant::ALL_COMPONENTS);
        auto error    = Utils::sampleErrorIidPauliNoise(N, per);
        auto syndrome = code.getSyndrome(error);
        decoder->decode(syndrome);
        auto const&       decodingResult = decoder->result;
        std::vector<bool> residualErr    = decodingResult.estimBoolVector;
        Utils::computeResidualErr(error, residualErr);
        auto                 success = code.isVectorStabilizer(residualErr);
        if (!success) {
            nrOfFailedRuns++;
        }
        decoder->reset();
        delete decoder;
    }
    //compute word error rate WER
    const auto blockErrRate    = static_cast<double>(nrOfFailedRuns) / static_cast<double>(nrOfRunsPerRate);
    const auto wordErrRate     = blockErrRate / static_cast<double>(K);                                                   // rate of codewords for decoder does not give correct answer (fails or introduces logical operator)
    std::cout << "per:wer = " << per << ":" << wordErrRate << std::endl;
    std::cout.flush();
    flint_cleanup();
}

int main(int argc, char* argv[]) {
    std::string codeName = argv[1];
    double      per          = std::stod(argv[1]);
    //runtime(codeName);
    decodingPerformance(per);
}
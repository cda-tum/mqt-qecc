//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "ImprovedUFD.hpp"
#include "OriginalUFD.hpp"

#include <bitset>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <locale>
#include <random>
using json = nlohmann::json;
class UnionFindSimulation: public testing::TestWithParam<std::string> {
};

TEST(UnionFindSimulation, TestSimulator) {
    std::string rawOut = "./testRawFile", testOut = "./testStatFile";
    double      minErate = 0.01, maxErate = 0.03, stepSize = 0.01;
    std::size_t runsPerRate = 2, runsPerCode = 2;
    auto        code = SteaneXCode();
    try {
        ImprovedUFD decoder;
        decoder.setCode(code);
        DecodingSimulator::simulateWER(rawOut, testOut, minErate, maxErate, stepSize, runsPerRate, decoder);
        std::vector<double> erRates = {minErate, maxErate};
        DecodingSimulator::simulateRuntime(rawOut, testOut, erRates, runsPerCode, decoder);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
    }
    EXPECT_TRUE(true);
}

/**
 * Simulate WER for growing number of physical err rate
 * Can also be used for threshold simulations
 */
TEST(UnionFindSimulation, EmpiricalEvaluationDecodingPerformance) {
    /**
     * ***************** Comment out accordingly *****************
     */
    //****server
    const std::string outpath    = "/home/berent/ufpaper/simulations/decodingPerfSim/final/out/";
    const std::string inCodePath = "/home/berent/ufpaper/simulations/decodingPerfSim/final/source/code/hgp_(4,8)-[[5408,18,26]]_hx.txt";
    const std::size_t code_K     = 18;
    //**** local
    //    const std::string outpath            = "/home/luca/Documents/uf-simulations/testrun/";
    //    const std::string inCodePath         = "/home/luca/Documents/codeRepos/qunionfind/examples/hgp_(4,7)-[[900,36,10]]_hx.txt";
    //    const std::size_t code_K             = 36;
    // ***************** configure end *****************

    const std::string outFilePath        = outpath + "results";
    const std::string dataFilePathInterm = outpath + "raw-interm";
    const std::string dataFilePath       = outpath + "raw-final";
    std::cout << "writing output to " << outpath << std::endl;
    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream decodingResOutput(outFilePath + timestamp + ".json");
    std::ofstream rawDataOutput(dataFilePath + timestamp + ".json");
    std::ofstream rawIntermediateOut(dataFilePathInterm + timestamp + ".json");
    rawIntermediateOut.rdbuf()->pubsetbuf(0, 0);
    rawDataOutput.rdbuf()->pubsetbuf(0, 0);

    /**
     * ***************** Comment out accordingly *****************
     */
    //**** Paper eval
    double            physicalErrRate    = 0.0001;
    const double      stepSize           = 0.0002;
    const double      maxPhysicalErrRate = 0.1;
    const std::size_t nrOfRuns           = std::floor(maxPhysicalErrRate / physicalErrRate);
    const std::size_t nrOfRunsPerRate    = 4096; // todo how deep to go?
    //**** tests
    //    const double      normalizationConstant = 100; //
    //    double            physicalErrRate       = 1.0 / normalizationConstant;
    //    double            stepSize              = 1.0 / normalizationConstant;
    //    const double      maxPhysicalErrRate    = 0.1;
    //    const std::size_t nrOfRuns              = std::floor(maxPhysicalErrRate / physicalErrRate);
    //    std::size_t       nrOfRunsPerRate       = 1;
    // ***************** configure end *****************

    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;
    decodingResOutput << "{ \"runs\" : [ ";
    rawIntermediateOut << "{ ";

    auto       code = HGPcode(inCodePath, code_K);
    const auto K    = code.getK();
    const auto N    = code.getN();
    ImprovedUFD decoder;
    for (std::size_t i = 0; i < nrOfRuns && physicalErrRate <= maxPhysicalErrRate; i++) {
        std::size_t nrOfFailedRuns = 0U;
        decodingResOutput << R"({ "run": { "physicalErrRate":)" << physicalErrRate << ", \"data\": [ ";
        for (size_t j = 0; j < nrOfRunsPerRate; j++) {
            std::cout << "run nr " << j << std::endl;
            decoder.setCode(code);
            auto        error    = Utils::sampleErrorIidPauliNoise(N, physicalErrRate);
            auto        syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            auto const&       decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            Utils::computeResidualErr(error, residualErr);
            auto                   success = code.isVectorStabilizer(residualErr);
            DecodingRunInformation stats;
            stats.result = decodingResult;

            if (success) {
                stats.status = SUCCESS;
            } else {
                stats.status = FAILURE;
                nrOfFailedRuns++;
            }
            stats.physicalErrR = physicalErrRate;
            stats.error        = error;
            stats.syndrome     = syndrome;
            stats.codeSize     = N;
            decodingResOutput << stats.to_json().dump(2U);
            if (j != nrOfRunsPerRate - 1) {
                decodingResOutput << ", ";
            }
            decoder.reset();
        }
        std::cout << "computing wer " << std::endl;
        //compute word error rate WER
        const auto blockErrRate    = static_cast<double>(nrOfFailedRuns) / static_cast<double>(nrOfRunsPerRate);
        const auto wordErrRate     = blockErrRate / static_cast<double>(K);                                                   // rate of codewords for decoder does not give correct answer (fails or introduces logical operator)
        const auto& [it, inserted] = wordErrRatePerPhysicalErrRate.try_emplace(std::to_string(physicalErrRate), wordErrRate); // to string for json parsing
        rawIntermediateOut << R"( ")" << it->first << R"(" )"
                           << ":" + std::to_string(it->second);
        // only for json output
        if (i != nrOfRuns - 1) {
            decodingResOutput << "]}},";
            rawIntermediateOut << ",";
        } else {
            decodingResOutput << "]}}";
            rawIntermediateOut << "}";
        }
        std::cout << "stepping to next PER " << std::endl;
        physicalErrRate += stepSize;
    }
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    decodingResOutput.close();
    rawDataOutput.close();
    rawIntermediateOut.close();
}

/**
 * Simulate average runtime for codes with growing nr of N for several physical err rates (err rates should only increase slope of curve)
 */
TEST(UnionFindSimulation, EmpiricalEvaluationDecoderRuntime) {
    /**
     * ***************** Comment out accordingly *****************
     */
    //**** server:
    //const std::string codeN   = "toric_(nan,nan)-[[1058,2,23]]_hx.txt";
    const std::string outPath = "/home/berent/ufpaper/simulations/montecarlo/final/out/";
    const std::string inPath  = "/home/berent/ufpaper/simulations/montecarlo/final/in/toricCodes/";
    //**** local:
   // const std::string outPath = "/home/luca/Documents/uf-simulations/runtime/repaired/";
    //const std::string inPath  = "/home/luca/Documents/codeRepos/qecc/examples/toricCodes2/";
    // ***************** config end *****************

    const std::string outFile         = outPath + "results_";
    const std::string runningDataFile = outPath + "raw-running_";
    const std::string finalDataFile   = outPath + "raw-final_";
    std::cout << "writing results to " << outPath << std::endl;
    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream dataOutStream(outFile + timestamp + ".json");
    std::ofstream intermediateRawOut(runningDataFile + timestamp + ".json"); // appends data after each step for intermediate results
    std::ofstream finalRawOut(finalDataFile + timestamp + ".json");          // single, final data dump at end
    intermediateRawOut.rdbuf()->pubsetbuf(0, 0);
    finalRawOut.rdbuf()->pubsetbuf(0, 0);

    /**
     * ***************** Basic parameters, comment out accordingly *****************
     */
    // Basic Parameter setup
    //**** paper eval:
    //const double      physErrRates[] = {0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05};
    const std::size_t nrOfTrials = 100'000;
    //****tests:
    const double physErrRate = 0.02;
    //const std::size_t nrOfTrials     = 1'0;
    // ***************** configure end *****************

    std::size_t                                                                    avgDecodingTimeAcc = 0U;
    double                                                                         avgDecTime         = 0.0;
    std::map<std::string, std::map<std::string, double, std::less<>>, std::less<>> dataPerRate;
    std::map<std::string, double, std::less<>>                                     tmp;
    std::vector<std::string>                                                       codePaths;
    std::cout << "reading codes " << std::endl;
    std::map<std::string, std::size_t, std::less<>> avgSampleRuns;
    DecodingRunInformation                          info;
    for (const auto& file: std::filesystem::directory_iterator(inPath)) {
        codePaths.emplace_back(file.path());
    }
    std::map<std::string, double, std::less<>> avgTimePerSizeData;
    std::cout << "Simulating physical err rate " << physErrRate << std::endl;
    try {
        for (const auto& currPath: codePaths) {
            std::cout << "next code : " << currPath << std::endl;
            avgDecodingTimeAcc = 0U;
            auto       code    = Code(currPath); // construct new for each trial
            const auto codeN   = code.getN();
            auto       decoder = ImprovedUFD();
            for (std::size_t j = 0; j < 10; j++) {              // 5 runs to compute average
                for (std::size_t i = 0; i < nrOfTrials; i++) { // nr of monte carlo samples
                    decoder.setCode(code);
                    std::cout << "run nr " << i << std::endl;
                    auto error = Utils::sampleErrorIidPauliNoise(codeN, physErrRate);
                    std::cout << "computing syndrome" << std::endl;
                    auto syndrome = code.getSyndrome(error);
                    std::cout << "starting decoding" << std::endl;
                    decoder.decode(syndrome);
                    std::cout << "decoding done " << std::endl;
                    auto const& decodingResult = decoder.result;
                    info.result                = decodingResult;
                    info.physicalErrR          = physErrRate;
                    info.codeSize              = codeN;
                    info.syndrome              = syndrome;
                    info.error                 = error;
                    std::cout << "summing times" << std::endl;
                    avgDecodingTimeAcc = avgDecodingTimeAcc + decodingResult.decodingTime;
                    if (avgDecodingTimeAcc > std::numeric_limits<std::size_t>::max()) {
                        throw QeccException("Accumulator too large");
                    }
                    nlohmann::json json = info.to_json();
                    dataOutStream << json.dump(2U) << ",";
                    dataOutStream.flush();
                    json = {};
                    info = {};
                    decoder.reset();
                }
                avgSampleRuns.try_emplace(std::to_string(j), avgDecodingTimeAcc);
                std::cout << currPath << " \"run\": " << j << ", \"time\":" << avgDecodingTimeAcc << std::endl;
            }
            json avgData = avgSampleRuns;
            std::cout << "interm:" << avgData.dump(2U) << std::endl;
            intermediateRawOut << "\"" << currPath << "\""
                               << "{";
            intermediateRawOut << avgData.dump(2U);
            intermediateRawOut << "},";
            intermediateRawOut.flush();
            avgSampleRuns = {};
        }
    } catch (std::exception& e) {
        std::cerr << "Exception occurred " << e.what() << std::endl;
    }
    flint_cleanup();
    finalRawOut.close();
    intermediateRawOut.close();
    dataOutStream.close();
}

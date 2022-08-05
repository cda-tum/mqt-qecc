/**
* Simulate average runtime for codes with growing nr of N for several physical err rates (err rates should only increase slope of curve)
*/
#include "Code.hpp"
#include "DecodingRunInformation.hpp"
#include "ImprovedUFD.hpp"
#include "OriginalUFD.hpp"
#include "Utils.hpp"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void runtime(const std::string& outpath, const std::string& codePath, const double per) {
    auto conden = codePath.substr(codePath.find("[[")+2, codePath.find("]]"));
    conden = conden.substr(0, conden.find(','));
    const std::string                               finalDataFile = outpath+conden+".txt";

    auto decoder = ImprovedUFD();
    auto       code  = Code(codePath);
    const auto codeN = code.getN();
    decoder.setCode(code);
    auto error    = Utils::sampleErrorIidPauliNoise(codeN, per);
    auto syndrome = code.getSyndrome(error);
    decoder.decode(syndrome);
    auto const& decodingResult = decoder.result;

    std::cout << codeN << ":" << decodingResult.decodingTime << std::endl;
    std::ofstream outfile;
    outfile.open(finalDataFile, std::ios_base::app);
    std::cout << "writing to " << finalDataFile << std::endl;
    outfile << decodingResult.decodingTime << std::endl;
    decoder.reset();
    outfile.flush();
    outfile.close();
}

void decodingPerformance(const double per) {
    std::cout << "performance" << std::endl;
    /**
     * ***************** Comment out accordingly *****************
     */
    //****server
    const std::string rootPath   = "/home/berent/ufpaper/simulations/decodingPerfSim/final/";
    const std::string outpath    = rootPath + "out/";
    const std::string inCodePath = rootPath + "source/code/lp_(4,8)-[[416,18,nan]]_hx.txt";
    const std::size_t code_K     = 18;
    //**** local
    //const std::string outpath    = "/home/luca/Documents/uf-simulations/final/"; //TODO adapt
    //const std::string inCodePath = "/home/luca/Documents/codeRepos/qecc/examples/lp_(4,8)-[[416,18,nan]]_hx.txt"; // TODO adapt
    //const std::size_t code_K     = 18;
    // ***************** configure end *****************

    const std::string outFilePath        = outpath + "results-" + std::to_string(per);
    const std::string dataFilePathInterm = outpath + "interm-" + std::to_string(per);
    const std::string dataFilePath       = outpath + "final-" + std::to_string(per);
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
    double            physicalErrRate    = per;
    const double      stepSize           = 0.00002;
    const double      maxPhysicalErrRate = 0.1;
    const std::size_t nrOfRunsPerRate    = 100'000;
    //**** tests
    // double            physicalErrRate    = 0.06;
    // double            stepSize           = 0.02;
    // const double      maxPhysicalErrRate = 0.1;
    // const std::size_t nrOfRuns           = std::floor(maxPhysicalErrRate / stepSize); // to avoid float type loop increment
    // std::size_t       nrOfRunsPerRate    = 2;
    // ***************** configure end *****************

    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;
    //    decodingResOutput << "{ \"runs\" : [ ";
    //    rawIntermediateOut << "{ ";

    auto        code = HGPcode(inCodePath, code_K);
    const auto  K    = code.getK();
    const auto  N    = code.getN();
    OriginalUFD decoder;
    //decoder.setGrowth(GrowthVariant::SINGLE_SMALLEST);
    decoder.setCode(code);
    // for (std::size_t i = 0; i < nrOfRuns && physicalErrRate <= maxPhysicalErrRate; i++) {
    std::size_t nrOfFailedRuns = 0U;
    //decodingResOutput << R"({ "run": { "physicalErrRate":)" << physicalErrRate << ", \"data\": [ ";
    for (size_t j = 0; j < nrOfRunsPerRate; j++) {
        auto error    = Utils::sampleErrorIidPauliNoise(N, physicalErrRate);
        auto syndrome = code.getSyndrome(error);
        decoder.decode(syndrome);
        auto const&       decodingResult = decoder.result;
        std::vector<bool> residualErr    = decodingResult.estimBoolVector;
        Utils::computeResidualErr(error, residualErr);
        auto                 success = code.isVectorStabilizer(residualErr);
        DecodingResultStatus status;
        if (success) {
            status = SUCCESS;
        } else {
            status = FAILURE;
            nrOfFailedRuns++;
        }
        //        json resj = DecodingRunInformation(physicalErrRate,
        //                                           code_K,
        //                                           error,
        //                                           syndrome,
        //                                           status,
        //                                           decodingResult)
        //                            .to_json();
        //        //decodingResOutput << resj.dump(2U);
        //            if (j != nrOfRunsPerRate - 1) {
        //                decodingResOutput << ", ";
        //            }
    }
    //compute word error rate WER
    const auto blockErrRate    = static_cast<double>(nrOfFailedRuns) / static_cast<double>(nrOfRunsPerRate);
    nrOfFailedRuns             = 0;
    const auto wordErrRate     = blockErrRate / static_cast<double>(K);                                                   // rate of codewords for decoder does not give correct answer (fails or introduces logical operator)
    const auto& [it, inserted] = wordErrRatePerPhysicalErrRate.try_emplace(std::to_string(physicalErrRate), wordErrRate); // to string for json parsing
    std::cout << "per:wer = " << it->first << ":" << it->second << std::endl;
    //        rawIntermediateOut << R"( ")" << it->first << R"(" )"
    //                           << ":" + std::to_string(it->second);
    //        rawIntermediateOut.flush();
    // only for json output
    //        if (i != nrOfRuns - 1) {
    //            decodingResOutput << "]}},";
    //            rawIntermediateOut << ",";
    //        } else {
    //            decodingResOutput << "]}}";
    //            rawIntermediateOut << "}";
    //        }
    //physicalErrRate += stepSize;
    //}
    decoder.reset();
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    decodingResOutput.close();
    rawDataOutput.close();
    rawIntermediateOut.close();
}

int main(int argc, char* argv[]) {
    std::string codeFileName = argv[3];
    std::string outpath      = argv[2];
    double      per          = std::stod(argv[1]);
    runtime(outpath, codeFileName, per);
    //double per = std::stod(argv[1]);
    //decodingPerformance(per);
}
#include "DecodingSimulator.hpp"

#include "DecodingRunInformation.hpp"
#include "OriginalUFD.hpp"
//
// Created by lucas on 21/06/22.
//
void DecodingSimulator::simulateWER(const std::string& rawDataOutputFilepath,
                                    const std::string& statsOutputFilepath,
                                    const double       minPhysicalErrRate,
                                    const double       maxPhysicalErrRate,
                                    const double       physErrRateStepSize,
                                    const size_t       nrRunsPerRate,
                                    Decoder&           decoder) {


    // Basic Parameter setup
    double              physicalErrRate = minPhysicalErrRate;
    double              stepSize        = physErrRateStepSize;
    const double        maxPhErrRate    = maxPhysicalErrRate;
    const auto          nrOfRuns        = static_cast<std::size_t>(std::floor(maxPhErrRate / minPhysicalErrRate));
    std::vector<double> pers(nrRunsPerRate + 1);
    double              perCnt = minPhysicalErrRate;

    for (size_t i = 0; i <= nrRunsPerRate; i++) {
        pers.emplace_back(perCnt);
        perCnt += stepSize;
    }


    const auto& code = decoder.getCode();
    const auto  K    = code->getK();
    const auto  N    = code->getN();

}

std::string generateOutFileName(const std::string& filepath) {
    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto timestamp = oss.str();
    return filepath + "-" + timestamp + ".json";
}


void simulatePer(const Decoder& decoder, const double per,const std::string& rawDataOutputFilepath,
                 const std::string& statsOutputFilepath){

   /* auto          jsonFileName = generateOutFileName(statsOutputFilepath);
    auto          dataFileName = generateOutFileName(rawDataOutputFilepath);
    std::ofstream statisticsOutstr(jsonFileName);
    std::ofstream rawDataOutput(dataFileName);

    std::cout << "Writing stats output to " << jsonFileName << std::endl;
    std::cout << "Writing raw data to " << dataFileName << std::endl;

    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;
    statisticsOutstr << "{ \"runs\" : [ ";
    statisticsOutstr << R"({ "run": { "physicalErrRate":)" << currPer << ", \"data\": [ ";

    auto nrOfFailedRuns = 0U;

    for (size_t j = 0; j < nrRunsPerRate; j++) {
        // reset decoder results
        decoder.reset();
        const auto error    = Utils::sampleErrorIidPauliNoise(N, currPer);
        auto       syndrome = code->getSyndrome(error);
        decoder.decode(syndrome);
        auto const&       decodingResult = decoder.result;
        std::vector<bool> residualErr    = decodingResult.estimBoolVector;
        Utils::computeResidualErr(error, residualErr);
        const auto success = code->isVectorStabilizer(residualErr);

        DecodingRunInformation stats;
        stats.result = decoder.result;
        if (success) {
            stats.status = SUCCESS;
            Utils::printGF2vector(residualErr);
        } else {
            stats.status = FAILURE;
            nrOfFailedRuns++;
        }
        statisticsOutstr << stats.to_json().dump(2U);
        if (j != nrRunsPerRate - 1) {
            statisticsOutstr << ", ";
        }
    }
    //compute word error rate WER
    const auto blockErrRate = static_cast<double>(nrOfFailedRuns) / static_cast<double>(nrRunsPerRate);
    const auto wordErrRate  = blockErrRate / static_cast<double>(K);                         // rate of codewords re decoder does not give correct answer (fails or introduces logical operator)
    wordErrRatePerPhysicalErrRate.try_emplace(std::to_string(physicalErrRate), wordErrRate); // to string for json parsing
    physicalErrRate += stepSize;

    statisticsOutstr << "}";
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    statisticsOutstr.close();
    rawDataOutput.close();*/
}


void DecodingSimulator::simulateAverageRuntime(const std::string& rawDataOutputFilepath,
                                               const std::string& decodingInfoOutfilePath,
                                               const double&      physicalErrRate,
                                               const std::size_t  nrRuns,
                                               const std::string& codesPath,
                                               const std::size_t  nrSamples) {
    std::cout << "writing raw data to " << rawDataOutputFilepath << std::endl;
    std::cout << "writing statistics to " << decodingInfoOutfilePath << std::endl;
    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream dataOutStream(generateOutFileName(decodingInfoOutfilePath));
    std::ofstream finalRawOut(generateOutFileName(rawDataOutputFilepath)); // single, final data dump at end
    finalRawOut.rdbuf()->pubsetbuf(0, 0);

    std::size_t                                                                    avgDecodingTimeAcc = 0U;
    std::map<std::string, std::map<std::string, double, std::less<>>, std::less<>> dataPerRate;
    std::map<std::string, double, std::less<>>                                     tmp;
    std::vector<std::string>                                                       codePaths{};
    std::cout << "reading codes " << std::endl;
    std::map<std::string, std::size_t, std::less<>>                        avgSampleRuns;
    std::map<std::string, std::map<std::string, std::size_t, std::less<>>> avgSampleRunsPerCode;
    auto                                                                   decoder = OriginalUFD();
    DecodingRunInformation                                                 info;
    for (const auto& file: std::filesystem::directory_iterator(codesPath)) {
        codePaths.emplace_back(file.path());
    }
    std::map<std::string, double, std::less<>> avgTimePerSizeData;
    try {
        for (const auto& currPath: codePaths) {
            avgDecodingTimeAcc = 0U;
            auto       code    = Code(currPath);
            const auto codeN   = code.getN();
            for (std::size_t j = 0; j < nrRuns; j++) {
                for (std::size_t i = 0; i < nrSamples; i++) {
                    decoder.setCode(code);
                    auto error    = Utils::sampleErrorIidPauliNoise(codeN, physicalErrRate);
                    auto syndrome = code.getSyndrome(error);
                    decoder.decode(syndrome);
                    auto const& decodingResult = decoder.result;
                    decoder.reset();
                    if (!decodingInfoOutfilePath.empty()) {
                        info.result       = decodingResult;
                        info.physicalErrR = physicalErrRate;
                        info.codeSize     = codeN;
                        info.syndrome     = syndrome;
                        info.error        = error;
                        info.print();
                    }
                    avgDecodingTimeAcc = avgDecodingTimeAcc + decodingResult.decodingTime;
                    if (avgDecodingTimeAcc > std::numeric_limits<std::size_t>::max()) {
                        throw QeccException("Accumulator too large");
                    }
                    decoder.reset();
                }
                auto average = avgDecodingTimeAcc / nrSamples;
                avgSampleRuns.try_emplace(std::to_string(j), average);
            }
            avgSampleRunsPerCode.try_emplace(currPath, avgSampleRuns);
            avgSampleRuns = {};
        }
    } catch (std::exception& e) {
        std::cerr << "Exception occurred " << e.what() << std::endl;
    }
    json j = avgSampleRunsPerCode;
    finalRawOut << j.dump(2U);
    flint_cleanup();
    finalRawOut.close();
    dataOutStream.close();
}

#include "DecodingSimulator.hpp"

#include "DecodingRunInformation.hpp"
//
// Created by lucas on 21/06/22.
//
void DecodingSimulator::simulateWER(const std::string& rawDataOutputFilepath,
                                    const std::string& statsOutputFilepath,
                                    const double       minPhysicalErrRate,
                                    const double       maxPhysicalErrRate,
                                    const double       physErrRateStepSize,
                                    const size_t       nrOfRunsPerErrRate,
                                    Decoder&           decoder) {
    auto          jsonFileName = generateOutFileName(statsOutputFilepath);
    auto          dataFileName = generateOutFileName(rawDataOutputFilepath);
    std::ofstream statisticsOutstr(jsonFileName);
    std::ofstream rawDataOutput(dataFileName);

    std::cout << "Writing stats output to " << jsonFileName << std::endl;
    std::cout << "Writing raw data to " << dataFileName << std::endl;

    // Basic Parameter setup
    double       physicalErrRate = minPhysicalErrRate;
    double       stepSize        = physErrRateStepSize;
    const double maxPhErrRate    = maxPhysicalErrRate;
    const auto   nrOfRuns        = static_cast<std::size_t>(std::floor(maxPhErrRate / minPhysicalErrRate));
    std::size_t  nrRunsPerRate   = nrOfRunsPerErrRate;

    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;
    statisticsOutstr << "{ \"runs\" : [ ";

    const auto& code = decoder.getCode();
    const auto  K    = code->getK();
    const auto  N    = code->getN();

    for (std::size_t i = 0; i < nrOfRuns && physicalErrRate <= maxPhErrRate; i++) {
        auto nrOfFailedRuns = 0U;
        statisticsOutstr << R"({ "run": { "physicalErrRate":)" << physicalErrRate << ", \"data\": [ ";

        for (size_t j = 0; j < nrRunsPerRate; j++) {
            // reset decoder results
            decoder.reset();
            const auto error    = Utils::sampleErrorIidPauliNoise(N, physicalErrRate);
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
        // only for json output
        if (i != nrOfRuns - 1) {
            statisticsOutstr << "]}},";
        } else {
            statisticsOutstr << "]}}";
        }
        physicalErrRate += stepSize;
    }
    statisticsOutstr << "}";
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    statisticsOutstr.close();
    rawDataOutput.close();
}

std::string DecodingSimulator::generateOutFileName(const std::string& filepath) {
    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto timestamp = oss.str();
    return filepath + "-" + timestamp + ".json";
}

void DecodingSimulator::simulateRuntime(const std::string&         rawDataOutputFilepath,
                                        const std::string&         decodingStatisticsOutputFilepath,
                                        const std::vector<double>& physicalErrRates,
                                        const std::size_t          nrRuns,
                                        Decoder&                   decoder) {
    auto          jsonFileName = generateOutFileName(decodingStatisticsOutputFilepath);
    auto          dataFileName = generateOutFileName(rawDataOutputFilepath);
    std::ofstream statisticsOutstr(jsonFileName);
    std::ofstream rawDataOutput(dataFileName);
    std::cout << "Writing stats output to " << jsonFileName << std::endl;
    std::cout << "Writing raw data to " << dataFileName << std::endl;

    // Basic Parameter setup
    const std::size_t                          nrOfTrials = nrRuns;
    double                                     avgDecTime = 0.0;
    std::map<std::string, double, std::less<>> avgDecodingTimePerSize;

    const auto& code = decoder.getCode();
    const auto  N    = code->getN();

    for (auto physErrRate: physicalErrRates) {
        auto avgDecodingTimeAcc = 0U;
        for (size_t i = 0; i < nrOfTrials; i++) {
            decoder.reset();
            const auto error    = Utils::sampleErrorIidPauliNoise(N, physErrRate);
            auto       syndrome = code->getSyndrome(error);
            decoder.decode(syndrome);
            const auto& decodingResult = decoder.result;
            avgDecodingTimeAcc         = avgDecodingTimeAcc + decodingResult.decodingTime;
        }
        avgDecTime = static_cast<double>(avgDecodingTimeAcc) / static_cast<double>(nrOfTrials);
        avgDecodingTimePerSize.try_emplace(std::to_string(N), avgDecTime);
    }
    json dataj = avgDecodingTimePerSize;
    rawDataOutput << dataj.dump(2U);
}

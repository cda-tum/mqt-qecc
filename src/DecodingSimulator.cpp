#include "DecodingSimulator.hpp"

#include "DecodingRunInformation.hpp"
//
// Created by lucas on 21/06/22.
//
void DecodingSimulator::simulateWER(std::string&   rawDataOutputFilepath,
                                    std::string&   statsOutputFilepath,
                                    const double   minPhysicalErrRate,
                                    const double   maxPhysicalErrRate,
                                    const double   physErrRateStepSize,
                                    const size_t   nrOfRunsPerErrRate,
                                    const Decoder& inDecoder) {
    auto          jsonFileName = generateOutFileName(statsOutputFilepath);
    auto          dataFileName = generateOutFileName(rawDataOutputFilepath);
    std::ofstream statisticsOutstr(jsonFileName);
    std::ofstream rawDataOutput(dataFileName);

    std::cout << "Writing stats output to " << jsonFileName << std::endl;
    std::cout << "Writing raw data to " << dataFileName << std::endl;

    // Basic Parameter setup
    double                        physicalErrRate = minPhysicalErrRate;
    double                        stepSize        = physErrRateStepSize;
    const double                  maxPhErrRate    = maxPhysicalErrRate;
    const size_t                  nrOfRuns        = std::floor(maxPhErrRate / minPhysicalErrRate);
    std::size_t                   nrRunsPerRate   = nrOfRunsPerErrRate; // todo how deep to go?
    std::size_t                   nrOfFailedRuns  = 0U;
    double                        blockErrRate    = 0.0;
    double                        wordErrRate     = 0.0;
    std::size_t                   K               = 0.0;
    std::map<std::string, double> wordErrRatePerPhysicalErrRate;
    statisticsOutstr << "{ \"runs\" : [ ";
    for (std::size_t i = 0; i < nrOfRuns && physicalErrRate <= maxPhErrRate; i++) {
        nrOfFailedRuns = 0U;
        blockErrRate   = 0.0;
        wordErrRate    = 0.0;
        statisticsOutstr << R"({ "run": { "physicalErrRate":)" << physicalErrRate << ", \"data\": [ ";

        for (size_t j = 0; j < nrRunsPerRate; j++) {
            Code code(inDecoder.getCode().Hz); // construct objects new to ensure clean state
            K = code.getK();
            Decoder decoder(inDecoder);
            auto    error    = Utils::sampleErrorIidPauliNoise(code.getN(), physicalErrRate);
            auto    syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            auto              decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            Utils::computeResidualErr(error, residualErr);
            auto success = code.isVectorStabilizer(residualErr);

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
        blockErrRate = (double)nrOfFailedRuns / (double)nrRunsPerRate;
        wordErrRate  = blockErrRate / (double)K;                                                            // rate of codewords re decoder does not give correct answer (fails or introduces logical operator)
        wordErrRatePerPhysicalErrRate.insert(std::make_pair(std::to_string(physicalErrRate), wordErrRate)); // to string for json parsing
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
                                        Code&                      code) {
    auto          jsonFileName = generateOutFileName(decodingStatisticsOutputFilepath);
    auto          dataFileName = generateOutFileName(rawDataOutputFilepath);
    std::ofstream statisticsOutstr(jsonFileName);
    std::ofstream rawDataOutput(dataFileName);
    std::cout << "Writing stats output to " << jsonFileName << std::endl;
    std::cout << "Writing raw data to " << dataFileName << std::endl;

    // Basic Parameter setup
    std::size_t                   avgDecodingTimeAcc = 0U;
    const std::size_t             nrOfTrials         = nrRuns;
    double                        avgDecTime         = 0.0;
    std::map<std::string, double> avgDecodingTimePerSize;
    for (auto physErrRate: physicalErrRates) {
        avgDecodingTimeAcc = 0U;
        for (size_t i = 0; i < nrOfTrials; i++) {
            auto        c = Code(code.Hz); // construct new for each trial
            ImprovedUFD decoder(c);
            auto        error    = Utils::sampleErrorIidPauliNoise(c.getN(), physErrRate);
            auto        syndrome = c.getSyndrome(error);
            decoder.decode(syndrome);
            auto decodingResult = decoder.result;
            avgDecodingTimeAcc  = avgDecodingTimeAcc + decodingResult.decodingTime;
        }
        avgDecTime = (double)avgDecodingTimeAcc / (double)nrOfTrials;
        avgDecodingTimePerSize.insert(std::make_pair<>(std::to_string(code.getN()), avgDecTime));
    }
    json dataj = avgDecodingTimePerSize;
    rawDataOutput << dataj.dump(2U);
}

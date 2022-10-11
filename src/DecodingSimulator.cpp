//
// Created by luca on 09/08/22.
//
#include "DecodingSimulator.hpp"

#include "DecodingRunInformation.hpp"
#include "UFDecoder.hpp"

std::string generateOutFileName(const std::string& filepath) {
    std::time_t now     = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char        buf[16] = {0};
    std::strftime(buf, sizeof(buf), "%Y-%m-%d", std::localtime(&now));

    std::ostringstream oss;
    oss << std::string(buf);
    auto timestamp = oss.str();
    return filepath + "-" + timestamp + ".json";
}

void DecodingSimulator::simulateWER(const std::string& rawDataOutputFilepath,
                                    const std::string& statsOutputFilepath,
                                    double             minPhysicalErrRate,
                                    double             maxPhysicalErrRate,
                                    std::size_t        nrRunsPerRate,
                                    Code&              code,
                                    const double       perStepSize,
                                    const DecoderType& decoderType) {
    bool                                       rawOut   = !rawDataOutputFilepath.empty();
    bool                                       statsOut = !statsOutputFilepath.empty();
    std::ofstream                              statisticsOutstr;
    std::ofstream                              rawDataOutput;
    std::map<std::string, double, std::less<>> wordErrRatePerPhysicalErrRate;

    if (rawOut) {
        auto dataFileName = generateOutFileName(rawDataOutputFilepath);
        rawDataOutput.open(dataFileName);
        std::cout << "Writing raw data to " << dataFileName << std::endl;
    }
    if (statsOut) {
        auto jsonFileName = generateOutFileName(statsOutputFilepath);
        statisticsOutstr.open(jsonFileName);
        std::cout << "Writing stats output to " << jsonFileName << std::endl;
        statisticsOutstr << "{ \"runs\" : [ ";
        statisticsOutstr << R"({ "run": { "physicalErrRate":)" << minPhysicalErrRate << ", \"data\": [ ";
    }

    auto nrOfFailedRuns = 0U;
    auto currPer        = minPhysicalErrRate;
    while (currPer < maxPhysicalErrRate) {
        nrOfFailedRuns = 0;
        for (std::size_t j = 0; j < nrRunsPerRate; j++) {
            Decoder* decoder;
            if (decoderType == DecoderType::UF_DECODER) {
                decoder = new UFDecoder();
            } else if (decoderType == DecoderType::UF_HEURISTIC) {
                decoder = new UFHeuristic();
            } else {
                throw QeccException("Invalid DecoderType, cannot simulate");
            }
            decoder->setCode(code);
            const auto error    = Utils::sampleErrorIidPauliNoise(code.getN(), currPer);
            const auto syndrome = decoder->getCode()->getXSyndrome(error);
            decoder->decode(syndrome);
            const auto& decodingResult = decoder->result;
            auto        residualErr    = decodingResult.estimBoolVector;
            Utils::computeResidualErr(error, residualErr);
            const auto success = decoder->getCode()->isXStabilizer(residualErr);

            DecodingRunInformation stats;
            stats.result = decoder->result;
            if (success) {
                stats.status = SUCCESS;
            } else {
                stats.status = FAILURE;
                nrOfFailedRuns++;
            }
            if (statsOut) {
                statisticsOutstr << stats.to_json().dump(2U);
                if (j != nrRunsPerRate - 1) {
                    statisticsOutstr << ", ";
                }
            }
        }
        // compute word error rate WER
        const auto blockErrRate = static_cast<double>(nrOfFailedRuns) / static_cast<double>(nrRunsPerRate);
        const auto wordErrRate  = blockErrRate / static_cast<double>(code.getK());       // rate of codewords re decoder does not give correct answer (fails or introduces logical operator)
        wordErrRatePerPhysicalErrRate.try_emplace(std::to_string(currPer), wordErrRate); // to string for json parsing

        currPer += perStepSize;
    }

    statisticsOutstr << "}";
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    statisticsOutstr.close();
    rawDataOutput.close();
    flint_cleanup();
}

void DecodingSimulator::simulateAverageRuntime(const std::string& rawDataOutputFilepath,
                                               const std::string& decodingInfoOutfilePath,
                                               const double&      physicalErrRate,
                                               const std::size_t  nrRuns,
                                               const std::string& codesPath,
                                               const std::size_t  nrSamples,
                                               const DecoderType& decoderType) {
    bool          rawOut  = !rawDataOutputFilepath.empty();
    bool          infoOut = !decodingInfoOutfilePath.empty();
    std::ofstream finalRawOut;
    std::ofstream dataOutStream;

    if (rawOut) {
        std::cout << "writing raw data to " << rawDataOutputFilepath << std::endl;
        finalRawOut.open(generateOutFileName(rawDataOutputFilepath));
        finalRawOut.rdbuf()->pubsetbuf(0, 0);
    }

    if (infoOut) {
        std::cout << "writing statistics to " << decodingInfoOutfilePath << std::endl;
        dataOutStream.open(generateOutFileName(decodingInfoOutfilePath));
        dataOutStream.rdbuf()->pubsetbuf(0, 0);
    }

    std::size_t                                                                    avgDecodingTimeAcc = 0U;
    std::map<std::string, std::map<std::string, double, std::less<>>, std::less<>> dataPerRate;
    std::map<std::string, double, std::less<>>                                     tmp;
    std::vector<std::string>                                                       codePaths{};
    std::map<std::string, std::size_t, std::less<>>                                avgSampleRuns;
    std::map<std::string, std::map<std::string, std::size_t, std::less<>>>         avgSampleRunsPerCode;

    DecodingRunInformation info;
    for (const auto& file : std::filesystem::directory_iterator(codesPath)) {
        codePaths.emplace_back(file.path());
    }
    std::map<std::string, double, std::less<>> avgTimePerSizeData;
    try {
        for (const auto& currPath : codePaths) {
            avgDecodingTimeAcc = 0U;
            auto       code    = Code(currPath);
            const auto codeN   = code.getN();
            for (std::size_t j = 0; j < nrRuns; j++) {
                for (std::size_t i = 0; i < nrSamples; i++) {
                    Decoder* decoder;
                    if (decoderType == DecoderType::UF_DECODER) {
                        decoder = new UFDecoder();
                    } else if (decoderType == DecoderType::UF_HEURISTIC) {
                        decoder = new UFHeuristic();
                    } else {
                        throw QeccException("Invalid DecoderType, cannot simulate");
                    }
                    decoder->setCode(code);
                    auto error    = Utils::sampleErrorIidPauliNoise(codeN, physicalErrRate);
                    auto syndrome = code.getXSyndrome(error);
                    decoder->decode(syndrome);
                    auto const& decodingResult = decoder->result;
                    if (infoOut) {
                        info.result       = decodingResult;
                        info.physicalErrR = physicalErrRate;
                        info.codeSize     = codeN;
                        info.syndrome     = syndrome;
                        info.error        = error;
                        info.print();
                    }
                    avgDecodingTimeAcc = avgDecodingTimeAcc + decodingResult.decodingTime;
                    decoder->reset();
                    delete decoder;
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
    if (rawOut) {
        json j = avgSampleRunsPerCode;
        finalRawOut << j.dump(2U);
        finalRawOut.close();
    }
    flint_cleanup();
    dataOutStream.close();
}

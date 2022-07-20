/**
* Simulate average runtime for codes with growing nr of N for several physical err rates (err rates should only increase slope of curve)
*/

int main() {
    //**** server:
    //const std::string codeN   = "toric_(nan,nan)-[[1058,2,23]]_hx.txt";
    const std::string outPath = "/home/berent/ufpaper/simulations/montecarlo/final/out/";
    const std::string inPath  = "/home/berent/ufpaper/simulations/montecarlo/final/in/toricCodes/";
    //**** local:
    //const std::string outPath = "/home/luca/Documents/uf-simulations/runtime/repaired/";
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
    const std::size_t nrOfTrials = 5;
    const std::size_t nrSamples  = 2;
    //****tests:
    const std::vector<double> physErrRates = {0.03, 0.06, 0.09, 0.1};
    //const std::size_t nrOfTrials     = 1'0;
    // ***************** configure end *****************

    std::size_t                                                                    avgDecodingTimeAcc = 0U;
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
    std::map<std::string, double, std::less<>> finalRawData;
    for (auto r: physErrRates) {
        finalRawOut << "{ \"" << r << "\":" << std::endl;
        std::cout << "Simulating physical err rate " << r << std::endl;
        try {
            for (const auto& currPath: codePaths) {
                std::cout << "next code : " << currPath << std::endl;
                avgDecodingTimeAcc = 0U;
                auto       code    = Code(currPath);
                const auto codeN   = code.getN();
                auto       decoder = ImprovedUFD();
                for (std::size_t j = 0; j < nrSamples; j++) {      // #sample runs to compute average
                    for (std::size_t i = 0; i < nrOfTrials; i++) { // nr of monte carlo samples
                        decoder.setCode(code);
                        auto error    = Utils::sampleErrorIidPauliNoise(codeN, r);
                        auto syndrome = code.getSyndrome(error);
                        decoder.decode(syndrome);
                        auto const& decodingResult = decoder.result;
                        info.result                = decodingResult;
                        info.physicalErrR          = r;
                        info.codeSize              = codeN;
                        info.syndrome              = syndrome;
                        info.error                 = error;
                        avgDecodingTimeAcc         = avgDecodingTimeAcc + decodingResult.decodingTime;
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
                intermediateRawOut.flush();
                avgSampleRuns = {};
                finalRawData.try_emplace(currPath, avgDecodingTimeAcc / nrSamples);
                avgDecodingTimeAcc = 0U;
            }
        } catch (std::exception& e) {
            std::cerr << "Exception occurred " << e.what() << std::endl;
            throw QeccException("Simulation failed with exception: ");
        }
        json fr = finalRawData;
        finalRawOut << fr.dump(2U) << "},";
        finalRawOut.flush();
    }
    finalRawOut.flush();
    finalRawOut.close();
    intermediateRawOut.close();
    dataOutStream.close();
    flint_cleanup();
}
#pragma once

#include "QeccException.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using json = nlohmann::basic_json<>;

/**
 * Used as return object in @Code
 */
struct CodeProperties {
    std::size_t n;
    std::size_t k;
    std::size_t d;
};

struct ParityCheckMatrix {
    std::unique_ptr<gf2Mat>                                   pcm;
    std::unordered_map<std::size_t, std::vector<std::size_t>> nbrCache;

    ParityCheckMatrix(const ParityCheckMatrix& m)          = delete;
    ParityCheckMatrix& operator=(const ParityCheckMatrix&) = delete;

    explicit ParityCheckMatrix(gf2Mat& mat) : pcm(std::make_unique<gf2Mat>(mat)) {}

    explicit ParityCheckMatrix(const std::string& filePath) {
        if (filePath.empty()) {
            throw QeccException("[PCM::ctor] - Cannot open pcm, filepath empty");
        }
        std::string   line;
        int           word = 0;
        std::ifstream inFile;
        gf2Mat        result;
        try {
            inFile.open(filePath);
            while (getline(inFile, line, '\n')) {
                gf2Vec             tempVec;
                std::istringstream instream(line);
                while (instream >> word) {
                    tempVec.push_back(static_cast<bool>(word));
                }
                result.emplace_back(tempVec);
            }
            pcm = std::make_unique<gf2Mat>(result);
        } catch (const std::exception& e) {
            std::cerr << "[PCM::ctor] - error opening file " << filePath << '\n';
            throw QeccException(e.what());
        }
        inFile.close();
    }

    /**
     * If H is nxm we have n checks and m bit nodes.
     * Indices of bit nodes range from 0 to m-1, and indices of check nodes from m to n-1
     * @param nodeIdx
     * @return a list of node indices of adjacent nodes
     */
    std::vector<std::size_t> getNbrs(const std::size_t& nodeIdx) {
        if (auto it = nbrCache.find(nodeIdx); it != nbrCache.end()) {
            return it->second;
        }
        if (pcm->empty() || pcm->front().empty()) {
            std::cerr << "error getting nbrs for node " << nodeIdx << '\n';
            throw QeccException("Cannot return neighbours, pcm empty");
        }
        const auto               nrChecks = pcm->size();
        const auto               nrBits   = pcm->front().size();
        std::vector<std::size_t> res;
        if (nodeIdx < nrBits) {
            for (std::size_t i = 0; i < nrChecks; i++) {
                if (pcm->at(i).at(nodeIdx)) {
                    res.emplace_back(nrBits + i);
                }
            }
        } else {
            for (std::size_t i = 0; i < nrBits; i++) {
                if (pcm->at(nodeIdx - nrBits).at(i)) {
                    res.emplace_back(i);
                }
            }
        }
        const auto& ins = nbrCache.try_emplace(nodeIdx, res);
        return ins.first->second;
    }
    [[nodiscard]] json to_json() const { // NOLINT(readability-identifier-naming)
        return json{
                {"pcm", *this->pcm}};
    }

    [[nodiscard]] std::string toString() const {
        return this->to_json().dump(2U);
    }
};

class Code {
private:
    std::unique_ptr<ParityCheckMatrix> hX;
    std::unique_ptr<ParityCheckMatrix> hZ;

public:
    std::size_t n = 0U;
    std::size_t k = 0U;
    std::size_t d = 0U;

    [[nodiscard]] const std::unique_ptr<ParityCheckMatrix>& gethX() const {
        return hX;
    }

    Code() = default;
    [[nodiscard]] const std::unique_ptr<ParityCheckMatrix>& gethZ() const {
        return hZ;
    }

    gf2Mat getHxMat() {
        return *this->hX->pcm;
    }

    gf2Mat getHzMat() {
        return *this->hZ->pcm;
    }

    void setHx(std::vector<std::vector<bool>>& hx) {
        hX = std::make_unique<ParityCheckMatrix>(hx);
    }

    void setHz(std::vector<std::vector<bool>>& hz) {
        hZ = std::make_unique<ParityCheckMatrix>(hz);
    }
    /*
     * Takes matrix hZ over GF(2) and constructs respective code for X errors with Z checks represented by hZ
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hz) : hZ(std::make_unique<ParityCheckMatrix>(hz)), n(hZ->pcm->front().size()) {
    }

    /*
     * Takes two pcms over GF(2) and constructs respective code
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hx, std::vector<std::vector<bool>>& hz) : hX(std::make_unique<ParityCheckMatrix>(hx)),
                                                                                            hZ(std::make_unique<ParityCheckMatrix>(hz)),
                                                                                            n(hZ->pcm->front().size()) {
    }

    /**
     * Constructs the X check part of a code given
     * @param pathToPcm
     */
    explicit Code(const std::string& pathToPcm) : hZ(std::make_unique<ParityCheckMatrix>(pathToPcm)) {
        if (hZ->pcm->empty() || hZ->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, hZ empty");
        }
        n = hZ->pcm->front().size();
    }

    explicit Code(const std::string& pathTohX, const std::string& pathTohZ) : hX(std::make_unique<ParityCheckMatrix>(pathTohX)), hZ(std::make_unique<ParityCheckMatrix>(pathTohZ)) {
        if (hZ->pcm->empty() || hZ->pcm->front().empty() || hX->pcm->empty() || hX->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, hX or hZ empty");
        }
        n = hZ->pcm->front().size();
        // todo hXhZ^T=0
        if (!hX->pcm || !hZ->pcm || hX->pcm->front().size() != hZ->pcm->front().size()) {
            throw QeccException("[Code::ctor] - hX and hZ dimensions do not match");
        }
    }

    [[nodiscard]] std::size_t getN() const {
        return n;
    }

    [[nodiscard]] std::size_t getK() const {
        return k;
    }

    /**
     * Returns the syndrome, given a X-error represented as binary vector
     * @param err
     * @return
     */
    [[nodiscard]] gf2Vec getXSyndrome(const gf2Vec& err) const {
        if (err.empty()) {
            throw QeccException("Cannot compute syndrome, err empty");
        }
        if (err.size() > this->getN()) {
            std::vector<bool> xerr;
            xerr.reserve(getN());
            std::vector<bool> zerr;
            zerr.reserve(getN());
            std::move(err.begin(), err.begin() + static_cast<std::int64_t>(getN()), std::back_inserter(xerr));
            std::move(err.begin() + static_cast<std::int64_t>(getN()), err.end(), std::back_inserter(zerr));
            return getSyndrome(xerr, zerr);
        }
        // per default X errs only
        gf2Vec syndr((*hZ->pcm).size(), false);
        Utils::rectMatrixMultiply(*hZ->pcm, err, syndr);
        return syndr;
    }

    /**
     * Returns the syndrome, given an error represented as two binary vector, one holding the X and one the Z part
     * @param err
     * @return
     */
    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& xerr, const gf2Vec& zerr) const {
        if (xerr.empty() || zerr.empty() || xerr.size() != getN() || zerr.size() != getN()) {
            throw QeccException("Cannot compute syndrome, err empty or wrong size");
        }

        gf2Vec xsyndr((*hZ->pcm).size(), false);
        Utils::rectMatrixMultiply(*hZ->pcm, xerr, xsyndr);

        gf2Vec zsyndr((*hX->pcm).size(), false);
        Utils::rectMatrixMultiply(*hX->pcm, zerr, zsyndr);
        gf2Vec res;
        res.reserve(xsyndr.size() + zsyndr.size());
        std::move(xsyndr.begin(), xsyndr.end(), std::back_inserter(res));
        std::move(zsyndr.begin(), zsyndr.end(), std::back_inserter(res));
        return res;
    }

    /**
     * Checks if the given vector is a X stabilizer of the code
     * @param est
     * @return
     */
    [[nodiscard]] bool isXStabilizer(const gf2Vec& est) const {
        if (!hX) {
            throw QeccException("hX not set, cannot check if vector is a stabilizer");
        }
        return Utils::isVectorInRowspace(*hX->pcm, est);
    }

    /**
     * Determines if the given vector represented as two components, X and Z
     * Is a stabilizer
     * @param xest
     * @param zest
     * @return
     */
    [[nodiscard]] bool isStabilizer(const gf2Vec& xest, const gf2Vec& zest) const {
        return Utils::isVectorInRowspace(*hX->pcm, xest) && Utils::isVectorInRowspace(*hZ->pcm, zest);
    }

    /**
     * Determines if the given vector which is assumed to hold a sympeltic representation of two components, X and Z
     * is a stabilizer
     * @param Xest
     * @param Zest
     * @return
     */
    [[nodiscard]] bool isStabilizer(const gf2Vec& est) const {
        if (std::all_of(est.begin(), est.end(), [](int i) { return !static_cast<bool>(i); })) { // trivial case, all 0 vector
            return true;
        }
        if (est.size() > getN()) {
            std::vector<bool> xEst;
            xEst.reserve(getN());
            std::vector<bool> zEst;
            zEst.reserve(getN());
            std::move(est.begin(), est.begin() + static_cast<std::int64_t>(est.size()) / 2, std::back_inserter(xEst));
            std::move(est.begin() + static_cast<std::int64_t>(est.size()) / 2, est.end(), std::back_inserter(zEst));
            return Utils::isVectorInRowspace(*hX->pcm, xEst) && Utils::isVectorInRowspace(*hZ->pcm, zEst);
        }
        return Utils::isVectorInRowspace(*hX->pcm, est);
    }

    [[nodiscard]] CodeProperties getProperties() const {
        CodeProperties res{n, getK(), d};
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, const Code& c) {
        auto   nrChecks = c.hZ->pcm->size();
        auto   nrData   = c.hZ->pcm->front().size();
        auto   dim      = nrChecks + nrData;
        gf2Mat res(dim);
        os << "hZ:\n";
        for (size_t i = 0; i < dim; i++) {
            gf2Vec row(dim);
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    row.at(nrData + j) = c.hZ->pcm->at(j).at(i);
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    row.at(j) = c.hZ->pcm->at(i - nrData).at(j);
                }
            }
            res.at(i) = row;
        }
        if (c.gethX()) {
            os << Utils::getStringFrom(res) << "hX:\n";
            for (size_t i = 0; i < dim; i++) {
                gf2Vec row(dim);
                if (i < dim - nrChecks) {
                    for (size_t j = 0; j < nrChecks; j++) {
                        row.at(nrData + j) = c.hX->pcm->at(j).at(i);
                    }
                } else {
                    for (size_t j = 0; j < nrData; j++) {
                        row.at(j) = c.hX->pcm->at(i - nrData).at(j);
                    }
                }
                res.at(i) = row;
            }
        }
        return os;
    }
    [[nodiscard]] json to_json() const { // NOLINT(readability-identifier-naming)
        return json{
                {"hZ", hZ ? hZ->to_json() : ""},
                {"hX", hX ? hX->to_json() : ""},
                {"n", n},
                {"k", k},
                {"d", d}};
    }

    void from_json(const json& j) { // NOLINT(readability-identifier-naming)
        j.at("n").get_to(n);
        j.at("k").get_to(k);
        j.at("d").get_to(d);
    }

    [[nodiscard]] std::string toString() const {
        return this->to_json().dump(2U);
    }
};

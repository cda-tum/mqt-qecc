//
// Created by lucas on 26/04/2022.
//

#ifndef QUNIONFIND_CODE_HPP
#define QUNIONFIND_CODE_HPP
#include "QeccException.hpp"
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
using json = nlohmann::json;

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
    std::unordered_map<std::size_t, std::vector<std::size_t>> nbrCache{};

    ParityCheckMatrix(const ParityCheckMatrix& m)          = delete;
    ParityCheckMatrix& operator=(const ParityCheckMatrix&) = delete;

    explicit ParityCheckMatrix(gf2Mat& pcm) : pcm(std::make_unique<gf2Mat>(pcm)) {}

    explicit ParityCheckMatrix(const std::string& filePath) {
        if (filePath.empty()) {
            throw QeccException("[PCM::ctor] - Cannot open pcm, filepath empty");
        }
        std::string   line;
        int           word;
        std::ifstream inFile;
        gf2Mat        result;
        try {
            inFile.open(filePath);
            while (getline(inFile, line, '\n')) {
                gf2Vec             tempVec;
                std::istringstream instream(line);
                while (instream >> word) {
                    tempVec.push_back(word);
                }
                result.emplace_back(tempVec);
            }
            pcm = std::make_unique<gf2Mat>(result);
        } catch (const std::exception& e) {
            std::cerr << "[PCM::ctor] - error opening file " << filePath << std::endl;
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
        } else {
            if (pcm->empty() || pcm->front().empty()) {
                std::cerr << "error getting nbrs for node " << nodeIdx << std::endl;
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
            const auto [nbrIt, inserted] = nbrCache.try_emplace(nodeIdx, res);
            return nbrIt->second;
        }
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
    std::unique_ptr<ParityCheckMatrix> Hx;
    std::unique_ptr<ParityCheckMatrix> Hz;

public:
    std::size_t n = 0U;
    std::size_t k = 0U;
    std::size_t d = 0U;

    [[nodiscard]] const std::unique_ptr<ParityCheckMatrix>& getHx() const {
        return Hx;
    }

    Code() = default;
    [[nodiscard]] const std::unique_ptr<ParityCheckMatrix>& getHz() const {
        return Hz;
    }

    gf2Mat getHxMat() {
        return *this->Hx->pcm;
    }

    gf2Mat getHzMat() {
        return *this->Hz->pcm;
    }

    void setHx(std::vector<std::vector<bool>>& hx) {
        Hx = std::make_unique<ParityCheckMatrix>(hx);
    }

    void setHz(std::vector<std::vector<bool>>& hz) {
        Hz = std::make_unique<ParityCheckMatrix>(hz);
    }
    /*
     * Takes matrix Hz over GF(2) and constructs respective code for X errors with Z checks represented by Hz
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hz) : Hz(std::make_unique<ParityCheckMatrix>(hz)) {
        n = Hz->pcm->front().size();
    }

    /*
     * Takes two pcms over GF(2) and constructs respective code
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hx, std::vector<std::vector<bool>>& hz) : Hx(std::make_unique<ParityCheckMatrix>(hx)), Hz(std::make_unique<ParityCheckMatrix>(hz)) {
        n = Hz->pcm->front().size();
    }

    /**
     * Constructs the X check part of a code given
     * @param pathToPcm
     */
    explicit Code(const std::string& pathToPcm) : Hz(std::make_unique<ParityCheckMatrix>(pathToPcm)) {
        if (Hz->pcm->empty() || Hz->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, Hz empty");
        }
        n = Hz->pcm->front().size();
    }

    explicit Code(const std::string& pathToHx, const std::string& pathToHz) : Hx(std::make_unique<ParityCheckMatrix>(pathToHx)), Hz(std::make_unique<ParityCheckMatrix>(pathToHz)) {
        if (Hz->pcm->empty() || Hz->pcm->front().empty() || Hx->pcm->empty() || Hx->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, Hx or Hz empty");
        }
        n = Hz->pcm->front().size();
        // todo HxHz^T=0
        if (!Hx->pcm || !Hz->pcm || Hx->pcm->front().size() != Hz->pcm->front().size()) {
            throw QeccException("[Code::ctor] - Hx and Hz dimensions do not match");
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
            throw QeccException("Cannot compute syndrome, err empy");
        } else if (err.size() > this->getN()) {
            std::vector<bool> xerr;
            xerr.reserve(getN());
            std::vector<bool> zerr;
            zerr.reserve(getN());
            std::move(err.begin(), err.begin() + getN(), std::back_inserter(xerr));
            std::move(err.begin() + getN(), err.end(), std::back_inserter(zerr));
            return getSyndrome(xerr, zerr);
        } else { // per defalut X errs only
            gf2Vec syndr((*Hz->pcm).size(), false);
            Utils::rectMatrixMultiply(*Hz->pcm, err, syndr);
            return syndr;
        }
    }

    /**
     * Returns the syndrome, given an error represented as two binary vector, one holding the X and one the Z part
     * @param err
     * @return
     */
    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& xerr, const gf2Vec& zerr) const {
        if (xerr.empty() || zerr.empty() || xerr.size() != getN() || zerr.size() != getN()) {
            throw QeccException("Cannot compute syndrome, err empy or wrong size");
        }

        gf2Vec xsyndr((*Hz->pcm).size(), false);
        Utils::rectMatrixMultiply(*Hz->pcm, xerr, xsyndr);

        gf2Vec zsyndr((*Hx->pcm).size(), false);
        Utils::rectMatrixMultiply(*Hx->pcm, zerr, zsyndr);
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
        if (!Hx) {
            throw QeccException("Hx not set, cannot check if vector is a stabilizer");
        }
        return Utils::isVectorInRowspace(*Hx->pcm, est);
    }

    /**
     * Determines if the given vector represented as two components, X and Z
     * Is a stabilizer
     * @param xest
     * @param zest
     * @return
     */
    [[nodiscard]] bool isStabilizer(const gf2Vec& xest, const gf2Vec& zest) const {
        return Utils::isVectorInRowspace(*Hx->pcm, xest) && Utils::isVectorInRowspace(*Hz->pcm, zest);
    }

    /**
     * Determines if the given vector which is assumed to hold a sympeltic representation of two components, X and Z
     * is a stabilizer
     * @param Xest
     * @param Zest
     * @return
     */
    [[nodiscard]] bool isStabilizer(const gf2Vec& est) const {
        if (std::all_of(est.begin(), est.end(), [](int i) { return !i; })) { // trivial case, all 0 vector
            return true;
        }
        if (est.size() > getN()) {
            std::vector<bool> xEst;
            xEst.reserve(getN());
            std::vector<bool> zEst;
            zEst.reserve(getN());
            std::move(est.begin(), est.begin() + (est.size()) / 2, std::back_inserter(xEst));
            std::move(est.begin() + (est.size()) / 2, est.end(), std::back_inserter(zEst));
            return Utils::isVectorInRowspace(*Hx->pcm, xEst) && Utils::isVectorInRowspace(*Hz->pcm, zEst);
        } else {
            return Utils::isVectorInRowspace(*Hx->pcm, est);
        }
    }

    [[nodiscard]] CodeProperties getProperties() const {
        CodeProperties res{.n = n, .k = getK(), .d = d};
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, const Code& c) {
        auto   nrChecks = c.Hz->pcm->size();
        auto   nrData   = c.Hz->pcm->front().size();
        auto   dim      = nrChecks + nrData;
        gf2Mat res(dim);
        os << "Hz: " << std::endl;
        for (size_t i = 0; i < dim; i++) {
            gf2Vec row(dim);
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    row.at(nrData + j) = c.Hz->pcm->at(j).at(i);
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    row.at(j) = c.Hz->pcm->at(i - nrData).at(j);
                }
            }
            res.at(i) = row;
        }
        if (c.getHx()) {
            os << Utils::getStringFrom(res) << "Hx: " << std::endl;
            for (size_t i = 0; i < dim; i++) {
                gf2Vec row(dim);
                if (i < dim - nrChecks) {
                    for (size_t j = 0; j < nrChecks; j++) {
                        row.at(nrData + j) = c.Hx->pcm->at(j).at(i);
                    }
                } else {
                    for (size_t j = 0; j < nrData; j++) {
                        row.at(j) = c.Hx->pcm->at(i - nrData).at(j);
                    }
                }
                res.at(i) = row;
            }
        }
        return os;
    }
    [[nodiscard]] json to_json() const { // NOLINT(readability-identifier-naming)
        return json{
                {"Hz", Hz ? Hz->to_json() : ""},
                {"Hx", Hx ? Hx->to_json() : ""},
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
#endif // QUNIONFIND_CODE_HPP

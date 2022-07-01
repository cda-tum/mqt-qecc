//
// Created by lucas on 26/04/2022.
//

#ifndef QUNIONFIND_CODE_HPP
#define QUNIONFIND_CODE_HPP
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <unordered_map>
#include <utility>
#include <vector>

struct CodeProperties{
    std::size_t n;
    std::size_t k;
    std::size_t d;
};

struct ParityCheckMatrix {
    explicit ParityCheckMatrix(gf2Mat pcm):
        pcm(std::move(pcm)) {}
    const gf2Mat pcm;

    std::unordered_map<std::size_t, std::vector<std::size_t>> nbrCache;

    /**
     * If H is nxm we have n checks and m bit nodes.
     * Indices of bit nodes range from 0 to m-1, and indices of check nodes from m to n-1
     * @param nodeIdx
     * @return a list of node indices of adjacent nodes
     */
    std::vector<std::size_t> getNbrs(const std::size_t& nodeIdx) {
        if (nbrCache.contains(nodeIdx)) {
            return nbrCache.at(nodeIdx);
        } else {
            auto                     nrChecks = pcm.size();
            auto                     nrBits   = pcm.at(0).size();
            std::vector<std::size_t> res;
            if (nodeIdx < nrBits) {
                for (std::size_t i = 0; i < nrChecks; i++) {
                    if (pcm.at(i).at(nodeIdx)) {
                        res.emplace_back(nrBits + i);
                    }
                }
            } else {
                for (std::size_t i = 0; i < nrBits; i++) {
                    if (pcm.at(nodeIdx - nrBits).at(i)) {
                        res.emplace_back(i);
                    }
                }
            }
            nbrCache.insert(std::make_pair(nodeIdx, res));
            return res;
        }
    }
};
/**
 * A graph representation for convenience using adjacency lists
 */
struct TannerGraph {
    std::vector<std::vector<std::size_t>> adjList;

    std::vector<std::size_t> getNeighbours(const std::size_t& node) {
        return adjList.at(node);
    }

    friend std::ostream& operator<<(std::ostream& os, TannerGraph const& c) {
        for (const auto& i: c.adjList) {
            os << ": | ";
            for (auto j: i) {
                os << j << " ";
            }
            os << "|";
            os << std::endl;
        }
        return os;
    }
};
/**
 * Considers X errors with Z checks only.
 * Z errors are corrected completely analogously in a symmetric way for Z and X.
 */
class Code {
public:
    TannerGraph       tannerGraph;
    ParityCheckMatrix Hz;
    std::size_t       K = 0U;
    std::size_t       N = 0U;

    /*
     * Takes matrix Hz over GF(2) and constructs respective code for X errors with Z checks represented by Hz
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(ParityCheckMatrix hz):
        Hz(std::move(hz)) {
        auto                                  nrChecks = Hz.pcm.size();
        auto                                  nrData   = Hz.pcm.at(0).size();
        std::size_t                           dim      = nrChecks + nrData;
        std::vector<std::vector<std::size_t>> adjMat(dim);

        for (size_t i = 0; i < dim; i++) {
            std::vector<std::size_t> row;
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    auto val = Hz.pcm.at(j).at(i);
                    if (val) {
                        row.emplace_back(nrData + j);
                    }
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    auto val = Hz.pcm.at(i - nrData).at(j);
                    if (val) {
                        row.emplace_back(j);
                    }
                }
            }
            adjMat.at(i) = row;
        }
        tannerGraph.adjList = adjMat;
        N                   = nrData;
    }

    [[nodiscard]] std::size_t getN() const {
        return N;
    }

    [[nodiscard]] std::size_t getK() const {
        return K;
    }

    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& err) const {
        gf2Mat errMat(err.size());
        for (size_t i = 0; i < err.size(); i++) {
            errMat.at(i) = gf2Vec{err.at(i)}; //transpose
        }
        auto res = Utils::rectMatrixMultiply(Hz.pcm, errMat);
        if (!res.empty()) {
            gf2Vec rres(Hz.pcm.size());
            for (size_t i = 0; i < rres.size(); i++) {
                rres.at(i) = res.at(i).at(0); // transpose back
            }
            return rres;
        } else {
            return gf2Vec{};
        }
    }

    [[nodiscard]] bool isVectorStabilizer(const gf2Vec& est) const {
        return Utils::isVectorInRowspace(Hz.pcm, est);
    }

    CodeProperties getProperties() const{
        CodeProperties res;
        res.n = this->N;
        res.k = this->getK();
        //todo res.d
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, const Code& c) {
        auto   nrChecks = c.Hz.pcm.size();
        auto   nrData   = c.Hz.pcm.at(0).size();
        auto   dim      = nrChecks + nrData;
        gf2Mat res(dim);

        for (size_t i = 0; i < dim; i++) {
            gf2Vec row(dim);
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    row.at(nrData + j) = c.Hz.pcm.at(j).at(i);
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    row.at(j) = c.Hz.pcm.at(i - nrData).at(j);
                }
            }
            res.at(i) = row;
        }
        return os << Utils::getStringFrom(res);
    }
};
#endif //QUNIONFIND_CODE_HPP

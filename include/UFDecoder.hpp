#pragma once

#include "Code.hpp"
#include "Decoder.hpp"

#include <cstddef>
#include <memory>
#include <unordered_set>
#include <vector>

class UFDecoder : public Decoder {
public:
    using Decoder::Decoder;
    void decode(const std::vector<bool>& syndrome) override;
    void reset() override;

private:
    void                                                       doDecode(const std::vector<bool>& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm);
    [[nodiscard]] bool                                         isValidComponent(const std::unordered_set<std::size_t>& nodeSet, const std::unordered_set<std::size_t>& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm) const;
    bool                                                       containsInvalidComponents(const std::unordered_set<std::size_t>& nodeSet, const std::unordered_set<std::size_t>& syndrome,
                                                                                         std::vector<std::unordered_set<std::size_t>>& invalidComps, const std::unique_ptr<ParityCheckMatrix>& pcm) const;
    [[nodiscard]] std::vector<std::size_t>                     computeInteriorBitNodes(const std::unordered_set<std::size_t>& nodeSet) const;
    [[nodiscard]] std::unordered_set<std::size_t>              getEstimateForComponent(const std::unordered_set<std::size_t>& nodeSet, const std::unordered_set<std::size_t>& syndrome,
                                                                                       const std::unique_ptr<ParityCheckMatrix>& pcm) const;
    void                                                       standardGrowth(std::unordered_set<std::size_t>& comps);
    void                                                       singleClusterSmallestFirstGrowth(std::unordered_set<std::size_t>& nodeSet);
    void                                                       singleClusterRandomFirstGrowth(std::unordered_set<std::size_t>& nodeSet);
    void                                                       singleQubitRandomFirstGrowth(std::unordered_set<std::size_t>& comps);
    [[nodiscard]] std::vector<std::unordered_set<std::size_t>> getConnectedComps(const std::unordered_set<std::size_t>& nodes) const;
};

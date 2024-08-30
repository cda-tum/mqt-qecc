#pragma once

#include "Code.hpp"
#include "Decoder.hpp"
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace std {
template <>
struct hash<std::pair<std::size_t, std::size_t>> {
    std::size_t operator()(const std::pair<std::size_t, std::size_t>& k) const {
        std::size_t h = 17;
        h             = h * 31 + std::hash<std::size_t>()(k.first);
        h             = h * 31 + std::hash<std::size_t>()(k.second);
        return h;
    }
};

} // namespace std

class UFHeuristic : public Decoder {
public:
    using Decoder::Decoder;
    void decode(const gf2Vec& syndrome) override;
    void reset() override;

private:
    // do not call.at only getNodeFromIdx()
    std::unordered_map<std::size_t, std::unique_ptr<TreeNode>> nodeMap;
    TreeNode*                                                  getNodeFromIdx(std::size_t idx);
    void                                                       standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                              std::unordered_map<std::size_t, bool>& presentMap, const std::unordered_set<std::size_t>& components, const std::unique_ptr<ParityCheckMatrix>& pcm);
    void                                                       singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                              std::unordered_map<std::size_t, bool>& presentMap, const std::unordered_set<std::size_t>& components, const std::unique_ptr<ParityCheckMatrix>& pcm);
    void                                                       singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                                std::unordered_map<std::size_t, bool>& presentMap, const std::unordered_set<std::size_t>& components, const std::unique_ptr<ParityCheckMatrix>& pcm);
    bool                                                       isValidComponent(const std::size_t& compId, const std::unique_ptr<ParityCheckMatrix>& pcm);
    std::vector<std::size_t>                                   erasureDecoder(std::unordered_set<std::size_t>& erasure, std::unordered_set<std::size_t>& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm);
    void                                                       extractValidComponents(std::unordered_set<std::size_t>& invalidComponents, std::unordered_set<std::size_t>& validComponents, const std::unique_ptr<ParityCheckMatrix>& pcm);
    std::unordered_set<std::size_t>                            computeInitTreeComponents(const gf2Vec& syndrome);
    void                                                       doDecoding(const gf2Vec& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm);
};

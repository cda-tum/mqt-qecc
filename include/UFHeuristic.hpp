//
// Created by lucas on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUFD_HPP
#define QUNIONFIND_IMPROVEDUFD_HPP
#include "Decoder.hpp"

#include <unordered_set>
namespace std {
    template<>
    struct hash<std::pair<std::size_t, std::size_t>> {
        std::size_t operator()(const std::pair<std::size_t, std::size_t>& k) const {
            std::size_t hash = 17;
            hash             = hash * 31 + std::hash<std::size_t>()(k.first);
            hash             = hash * 31 + std::hash<std::size_t>()(k.second);
            return hash;
        }
    };

} // namespace std

class UFHeuristic: public Decoder {
public:
    using Decoder::Decoder;
    void decode(const gf2Vec & syndrome) override;
    void reset() override;

private:
    // do not call.at only getNodeFromIdx()
    std::unordered_map<std::size_t, std::unique_ptr<TreeNode>> nodeMap{};
    TreeNode*                                    getNodeFromIdx(std::size_t idx);
    void                                                       standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                              std::unordered_map<std::size_t, bool>& presentMap, const std::vector<std::size_t>& components);
    void                                                       singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                              std::unordered_map<std::size_t, bool>& presentMap, const std::vector<std::size_t>& components);
    void                                                       singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                                std::unordered_map<std::size_t, bool>& presentMap, const std::vector<std::size_t>& components);
    bool                                                       isValidComponent(const std::size_t& compId);
    std::vector<std::size_t>                            erasureDecoder(std::vector<std::size_t>& erasure, std::vector<std::size_t>& syndrome);
    void                                                       extractValidComponents(std::vector<std::size_t>& invalidComponents, std::vector<std::size_t>& erasure);
    std::vector<std::size_t>                            computeInitTreeComponents(const gf2Vec& syndrome);
};
#endif //QUNIONFIND_IMPROVEDUFD_HPP

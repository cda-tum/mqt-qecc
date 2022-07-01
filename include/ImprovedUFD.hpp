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

template<>
struct std::hash<TreeNode> {
    size_t operator()(const TreeNode& k) const {
        // Compute individual hash values for first, second and third
        // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res        = res * 31 + hash<std::size_t>()(k.vertexIdx);
        return res;
    }
};

class ImprovedUFD: virtual public Decoder {
public:
    explicit ImprovedUFD(const Code& code):
        Decoder(code){};
    void decode(std::vector<bool>& syndrome) override;

private:
    std::unordered_map<std::size_t, std::shared_ptr<TreeNode>>    nodeMap;
    std::shared_ptr<TreeNode> getNodeFromIdx(const std::size_t idx);
    void                                          standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                 std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components);
    void                                          singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                 std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components);
    void                                          singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                                   std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components);
    bool                                          isValidComponent(const std::shared_ptr<TreeNode>& component);
    std::unordered_set<std::size_t>               erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::unordered_set<std::shared_ptr<TreeNode>>& syndrome);
    void                                          extractValidComponents(std::unordered_set<std::shared_ptr<TreeNode>>& invalidComponents, std::vector<std::shared_ptr<TreeNode>>& erasure);
    std::unordered_set<size_t>                    peeling(std::vector<std::shared_ptr<TreeNode>>& erasure, std::unordered_set<std::shared_ptr<TreeNode>>& syndrome);
    std::unordered_set<std::shared_ptr<TreeNode>> computeInitTreeComponents(const std::vector<bool>& syndrome);
};
#endif //QUNIONFIND_IMPROVEDUFD_HPP

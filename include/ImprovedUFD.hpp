//
// Created by luca on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUFD_HPP
#define QUNIONFIND_IMPROVEDUFD_HPP
#include "Decoder.hpp"
class ImprovedUFD: virtual public Decoder {
public:
    explicit ImprovedUFD(Code& code):
        Decoder(code){};
    void decode(std::vector<bool>& syndrome) override;

private:
    void                                standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                       std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    void                                singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                       std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    void                                singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                                         std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    bool                                isValidComponent(const std::shared_ptr<TreeNode>& component);
    std::vector<std::size_t>            erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
    void                                extractValidComponents(std::set<std::shared_ptr<TreeNode>>& components, std::vector<std::shared_ptr<TreeNode>>& erasure);
    std::vector<size_t>                 peelingDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
    std::set<std::shared_ptr<TreeNode>> computeInitTreeComponents(const std::vector<bool>& syndrome);
};
#endif //QUNIONFIND_IMPROVEDUFD_HPP

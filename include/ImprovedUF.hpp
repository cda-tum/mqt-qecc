//
// Created by luca on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUF_HPP
#define QUNIONFIND_IMPROVEDUF_HPP
#include "Decoder.hpp"
class ImprovedUF: virtual public Decoder {
public:
    explicit ImprovedUF(Code& code):
        Decoder(code) {};
    void decode(std::set<std::shared_ptr<TreeNode>> &syndrome) override;
private:
    void                     standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                            std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    void                     singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                            std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    void                     singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                              std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    bool                     isValidComponent(const std::shared_ptr<TreeNode>& component);
    std::vector<std::size_t> erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
    void                     extractValidComponents(std::set<std::shared_ptr<TreeNode>>& components, std::vector<std::shared_ptr<TreeNode>>& erasure);
    std::vector<size_t>      peelingDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
};
#endif //QUNIONFIND_IMPROVEDUF_HPP

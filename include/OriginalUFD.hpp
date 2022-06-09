//
// Created by luca on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUF_HPP
#define QUNIONFIND_IMPROVEDUF_HPP
#include "Decoder.hpp"
class OriginalUFD: virtual public Decoder {
public:
    explicit OriginalUFD(Code& code):
        Decoder(code) {};
    void decode(std::vector<bool> &syndrome) override;
private:
    void                     standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                            std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components);
    bool                     isValidComponent(const std::shared_ptr<TreeNode>& component);
    bool containsInvalidComponents(const std::set<std::shared_ptr<TreeNode>>& components);
};
#endif //QUNIONFIND_IMPROVEDUF_HPP

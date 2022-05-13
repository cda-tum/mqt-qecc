//
// Created by lucas on 21/04/2022.
//

#include "Code.hpp"
#include "TreeNode.hpp"

#include <vector>

#ifndef QUNIONFIND_DECODER_HPP
    #define QUNIONFIND_DECODER_HPP

#endif //QUNIONFIND_DECODER_HPP

class Decoder {
public:
    std::vector<TreeNode> decode(std::set<TreeNode>& syndrome);

private:
    Code code;

    bool                  isValidComponent(const TreeNode& component);
    std::vector<TreeNode> erasureDecoder(std::vector<TreeNode>& erasuree);

    bool                     existsInvalidComponent(std::vector<TreeNode>& vector1);
    std::vector<std::size_t> peelingDecoder(std::vector<TreeNode>& erasure, std::set<TreeNode>& syndrome);
    void                     extractValidComponents(std::set<TreeNode>& components, std::vector<TreeNode>& erasure);
};
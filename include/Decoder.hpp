//
// Created by lucas on 21/04/2022.
//

#include "Code.hpp"
#include "TreeNode.hpp"
#include "Codes.hpp"

#include <utility>
#include <vector>

#ifndef QUNIONFIND_DECODER_HPP
#define QUNIONFIND_DECODER_HPP


class Decoder {
public:
    explicit Decoder(Code& code):
        code(code) {}
    std::vector<std::size_t> decode(std::set<std::shared_ptr<TreeNode>>& syndrome); // todo set or vector more efficient?

private:
    Code code;
    bool                  isValidComponent(const std::shared_ptr<TreeNode>& component);
    std::vector<std::size_t> erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
    void                     extractValidComponents(std::set<std::shared_ptr<TreeNode>>& components, std::vector<std::shared_ptr<TreeNode>>& erasure);
    std::vector<size_t>                    peelingDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
};
#endif //QUNIONFIND_DECODER_HPP

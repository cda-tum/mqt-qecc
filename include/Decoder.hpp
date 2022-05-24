//
// Created by lucas on 21/04/2022.
//

#include "Code.hpp"
#include "Codes.hpp"
#include "TreeNode.hpp"

#include <chrono>
#include <utility>
#include <vector>
#ifndef QUNIONFIND_DECODER_HPP
#define QUNIONFIND_DECODER_HPP

struct DecodingResult{
    std::size_t decodingTime;
    std::vector<std::size_t> estimIdxVector;
    std::vector<bool>        estimBoolVector;
};

class Decoder {
public:
    DecodingResult result;
    explicit Decoder(Code& code):
        code(code) {}
    void decode(std::set<std::shared_ptr<TreeNode>>& syndrome); // todo set or vector more efficient?

private:
    void standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                        std::map<std::size_t, bool>&                      presentMap, const std::set<std::shared_ptr<TreeNode>>& components );
    void singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                          std::map<std::size_t, bool>&                      presentMap,const std::set<std::shared_ptr<TreeNode>>& components);
    Code code;
    bool                  isValidComponent(const std::shared_ptr<TreeNode>& component);
    std::vector<std::size_t> erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
    void                     extractValidComponents(std::set<std::shared_ptr<TreeNode>>& components, std::vector<std::shared_ptr<TreeNode>>& erasure);
    std::vector<size_t>                    peelingDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::set<std::shared_ptr<TreeNode>>& syndrome);
};
#endif //QUNIONFIND_DECODER_HPP

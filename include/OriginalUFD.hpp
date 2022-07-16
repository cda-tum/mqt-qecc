//
// Created by lucas on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUF_HPP
#define QUNIONFIND_IMPROVEDUF_HPP
#include "Decoder.hpp"
class OriginalUFD: public Decoder {
public:
    using Decoder::Decoder;
    void decode(std::vector<bool>& syndrome) override;

private:
    bool                     isValidComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const;
    bool                     containsInvalidComponents(const std::vector<std::set<std::size_t>>& components, const std::vector<bool>& syndrome) const;
    std::vector<std::size_t> computeInteriorBitNodes(const std::set<std::size_t>& component) const;
    std::set<std::size_t>    getEstimateForComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const;
};
#endif //QUNIONFIND_IMPROVEDUF_HPP

//
// Created by lucas on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUF_HPP
#define QUNIONFIND_IMPROVEDUF_HPP
#include "Decoder.hpp"
class OriginalUFD: virtual public Decoder {
public:
    explicit OriginalUFD(Code& code):
        Decoder(code){};
    void decode(std::vector<bool>& syndrome) override;

private:
    bool                     isValidComponent(std::set<std::size_t>& component, const std::vector<bool>& syndrome);
    bool                     containsInvalidComponents(std::vector<std::set<std::size_t>>& components, const std::vector<bool>& syndrome);
    std::vector<std::size_t> computeInteriorBitNodes(std::set<std::size_t>& component);
    std::set<std::size_t>    getEstimateForComponent(std::set<std::size_t>& component, const std::vector<bool>& syndrome);
};
#endif //QUNIONFIND_IMPROVEDUF_HPP

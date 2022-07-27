//
// Created by lucas on 09/06/22.
//

#ifndef QUNIONFIND_IMPROVEDUF_HPP
#define QUNIONFIND_IMPROVEDUF_HPP
#include "Decoder.hpp"
class OriginalUFD: public Decoder {
public:
    using Decoder::Decoder;
    void decode(const std::vector<bool>& syndrome) override;
    void reset() override;
private:
    [[nodiscard]] bool                     isValidComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const;
    bool                     containsInvalidComponents(const std::vector<std::set<std::size_t>>& components, const std::vector<bool>& syndrome,
                                                       std::vector<std::set<std::size_t>>& invalidComps) const;
    [[nodiscard]] std::vector<std::size_t> computeInteriorBitNodes(const std::set<std::size_t>& component) const;
    [[nodiscard]] std::set<std::size_t>    getEstimateForComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const;
    void standardGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd);
    void singleClusterSmallestFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd);
    void singleClusterRandomFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd);
    void singleQubitRandomFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd);
};
#endif //QUNIONFIND_IMPROVEDUF_HPP

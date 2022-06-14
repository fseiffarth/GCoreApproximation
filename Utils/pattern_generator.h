//
// Created by anonymous on 10.05.2021.
//

#ifndef CLOSURES_PATTERN_GENERATOR_H
#define CLOSURES_PATTERN_GENERATOR_H

#include "../Data/GraphData.h"
#include "../Utils/Enums.h"

class PatternGenerator {
public:
    template <typename T>
    static std::vector<GraphData> generatePatterns(GraphData& graphData, T& closure, PatternType patternType, int number, std::mt19937_64& gen, bool basedOnTraining = false);
};

#include "pattern_generator.txx"

#endif //CLOSURES_PATTERN_GENERATOR_H

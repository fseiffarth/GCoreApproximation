//
// Created by florian on 15.08.23.
//

#ifndef GCOREAPPROXIMATION_COREGROWALGORITHM_H
#define GCOREAPPROXIMATION_COREGROWALGORITHM_H

#include "typedefs.h"
#include "../Data/GraphData.h"

struct CoreGrowAlgorithmInputParameters {
    int num_runs = 100;
    int grow_steps = 15;
    double core_percentage = 0.9;
    int seed;
    bool print;
    bool save;
};

struct CoreGrowAlgorithmOutputParameters {
    std::vector<NodeId> core_nodes;
    double runtime;
};

class CoreGrowAlgorithm { ;

    GraphData& graphData;
public:
    void Run(CoreGrowAlgorithmOutputParameters& outputParameters, const CoreGrowAlgorithmInputParameters& inputParameters);

    explicit CoreGrowAlgorithm(GraphData &graphData) : graphData(graphData) {}
};


#endif //GCOREAPPROXIMATION_COREGROWALGORITHM_H

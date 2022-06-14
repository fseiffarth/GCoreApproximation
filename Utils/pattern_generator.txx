//
// Created by anonymous on 10.05.2021.
//

#include "pattern_generator.h"
#include "subgraphs.h"
#include "GraphFunctions.h"
#include "OuterplanarSubgraphDFS.h"

template <typename T>
std::vector<GraphData> PatternGenerator::generatePatterns(GraphData &graphData, T& closure, PatternType patternType, int number,
                                   std::mt19937_64 &gen, bool basedOnTraining) {
    std::vector<GraphData> patterns;
    GraphData subgraph = GraphData(new TUNGraph());
    std::vector<int> neighbors;
    GraphFunctions::generateNeighborVector(graphData.get_graph(), neighbors);
    if (basedOnTraining){
        for (int i = 0; i < number; ++i) {
            switch (patternType) {
                case PatternType::BFS_TREE:
                    subgraph.set_graph(GraphFunctions::bfsSubtree(graphData, neighbors, closure, gen, basedOnTraining));
                    patterns.emplace_back(GraphData(subgraph.get_graph(), "Pattern" + std::to_string(i), graphData.labels()));
                    patterns.back().graphType = GraphType::TREE;
                    break;
                case PatternType::OUTERPLANAR:
                    subgraph.set_graph(GraphFunctions::maximalOuterplanarSubgraphGreedy(graphData, closure, gen,
                                                                                basedOnTraining));
                    patterns.emplace_back(GraphData(subgraph.get_graph(), "Pattern" + std::to_string(i), graphData.labels()));
                    patterns.back().graphType = GraphType::OUTERPLANAR;
                    break;
                default:
                    break;
            }
        }
    }
    else {
        OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graphData.get_graph());
        for (int i = 0; i < number; ++i) {
            switch (patternType) {
                case PatternType::BFS_TREE:
                    GraphFunctions::bfsSubtree(graphData.get_graph(), subgraph, neighbors, gen);
                    patterns.emplace_back(GraphData(subgraph.get_graph(), "Pattern" + std::to_string(i), graphData.labels()));
                    patterns.back().graphType = GraphType::TREE;
                    break;
                case PatternType::OUTERPLANAR:
                    subgraphGeneration.generate(subgraph, gen, false);
                    //_subgraph = GraphFunctions::maximalOuterplanarSubgraphGreedy(graphData.graph(), _gen);
                    patterns.emplace_back(GraphData(subgraph.get_graph(), "Pattern" + std::to_string(i), graphData.labels()));
                    patterns.back().graphType = GraphType::OUTERPLANAR;
                    break;
                default:
                    break;
            }
        }
    }
    return patterns;
}

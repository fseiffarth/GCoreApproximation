//
// Created by anonymous on 30.05.2021.
//

#ifndef CLOSURES_GRAPHFUNCTIONS_TXX
#define CLOSURES_GRAPHFUNCTIONS_TXX

#include "../ClosureOperators/BaseOperator.h"
#include "GraphFunctions.h"

template<typename T>
PUNGraph GraphFunctions::bfsSubtree(GraphData &graph, std::vector<NodeId>& neighborIds, T &closure, std::mt19937_64 &gen, bool basedOnTraining) {
    PUNGraph subtree = new TUNGraph();
    if (basedOnTraining){
        std::vector<Nodes> closedSetElements;
        for (auto set : graph.trainingSet) {
            ClosureParameters output = ClosureParameters(set);
            closure.closure(graph, output);
            closedSetElements.emplace_back(output.closed_set.begin(), output.closed_set.end());
        }
        PUNGraph reducedGraph = GraphFunctions::GetReducedGraph(graph.get_graph(), closedSetElements);
        bfsSubtree(reducedGraph,subtree, neighborIds, gen);
        return GraphFunctions::RebuildReducedComponents(graph.get_graph(), subtree, closedSetElements);
    }
    else{
        bfsSubtree(graph.get_graph(), subtree, neighborIds, gen);
        return subtree;
    }
}

template<typename T>
PUNGraph
GraphFunctions::maximalOuterplanarSubgraphGreedy(GraphData &graph, T &closure, std::mt19937_64 &gen, bool basedOnTraining) {
    if (basedOnTraining){
        std::vector<Nodes> closedSetElements;
        for (auto set : graph.trainingSet) {
            ClosureParameters output = ClosureParameters(set);
            closure.closure(graph, output);
            closedSetElements.emplace_back(output.closed_set.begin(), output.closed_set.end());
        }
        PUNGraph reducedGraph = GraphFunctions::GetReducedGraph(graph.get_graph(), closedSetElements);
        PUNGraph reducedOuterplanar = maximalOuterplanarSubgraphGreedy(reducedGraph, gen);
        return GraphFunctions::RebuildReducedComponents(graph.get_graph(), reducedOuterplanar, closedSetElements);
    }
    else{
        return maximalOuterplanarSubgraphGreedy(graph, gen);
    }
}

template<typename T>
PUNGraph
GraphFunctions::outerplanarSubgraphLinear(GraphData &graph, T &closure, std::mt19937_64 &gen, bool basedOnTraining) {
    if (basedOnTraining){
        std::vector<Nodes> closedSetElements;
        for (auto set : graph.trainingSet) {
            ClosureParameters output = ClosureParameters(set);
            closure.closure(graph, output);
            closedSetElements.emplace_back(output.closed_set.begin(), output.closed_set.end());
        }
        PUNGraph reducedGraph = GraphFunctions::GetReducedGraph(graph.get_graph(), closedSetElements);
        PUNGraph reducedOuterplanar = outerplanarSubgraphLinear(reducedGraph, gen);
        return GraphFunctions::RebuildReducedComponents(graph.get_graph(), reducedOuterplanar, closedSetElements);
    }
    else{
        return outerplanarSubgraphLinear(graph, gen);
    }
}

template<typename T>
PUNGraph GraphFunctions::maximalOuterplanarSubgraphLinear(GraphData &graph, T &closure, std::mt19937_64 &gen,
                                                          bool basedOnTraining) {
    if (basedOnTraining){
        std::vector<Nodes> closedSetElements;
        for (auto set : graph.trainingSet) {
            ClosureParameters output = ClosureParameters(set);
            closure.closure(graph, output);
            closedSetElements.emplace_back(output.closed_set.begin(), output.closed_set.end());
        }
        PUNGraph reducedGraph = GraphFunctions::GetReducedGraph(graph.get_graph(), closedSetElements);
        PUNGraph reducedOuterplanar = maximalOuterplanarSubgraphLinear(reducedGraph, gen);
        return GraphFunctions::RebuildReducedComponents(graph.get_graph(), reducedOuterplanar, closedSetElements);
    }
    else{
        return maximalOuterplanarSubgraphLinear(graph, gen);
    }
}



#endif //CLOSURES_GRAPHFUNCTIONS_TXX

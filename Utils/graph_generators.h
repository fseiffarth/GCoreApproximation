//
// Created by anonymous on 07.05.2021.
//

#ifndef DS2021_GRAPH_GENERATORS_H
#define DS2021_GRAPH_GENERATORS_H

#include "typedefs.h"
#include "../Data/GraphData.h"
#include "../ClosureOperators/GraphClosures.h"

enum class GenerationType{
    TWO_COMPONENTS,
    ONE_COMPONENT,
};

class GraphGeneration {
public:
    static GraphData generateGraph(int size, double density, std::mt19937_64& gen);
    static GraphData generateERGraph();
    static GraphData generateSeparableComponents(GraphClosureSP& graphClosure, int sizeA, double densityA, int sizeB, double densityB, int connectivity, bool& valid, bool random_connections, std::mt19937_64& gen);
    static GraphData generateOneComponent(GraphClosureSP &graphClosure, int sizeA, int sizeB, double density, bool &valid, std::mt19937_64 &gen);

private:
    static bool checkEdgePreservesSeparability(GraphClosureSP& graphClosure, GraphData& mergedGraph, const GraphData& componentA, const GraphData& componentB, ClosureParameters& closureOutputA, ClosureParameters& closureOutputB, std::set<NodeId>& nodeSetA, std::set<NodeId>& nodeSetB, Nodes& newEndpointsA, Nodes& newEndpointsB, int source, int destination);
    static void sampleEdge(int &source, int &destination, int sizeA, int sizeB, std::mt19937_64& gen);
};

#endif //DS2021_GRAPH_GENERATORS_H

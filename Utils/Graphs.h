//
// Created by anonymous on 10.05.2021.
//

#ifndef CLOSURES_GRAPHS_H
#define CLOSURES_GRAPHS_H

#include "typedefs.h"
#include "GraphStructs.h"

class Graphs {
public:
    static PUNGraph circle(int size);
    static PUNGraph triangle();
    static PUNGraph complete(int a);
    static PUNGraph completeBipartite(int a, int b);
    static PUNGraph path(int size);
    static PUNGraph mergeGraphs(const std::vector<PUNGraph>& graphs);
    static PUNGraph componentGraph(const std::vector<PUNGraph>& graphs, const std::vector<std::pair<int, int>>& nodePairs);
    static PUNGraph connectGraphs(const std::pair<PUNGraph, PUNGraph> &graphs, const std::vector<std::pair<int, int>> &nodePairs);
};


#endif //CLOSURES_GRAPHS_H

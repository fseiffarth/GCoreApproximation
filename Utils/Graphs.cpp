//
// Created by anonymous on 10.05.2021.
//

#include "Graphs.h"

PUNGraph Graphs::circle(int size) {
    return TSnap::GenCircle<PUNGraph>(size);
}

PUNGraph Graphs::completeBipartite(int a, int b) {
    PUNGraph graph = new TUNGraph(a, b);
    for (int i = 0; i < a + b; ++i) {
        graph->AddNode();
    }
    for (int i = 0; i < a; ++i) {
        for (int j = 0; j < b; ++j) {
            graph->AddEdge(i, a + j);
        }
    }
    return graph;
}

PUNGraph Graphs::complete(int size) {
    return TSnap::GenFull<PUNGraph>(size);
}

PUNGraph Graphs::triangle() {
    return circle(3);
}

PUNGraph Graphs::path(int size) {
    PUNGraph graph = new TUNGraph(size, size-1);
    for (int i = 0; i < size; ++i) {
        graph->AddNode();
    }
    for (int i = 0; i < size - 1; ++i) {
        graph->AddEdge(i, i+1);
    }
    return graph;
}

PUNGraph Graphs::connectGraphs(const std::pair<PUNGraph, PUNGraph>& graphs, const std::vector<std::pair<int, int>>& nodePairs) {
    PUNGraph graph = mergeGraphs({graphs.first, graphs.second});
    int i = graphs.first->GetNodes();
    for (auto const& nodePair : nodePairs) {
        graph->AddEdge(nodePair.first, nodePair.second + i);
    }
    return graph;
}

PUNGraph Graphs::componentGraph(const std::vector<PUNGraph>& graphs, const std::vector<std::pair<int, int>>& nodePairs) {
    PUNGraph graph = mergeGraphs(graphs);
    if (nodePairs.size() - 1 == graphs.size()) {
        int Id = graphs[0]->GetNodes();
        int Counter = 1;
        for (auto const & nodePair : nodePairs) {
            graph->AddEdge(nodePair.first, nodePair.second + Id);
            Id += graphs[Counter]->GetNodes();
            ++Counter;
        }
    }
    return graph;
}

PUNGraph Graphs::mergeGraphs(const std::vector<PUNGraph>& graphs) {
    PUNGraph graph = new TUNGraph();
    int Id = 0;
    for (PUNGraph g : graphs) {
        for (TUNGraph::TNodeI node =  g->BegNI(); node != g->EndNI(); node++) {
            graph->AddNode();
        }
        for (TUNGraph::TEdgeI edge =  g->BegEI(); edge != g->EndEI(); edge++) {
            graph->AddEdge(edge.GetSrcNId() + Id, edge.GetDstNId() + Id);
        }
        Id += g->GetNodes();
    }
    return graph;
}

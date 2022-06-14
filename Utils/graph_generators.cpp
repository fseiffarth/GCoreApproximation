//
// Created by anonymous on 07.05.2021.
//

#include "graph_generators.h"
#include "drawing.h"
#include "../ClosureOperators/GraphClosures.h"
#include "GraphFunctions.h"

void GraphGeneration::sampleEdge(int &source, int &destination, int sizeA, int sizeB, std::mt19937_64& gen) {
    source = std::uniform_int_distribution<int>(0, sizeA - 1)(gen);
    destination = std::uniform_int_distribution<int>(0, sizeB - 1)(gen);
}

GraphData GraphGeneration::generateGraph(int size, double density, std::mt19937_64& gen) {
    int edgeNum = static_cast<int>((size - 1) * density);
    PUNGraph punGraph = GraphFunctions::randomTree(size, gen);
    edgeNum -= (size - 1);
    while(edgeNum > 0){
        int source = 0, destination = 0;
        while (source == destination || punGraph->IsEdge(source, destination)) {
            sampleEdge(source, destination, punGraph->GetNodes(),punGraph->GetNodes(), gen);
        }
        punGraph->AddEdge(source, destination);
        --edgeNum;
    }
    GraphData graph = GraphData(punGraph, "");
    return graph;
}

bool GraphGeneration::checkEdgePreservesSeparability(GraphClosureSP &graphClosure, GraphData &mergedGraph,
                                                     const GraphData &componentA, const GraphData &componentB, ClosureParameters& closureOutputA, ClosureParameters& closureOutputB, std::set<NodeId>& nodeSetA, std::set<NodeId>& nodeSetB,
                                                     Nodes &newEndpointsA, Nodes &newEndpointsB, int newNodeA,
                                                     int newNodeB) {
    mergedGraph.get_graph()->AddEdge(newNodeA, newNodeB);
    for (NodeId nodeId : newEndpointsA) {
        closureOutputA.reset();
        closureOutputA.closed_set.clear();
        nodeSetA = std::set<NodeId>{nodeId, newNodeA};
        graphClosure.closure(mergedGraph, closureOutputA);
        if (closureOutputA.output_intersects_forbidden){
            mergedGraph.get_graph()->DelEdge(newNodeA, newNodeB);
            return false;
        }
    }
    for (NodeId nodeId : newEndpointsB) {
        closureOutputB.reset();
        closureOutputB.closed_set.clear();
        nodeSetB = std::set<NodeId>{nodeId, newNodeB};
        graphClosure.closure(mergedGraph, closureOutputB);
        if (closureOutputB.output_intersects_forbidden){
            mergedGraph.get_graph()->DelEdge(newNodeA, newNodeB);
            return false;
        }
    }
    mergedGraph.get_graph()->DelEdge(newNodeA, newNodeB);
    return true;
}

GraphData GraphGeneration::generateSeparableComponents(GraphClosureSP& graphClosure, int sizeA, double densityA, int sizeB, double densityB, int connectivity, bool& valid, bool random_connections, std::mt19937_64& gen) {
    GraphData componentA = generateGraph(sizeA, densityA, gen);
    GraphData componentB = generateGraph(sizeB, densityB, gen);
    GraphData graph = componentA + componentB;
    Nodes newNodesA, newNodesB;
    std::set<NodeId> nodeSetA, nodeSetB;
    ClosureParameters closureOutputA(nodeSetA);
    ClosureParameters closureOutputB(nodeSetB);
    std::set<NodeId> forbiddenA;
    std::set<NodeId> forbiddenB;
    for (int i = 0; i < componentB.size(); ++i) {
        forbiddenA.insert((int) componentA.size() + i);
    }
    for (int i = 0; i < (int) componentA.size(); ++i) {
        forbiddenB.insert(i);
    }
    closureOutputA.forbidden_elements = {&forbiddenA};
    closureOutputB.forbidden_elements = {&forbiddenB};
    std::vector<std::pair<NodeId, NodeId>> edges;
    for (NodeId i = 0; i < sizeA; ++i) {
        for (NodeId j = sizeA; j < sizeA + sizeB; ++j) {
            edges.emplace_back(std::pair<NodeId, NodeId>{i, j});
        }
    }
    std::shuffle(edges.begin(), edges.end(), gen);
    valid = true;
    for (auto const & [src, dst] : edges) {
        if (connectivity == 0){
            valid = true;
            break;
        }
        if (!graph.get_graph()->IsEdge(src, dst) && (random_connections || checkEdgePreservesSeparability(graphClosure, graph, componentA, componentB, closureOutputA, closureOutputB, nodeSetA, nodeSetB, newNodesA, newNodesB, src, dst))){
            graph.get_graph()->AddEdge(src, dst);
            newNodesA.emplace_back(src);
            newNodesB.emplace_back(dst);
            --connectivity;
        }
    }
    if (connectivity != 0){
        valid = false;
    }
    Labels labels;
    for (int i = 0; i < componentA.get_graph()->GetNodes(); ++i) {
        labels.emplace_back(0);
    }
    for (int i = 0; i < componentB.get_graph()->GetNodes(); ++i) {
        labels.emplace_back(1);
    }
    GraphData labeledGraph = GraphData(graph.get_graph(), "", labels);
    return labeledGraph;
}

GraphData GraphGeneration::generateOneComponent(GraphClosureSP &graphClosure, int sizeA, int sizeB, double density, bool &valid, std::mt19937_64 &gen) {
    valid = true;
    int edgeNum = static_cast<int>((sizeA + sizeB - 1) * density);
    PUNGraph graph = GraphFunctions::randomTree(sizeA, gen);
    for (int i = 0; i < sizeB; ++i) {
        graph->AddNode();
    }
    GraphData oneComponent = GraphData(graph);
    std::set<NodeId> inputSet;
    std::set<NodeId> forbiddenSet;
    for (int i = 0; i < sizeA; ++i) {
        inputSet.insert(i);
    }
    ClosureParameters closureParameters(inputSet);
    for (int i = sizeA; i < sizeA + sizeB; ++i) {
        forbiddenSet.insert(i);
    }
    closureParameters.forbidden_elements = {&forbiddenSet};
    std::vector<int> endPoints = std::vector<int>(sizeA + sizeB);
    std::iota(endPoints.begin(), endPoints.end(), 0);
    //Add nodes of other class plus an edge
    for (int i = 0; i < sizeB; ++i) {
        int newNodeId = sizeA + i;
        while(true) {
            int endPointId = std::uniform_int_distribution<int>(0, newNodeId - 1)(gen);
            oneComponent.get_graph()->AddEdge(endPointId, newNodeId);
            closureParameters.reset();
            closureParameters.closed_set.clear();
            graphClosure.closure(oneComponent, closureParameters);
            if (!closureParameters.output_intersects_forbidden){
                break;
            }
            else{
                oneComponent.get_graph()->DelEdge(endPointId, newNodeId);
            }
        }
    }
    edgeNum -= (sizeA + sizeB - 1);
    std::vector<std::pair<NodeId, NodeId>> edges;
    for (int i = 0; i < sizeA + sizeB; ++i) {
        for (int j = i + 1; j < sizeA + sizeB; ++j) {
            edges.emplace_back(std::pair<NodeId, NodeId>{i, j});
        }
    }
    std::shuffle(edges.begin(), edges.end(), gen);
    //Add _valid random possibleEdges to get density
    for (auto& edge : edges) {
        if (edgeNum == 0){
            break;
        }
        if (!oneComponent.get_graph()->IsEdge(edge.first, edge.second)){
            oneComponent.get_graph()->AddEdge(edge.first, edge.second);
            closureParameters.reset();
            closureParameters.closed_set.clear();
            graphClosure.closure(oneComponent, closureParameters);
            if (closureParameters.output_intersects_forbidden){
                oneComponent.get_graph()->DelEdge(edge.first, edge.second);
                continue;
            }
            else{
                --edgeNum;
            }
        }
    }
    if (edges.empty() && edgeNum > 0){
        valid = false;
        return oneComponent;
    }
    for (int i = 0; i < sizeA; ++i) {
        oneComponent.labels().emplace_back(0);
    }
    for (int i = 0; i < sizeB; ++i) {
        oneComponent.labels().emplace_back(1);
    }
    return oneComponent;
}

GraphData GraphGeneration::generateERGraph() {
    return GraphData();
}



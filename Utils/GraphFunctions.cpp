//
// Created by anonymous on 10.05.2021.
//

#include <list>
#include <deque>
#include "GraphFunctions.h"
#include "GraphStructs.h"
#include "directed_graph.h"
#include "StaticFunctions.h"
#include "../Data/OuterplanarGraphData.h"
#include "Enums.h"

void GraphFunctions::mergeGraphs(PUNGraph& graph, const PUNGraph &other) {
    int initialSize = graph->GetNodes();
    for(TUNGraph::TNodeI node = other->BegNI(); node != other->EndNI(); node++){
        graph->AddNode();
    }
    for(TUNGraph::TEdgeI edge = other->BegEI(); edge != other->EndEI(); edge++){
        graph->AddEdge(initialSize + edge.GetSrcNId(), initialSize + edge.GetDstNId());
    }
}

void GraphFunctions::addEdgesToGraph(PUNGraph& graph, const PUNGraph &other) {
    for(TUNGraph::TEdgeI edge = other->BegEI(); edge != other->EndEI(); edge++){
        int src = edge.GetSrcNId();
        int dst = edge.GetDstNId();
        if (graph->IsNode(src) && graph->IsNode(dst)) {
            graph->AddEdge(src, dst);
        }
    }
}

void GraphFunctions::bfsTree(const PUNGraph& undirectedGraph, const PUNGraph& subTree, NodeId rootNodeId, std::vector<NodePair>& missingEdges, std::mt19937_64& gen) {
    subTree->Clr();
    int maxDegree = 0;
    for (auto node = undirectedGraph->BegNI(); node != undirectedGraph->EndNI(); node++) {
        subTree->AddNode();
        maxDegree = std::max(maxDegree, node.GetDeg());
    }
    std::vector<int> neighborIds = std::vector<int>(maxDegree);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
    std::vector<std::pair<int, int>> swappedIds;
    std::vector<bool> VisitedNodes = std::vector<bool>(undirectedGraph->GetNodes(), false);
    VisitedNodes[rootNodeId] = true;
    std::deque<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = undirectedGraph->GetNI(NextNodeId).GetDeg();
        int max = degree;
        for (NodeId i = 0; i < degree; ++i) {
            --max;
            //Get random neighbor
            int idx = std::uniform_int_distribution<int>(0, max)(gen);
            NodeId NeighborNodeId = undirectedGraph->GetNI(NextNodeId).GetNbrNId(neighborIds[idx]);
            std::swap(neighborIds[idx], neighborIds[max]);
            swappedIds.emplace_back(std::pair<int, int>{idx, max});
            if (!VisitedNodes[NeighborNodeId]){
                VisitedNodes[NeighborNodeId] = true;
                CurrentNodes.push_front(NeighborNodeId);
                subTree->AddEdge(NextNodeId, NeighborNodeId);
            }
            else{
                missingEdges.emplace_back(NodePair(NextNodeId, NeighborNodeId));
            }
        }
        int forLoopSize = (int) swappedIds.size();
        for (int i = 0; i < forLoopSize; ++i) {
            std::swap(neighborIds[swappedIds.back().first], neighborIds[swappedIds.back().second]);
            swappedIds.pop_back();
        }
    }
}

void GraphFunctions::getNodesInBall(const PUNGraph& undirectedGraph, TIntV& reachedNodes, std::mt19937_64& gen, int size,NodeId rootNodeId, int hop) {
    if (rootNodeId == -1){
        rootNodeId = std::uniform_int_distribution<int>(0, undirectedGraph->GetNodes() - 1)(gen);
    }
    std::vector<std::vector<NodeId>> hopNodes = std::vector<std::vector<NodeId>>{std::vector<NodeId>{rootNodeId}};
    std::vector<int> hop_sizes = std::vector<int>{1};
    std::vector<int> hop_distance = std::vector<int>(undirectedGraph->GetNodes(), -1);
    hop_distance[rootNodeId] = 0;
    std::deque<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    int visited = 1;
    while (!CurrentNodes.empty()){
        if (size != -1 && size <= visited){
            break;
        }
        NodeId CurrentNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = undirectedGraph->GetNI(CurrentNodeId).GetDeg();
        for (NodeId i = 0; i < degree; ++i) {
            //Get random neighbor
            NodeId NeighborNodeId = undirectedGraph->GetNI(CurrentNodeId).GetNbrNId(i);
            if (hop_distance[NeighborNodeId] == -1){
                hop_distance[NeighborNodeId] = hop_distance[CurrentNodeId] + 1;
                if (hop_distance[NeighborNodeId] >= hopNodes.size()){
                    hopNodes.emplace_back(std::vector<NodeId>{NeighborNodeId});
                    hop_sizes.emplace_back(1);
                }
                else{
                    hopNodes[hop_distance[NeighborNodeId]].emplace_back(NeighborNodeId);
                    ++hop_sizes[hop_distance[NeighborNodeId]];
                }
                CurrentNodes.push_front(NeighborNodeId);
                ++visited;
            }
        }
    }
    int reachedNodesSize = 0;
    int current_hop = 0;
    reachedNodes.Clr();
    if (size != -1 && hop != -1) {
        while (current_hop < hop_sizes.size() && (reachedNodes.Len() + hop_sizes[current_hop] < size || current_hop != hop + 1)) {
            for (auto nodeId : hopNodes[current_hop]) {
                reachedNodes.Add(nodeId);
            }
            ++current_hop;
        }
    }
    else if(hop != -1){
        while (current_hop < hop_sizes.size() &&current_hop != hop + 1) {
            for (auto nodeId : hopNodes[current_hop]) {
                reachedNodes.Add(nodeId);
            }
            ++current_hop;
        }
    }
    else if(size != -1){
        while (current_hop < hop_sizes.size() && reachedNodes.Len() + hop_sizes[current_hop] < size) {
            for (auto nodeId : hopNodes[current_hop]) {
                reachedNodes.Add(nodeId);
            }
            ++current_hop;
        }
    }
    else{
        while (current_hop < hop_sizes.size()){
            for (auto nodeId : hopNodes[current_hop]) {
                reachedNodes.Add(nodeId);
            }
            ++current_hop;
        }
    }
}


void GraphFunctions::bfsSubtree(const PUNGraph &graph, std::vector<NodePair>& edges, NodeId rootNodeId, std::mt19937_64 &gen) {
    int maxDegree = 0;
    for (auto node = graph->BegNI(); node != graph->EndNI(); node++) {
        int degree = node.GetDeg();
        maxDegree = std::max(maxDegree, degree);
    }
    std::vector<int> neighborIds = std::vector<int>(maxDegree);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
    std::vector<std::pair<int, int>> swappedIds;
    std::vector<bool> VisitedNodes = std::vector<bool>(graph->GetNodes(), false);
    VisitedNodes[rootNodeId] = true;
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = graph->GetNI(NextNodeId).GetDeg();
        int max = degree;
        for (NodeId i = 0; i < degree; ++i) {
            --max;
            //Get random neighbor
            int idx = std::uniform_int_distribution<int>(0, max)(gen);
            NodeId NeighborNodeId = graph->GetNI(NextNodeId).GetNbrNId(neighborIds[idx]);
            std::swap(neighborIds[idx], neighborIds[max]);
            swappedIds.emplace_back(std::pair<int, int>{idx, max});
            if (!VisitedNodes[NeighborNodeId]){
                VisitedNodes[NeighborNodeId] = true;
                CurrentNodes.insert(CurrentNodes.begin(), NeighborNodeId);
                edges.emplace_back(NodePair(NextNodeId, NeighborNodeId));
            }
        }
        int forLoopSize = (int) swappedIds.size();
        for (int i = 0; i < forLoopSize; ++i) {
            std::swap(neighborIds[swappedIds.back().first], neighborIds[swappedIds.back().second]);
            swappedIds.pop_back();
        }
    }
}

void GraphFunctions::dfsTree(const PUNGraph &undirectedGraph, std::vector<NodePair>& edges, NodeId rootNodeId, std::mt19937_64& gen) {

}

void GraphFunctions::dfsRandomSubtree(const PUNGraph &undirectedGraph, std::vector<NodePair> &edges,
                                      std::mt19937_64 &gen) {
    NodeId randNode = std::uniform_int_distribution<int>(0, undirectedGraph->GetNodes() - 1)(gen);
    dfsTree(undirectedGraph, edges, randNode, gen);
}

bool GraphFunctions::IsOuterPlanarIncremental(const PUNGraph &graph, const NodePair& newEdge, std::vector<int>& nodeDegrees, std::vector<int>& degree2Nodes) {
    if (graph->GetEdges() > 2 * graph->GetNodes() - 3) {
        return false;
    } else {
        std::unordered_map<NodePair, int, hashNodePair> triangulationCount;
        std::vector<NodePair> pairs;
        std::vector<NodePair> edges;

        for (TUNGraph::TEdgeI edge = graph->BegEI(); edge != graph->EndEI(); edge++) {
            int src = edge.GetSrcNId();
            int dest = edge.GetDstNId();
            triangulationCount[NodePair(src, dest)] = 0;
        }
        while (!degree2Nodes.empty()) {
            int currentNodeId = degree2Nodes.back();
            degree2Nodes.pop_back();
        }
        return true;
    }
}

bool GraphFunctions::IsMaximalOuterplanarSubgraph(const PUNGraph &graph, const PUNGraph & subgraph, std::vector<NodePair>& missingEdges) {
    bool outerplanar = IsOuterPlanar(subgraph);
    if (!outerplanar){
        return false;
    }
    bool additionalEdge = false;
    missingEdges.clear();

    PUNGraph check_graph = new TUNGraph();
    for (int i = 0; i < subgraph->GetNodes(); ++i) {
        check_graph->AddNode();
    }
    for (auto EdgeIt = subgraph->BegEI(); EdgeIt != subgraph->EndEI(); EdgeIt++) {
        check_graph->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
    }

    if (!TSnap::IsConnected(check_graph)){
        throw std::length_error("Graph and subgraph do not have the same number of nodes");
    }

    for (auto edge = graph->BegEI(); edge != graph->EndEI(); edge++) {
        NodeId src = edge.GetSrcNId();
        NodeId dst = edge.GetDstNId();
        if (!check_graph->IsEdge(src,dst)){
            check_graph->AddEdge(src, dst);
            additionalEdge = IsOuterPlanar(check_graph, src, dst);
            if (additionalEdge){
                missingEdges.emplace_back(NodePair(edge.GetSrcNId(), edge.GetDstNId()));
            }
            else{
                check_graph->DelEdge(edge.GetSrcNId(), edge.GetDstNId());
            }
        }
    }
    if (missingEdges.empty()){
        return true;
    }
    return false;
}

bool GraphFunctions::IsOuterPlanar(const PUNGraph &graph, NodeId src, NodeId dst) {
    TCnComV components;
    TSnap::GetBiCon(graph, components);
    for (int i = 0; i < components.Len(); ++i) {
        if (components[i].Len() > 2) {
            if (src == -1 || dst == -1 || (components[i].IsNIdIn(src) && components[i].IsNIdIn(dst))) {
                PUNGraph graphComponent = TSnap::GetSubGraph(graph, components[i](), true);
                if (graphComponent->GetEdges() > 2 * graphComponent->GetNodes() - 3) {
                    return false;
                } else {
                    PUNGraph algorithmGraph = new TUNGraph();
                    Nodes degree2Nodes;
                    std::unordered_map<NodePair, int, hashNodePair> triangulationCount;
                    std::vector<NodePair> pairs;
                    std::vector<NodePair> edges;
                    //Preprocessing on get_node degrees
                    for (TUNGraph::TNodeI node = graphComponent->BegNI(); node != graphComponent->EndNI(); node++) {
                        algorithmGraph->AddNode();
                        int id = node.GetId();
                        int degree = node.GetDeg();
                        if (degree == 2) {
                            degree2Nodes.emplace_back(id);
                        }
                    }
                    for (TUNGraph::TEdgeI edge = graphComponent->BegEI(); edge != graphComponent->EndEI(); edge++) {
                        int a = edge.GetSrcNId();
                        int b = edge.GetDstNId();
                        algorithmGraph->AddEdge(a, b);
                        triangulationCount[NodePair(a, b)] = 0;
                    }
                    while (!degree2Nodes.empty()) {
                        int currentNodeId = degree2Nodes.back();
                        degree2Nodes.pop_back();
                        TUNGraph::TNodeI currentNode = algorithmGraph->GetNI(currentNodeId);
                        if (currentNode.GetDeg() == 2) {
                            int near = currentNode.GetNbrNId(0);
                            TUNGraph::TNodeI nearNode = algorithmGraph->GetNI(near);
                            int next = currentNode.GetNbrNId(1);
                            TUNGraph::TNodeI nextNode = algorithmGraph->GetNI(next);
                            NodePair nodePair = NodePair(near, next);
                            algorithmGraph->DelEdge(currentNodeId, near);
                            algorithmGraph->DelEdge(currentNodeId, next);
                            if (!algorithmGraph->IsEdge(near, next)) {
                                algorithmGraph->AddEdge(near, next);
                                triangulationCount[nodePair] = std::max(1, std::max(
                                        triangulationCount[NodePair(currentNodeId, near)],
                                        triangulationCount[NodePair(currentNodeId, next)]));
                                if (triangulationCount[nodePair] > 2) {
                                    return false;
                                }
                            } else {
                                triangulationCount[nodePair] += 1;
                                if (triangulationCount[nodePair] > 2 ||
                                    triangulationCount[NodePair(currentNodeId, near)] == 2 ||
                                    triangulationCount[NodePair(currentNodeId, next)] == 2) {
                                    return false;
                                }
                                if (nearNode.GetDeg() == 2) {
                                    degree2Nodes.emplace_back(near);
                                }
                                if (nextNode.GetDeg() == 2) {
                                    degree2Nodes.emplace_back(next);
                                }
                            }
                        }
                    }
                    if (algorithmGraph->GetEdges() > 1) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}



PUNGraph GraphFunctions::GetReducedGraph(const PUNGraph &graph, std::vector<Nodes> &reducedComponents, bool renumberNodes) {
    PUNGraph reducedGraph = new TUNGraph();
    std::map<int, NodeId> componentToNodeId;
    for (auto nodeIterator=graph->BegNI(); nodeIterator != graph->EndNI(); nodeIterator++) {
        int idx;
        if (!FoundInComponent(nodeIterator.GetId(), reducedComponents, idx)){
            reducedGraph->AddNode(nodeIterator.GetId());
        }
    }
    for (size_t i = 0; i< reducedComponents.size(); ++i) {
        int newNode = reducedGraph->AddNode();
        componentToNodeId[i] = newNode;
    }
    for (auto edgeIt=graph->BegEI(); edgeIt != graph->EndEI(); edgeIt++) {
        int src = edgeIt.GetSrcNId();
        int dst = edgeIt.GetDstNId();
        int componentIdSrc;
        int componentIdDst;
        bool srcInComponent = FoundInComponent(src, reducedComponents, componentIdSrc);
        bool dstInComponent = FoundInComponent(dst, reducedComponents, componentIdDst);

        if (!srcInComponent && !dstInComponent){
            reducedGraph->AddEdge(src, dst);
        }
        else if (!srcInComponent){
            reducedGraph->AddEdge(src, componentToNodeId[componentIdDst]);
        }
        else if (!dstInComponent){
            reducedGraph->AddEdge(componentToNodeId[componentIdSrc], dst);
        }
    }
    return reducedGraph;
}

PUNGraph GraphFunctions::ResetGraphIds(const PUNGraph &graph) {
    PUNGraph normalizedGraph = new TUNGraph();
    std::map<NodeId, NodeId> mapPair;
    for (auto nodeIterator=graph->BegNI(); nodeIterator != graph->EndNI(); nodeIterator++) {
        int oldId = nodeIterator.GetId();
        int newId = normalizedGraph->AddNode();
        mapPair[oldId] = newId;
    }
    for (auto edgeIt=graph->BegEI(); edgeIt != graph->EndEI(); edgeIt++) {
        normalizedGraph->AddEdge(mapPair[edgeIt.GetSrcNId()], mapPair[edgeIt.GetDstNId()]);
    }
    return normalizedGraph;
}

void GraphFunctions::NormalizeGraphIds(PUNGraph &graph) {
    graph = ResetGraphIds(graph);
}

bool GraphFunctions::FoundInComponent(NodeId nodeId, std::vector<Nodes> &components, int& componentIdx) {
    int Counter = 0;
    for (auto& component : components) {
        if (std::find(component.begin(), component.end(), nodeId) != component.end()){
            componentIdx = Counter;
            return true;
        }
        ++Counter;
    }
    return false;
}

PUNGraph GraphFunctions::RebuildReducedComponents(const PUNGraph &fullGraph, const PUNGraph &reducedGraph,
                                                  std::vector<Nodes> &reducedComponents, bool renumberNodes) {
    PUNGraph rebuildGraph = new TUNGraph();
    for (auto nodeIterator=fullGraph->BegNI(); nodeIterator != fullGraph->EndNI(); nodeIterator++) {
        reducedGraph->AddNode();
    }
    for (auto edgeIt=reducedGraph->BegEI(); edgeIt != reducedGraph->EndEI(); edgeIt++) {
        int src = edgeIt.GetSrcNId();
        int dst = edgeIt.GetDstNId();
        if (src < reducedGraph->GetNodes() - reducedComponents.size() && dst < reducedGraph->GetNodes() - reducedComponents.size()){
            rebuildGraph->AddEdge(src, dst);
        }
    }
    for (auto edgeIt=fullGraph->BegEI(); edgeIt != fullGraph->EndEI(); edgeIt++) {
        int src = edgeIt.GetSrcNId();
        int dst = edgeIt.GetDstNId();
        int srcComponentIdx;
        bool srcInComponent = FoundInComponent(src, reducedComponents, srcComponentIdx);
        int dstComponentIdx;
        bool dstInComponent = FoundInComponent(dst, reducedComponents, dstComponentIdx);
        if (srcInComponent && dstComponentIdx){
            rebuildGraph->AddEdge(src, dst);
        }
        else if (srcInComponent){
            if (reducedGraph->IsEdge(static_cast<int>(reducedGraph->GetNodes()) - reducedComponents.size() + srcComponentIdx, dst)){
                rebuildGraph->AddEdge(src, dst);
            }
        }
        else if (dstInComponent){
            if (reducedGraph->IsEdge(src, static_cast<int>(reducedGraph->GetNodes()) - reducedComponents.size() + dstComponentIdx)){
                rebuildGraph->AddEdge(src, dst);
            }
        }
    }
    return rebuildGraph;
}

PUNGraph GraphFunctions::randomTree(int size, int seed) {
    std::mt19937_64 gen(seed);
    PUNGraph randomTree = new TUNGraph();
    //Create Tree
    for (int i = 0; i < size; ++i) {
        int newNodeId = randomTree->AddNode();
        if (randomTree->GetNodes() > 1){
            int randomNode = std::uniform_int_distribution<int>(0, randomTree->GetNodes() - 2)(gen);
            randomTree->AddEdge(randomNode, newNodeId);
        }
    }
    return randomTree;
}

PUNGraph GraphFunctions::randomTree(int size, std::mt19937_64 & gen) {
    PUNGraph randomTree = new TUNGraph();
    //Create Tree
    for (int i = 0; i < size; ++i) {
        int newNodeId = randomTree->AddNode();
        if (randomTree->GetNodes() > 1){
            int randomNode = std::uniform_int_distribution<int>(0, randomTree->GetNodes() - 2)(gen);
            randomTree->AddEdge(randomNode, newNodeId);
        }
    }
    return randomTree;
}

PUNGraph GraphFunctions::randomOuterplanar(int size, int seed) {
    PUNGraph outerplanar = randomTree(size, seed);
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            if (!outerplanar->IsEdge(i, j)){
                outerplanar->AddEdge(i, j);
                if (!GraphFunctions::IsOuterPlanar(outerplanar)){
                    outerplanar->DelEdge(i, j);
                }
            }
        }
    }
    return outerplanar;
}

void GraphFunctions::bfsSubtree(const PUNGraph &undirectedGraph, GraphData &subTree, std::vector<NodeId> &neighborIds,
                                std::mt19937_64 &gen, int rootNodeId) {
    if (rootNodeId == -1) {
        rootNodeId = std::uniform_int_distribution<int>(0, undirectedGraph->GetNodes() - 1)(gen);
    }
    std::deque<std::pair<int, int>> swappedIds;
    for (int i=0; i<undirectedGraph->GetNodes(); ++i) {
        subTree.graph()->AddNode();
    }
    subTree.graphType = GraphType::TREE;
    ++subTree.Id;
    subTree.ContainmentList[rootNodeId] = subTree.Id;
    std::deque<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = undirectedGraph->GetNI(NextNodeId).GetDeg();
        for (NodeId i = 0; i < degree; ++i) {
            //Get random neighbor
            int idx = std::uniform_int_distribution<int>(i, degree - 1)(gen);
            NodeId NeighborNodeId = undirectedGraph->GetNI(NextNodeId).GetNbrNId(neighborIds[idx]);
            std::swap(neighborIds[idx], neighborIds[i]);
            swappedIds.push_front(std::pair<int, int>{idx, i});
            if (subTree.ContainmentList[NeighborNodeId] != subTree.Id){
                CurrentNodes.push_front(NeighborNodeId);
                subTree.graph()->AddEdge(NextNodeId, NeighborNodeId);
                subTree.ContainmentList[NeighborNodeId] = subTree.Id;
            }
        }
        for (auto const & [a, b] : swappedIds) {
            std::swap(neighborIds[b], neighborIds[a]);
        }
        swappedIds.clear();
    }
}

void GraphFunctions::dfsSubtree(const PUNGraph &undirectedGraph, GraphData &subTree, std::vector<NodeId> &neighborIds,
                                std::mt19937_64 &gen, int rootNodeId) {
    if (rootNodeId == -1) {
        rootNodeId = std::uniform_int_distribution<int>(0, undirectedGraph->GetNodes() - 1)(gen);
    }
    std::deque<std::pair<int, int>> swappedIds;
    for (int i=0; i<undirectedGraph->GetNodes(); ++i) {
        subTree.graph()->AddNode();
    }
    subTree.graphType = GraphType::TREE;
    ++subTree.Id;
    subTree.ContainmentList[rootNodeId] = subTree.Id;
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = undirectedGraph->GetNI(NextNodeId).GetDeg();
        for (NodeId i = 0; i < degree; ++i) {
            //Get random neighbor
            int idx = std::uniform_int_distribution<int>(i, degree - 1)(gen);
            NodeId NeighborNodeId = undirectedGraph->GetNI(NextNodeId).GetNbrNId(neighborIds[idx]);
            std::swap(neighborIds[idx], neighborIds[i]);
            swappedIds.push_front(std::pair<int, int>{idx, i});
            if (subTree.ContainmentList[NeighborNodeId] != subTree.Id){
                CurrentNodes.push_back(NeighborNodeId);
                subTree.graph()->AddEdge(NextNodeId, NeighborNodeId);
                subTree.ContainmentList[NeighborNodeId] = subTree.Id;
            }
        }
        for (auto const & [a, b] : swappedIds) {
            std::swap(neighborIds[b], neighborIds[a]);
        }
        swappedIds.clear();
    }
}


void
GraphFunctions::bfsSubtree(const PUNGraph &undirectedGraph, PUNGraph &subtree, std::vector<NodeId> &neighborIds,
                           std::mt19937_64 &gen, int rootNodeId) {
    if (rootNodeId == -1) {
        rootNodeId = std::uniform_int_distribution<int>(0, undirectedGraph->GetNodes() - 1)(gen);
    }
    std::deque<std::pair<int, int>> swappedIds;
    for (auto NodeIt = undirectedGraph->BegNI(); NodeIt != undirectedGraph->EndNI(); NodeIt++) {
        subtree->AddNode();
    }
    std::deque<NodeId> CurrentNodes;
    std::vector<bool> visitedNodes = std::vector<bool>(undirectedGraph->GetNodes(), false);
    visitedNodes[rootNodeId] = true;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        int degree = undirectedGraph->GetNI(NextNodeId).GetDeg();
        for (NodeId i = 0; i < degree; ++i) {
            //Get random neighbor
            int idx = std::uniform_int_distribution<int>(i, degree - 1)(gen);
            NodeId NeighborNodeId = undirectedGraph->GetNI(NextNodeId).GetNbrNId(neighborIds[idx]);
            std::swap(neighborIds[idx], neighborIds[i]);
            swappedIds.push_front(std::pair<int, int>{idx, i});
            if (!visitedNodes[NeighborNodeId]){
                CurrentNodes.push_front(NeighborNodeId);
                subtree->AddEdge(NextNodeId, NeighborNodeId);
                visitedNodes[NeighborNodeId] = true;
            }
        }
        for (auto const & [a, b] : swappedIds) {
            std::swap(neighborIds[b], neighborIds[a]);
        }
        swappedIds.clear();
    }
}


PUNGraph GraphFunctions::maximalOuterplanarSubgraphGreedy(const PUNGraph &graph, std::mt19937_64& gen) {
    PUNGraph outerplanarSubgraph = new TUNGraph();
    for (int i = 0; i < graph->GetNodes(); ++i) {
        outerplanarSubgraph->AddNode();
    }
    std::vector<TUNGraph::TEdgeI> edges;
    for (TUNGraph::TEdgeI edge = graph->BegEI(); edge != graph->EndEI(); edge++) {
        edges.emplace_back(edge);
    }
    for (int i = 0; i < graph->GetEdges(); ++i) {
        int id = std::uniform_int_distribution<int>(0, static_cast<int>(edges.size()) - 1)(gen);
        TUNGraph::TEdgeI edge = edges[id];
        std::swap(edges[id], edges[edges.size() - 1]);
        edges.pop_back();
        int src = edge.GetSrcNId();
        int dest = edge.GetDstNId();
        outerplanarSubgraph->AddEdge(src, dest);
        if (!GraphFunctions::IsOuterPlanar(outerplanarSubgraph)){
            outerplanarSubgraph->DelEdge(edge.GetSrcNId(), edge.GetDstNId());
        }
    }
    return outerplanarSubgraph;
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphBFSGreedy(const PUNGraph &graph, std::mt19937_64& gen) {
    std::vector<NodePair> edges;
    NodeId rootNodeId = std::uniform_int_distribution<int>(0, graph->GetNodes() - 1)(gen);
    PUNGraph outerplanarSubgraph = new TUNGraph();
    GraphFunctions::bfsTree(graph, outerplanarSubgraph, rootNodeId, edges, gen);
    for (int i = 0; i < edges.size();++i) {
        int id = std::uniform_int_distribution<int>(0, static_cast<int>(edges.size()) - 1)(gen);
        NodePair edge = edges[id];
        std::swap(edges[id], edges.back());
        edges.pop_back();
        int src = edge.first();
        int dest = edge.second();
        outerplanarSubgraph->AddEdge(src, dest);
        if (!GraphFunctions::IsOuterPlanar(outerplanarSubgraph)){
            outerplanarSubgraph->DelEdge(edge.first(), edge.second());
        }
    }
    return outerplanarSubgraph;
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphBFSGreedyIncremental(const PUNGraph &graph, std::mt19937_64 &gen) {
    std::vector<NodePair> edges;
    NodeId rootNodeId = std::uniform_int_distribution<int>(0, graph->GetNodes() - 1)(gen);
    PUNGraph outerplanarSubgraph = new TUNGraph();
    GraphFunctions::bfsTree(graph, outerplanarSubgraph, rootNodeId, edges, gen);
    std::vector<int> nodeDegrees;
    std::vector<NodeId> degree2Nodes;
    for(auto node = outerplanarSubgraph->BegNI(); node != outerplanarSubgraph->EndNI(); node++){
        int degree = node.GetDeg();
        nodeDegrees.emplace_back(degree);
        if (degree == 2){
            degree2Nodes.emplace_back(node.GetId());
        }
    }
    for (int i = 0; i < edges.size();++i) {
        int id = std::uniform_int_distribution<int>(0, static_cast<int>(edges.size()) - 1)(gen);
        NodePair edge = edges[id];
        std::swap(edges[id], edges.back());
        edges.pop_back();
        int src = edge.first();
        int dest = edge.second();
        outerplanarSubgraph->AddEdge(src, dest);
        if (!GraphFunctions::IsOuterPlanarIncremental(outerplanarSubgraph, NodePair(src, dest), nodeDegrees, degree2Nodes)){
            outerplanarSubgraph->DelEdge(edge.first(), edge.second());
        }
    }
    return outerplanarSubgraph;
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphDFSGreedy(const PUNGraph &graph, std::mt19937_64& gen) {
    std::vector<NodeId> visitedNodes;
    PUNGraph outerplanarSubgraph = GraphFunctions::dfsRandomTree(graph, visitedNodes, gen);
    std::vector<TUNGraph::TEdgeI> edges;
    for (TUNGraph::TEdgeI edge = graph->BegEI(); edge != graph->EndEI(); edge++) {
        if (outerplanarSubgraph->IsEdge(edge.GetSrcNId(), edge.GetDstNId())) {
            edges.emplace_back(edge);
        }
    }
    for (int i = 0; i < graph->GetEdges(); ++i) {
        int id = std::uniform_int_distribution<int>(0, static_cast<int>(edges.size()) - 1)(gen);
        TUNGraph::TEdgeI edge = edges[id];
        std::swap(edges[id], edges[edges.size() - 1]);
        edges.pop_back();
        int src = edge.GetSrcNId();
        int dest = edge.GetDstNId();
        outerplanarSubgraph->AddEdge(src, dest);
        if (!GraphFunctions::IsOuterPlanar(outerplanarSubgraph)){
            outerplanarSubgraph->DelEdge(edge.GetSrcNId(), edge.GetDstNId());
        }
    }
    return outerplanarSubgraph;
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphGreedy(const GraphData &graph, std::mt19937_64 &gen) {
    return maximalOuterplanarSubgraphGreedy(graph.get_graph(), gen);
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphLinear(const GraphData &graph, std::mt19937_64 &gen) {
    return maximalOuterplanarSubgraphLinear(graph.get_graph(), gen);
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphLinear(const PUNGraph &graph, std::mt19937_64 &gen) {
    PUNGraph outerplanarGraph = new TUNGraph();
    for (int i = 0; i < graph->GetNodes(); ++i) {
        outerplanarGraph->AddNode();
    }
    TCnComV components;
    TSnap::GetBiCon(graph, components);
    for (int i = 0; i < components.Len(); ++i) {
        biconnectedComponentMaximalOuterPlanarSubgraph(graph, components[i], outerplanarGraph, gen);
    }
    return outerplanarGraph;
}


///
/// \param algorithmGraph output algorithmGraph which should be outerplanar
/// \param componentGraph biconnected componentGraph of a algorithmGraph
/// \param gen random generator
void GraphFunctions::biconnectedComponentMaximalOuterPlanarSubgraph(const PUNGraph& graph, TCnCom& component, PUNGraph& outerplanarGraph, std::mt19937_64& gen) {
    PUNGraph componentGraph;
    componentGraph = TSnap::GetSubGraph(graph, component());
    //Component is a tree
    if (componentGraph->GetEdges() == componentGraph->GetNodes() - 1){
        //Add componentGraph possibleEdges to algorithmGraph
        for (auto edge = componentGraph->BegEI(); edge != componentGraph->EndEI(); edge++) {
            NodeId src = edge.GetSrcNId();
            NodeId dst = edge.GetDstNId();
            outerplanarGraph->AddEdge(src, dst);
        }
        return;
    }
        //Not a tree
    else {
        Nodes degree2Nodes;
        std::vector<std::vector<NodePair>> dependentEdges;
        std::unordered_map<NodePair, EdgeInfo, hashNodePair> algorithmEdges;
        std::vector<NodeId> currentNeighbors = std::vector<NodeId>(2);
        std::vector<NodePair> removableEdges;
        //std::unordered_map<NodePair, bool, hashNodePair> resultingEdges;
        std::unordered_map<int, int> nodeDegrees;
        int edgeNum = componentGraph->GetEdges();

        //Add nodes to componentGraph graphs
        for (auto node = componentGraph->BegNI(); node != componentGraph->EndNI(); node++) {
            int nodeId = node.GetId();
            nodeDegrees[nodeId] = node.GetDeg();
            if (nodeDegrees[nodeId] == 2){
                degree2Nodes.emplace_back(nodeId);
            }
        }
        //Add possibleEdges to componentGraph graph
        for (auto edge = componentGraph->BegEI(); edge != componentGraph->EndEI(); edge++) {
            int src = edge.GetSrcNId();
            int dst = edge.GetDstNId();
            algorithmEdges[NodePair(src, dst)] = EdgeInfo(0, true);
        }
        //Run DFS to mark undeletable possibleEdges
        std::vector<NodePair> bfsEdges;
        GraphFunctions::bfsRandomSubtree(componentGraph, bfsEdges, gen);
        for (auto const& nodePair :bfsEdges) {
            algorithmEdges[nodePair].set_bfs_edge();
        }
        //Find all removable possibleEdges
        for (auto const& [nodePair, edgeInfo]: algorithmEdges) {
            if (!edgeInfo.is_bfs_edge()) {
                removableEdges.emplace_back(nodePair);
            }
        }
        //Start algorithm
        int currentNodeId;
        while (edgeNum > 1) {
            if (!degree2Nodes.empty()) {
                NodeId n1, n2;
                //Get random get_node of out_degree 2
                int randIdx = std::uniform_int_distribution<int>(0, ((int) degree2Nodes.size()) - 1)(gen);
                std::swap(degree2Nodes[randIdx], degree2Nodes.back());
                currentNodeId = degree2Nodes.back();
                degree2Nodes.pop_back();
                if (nodeDegrees[currentNodeId] == 2) {
                    getValidNeighbor(componentGraph, algorithmEdges, currentNodeId, 0, n1);
                    getValidNeighbor(componentGraph, algorithmEdges, currentNodeId, 1, n2);
                    NodePair triangulationPair = NodePair(n1, n2);
                    currentNeighbors[0] = n1, currentNeighbors[1] = n2;
                    std::shuffle(currentNeighbors.begin(), currentNeighbors.end(), gen);
                    NodePair firstNeighborPair = NodePair(currentNodeId, currentNeighbors[0]), secondNeighborPair = NodePair(currentNodeId, currentNeighbors[1]);
                    EdgeInfo &firstNeighborEdge = algorithmEdges[firstNeighborPair], &secondNeighborEdge = algorithmEdges[secondNeighborPair];
                    int countFirstNeighbor = firstNeighborEdge.count(), countSecondNeighbor = secondNeighborEdge.count();
                    //Delete possibleEdges implicitly
                    firstNeighborEdge.set_valid(false), secondNeighborEdge.set_valid(false);
                    edgeNum -= 2;
                    nodeDegrees[currentNodeId] = 0;
                    if (algorithmEdges.find(triangulationPair) != algorithmEdges.end() &&
                        algorithmEdges[triangulationPair].is_valid()) {
                        --nodeDegrees[n1], --nodeDegrees[n2];
                        EdgeInfo &triangleEdge = algorithmEdges[triangulationPair];
                        int countTriangle = triangleEdge.count();
                        //Check if possibleEdges need to be deleted from original graph
                        if (countFirstNeighbor == 2 || countTriangle == 2 || countSecondNeighbor == 2){
                            deleteAndReduceCount(algorithmEdges, firstNeighborEdge, secondNeighborEdge, triangleEdge);
                        }
                        else {
                            ++triangleEdge.count();
                            EdgeInfo* currentEdgeInfo = nullptr;
                            NodePair currentNeighborPair;
                            if (!firstNeighborEdge.is_bfs_edge()){
                                currentEdgeInfo = &firstNeighborEdge;
                                currentNeighborPair = firstNeighborPair;
                            }
                            else if(!secondNeighborEdge.is_bfs_edge()){
                                currentEdgeInfo = &secondNeighborEdge;
                                currentNeighborPair = secondNeighborPair;
                            }
                            if (currentEdgeInfo != nullptr){
                                if (currentEdgeInfo->dependentEdgesPointer.first == nullptr){
                                    dependentEdges.emplace_back(std::vector<NodePair>{currentNeighborPair});
                                    currentEdgeInfo->dependentEdgesPointer.first = &dependentEdges.back();
                                }
                                else{
                                    currentEdgeInfo->dependentEdgesPointer.first->emplace_back(currentNeighborPair);
                                }
                                if (triangleEdge.dependentEdgesPointer.first == nullptr){
                                    triangleEdge.dependentEdgesPointer.first = currentEdgeInfo->dependentEdgesPointer.first;
                                }
                                else{
                                    triangleEdge.dependentEdgesPointer.second = currentEdgeInfo->dependentEdgesPointer.first;
                                }
                            }
                        }
                        if (nodeDegrees[n1] == 2) {
                            degree2Nodes.emplace_back(n1);
                        }
                        if (nodeDegrees[n2] == 2) {
                            degree2Nodes.emplace_back(n2);
                        }
                        if (nodeDegrees[n1] == 1 || nodeDegrees[n2] == 1) {
                            triangleEdge.set_valid(-1);
                            edgeNum -= 1;
                        }
                    } else {
                        componentGraph->AddEdge(n1, n2);
                        removableEdges.emplace_back(triangulationPair);
                        edgeNum += 1;
                        int count = std::max(1, std::max(firstNeighborEdge.count(), secondNeighborEdge.count()));
                        algorithmEdges[triangulationPair] = EdgeInfo(count,
                                                                     true, true, true);

                        if (!firstNeighborEdge.is_bfs_edge()) {
                            algorithmEdges[triangulationPair].dependentEdgesPointer = firstNeighborEdge.dependentEdgesPointer;
                            if (!firstNeighborEdge.is_triangulationEdge()) {
                                if (algorithmEdges[triangulationPair].dependentEdgesPointer.first == nullptr){
                                    dependentEdges.emplace_back(std::vector<NodePair>{firstNeighborPair});
                                    algorithmEdges[triangulationPair].dependentEdgesPointer.first = &dependentEdges[dependentEdges.size() - 1];
                                }
                                else {
                                    algorithmEdges[triangulationPair].dependentEdgesPointer.first->emplace_back(
                                            firstNeighborPair);
                                }
                            }
                        }
                        else if(!secondNeighborEdge.is_bfs_edge()){
                            algorithmEdges[triangulationPair].dependentEdgesPointer = secondNeighborEdge.dependentEdgesPointer;
                            if (!secondNeighborEdge.is_triangulationEdge()) {
                                if (algorithmEdges[triangulationPair].dependentEdgesPointer.first == nullptr){
                                    dependentEdges.emplace_back(std::vector<NodePair>{secondNeighborPair});
                                    algorithmEdges[triangulationPair].dependentEdgesPointer.first = &dependentEdges[dependentEdges.size() - 1];
                                }
                                else {
                                    algorithmEdges[triangulationPair].dependentEdgesPointer.first->emplace_back(
                                            secondNeighborPair);
                                }
                            }
                        }
                    }
                }
            } else {
                while(!removableEdges.empty()){
                    int idx = std::uniform_int_distribution<int>(0, (int) removableEdges.size() - 1)(gen);
                    NodePair nodePair = removableEdges[idx];
                    EdgeInfo& edgeInfo = algorithmEdges[nodePair];
                    if (edgeInfo.is_valid() != false){
                        std::swap(removableEdges[idx], removableEdges.back());
                        removableEdges.pop_back();
                        auto src = nodePair.first(), dst = nodePair.second();
                        deleteEdges(algorithmEdges, edgeInfo);
                        edgeNum-=1;
                        --nodeDegrees[src], --nodeDegrees[dst];
                        if (nodeDegrees[src] == 2){
                            degree2Nodes.emplace_back(src);
                        }
                        if (nodeDegrees[dst] == 2){
                            degree2Nodes.emplace_back(dst);
                        }
                        break;
                    }
                    std::swap(removableEdges[idx], removableEdges.back());
                    removableEdges.pop_back();
                }
            }
        }
        //Add all _valid possibleEdges to the output graph
        for (auto const& [key, edgeInfo] : algorithmEdges) {
            if (!edgeInfo.is_deleted()){
                outerplanarGraph->AddEdge(key.first(), key.second());
            }
        }
    }
}

void GraphFunctions::deleteEdges(std::unordered_map<NodePair, bool, hashNodePair>& resultingEdges, const std::vector<NodePair> &edges) {
    for (const auto& edge : edges) {
        resultingEdges[edge] = false;
    }
}

void GraphFunctions::GetBiconnectedComponents(const PUNGraph& graph, std::vector<PUNGraph>& graphComponents){
    TCnComV components;
    TSnap::GetBiCon(graph, components);
    for (int i = 0; i < components.Len(); ++i) {
        graphComponents.emplace_back(TSnap::GetSubGraph(graph, components[i]()));
    }
}

void GraphFunctions::GetBiconnectedOuterplanarFaces(const PUNGraph &component, OuterplanarComponent& outerplanarComponent) {
    std::vector<std::vector<NodeId>> neighbors = std::vector<std::vector<NodeId>>(component->GetNodes());
    std::vector<PUNGraph> currentFace;
    int face_num = 0;
    Nodes degree2Nodes;
    std::vector<int> degrees = std::vector<int>(component->GetNodes(), 0);
    //Preprocessing on get_node degrees
    for (TUNGraph::TNodeI node = component->BegNI(); node != component->EndNI(); node++) {
        int id = node.GetId();
        int degree = node.GetDeg();
        if (degree == 2) {
            degree2Nodes.emplace_back(id);
        }
        neighbors[id] = {node.GetNbrNId(0), node.GetNbrNId(1)};
        degrees[id] = degree;
    }
    outerplanarComponent.faces.emplace_back(new TUNGraph());
    currentFace.emplace_back(outerplanarComponent.faces.back());
    while (!degree2Nodes.empty() && component->GetEdges() > 2){
        int currentNodeId = degree2Nodes.back();
        degree2Nodes.pop_back();
        if (degrees[currentNodeId] > 0) {
            degrees[currentNodeId] = 0;
            int n1 = neighbors[currentNodeId][0];
            int n2 = neighbors[currentNodeId][1];
            TUNGraph::TNodeI nearNode = component->GetNI(n1);
            TUNGraph::TNodeI nextNode = component->GetNI(n2);

            int near_degree = degrees[n1];
            int next_degree = degrees[n2];
            if (component->IsEdge(n1, n2)) {
                currentFace.back()->AddNode(n1);
                currentFace.back()->AddNode(n2);
                currentFace.back()->AddEdge(n1, n2);
                ++face_num;
                currentFace.pop_back();
                //Update degreeTwo nodes
                --degrees[n1];
                --degrees[n2];
            } else {
                currentFace.back()->AddNode(currentNodeId);
                currentFace.back()->AddNode(n1);
                currentFace.back()->AddNode(n2);
                currentFace.back()->AddEdge(currentNodeId, n1);
                currentFace.back()->AddEdge(currentNodeId, n2);
            }
            if (degrees[n1] == 2 || degrees[n2] == 2) {
                if (degrees[n1] == 2) {
                    degree2Nodes.emplace_back(n1);
                    for (int i = 0; i < nearNode.GetDeg(); ++i) {
                        if (nearNode.GetNbrNId(i) > 0) {
                            neighbors[n1].emplace_back(i);
                        }
                        if (neighbors[n1].size() == 2) {
                            break;
                        }
                    }
                } else {
                    degree2Nodes.emplace_back(n2);
                    for (int i = 0; i < nextNode.GetDeg(); ++i) {
                        if (nextNode.GetNbrNId(i) > 0) {
                            neighbors[n2].emplace_back(i);
                        }
                        if (neighbors[n2].size() == 2) {
                            break;
                        }
                    }
                }
            }
            else {
                outerplanarComponent.faces.emplace_back(new TUNGraph());
                currentFace.emplace_back(outerplanarComponent.faces.back());
            }
        }
    }

}

void GraphFunctions::GetBiconnectedOuterplanarFaceNum(const PUNGraph &component, int &face_num) {
    face_num = 0;
    Nodes degree2Nodes;
    //Preprocessing on get_node degrees
    for (TUNGraph::TNodeI node = component->BegNI(); node != component->EndNI(); node++) {
        int id = node.GetId();
        int degree = node.GetDeg();
        if (degree == 2) {
            degree2Nodes.emplace_back(id);
        }
    }
    while (!degree2Nodes.empty() && component->GetEdges() > 2){
        int currentNodeId = degree2Nodes.back();
        TUNGraph::TNodeI currentNode = component->GetNI(currentNodeId);
        degree2Nodes.pop_back();
        int n1 = currentNode.GetNbrNId(0);
        TUNGraph::TNodeI nearNode = component->GetNI(n1);
        int n2 = currentNode.GetNbrNId(1);
        TUNGraph::TNodeI nextNode = component->GetNI(n2);
        //Delete possibleEdges from component
        component->DelEdge(currentNodeId, n1);
        component->DelEdge(currentNodeId, n2);
        if (component->IsEdge(n1, n2)) {
            ++face_num;
            //Update degreeTwo nodes
            int deg1 = component->GetNI(n1).GetDeg();
            int deg2 = component->GetNI(n2).GetDeg();
            if (deg1 == 2){
                degree2Nodes.emplace_back(n1);
            }
            if (deg2 == 2){
                degree2Nodes.emplace_back(n2);
            }
        }
        else{
            component->AddEdge(n1, n2);
        }
    }
}

bool GraphFunctions::IsTree(const PUNGraph &graph) {
    return TSnap::IsConnected(graph) && graph->GetEdges() == graph->GetNodes() - 1;
}

void GraphFunctions::getValidNeighbor(PUNGraph &graph, std::unordered_map<NodePair, EdgeInfo, hashNodePair>& algorithmEdges, int nodeId, int neighborIdx, NodeId &neighborId) {
    neighborId = -1;
    auto node = graph->GetNI(nodeId);
    int validIdx = -1;
    for (int i = 0; i < node.GetDeg(); ++i) {
        neighborId = node.GetNbrNId(i);
        if(algorithmEdges[NodePair(nodeId, neighborId)].is_valid()){
            ++validIdx;
        }
        if (validIdx == neighborIdx){
            return;
        }
    }
}

std::vector<NodeId> GraphFunctions::nodesOnAllShortestPaths(const GraphData& graphData, NodeId A, NodeId B){
    std::vector<int> distances = std::vector<int>(graphData.get_graph()->GetNodes(), 0);

    return std::vector<NodeId>();
}

int GraphFunctions::deleteEdges(std::unordered_map<NodePair, EdgeInfo, hashNodePair>& algorithmEdges, EdgeInfo& currentEdge) {
    int numDeleted = 0;
    currentEdge.delete_edge();
    if (currentEdge.dependentEdgesPointer.first != nullptr){
        for(auto const & edge: *currentEdge.dependentEdgesPointer.first){
            algorithmEdges[edge].delete_edge();
        }
        ++numDeleted;
    }
    if (currentEdge.dependentEdgesPointer.second != nullptr){
        for(auto const & edge: *currentEdge.dependentEdgesPointer.first){
            algorithmEdges[edge].delete_edge();
        }
        ++numDeleted;
    }
    return numDeleted;
}

void GraphFunctions::bfsRandomSubtree(const PUNGraph &graph, std::vector<NodePair>& edges, std::mt19937_64 &gen) {
    NodeId randNodeId = std::uniform_int_distribution<int>(0, graph->GetNodes() - 1)(gen);
    bfsSubtree(graph, edges, randNodeId, gen);
}

void GraphFunctions::deleteAndReduceCount(std::unordered_map<NodePair, EdgeInfo, hashNodePair>& algorithmEdges,
                                          GraphFunctions::EdgeInfo &FirstEdge, GraphFunctions::EdgeInfo &SecondEdge,
                                          GraphFunctions::EdgeInfo &TriangleEdge) {

    bool canDeleteFirst = !FirstEdge.is_bfs_edge();
    bool canDeleteSecond = !SecondEdge.is_bfs_edge();
    bool canDeleteTriangle = !TriangleEdge.is_bfs_edge();

    int reduceFirstCount = (bool)FirstEdge.dependentEdgesPointer.first + (bool)FirstEdge.dependentEdgesPointer.second;
    int reduceSecondCount = (bool)SecondEdge.dependentEdgesPointer.first + (bool)SecondEdge.dependentEdgesPointer.second;
    int reduceTriangleCount = (bool)TriangleEdge.dependentEdgesPointer.first + (bool)TriangleEdge.dependentEdgesPointer.second;

    if (canDeleteFirst && reduceFirstCount >= FirstEdge.get_count()){
        deleteEdges(algorithmEdges, FirstEdge);
    }
    else if (canDeleteSecond && reduceSecondCount >= SecondEdge.get_count()){
        deleteEdges(algorithmEdges, SecondEdge);
    }
    else if (canDeleteTriangle && reduceTriangleCount >= TriangleEdge.get_count()){
        deleteEdges(algorithmEdges, TriangleEdge);
    }
}

PUNGraph GraphFunctions::dfsRandomTree(const PUNGraph &graph, std::vector<NodeId>& visitedNodes, std::mt19937_64 &gen) {
    PUNGraph dfsTree = new TUNGraph();
    for (auto node = graph->BegNI(); node != graph->EndNI(); node++) {
        dfsTree->AddNode();
    }
    return dfsTree;
}

PUNGraph GraphFunctions::maximalOuterplanarSubgraphBFSGreedy(const GraphData &graph, std::mt19937_64 &gen) {
    return maximalOuterplanarSubgraphBFSGreedy(graph.get_graph(), gen);
}

void GraphFunctions::print(const PUNGraph &graph, const std::string& name) {
    std::cout << name << std::endl;
    for(auto edge=graph->BegEI(); edge!=graph->EndEI(); edge++){
        std::cout << " " << edge.GetSrcNId() << "-" << edge.GetDstNId() << " ";
    }
    std::cout << std::endl;

}

void GraphFunctions::NormalizeIdsOfGraphs(const std::vector<std::string>& paths) {
    for(const auto& path : paths) {
        GraphData graph = GraphData(path);
        PUNGraph n_graph = ResetGraphIds(graph.get_graph());
        GraphData norm_graph = GraphData(n_graph);
        std::string new_path = std::filesystem::path(path).replace_extension();
        norm_graph.save(new_path);
    }
}

int GraphFunctions::getMaxDegree(const PUNGraph &graph) {
    int maxDegree = 0;
    for (auto NodeIt = graph->BegNI(); NodeIt != graph->EndNI(); NodeIt++) {
        maxDegree = std::max(NodeIt.GetDeg(), maxDegree);
    }
    return maxDegree;
}

int GraphFunctions::getMaxDegree(const SIMPLE_GRAPH& graph) {
    int maxDegree = 0;
    for (auto const & n : graph) {
        maxDegree = std::max(static_cast<int>(n.size()), maxDegree);
    }
    return maxDegree;
}

void GraphFunctions::generateNeighborVector(const PUNGraph& graph, std::vector<int>& neighborVector)
{
    int maxDegree = GraphFunctions::getMaxDegree(graph);
    neighborVector.clear();
    neighborVector.resize(maxDegree,0);
    std::iota(neighborVector.begin(), neighborVector.end(), 0);
}

void GraphFunctions::generateNeighborVector(const SIMPLE_GRAPH& graph, std::vector<int>& neighborVector)
{
    int maxDegree = GraphFunctions::getMaxDegree(graph);
    neighborVector.clear();
    neighborVector.resize(maxDegree,0);
    std::iota(neighborVector.begin(), neighborVector.end(), 0);
}

void GraphFunctions::GetLargestComponent(GraphData &data, bool renumber) {
    TCnComV ConComps;
    TSnap::GetWccs(data.get_graph(), ConComps);
    if (ConComps.Len() == 1) {
        std::cout << "Graph " << data.getName() << " is connected!" << std::endl;
    } else {
        std::cout << "Graph " << data.getName() << " is not connected!" << std::endl;
        std::cout << "Get biggest connected component!" << std::endl;

        TCnCom* biggestComponent = nullptr;
        int Components = ConComps.Len();
        for (int i = 0; i < Components; ++i) {
            if (biggestComponent == nullptr || biggestComponent->NIdV.Len() < ConComps[i].NIdV.Len()){
                biggestComponent = &ConComps[i];
            }
        }
        PUNGraph biggestComponentGraph = TSnap::GetSubGraph(data.get_graph(), biggestComponent->NIdV, renumber);
        GraphData Component = GraphData(GraphFunctions::ResetGraphIds(biggestComponentGraph));
        Component.setName(data.getName() + "_component");
        data = Component;
    }
}

void GraphFunctions::analyse_graph(const PUNGraph& graph, const std::string& name, bool component_analysis, FileEvaluation* fileEvaluation) {
    std::cout << std::endl;
    std::cout << name << ":" << std::endl;
    std::cout << "\tSize: " << graph->GetNodes() << std::endl;
    std::cout << "\tEdges: " << graph->GetEdges() << std::endl;
    std::map<int, int> degree_distribution;
    GraphFunctions::getDegreeDistribution(graph, degree_distribution);
    std::cout << "\tDegree distribution: " << StaticFunctions::printMap(degree_distribution) << std::endl;
    TCnComV conComponents;
    TSnap::GetWccs(graph, conComponents);
    std::cout << "\tConnected Components: " << conComponents.Len() << std::endl;
    if (fileEvaluation){
        //Graph values
        fileEvaluation->headerValueInsert({"Graph", "Nodes", "Edges", "Density"},
                                          {name, std::to_string(graph->GetNodes()), std::to_string(graph->GetEdges()),
                                           std::to_string((double) graph->GetEdges()/((double) graph->GetNodes()* (double) graph->GetNodes()))});
    }
    if (component_analysis) {
        // Component analysis
        for (int i = 0; i < conComponents.Len(); ++i) {
            TCnCom &com = conComponents[i];
            PUNGraph comGraph = TSnap::GetSubGraph(graph, com.NIdV, true);
            std::cout << "\tComponent analysis: " << std::endl;
            std::cout << "\t\t" << "Component " << i + 1 << std::endl;
            std::cout << "\t\t" << "Size:" << comGraph->GetNodes() << std::endl;
            std::cout << "\t\t" << "Edges:" << comGraph->GetEdges() << std::endl;
            std::cout << "\t\t" << "Diameter:" << TSnap::GetBfsFullDiam(comGraph, com.NIdV.Len()) << std::endl;

            GraphFunctions::getDegreeDistribution(comGraph, degree_distribution);
            std::cout << "\t\tComponent degree distribution" << std::endl;
            std::cout << StaticFunctions::printMap(degree_distribution) << std::endl;
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

PUNGraph GraphFunctions::RemoveSubgraph(const PUNGraph& baseGraph, const PUNGraph& subGraph) {
    PUNGraph newGraph = new TUNGraph();
    for (auto NodeIt = baseGraph->BegNI(); NodeIt != baseGraph->EndNI(); NodeIt++) {
        if (!subGraph->IsNode(NodeIt.GetId())){
            newGraph->AddNode(NodeIt.GetId());
        }
    }
    for (auto EdgeIt = baseGraph->BegEI(); EdgeIt != baseGraph->EndEI(); EdgeIt++) {
        if (newGraph->IsNode(EdgeIt.GetSrcNId()) && newGraph->IsNode(EdgeIt.GetDstNId())){
            newGraph->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
        }
    }
    return GraphFunctions::ResetGraphIds(newGraph);
}

void GraphFunctions::getDegreeDistribution(const PUNGraph& graph, std::map<int, int> &degree_distribution, std::vector<NodeId>* nodes) {
    degree_distribution.clear();
    if (nodes){
        for (auto const Node : *nodes) {
            auto NodeIt = graph->GetNI(Node);
            int degree = NodeIt.GetDeg();
            if (degree_distribution.find(degree) == degree_distribution.end()) {
                degree_distribution[degree] = 1;
            } else {
                degree_distribution[degree] += 1;
            }
        }
    }
    else {
        for (auto NodeIt = graph->BegNI(); NodeIt != graph->EndNI(); NodeIt++) {
            int degree = NodeIt.GetDeg();
            if (degree_distribution.find(degree) == degree_distribution.end()) {
                degree_distribution[degree] = 1;
            } else {
                degree_distribution[degree] += 1;
            }
        }
    }

}

void GraphFunctions::addRandomEdges(GraphData &data, int edgeNum, int seed) {
    std::mt19937_64 gen(seed);
    if (edgeNum > (data.nodes() * data.nodes() - 1) / 2 - data.edges()){
        throw std::range_error("Cannot add " + std::to_string(edgeNum) + " edges! Number of edges must be smaller than " + std::to_string((data.nodes() * data.nodes() - 1) / 2 - data.edges()));
    }
    for (int i = 0; i < edgeNum; ++i) {
        int src = std::uniform_int_distribution<int>(0, data.nodes() - 1)(gen);
        int dst = std::uniform_int_distribution<int>(0, data.nodes() - 1)(gen);
        if (src != dst && !data.get_graph()->IsEdge(src, dst)){
            data.graph()->AddEdge(src, dst);
        }
        else{
            --i;
        }
    }
}


void GraphFunctions::checkingOuterplanarity(const PUNGraph &graph, const PUNGraph &outerplanarSubgraph, int &notOuterplanarSubgraphs,
                                            int &notMaximalSubgraphs, std::vector<int> &nonOuterplanarSeeds,
                                            std::vector<int> &nonMaximalSeeds, std::vector<double> &algorithmMissingEdges,
                                            std::vector<double> &maximalEdges, int seed) {
//GraphFunctions::print(_subgraph);
    bool outerplanar = GraphFunctions::IsOuterPlanar(outerplanarSubgraph);
    std::vector<NodePair> missingEdges;
    bool maximalOuterplanar = GraphFunctions::IsMaximalOuterplanarSubgraph(graph, outerplanarSubgraph, missingEdges);
    std::string answer;
    std::string maximalAnswer;
    if (outerplanar) {
        answer = "Yes";
    } else {
        ++notOuterplanarSubgraphs;
        nonOuterplanarSeeds.emplace_back(seed);
        answer = "No";
    }
    if (maximalOuterplanar) {
        maximalAnswer = "Yes";

    } else {
        ++notMaximalSubgraphs;
        nonMaximalSeeds.emplace_back(seed);
        maximalAnswer = "No";
    }
    algorithmMissingEdges.emplace_back(missingEdges.size());
//    if (printLevel == PrintLevel::ALL) {
//        std::cout << "Outerplanar: " << answer << ", maximal: " << maximalAnswer << std::endl;
//
//        if (!maximalOuterplanar) {
//            std::cout << "Missing Edges ";
//            for (auto const &edge: missingEdges) {
//                std::cout << "(" << edge.first() << ", " << edge.second() << ")";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "//////////////////////////////////////" << std::endl;
//    }
    maximalEdges.emplace_back(outerplanarSubgraph->GetEdges() + missingEdges.size());
}

void GraphFunctions::GetSamples(GraphData &graph, PatternType type, std::vector<GraphData> &subgraphs,
                                 OuterplanarSubgraphDFS *outerplanarSubgraphDFS, std::vector<NodeId>& neighborIds, int samples, int seed, double& runtime, bool save_samples, bool p) {
    runtime = 0;
    subgraphs.clear();
    if (save_samples) {
#pragma omp critical
        {
            if (std::filesystem::exists("../out/Samples/")) {
                if (!std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                    std::filesystem::create_directory("../out/Samples/" + graph.getName());
                }
            }
            for (const auto &entry: std::filesystem::recursive_directory_iterator(
                    "../out/Samples/" + graph.getName() + "/")) {
                std::filesystem::remove(entry);
            }
        }
    }
    std::chrono::time_point<std::chrono::system_clock> start_generation = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < samples; ++j) {
        std::mt19937_64 generator(j + samples * seed);
        GraphData* subgraph_data = nullptr;
        GraphData data;
        switch (type) {
            case PatternType::BFS_TREE:
                if (save_samples && !SmallGraph(graph)) {
                    data = GraphData(new TUNGraph, (int) graph.size());
                    subgraph_data = &data;
                }
                else{
                    subgraphs.emplace_back(GraphData(new TUNGraph, (int) graph.size()));
                    subgraph_data = &subgraphs.back();
                }
                GraphFunctions::bfsSubtree(graph.get_graph(), *subgraph_data, neighborIds, generator);
                subgraph_data->graphType = GraphType::TREE;
                subgraph_data->setName("tree_sample_" + std::to_string(j));
                if (save_samples && !SmallGraph(graph)) {
                    if (std::filesystem::exists("../out/Samples/")) {
                        if (!std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                            std::filesystem::create_directory("../out/Samples/" + graph.getName());
                        }
                        subgraph_data->save_bin("../out/Samples/" + graph.getName() + "/");
                    }
                }
                if (p) {
                    std::cout << "Generated tree sample " << j + 1 << "/" << samples << " of " << graph.getName()
                              << std::endl;
                }
                break;
            case PatternType::OUTERPLANAR:
                if (save_samples && !SmallGraph(graph)) {
                    data = GraphData(new TUNGraph, (int) graph.size());
                    subgraph_data = &data;
                }
                else{
                    subgraphs.emplace_back(GraphData(new TUNGraph, (int) graph.size()));
                    subgraph_data = &subgraphs.back();
                }
                outerplanarSubgraphDFS->generate(subgraphs.back(), generator, false);
                subgraph_data->graphType = GraphType::OUTERPLANAR;
                subgraph_data->setName("outerplanar_sample_" + std::to_string(j));
                if (save_samples && !SmallGraph(graph)) {
                    if (std::filesystem::exists("../out/Samples/")) {
                        if (!std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                            std::filesystem::create_directory("../out/Samples/" + graph.getName());
                        }
                        subgraph_data->save_bin("../out/Samples/" + graph.getName() + "/");
                    }
                }
                if (p) {
                    std::cout << "Generated outerplanar sample " << j + 1 << "/" << samples << " of " << graph.getName()
                              << std::endl;
                }
                break;
        }
    }
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start_generation).count() /
               1000000.0);
}

void GraphFunctions::GetOuterplanarSamples(GraphData &graph, PatternType type, std::vector<OuterplanarGraphData> &subgraphs,
                                           OuterplanarSubgraphDFS *outerplanarSubgraphDFS, int samples, int seed, double &runtime, double &conversion_runtime, bool clear, bool save_samples, bool p) {
    runtime = 0;
    conversion_runtime = 0;
    subgraphs.clear();
    auto start_generation = std::chrono::high_resolution_clock::now();
#pragma omp critical
    {
        if (save_samples && !SmallGraph(graph)) {
            if (std::filesystem::exists("../out/Samples/")) {
                if (!std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                    std::filesystem::create_directory("../out/Samples/" + graph.getName());
                }
            }
            for (const auto &entry: std::filesystem::recursive_directory_iterator(
                    "../out/Samples/" + graph.getName() + "/")) {
                std::filesystem::remove(entry);
            }
        }
    }
    for (int j = 0; j < samples; ++j) {
        OuterplanarGraphData* subgraph_data = nullptr;
        OuterplanarGraphData data;
        std::mt19937_64 generator(j + samples * seed);
        if (save_samples && !SmallGraph(graph)){
            data = OuterplanarGraphData(new TUNGraph, (int) graph.size());
            subgraph_data = &data;
        }
        else {
            subgraphs.emplace_back(OuterplanarGraphData(new TUNGraph, (int) graph.size()));
            subgraph_data = &subgraphs.back();
        }
        outerplanarSubgraphDFS->generate(*subgraph_data, generator, false);
        auto conversion = std::chrono::high_resolution_clock::now();
        if (!save_samples || SmallGraph(graph)) {
            subgraph_data->set();
            if (clear) {
                subgraph_data->graph().Clr();
                subgraph_data->ContainmentList.clear();
                subgraph_data->DistanceList.clear();
                subgraph_data->neighborList.clear();
                subgraph_data->PredecessorsList.clear();
            }
        }
        conversion_runtime += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - conversion).count() /
                   1000000.0);
        subgraph_data->graphType = GraphType::OUTERPLANAR;
        subgraph_data->setName("outerplanar_sample_" + std::to_string(j));
        if (save_samples && !SmallGraph(graph)){
            if (std::filesystem::exists("../out/Samples/")) {
                if (!std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                    std::filesystem::create_directory("../out/Samples/" + graph.getName());
                }
                subgraph_data->save_bin("../out/Samples/" + graph.getName() + "/");
            }
        }
        std::cout << "Generated outerplanar sample " << j << "/" << samples << " of " << graph.getName() << std::endl;
    }
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start_generation).count() /
               1000000.0);
    runtime -= conversion_runtime;
}

void GraphFunctions::GetCoreGraphStats(GraphData &graph, TIntV &coreNodes, double runtime, const std::string& core_name, FileEvaluation* fileEvaluation) {
    //Statistics of core graph
    GraphData core_graph = GraphData(TSnap::GetSubGraph(graph.get_graph(), coreNodes));

    std::cout << std::endl;
    std::cout << "Core Subgraph: " << std::endl;
    std::cout << "\tSize: " << core_graph.size() << std::endl;
    std::cout << "\tEdges: " << core_graph.get_graph()->GetEdges() << std::endl;

    int overlapOutEdges = 0;

    for(auto NodeIt = graph.get_graph()->BegNI(); NodeIt != graph.get_graph()->EndNI(); NodeIt++){
        NodeId Id = NodeIt.GetId();
        if(!core_graph.get_graph()->IsNode(Id)) {
            int degree = NodeIt.GetDeg();
            for (int n = 0; n < degree; ++n) {
                int NeighborId = NodeIt.GetNbrNId(n);
                if (core_graph.get_graph()->IsNode(NeighborId)){
                    ++overlapOutEdges;
                }
            }
        }
    }

    std::map<int, int> degree_distribution;
    graph.getDegreeDistribution(degree_distribution, coreNodes);
    std::cout << "\tOverlap degree distribution: " << StaticFunctions::printMap(degree_distribution) << std::endl;
    std::cout << std::endl;
    std::cout << "\tOverlap out edges: " << overlapOutEdges << std::endl;

    if (fileEvaluation != nullptr){
        //Core Values
        fileEvaluation->headerValueInsert({core_name + " Nodes", core_name + " Relative Nodes", core_name + " Edges", core_name + " Relative Edges", core_name + " Out Edges"},
                                          {std::to_string(core_graph.size()), std::to_string((double) core_graph.size()/(double) graph.size()),
                                           std::to_string(core_graph.get_graph()->GetEdges()),
                                           std::to_string((double) core_graph.get_graph()->GetEdges()/graph.get_graph()->GetEdges()),
                                           std::to_string(overlapOutEdges)}, -2);
    }

}






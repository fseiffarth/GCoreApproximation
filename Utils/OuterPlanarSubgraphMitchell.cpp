//
// Created by anonymous on 29.06.21.
//

#include <list>
#include "OuterPlanarSubgraphMitchell.h"
#include "GraphFunctions.h"

OuterPlanarSubgraphMitchell::OuterPlanarSubgraphMitchell(const PUNGraph graph)  : OuterplanarSubgraph(graph) {
    GetBiconnectedComponents();
    Init();
}

AlgorithmEdge::AlgorithmEdge(const NodePair& nodePair, bool type)
        : _nodePair(nodePair), _type(type){}

bool AlgorithmEdge::original_edge() const {
    return _type;
}


void OuterPlanarSubgraphMitchell::getValidNeighbor(PUNGraph graph, std::unordered_map<NodePair, AlgorithmEdge, hashNodePair>& algEdges, int nodeId, int neighborIdx, NodeId &neighborId) {
    neighborId = -1;
    auto node = graph->GetNI(nodeId);
    int validIdx = -1;
    for (int i = 0; i < node.GetDeg(); ++i) {
        neighborId = node.GetNbrNId(i);
        if(algEdges[NodePair(nodeId, neighborId)].is_valid()){
            ++validIdx;
        }
        if (validIdx == neighborIdx){
            return;
        }
    }
}

void OuterPlanarSubgraphMitchell::Init() {
    int componentCounter = 0;
    nodeToComponentAndId = std::vector<std::list<int>>(this->_graph->GetNodes(), std::list<int>());
    for (auto const & component : biconnectedComponents) {
        biComponents.emplace_back(BiconnectedComponent());
        BiconnectedComponent& biCom = biComponents.back();
        biCom.nodeIds = std::vector<int>();
        biCom.nodeDegrees = std::vector<int>();
        biCom.algorithmEdges = std::unordered_map<NodePair, AlgorithmEdge, hashNodePair>();
        biCom.edges = std::unordered_map<NodePair, bool, hashNodePair>();
        biCom.degree2Nodes = std::vector<int>();
        biCom.componentId = componentCounter;
        int componentNodeId = 0;
        for (auto node = component->BegNI(); node != component->EndNI(); node++) {
            nodeToComponentAndId[node.GetId()].push_back(componentNodeId);
            biCom.nodeIds.emplace_back(node.GetId());
            biCom.nodeDegrees.emplace_back(node.GetDeg());
            if (node.GetDeg() == 2){
                biCom.degree2Nodes.emplace_back(node.GetId());
            }
            ++componentNodeId;
        }
        for (auto edge = component->BegEI(); edge != component->EndEI(); edge++) {
            NodePair nodePair = NodePair(edge.GetSrcNId(), edge.GetDstNId());
            biCom.algorithmEdges[nodePair] = AlgorithmEdge(nodePair);
            biCom.edges[nodePair] = true;
        }
        ++componentCounter;
    }
}

void OuterPlanarSubgraphMitchell::GetBiconnectedComponents() {
    TCnComV components;
    TSnap::GetBiCon(_graph, components);
    for (auto component = components.BegI(); component != components.EndI(); component++) {
        PUNGraph subgraph = new TUNGraph();
        subgraph = TSnap::GetSubGraph(_graph, (*component)());
        if (subgraph->GetEdges() == subgraph->GetNodes() - 1){
            //Add componentGraph possibleEdges to algorithmGraph
            for (auto edge = subgraph->BegEI(); edge != subgraph->EndEI(); edge++) {
                NodeId src = edge.GetSrcNId();
                NodeId dst = edge.GetDstNId();
                treeEdges.emplace_back(NodePair(src, dst));
            }
        }
        else{
            this->biconnectedComponents.emplace_back(subgraph);
        }
    }

}

void OuterPlanarSubgraphMitchell::generate(GraphData& subgraph, std::mt19937_64& gen, bool p) {
    this->_gen = gen;
    this->print = p;



    //Run spanning tree to mark non-deletable edges, TODO give spanning tree as parameter
    std::vector<NodePair> spanningTreeEdges;
    GraphFunctions::bfsRandomSubtree(_graph, spanningTreeEdges, _gen);
    std::unordered_map<NodePair, bool, hashNodePair> edgesInSpanningTree;
    if (this->print){
        std::cout << "Spanning Tree: " << std::endl;
    }
    for (auto const edge : spanningTreeEdges) {
        edgesInSpanningTree[edge] = true;
        if (this->print) {
            edge.print();
        }
    }
    if (this->print){
        std::cout << std::endl;
        std::cout  << std::endl << "\t" << "*************" << " Start Algorithm" << "******************" << std::endl;
    }
    for (BiconnectedComponent& biCom : biComponents) {
        for(auto & [nodePair, edgeInfo] : biCom.algorithmEdges){
            bool b = edgesInSpanningTree.find(nodePair) != edgesInSpanningTree.end();
            edgeInfo.setTreeEdge(b);
        }
        BiconnectedComponentSampling(subgraph.get_graph(), biCom);
        for (auto nodeId : biCom.nodeIds) {
            nodeToComponentAndId[nodeId].pop_front();
        }
    }
}

void OuterPlanarSubgraphMitchell::BiconnectedComponentSampling(const PUNGraph subgraph, BiconnectedComponent& biCom) {
    PUNGraph component = biconnectedComponents[biCom.componentId];
    int nodeNum = component->GetNodes();
    currentComponent = biCom;

    //printing of the algorithm
    if (this->print){
        std::cout  << std::endl << "\t" << "*************" << " New Component" << "******************" << std::endl;
    }

    removableEdges.clear();
    //Find all removable possibleEdges
    if (this->print){
        std::cout << "\t" << "Removable Edges:" << std::endl << "\t";
    }
    for (auto const& [nodePair, edgeInfo]: currentComponent.algorithmEdges) {
        if (!edgeInfo.isTreeEdge()) {
            removableEdges.emplace_back(nodePair);
            if (this->print){
                nodePair.print();
            }
        }
    }
    if (this->print){
        std::cout << std::endl;
    }

    for (auto edge = component->BegEI(); edge != component->EndEI(); edge++) {
        triangulationForest.AddNode(TriangulationForest::NodeType::EDGE, &currentComponent.algorithmEdges[NodePair(edge.GetSrcNId(), edge.GetDstNId())]);
    }
    for (auto node = component->BegNI(); node != component->EndNI(); node++) {
        triangulationForest.AddNode(TriangulationForest::NodeType::NODE, nullptr, node.GetId());
    }
    auto printing = PrintSteps(triangulationForest);
    while (nodeNum > 0){
        std::vector<NodeId>& degree2Nodes = currentComponent.degree2Nodes;
        if (this->print){
            std::cout << "\t" << "Next Step: " << std::endl;
            std::cout << "\t\t Node Degrees: ";
            for (int i = 0; i < currentComponent.nodeIds.size(); ++i) {
                std::cout << "(" << currentComponent.nodeIds[i] << ", " << currentComponent.nodeDegrees[i] << " deg)";
            }
            std::cout << std::endl;
            std::cout << "\t\t Degree 2 Nodes: ";
            for (auto nodeId : degree2Nodes) {
                std::cout << nodeId << ", ";
            }
            std::cout << std::endl;
        }
        if (!degree2Nodes.empty()){
            int rnd = std::uniform_int_distribution<int>(0, (int) degree2Nodes.size() - 1)(_gen);
            int v = degree2Nodes[rnd];
            std::swap(degree2Nodes[rnd], degree2Nodes.back());
            degree2Nodes.pop_back();
            if (currentComponent.nodeDegrees[nodeIdToComponentId(v)] == 2){
                currentComponent.nodeDegrees[nodeIdToComponentId(v)] = 0;
                if (this->print){
                    std::cout << "\t\t\t" << "Triangulate with get_node " << v << std::endl;
                }
                --nodeNum;
                NodeId x, y;
                getValidNeighbor(component, currentComponent.algorithmEdges, v, 0, x);
                getValidNeighbor(component, currentComponent.algorithmEdges, v, 1, y);
                NodePair edgeA = NodePair(v, x);
                NodePair edgeB = NodePair(v, y);
                AlgorithmEdge& algorithmEdgeA = currentComponent.algorithmEdges[edgeA];
                AlgorithmEdge& algorithmEdgeB = currentComponent.algorithmEdges[edgeB];
                NodePair triangulationEdge = NodePair(x, y);
                CollapseTriangle(component, algorithmEdgeA, algorithmEdgeB, triangulationEdge, printing);
                UpdateTriangulationForest(triangulationEdge, edgeA, edgeB, v);
            }
            else{
                --nodeNum;
                if (currentComponent.nodeDegrees[nodeIdToComponentId(v)] == 1){
                    NodeId x;
                    getValidNeighbor(component, currentComponent.algorithmEdges, v, 0, x);
                    --currentComponent.nodeDegrees[nodeIdToComponentId(x)];
                    if (currentComponent.nodeDegrees[nodeIdToComponentId(x)] == 0){
                        --nodeNum;
                    }
                }
                currentComponent.nodeDegrees[nodeIdToComponentId(v)] = 0;
                if (print){
                    printing.printRemovedNode(v);
                }
            }

        }
        else {
            DeleteRandomEdge(printing);
        }
        if (this->print) {
            printing.printTriangulationForest();
        }
    }
    //Add all undeleted edges to the _subgraph
    for (auto const& [edge, valid] : currentComponent.edges) {
        if (valid) {
            subgraph->AddEdge(edge.first(), edge.second());
        }
    }
}

void OuterPlanarSubgraphMitchell::DeleteRandomEdge(const PrintSteps& printing) {
    while(!removableEdges.empty()){
        //Get random edge from all the edges which can be removed
        int idx = std::uniform_int_distribution<int>(0, (int) removableEdges.size() - 1)(_gen);
        NodePair nodePair = removableEdges[idx];
        std::swap(removableEdges[idx], removableEdges.back());
        removableEdges.pop_back();

        //Get corresponding algorithm edge and check if it is valid
        AlgorithmEdge& algorithmEdge = currentComponent.algorithmEdges[nodePair];
        if (algorithmEdge.is_valid()) {
            auto src = nodePair.first(), dst = nodePair.second();
            auto componentSrc = nodeIdToComponentId(src), componentDst = nodeIdToComponentId(dst);
            int degree = triangulationForest.degree(algorithmEdge);
            if (degree == algorithmEdge.triangleCount()){
                bool validPath = triangulationForest.GetValidDeletePath(nodePair, _path1, _gen);
                if (validPath) {
                    if (degree == 2) {
                        if (triangulationForest.GetValidDeletePath(nodePair, _path2, _gen) && _path2.size() > 1) {
                            if (print){
                                printing.printEdge(nodePair);
                            }
                            DeleteEdges(_path1, printing, true);
                            DeleteEdges(_path2, printing, true);
                            --currentComponent.nodeDegrees[componentSrc];
                            --currentComponent.nodeDegrees[componentDst];
                            if (currentComponent.nodeDegrees[componentSrc] == 2) {
                                currentComponent.degree2Nodes.emplace_back(src);
                            }
                            if (currentComponent.nodeDegrees[componentDst] == 2) {
                                currentComponent.degree2Nodes.emplace_back(dst);
                            }
                            break;
                        }
                        else{
                            triangulationForest.AddPath(_path1);
                        }
                    } else {
                        if (print){
                            printing.printEdge(nodePair);
                        }
                        DeleteEdges(_path1, printing, true);
                        --currentComponent.nodeDegrees[componentSrc];
                        --currentComponent.nodeDegrees[componentDst];
                        if (currentComponent.nodeDegrees[componentSrc] == 2) {
                            currentComponent.degree2Nodes.emplace_back(src);
                        }
                        if (currentComponent.nodeDegrees[componentDst] == 2) {
                            currentComponent.degree2Nodes.emplace_back(dst);
                        }
                        break;
                    }
                }
            }
        }
    }
}

void OuterPlanarSubgraphMitchell::UpdateDegrees(NodeId id) {
    int componentId = nodeIdToComponentId(id);
    --currentComponent.nodeDegrees[componentId];
    if (currentComponent.nodeDegrees[componentId] == 2){
        currentComponent.degree2Nodes.emplace_back(id);
    }
}

void OuterPlanarSubgraphMitchell::CollapseTriangle(PUNGraph component, AlgorithmEdge& edgeA, AlgorithmEdge& edgeB, NodePair& triangulationPair, const PrintSteps& printing) {
    AlgorithmEdge* triangulationEdge = nullptr;
    edgeA.set_valid(false);
    edgeB.set_valid(false);
    if (currentComponent.algorithmEdges.find(triangulationPair) != currentComponent.algorithmEdges.end()){
        triangulationEdge = &currentComponent.algorithmEdges[triangulationPair];
    }

    if (edgeA.triangleCount() > 1){
        triangulationForest.GetValidPath(edgeA.edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        edgeA.set_count(edgeA.triangleCount() - 1);
    }
    if (edgeB.triangleCount() > 1){
        triangulationForest.GetValidPath(edgeB.edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        edgeB.set_count(edgeB.triangleCount() - 1);
    }
    if (triangulationEdge != nullptr && triangulationEdge->triangleCount() > 1){
        triangulationForest.GetValidPath(triangulationEdge->edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        triangulationEdge->set_count(triangulationEdge->triangleCount() - 1);
    }

    //Triangulated baseEdge is not contained in the graph
    if (triangulationEdge == nullptr || !triangulationEdge->is_valid()){
        if (triangulationEdge == nullptr){
            component->AddEdge(triangulationPair.first(), triangulationPair.second());
            currentComponent.algorithmEdges[triangulationPair] = AlgorithmEdge(triangulationPair, false);
            triangulationEdge = &currentComponent.algorithmEdges[triangulationPair];
            triangulationForest.AddNode(TriangulationForest::NodeType::EDGE, triangulationEdge);
        }
        triangulationEdge->set_valid(true);
        triangulationEdge->setTreeEdge(false);
        triangulationEdge->set_count(1);
        removableEdges.emplace_back(triangulationPair);
    }
    //Triangulation edge is in the graph
    else{
        triangulationEdge->set_count(triangulationEdge->triangleCount() + 1);
        UpdateDegrees(triangulationPair.first());
        UpdateDegrees(triangulationPair.second());
    }

    if (this->print) {
        printing.printTriangulationStep(Triangle(edgeA, edgeB, *triangulationEdge));
    }
}

void OuterPlanarSubgraphMitchell::UpdateTriangulationForest(const NodePair &triangulationEdge, const NodePair &edgeA,
                                                            const NodePair &edgeB, NodeId v) {

    //Add collapsed edges to collapse tree and connect to base edge
    triangulationForest.AddEdge(triangulationEdge, v);
    triangulationForest.AddEdge(v, edgeA);
    triangulationForest.AddEdge(v, edgeB);

}

void OuterPlanarSubgraphMitchell::DeleteEdges(std::vector<TriangulationForest::TFNode*> &path, const PrintSteps& printSteps, bool includingFirstEdge) {
    for (auto const node : path) {
        if (!includingFirstEdge){
            includingFirstEdge = true;
        }
        else {
            if (node->nodeType == TriangulationForest::NodeType::EDGE) {
                node->correspondingEdge->delete_edge();
                auto const it = currentComponent.edges.find(node->correspondingEdge->edge());
                if (it != currentComponent.edges.end()){
                    currentComponent.edges[node->correspondingEdge->edge()] = false;
                }
            }
        }
    }
    if (this->print) {
        printSteps.printDeleteStep(path);
    }
}


void TriangulationForest::AddNode(TriangulationForest::NodeType nodeType, AlgorithmEdge* edge, NodeId nodeId) {
    _nodes.emplace_back(TFNode(nodeType, edge, nodeId));
    if (nodeType == NodeType::EDGE){
        edgeMap[edge->edge()] = (int) _nodes.size() - 1;
    }
    else{
        nodeMap[nodeId] = (int) _nodes.size() - 1;
    }
}

void TriangulationForest::AddEdge(NodeId nodeId, const NodePair& nodePair) {
    TFNode* srcNode = GetNode(nodeId);
    srcNode->children.emplace_back(edgeMap[nodePair]);
}

void TriangulationForest::AddEdge(const NodePair &nodePair, NodeId nodeId) {
    TFNode* srcNode = GetNode(nodePair);
    srcNode->children.emplace_back(nodeMap[nodeId]);
}

TriangulationForest::TFNode* TriangulationForest::rand_neighbor(TFNode* node, std::mt19937_64 &gen) {
    TFNode* neighbor = nullptr;
    if (!node->children.empty()){
        int rndIdx = std::uniform_int_distribution<int>(0, (int) node->children.size() - 1)(gen);
        neighbor = _get_node(node->children[rndIdx]);
        node->children.erase(node->children.begin() + rndIdx);
    }
    return neighbor;
}

void TriangulationForest::GetValidPath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen) {
    path.clear();
    TFNode* pathRoot = GetNode(nodePair);
    path.emplace_back(GetNode(nodePair));
    bool pathUp = false;
    while (!path.empty()){
        TFNode* neighbor = rand_neighbor(path.back(), gen);
        if (neighbor == nullptr){
            if (!pathUp) {
                break;
            }
            else{
                path.pop_back();
            }
        }
        else{
            if (neighbor->correspondingEdge == nullptr || !neighbor->correspondingEdge->isTreeEdge()){
                path.emplace_back(neighbor);
                pathUp = false;
            }
            else{
                pathUp = true;
            }
        }
    }
    if (pathRoot->children.empty()){
        pathRoot->correspondingEdge->setTreeEdge();
    }
}

bool TriangulationForest::GetValidDeletePath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen) {
    GetValidPath(nodePair, path, gen);
    return (int) path.size() > 0 && path.back()->correspondingEdge->original_edge();
}

void TriangulationForest::print() const{
    for (auto const & node : _nodes) {
        for (auto const & children : node.children) {
            if (node.nodeType == NodeType::NODE){
                std::cout << node.nodeId;
            }
            else{
                std::cout << "(";
                node.correspondingEdge->edge().print();
                std::cout << ")";
            }
            std::cout << " - ";
            const TFNode& childNode = _nodes[children];
            if (childNode.nodeType == NodeType::NODE){
                std::cout << childNode.nodeId << " ";
            }
            else{
                std::cout << "(";
                childNode.correspondingEdge->edge().print();
                std::cout << ") ";
            }
        }
    }
    std::cout << std::endl;
}

void TriangulationForest::AddPath(std::vector<TFNode *> path) {
    for (int i = 0; i < path.size() - 1; ++i) {
        if (path[i+1]->nodeType == NodeType::NODE){
            path[i]->children.emplace_back(nodeMap[path[i+1]->nodeId]);
        }
        else{
            path[i]->children.emplace_back(edgeMap[path[i+1]->correspondingEdge->edge()]);
        }
    }
}


void OuterPlanarSubgraphMitchell::PrintSteps::printTriangulationStep(const Triangle& triangle) const{
    std::cout << "\t\t\t\t" << "Triangulate ";
    triangle.print();
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printDeleteStep(std::vector<TriangulationForest::TFNode*>& path) const {
    std::cout << "\t\t\t\t" << "Delete Edges ";
    for (auto const & node : path) {
        if (node->nodeType == TriangulationForest::NodeType::EDGE){
            std::cout << " ";
            node->correspondingEdge->edge().print();
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printRemovedNode(NodeId i) const {
    std::cout << "\t\t\t" << "Removed Node " << i;
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printTriangulationForest() const {
    std::cout << "\t\t\t\t" << "New Triangulation Forest ";
    _triangulationForest.print();
}

void OuterPlanarSubgraphMitchell::PrintSteps::printEdge(const NodePair & pair) const {
    std::cout << "\t\t\t\t" << "Delete Base Edge ";
    pair.print();
}

int OuterPlanarSubgraphMitchell::nodeIdToComponentId(NodeId nodeId) const {
    return nodeToComponentAndId[nodeId].front();
}

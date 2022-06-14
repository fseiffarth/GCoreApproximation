//
// Created by anonymous on 07.11.21.
//

#include "OuterplanarGraphData.h"
#include "../Utils/StaticFunctions.h"

void OuterplanarGraphData::set() {
    if (nodeToComponents.empty()){
        this->init_outerplanar();
    }
    bbTree.tree.graphType = GraphType::TREE;
    TCnComV components;
    TSnap::GetBiCon(get_graph(), components);
    //StaticFunctions::printComponents(components);
    for (int i = 0; i < components.Len(); ++i) {
        TCnCom& currentComponent = components[i];
        if (currentComponent.Len() > 2) {
            Components.emplace_back(OuterplanarComponent());
            Components.back().component = GraphData(TSnap::GetSubGraph(get_graph(), components[i].NIdV, true));
            Components.back().component.graphType = GraphType::OUTERPLANAR;
            Components.back().getInnerFaces();
        }
        for (int j = 0; j < currentComponent.Len(); ++j) {
            int currentNodeId = currentComponent[j];
            nodeToComponents[currentNodeId].emplace_back(i);
            maxCompSize[currentNodeId] = std::max(maxCompSize[currentNodeId], currentComponent.Len());
            if (currentComponent.Len() > 2) {
                Components.back().NodeIdToComponentNodeId[currentNodeId] = j;
                Components.back().ComponentNodeIdToNodeId[j] = currentNodeId;
            }
        }
    }
    for (int i = 0; i < nodeToComponents.size(); ++i) {
        auto & comps = nodeToComponents[i];
        if (comps.size() > 1 || max_component_size(i) == 2){
            int currentNode = bbTree.tree.graph()->AddNode();
            bbTree.nodeToBBNodes[i] = {currentNode};
            bbTree.BBNodesToNodeOrComponent[currentNode] = {i, -1};
        }
        else{
            bbTree.nodeToBBNodes[i] = {};
        }
    }
    int c = 0;
    for (int i = 0; i < components.Len(); ++i) {
        TCnCom& currentComponent = components[i];
        if (currentComponent.Len() > 2) {
            int component_node = bbTree.tree.graph()->AddNode();
            for (int j = 0; j < currentComponent.Len(); ++j) {
                bbTree.nodeToBBNodes[currentComponent[j]].emplace_back(component_node);
                bbTree.BBNodesToNodeOrComponent[component_node] = {-1, c};
                int currentNodeId = currentComponent[j];
                if (nodeToComponents[currentNodeId].size() > 1){
                    bbTree.tree.graph()->AddEdge(bbTree.nodeToBBNodes[currentNodeId][0], component_node);
                }
            }
        }
        else{
            --c;
            bbTree.tree.graph()->AddEdge(bbTree.nodeToBBNodes[components[i][0]][0], bbTree.nodeToBBNodes[components[i][1]][0]);
        }
        ++c;
    }
    //bbTree.tree.print();
    bbTree.tree.init();
}

void OuterplanarGraphData::get_components(NodeId nodeId, std::vector<NodeOrComponent> &nodeOrComponent) {
    nodeOrComponent.clear();
    for(auto const & bbnode : bbTree.nodeToBBNodes[nodeId]){
        nodeOrComponent.emplace_back(bbTree.BBNodesToNodeOrComponent[bbnode]);
    }
}

void OuterplanarGraphData::get_bb_tree_ids(NodeId nodeId, std::set<NodeId>& ids) {
    ids.insert(bbTree.nodeToBBNodes[nodeId].begin(), bbTree.nodeToBBNodes[nodeId].end());
}

void OuterplanarGraphData::init_outerplanar() {
    nodeToComponents = std::vector<std::vector<int>>(nodes());
    maxCompSize = std::vector<int>(nodes(), 0);
}

void OuterplanarComponent::getInnerFaces() {

}

BBTree::~BBTree() {
}

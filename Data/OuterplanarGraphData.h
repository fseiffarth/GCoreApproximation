//
// Created by anonymous on 07.11.21.
//

#ifndef CLOSURES_OUTERPLANARGRAPHDATA_H
#define CLOSURES_OUTERPLANARGRAPHDATA_H


#include "GraphData.h"

struct OuterplanarComponent{
    std::vector<PUNGraph> faces;
    std::vector<std::vector<int>> nodeToFaces;
    GraphData component = GraphData(new TUNGraph);
    std::map<NodeId, NodeId> NodeIdToComponentNodeId;
    std::map<NodeId, NodeId> ComponentNodeIdToNodeId;
    void getInnerFaces();
    NodeId componentId(NodeId nodeId){return NodeIdToComponentNodeId[nodeId];};
    NodeId nodeId(NodeId componentNodeId){return ComponentNodeIdToNodeId[componentNodeId];};
};



struct NodeOrComponent{
    int _nodeId = -1;
    int _componentId = -1;

    bool is_node(int& id) const{id = _nodeId; return _nodeId != -1;};
    bool is_component(int& id) const{id = _componentId; return _componentId != -1;};
};

struct BBTree{
    ~BBTree();
    GraphData tree = GraphData(new TUNGraph());
    std::unordered_map<NodeId, std::vector<int>> nodeToBBNodes;
    std::unordered_map<NodeId, NodeOrComponent> BBNodesToNodeOrComponent;
};

class OuterplanarGraphData : public GraphData {
public:
    OuterplanarGraphData() : GraphData(){};
    explicit OuterplanarGraphData(const std::string & path) : GraphData(path){
        nodeToComponents = std::vector<std::vector<int>>(this->nodes());
        maxCompSize = std::vector<int>(this->nodes(), 0);
        set();
    };
    OuterplanarGraphData(PUNGraph pGraph, int size) : GraphData(pGraph, size){
        nodeToComponents = std::vector<std::vector<int>>(pGraph->GetNodes());
        maxCompSize = std::vector<int>(pGraph->GetNodes(), 0);
        set();
    };
    explicit OuterplanarGraphData(PUNGraph pGraph) : GraphData(pGraph){
        nodeToComponents = std::vector<std::vector<int>>(pGraph->GetNodes());
        maxCompSize = std::vector<int>(pGraph->GetNodes(), 0);
        set();
    };
    explicit OuterplanarGraphData(GraphData& graphData) : GraphData(graphData){
        nodeToComponents = std::vector<std::vector<int>>(graphData.get_graph()->GetNodes());
        maxCompSize = std::vector<int>(graphData.get_graph()->GetNodes(), 0);
        set();
    };

    std::vector<OuterplanarComponent> Components;
    void set();
    GraphData& get_bbTree() {return bbTree.tree;};
    NodeOrComponent& get_bbNodeOrComponent(NodeId nodeId){return bbTree.BBNodesToNodeOrComponent[nodeId];};
    void get_components(NodeId nodeId, std::vector<NodeOrComponent>& nodeOrComponent);
    int get_component_num(NodeId nodeId){return (int) bbTree.nodeToBBNodes[nodeId].size();};
    int max_component_size(NodeId nodeId){return maxCompSize[nodeId];};
    void get_bb_tree_ids(NodeId nodeId, std::set<NodeId>& ids);

private:
    std::vector<std::vector<int>> nodeToComponents;
    std::vector<int> maxCompSize;
    BBTree bbTree;

    void init_outerplanar();
};


#endif //CLOSURES_OUTERPLANARGRAPHDATA_H

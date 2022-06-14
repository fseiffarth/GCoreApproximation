//
// Created by anonymous on 29.06.21.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPHMITCHELL_H
#define CLOSURES_OUTERPLANARSUBGRAPHMITCHELL_H


#include "GraphStructs.h"
#include "OuterplanarSubgraph.h"
#include <list>

struct AlgorithmEdge{
    bool original_edge() const;

public:
    AlgorithmEdge() = default;
    explicit AlgorithmEdge(const NodePair& nodePair, bool type = true);
    void setTreeEdge(bool b = true){ treeEdge = b;}
    [[nodiscard]] bool isTreeEdge() const{return treeEdge;};
    [[nodiscard]] bool is_valid() const{return _valid;};
    void set_valid(bool b){ _valid = b;};
    [[nodiscard]] int triangleCount() const{return _triangleCount;};
    [[nodiscard]] const NodePair& edge() const{return _nodePair;};
    void delete_edge(){_deleted = true; _valid = false; _triangleCount = 0;};
    [[nodiscard]] bool is_deleted() const {return _deleted;};
    void set_count(int count){_triangleCount = count;};

private:
    NodePair _nodePair{};
    int _triangleCount = 0;
    bool _type = true;
    bool _valid = true;
    bool treeEdge = false;
    bool _deleted = false;
};

struct BiconnectedComponent{
    BiconnectedComponent()= default;
    //Nodes and degrees
    std::vector<NodeId> nodeIds;
    std::vector<NodeId> degree2Nodes;
    std::vector<int> nodeDegrees;

    //edges
    std::unordered_map<NodePair, AlgorithmEdge, hashNodePair> algorithmEdges;
    std::unordered_map<NodePair, bool, hashNodePair> edges;
    int componentId = -1;
};

class TriangulationForest{
public:
    enum class NodeType{
        EDGE,
        NODE,
    };
    struct TFNode{
        explicit TFNode(NodeType nodeType, AlgorithmEdge* correspondingEdge = nullptr, NodeId nodeId = -1)
                : nodeType(nodeType), correspondingEdge(correspondingEdge), nodeId(nodeId) {}

    public:
        std::vector<NodeId> children;
        NodeType nodeType;
        AlgorithmEdge* correspondingEdge;
        NodeId nodeId;
    };

    void AddNode(NodeType nodeType, AlgorithmEdge* edge = nullptr, NodeId nodeId = -1);
    void AddEdge(const NodePair& nodePair, NodeId nodeId);
    void AddEdge(NodeId nodeId, const NodePair& nodePair);
    int degree(const NodePair& nodePair){return (int) Node(nodePair)->children.size();};
    int degree(const AlgorithmEdge& algorithmEdge){ return degree(algorithmEdge.edge());};
    void clear(){_nodes.clear(); edgeMap.clear(); nodeMap.clear();}
    TFNode* rand_neighbor(TFNode* node, std::mt19937_64& gen);
    void GetValidPath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen);
    bool GetValidDeletePath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen);

    const TFNode* Node(NodeId nodeId) {return _node(nodeMap[nodeId]);};
    TFNode* GetNode(NodeId nodeId) {return _get_node(nodeMap[nodeId]);};
    const TFNode* Node(const NodePair& nodePair) {return _node(edgeMap[nodePair]);};
    TFNode* GetNode(const NodePair& nodePair) {return _get_node(edgeMap[nodePair]);};
    void print() const;

    void AddPath(std::vector<TFNode *> path);

private:
    std::vector<TFNode> _nodes;
    std::unordered_map<NodePair, int, hashNodePair> edgeMap;
    std::unordered_map<NodeId, int> nodeMap;
    TFNode* _get_node(NodeId nodeId){return &_nodes[nodeId];};
    const TFNode* _node(NodeId nodeId){return &_nodes[nodeId];};
};

class Triangle{
public:
    Triangle(const AlgorithmEdge& edgeA,const AlgorithmEdge& edgeB,const AlgorithmEdge& edgeBase) : a(edgeA), b(edgeB), base(edgeBase){};
    const AlgorithmEdge& a;
    const AlgorithmEdge& b;
    const AlgorithmEdge& base;
    void print() const{
        NodeId head = a.edge().first();
        if (a.edge().first() == base.edge().first() || a.edge().first() == base.edge().second()){
            head = a.edge().second();
        }
        std::cout << base.edge().first() <<  "/" << head << "\\" << base.edge().second() << " weights " << a.triangleCount() << "/" << base.triangleCount() << "\\" << b.triangleCount();
    }
};


class OuterPlanarSubgraphMitchell : public OuterplanarSubgraph {
public:
    explicit OuterPlanarSubgraphMitchell(const PUNGraph graph);
    void generate(GraphData& subgraph, std::mt19937_64& gen, bool p) override;

    class PrintSteps {
    public:
        explicit PrintSteps(const TriangulationForest& triangulationForest) : _triangulationForest(triangulationForest){};
        NodePair currentDeleteEdge{};
        void printTriangulationStep(const Triangle& triangle) const;
        void printDeleteStep(std::vector<TriangulationForest::TFNode*>& path) const;
        void printRemovedNode(NodeId i) const;
        void printTriangulationForest() const;
        const TriangulationForest& _triangulationForest;
        void printEdge(const NodePair & pair) const;
    };
private:
    void getValidNeighbor(PUNGraph graph, std::unordered_map<NodePair, AlgorithmEdge, hashNodePair> &algEdges,
                     int nodeId,
                     int neighborIdx, NodeId &neighborId);
    int nodeIdToComponentId(NodeId nodeId) const;
    void GetBiconnectedComponents();
    void BiconnectedComponentSampling(const PUNGraph subgraph, BiconnectedComponent& biCom);
    void Init();
    void UpdateDegrees(NodeId id);
    void DeleteRandomEdge(const PrintSteps& printing);
    void CollapseTriangle(PUNGraph component, AlgorithmEdge& edgeA, AlgorithmEdge& edgeB, NodePair& triangulationPair, const PrintSteps& printing);
    void UpdateTriangulationForest(const NodePair &triangulationEdge, const NodePair &edgeA, const NodePair &edgeB, NodeId v);
    void DeleteEdges(std::vector<TriangulationForest::TFNode*>& path, const PrintSteps& printSteps, bool includingFirstEdge = false);


    std::vector<PUNGraph> biconnectedComponents;

    std::vector<std::list<int>> nodeToComponentAndId;
    std::vector<BiconnectedComponent> biComponents;
    BiconnectedComponent currentComponent;
    std::vector<NodePair> treeEdges;

    std::vector<NodePair> removableEdges;
    TriangulationForest triangulationForest;
    std::vector<TriangulationForest::TFNode*> _path1;
    std::vector<TriangulationForest::TFNode*> _path2;


};


#endif //CLOSURES_OUTERPLANARSUBGRAPHMITCHELL_H

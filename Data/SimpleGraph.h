//
// Created by anonymous on 20.09.21.
//

#ifndef CLOSURES_SIMPLEGRAPH_H
#define CLOSURES_SIMPLEGRAPH_H


#include "../Utils/typedefs.h"
#include "Data.h"
#include "GraphData.h"

class SimpleGraph : public Data<SIMPLE_GRAPH>{
public:
    explicit SimpleGraph(const std::string& graphPath);
    explicit SimpleGraph(int size);
    SimpleGraph(const std::string& graphPath, const std::string& labelPath);
    explicit SimpleGraph(const PUNGraph& Graph, const std::string& name = "");
    SimpleGraph(const PUNGraph& Graph, const std::string& name, Labels& labels);
    SimpleGraph(const PUNGraph &Graph, int size, const std::string& name = "");

    [[nodiscard]] const SIMPLE_GRAPH & get_graph() const;
    SIMPLE_GRAPH& graph();
    void clear_edges();
    [[nodiscard]] Nodes& node(NodeId Id);
    [[nodiscard]] const Nodes& get_node(NodeId Id) const;
    NodeId elem(NodeId Id) override;
    [[nodiscard]] int degree(NodeId Id) const;
    [[nodiscard]] NodeId neighbor(int NodeId, int NeighborIdx) const;
    NodeId random_neighbor(int NodeId, int minIdx, std::mt19937_64 &gen);
    [[nodiscard]] size_t size() const override;
    void AddEdge(NodeId src, NodeId dest, bool check = false);
    bool FindEdge(NodeId src, NodeId dest, bool directed = false);
    void save(const std::string& path) const;
    void save_edges(const std::string& path) const;
    [[nodiscard]] SimpleGraph mergeGraphs(const SimpleGraph& other) const;
    void mergeGraphs(const SimpleGraph& other);
    void operator+=(const SimpleGraph& other){ mergeGraphs(other);}
    SimpleGraph operator+(const SimpleGraph& other) const {return mergeGraphs(other);}

    struct NodeIterator{
        // Prefix increment
        void operator++() { ++_idx;};
        const Nodes& operator*() const { return _graph->get_node(_idx); }
        friend bool operator== (const NodeIterator& a, int b) { return a._idx == b; };
        friend bool operator!= (const NodeIterator& a, int b) { return a._idx != b; };

        const SimpleGraph* _graph{};
        int _idx{};
    };
    [[nodiscard]] NodeIterator BegNI() const { return NodeIterator{this, 0}; }
    [[nodiscard]] int EndNI() const   { return static_cast<int>(this->get_graph().size()); }
    
    struct NeighborIterator{
        // Prefix increment
        void operator++() { ++_idx;};
        NodeId operator*() const {
            if (_gen == nullptr) {
                return _graph->neighbor(_nodeId, _idx);
            }
            else{
                return _graph->random_neighbor(_nodeId, _idx, *_gen);
            }
        }
        friend bool operator== (const NeighborIterator& a, int b) { return a._idx == b; };
        friend bool operator!= (const NeighborIterator& a, int b) { return a._idx != b; };
    
        SimpleGraph* _graph{};
        NodeId _nodeId{};
        int _idx{};
        std::mt19937_64* _gen = nullptr;
    };
    NeighborIterator BegNeighI(NodeId nodeId, std::mt19937_64* gen = nullptr) { return NeighborIterator{this, nodeId, 0, gen}; }
    [[nodiscard]] int EndNeighI(NodeId nodeId) const   { return this->degree(nodeId); }
    
    
    std::vector<int> DistanceList;
    std::vector<bool> ContainmentList;
    std::vector<std::vector<int>> PredecessorsList;
    GraphType graphType = GraphType::GENERAL;
    std::vector<int> neighborList;
    int maxDegree;

    void AddNode();


};


#endif //CLOSURES_SIMPLEGRAPH_H

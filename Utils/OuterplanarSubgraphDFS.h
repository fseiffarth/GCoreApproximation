//
// Created by anonymous on 10.06.2021.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPHDFS_H
#define CLOSURES_OUTERPLANARSUBGRAPHDFS_H

#include "GraphStructs.h"
#include "OuterplanarSubgraph.h"
#include "../Data/GraphData.h"
#include "../Data/OuterplanarGraphData.h"


enum class Side {
    LEFT,
    RIGHT,
};

class NodeReachability{
public:
    NodeReachability(): _left(true), _right(true){};
    [[nodiscard]] bool left() const {return _left;};
    [[nodiscard]] bool right() const {return _right;};
    void set(bool left, bool right){ _left = left; _right = right;};
    void close_left(){ _left = false;};
    void open_left(){_left = true;};
    void close_right(){ _right = false;};
    void open_right(){_right = true;};
    void reset(){ _right = true, _left = true;};
    [[nodiscard]] int count() const{
        int count = 0;
        if (_left){++count;};
        if(_right){++count;};
        return count;
    }
private:
    bool _left = true;
    bool _right = true;

};

class NodeStruct{
public:
    NodeStruct() = default;
    explicit NodeStruct(int id) : node_id(id){};
    NodeId node_id = -1;
    bool visited = false;
    int dfs_depth = -1;
    NodeReachability nodeReachability = NodeReachability();
    NodeStruct* dfs_parent = nullptr;
    int last_R = 0;
    int last_L = 0;

    void reset(){
        dfs_parent = nullptr;
        visited = false;
        dfs_depth = -1;
        nodeReachability.reset();
        last_L = 0;
        last_R = 0;
        left_spanning = false;
        right_spanning = false;
    };

    bool operator <(const NodeStruct &b) const{
        return this->dfs_depth < b.dfs_depth;
    }
    bool operator >(const NodeStruct &b) const{
        return this->dfs_depth > b.dfs_depth;
    }
    bool operator ==(const NodeStruct &b) const{
        return this->dfs_depth == b.dfs_depth;
    }
    bool operator >(int i) const{
        return this->dfs_depth > i;
    }
    bool operator <(int i) const{
        return this->dfs_depth < i;
    }
    bool operator >=(int i) const{
        return this->dfs_depth >= i;
    }
    bool operator <=(int i) const{
        return this->dfs_depth <= i;
    }
    bool operator ==(int i) const{
        return this->dfs_depth == i;
    }

    NodeStruct *last_spanned_root = nullptr;
    bool left_spanning = false;
    bool right_spanning = false;
};

class DirectedEdgeStruct{
public:
    DirectedEdgeStruct(NodeId src, NodeId dst) : src(src), dst(dst){};
    DirectedEdgeStruct(const NodeStruct& src, const NodeStruct& dst) : src(src.node_id), dst(dst.node_id){};
    NodeId src;
    NodeId dst;
};


class OuterplanarSubgraphDFS : public OuterplanarSubgraph {
public:
    explicit OuterplanarSubgraphDFS(const PUNGraph graph);

    //Copy constructor
    OuterplanarSubgraphDFS(const OuterplanarSubgraphDFS& other);

    void generate(GraphData& subgraph, std::mt19937_64& gen, bool p) override;
    void generate(OuterplanarGraphData& subgraph, std::mt19937_64& gen, bool p);
    void get_next_node();
    bool CheckLeft(const NodeStruct& src, const NodeStruct& dst);
    bool CheckRight(const NodeStruct& src, const NodeStruct& dst);
    [[nodiscard]] bool crossesLeftDiagonal(const NodeStruct& dst) const;
    [[nodiscard]] bool crossesRightDiagonal(const NodeStruct& dst) const;
    [[nodiscard]] bool enclosingPointByLeftEdge(const NodeStruct &src,const NodeStruct &dst) const;
    [[nodiscard]] bool enclosingPointByRightEdge(const NodeStruct &src,const NodeStruct &dst) const;

    void GetCotreeEdges(NodeStruct& currentNode);

    void GetValidEdges(NodeStruct& currentNode, std::vector<DirectedEdgeStruct>& edges, Side side);
    void AddEdges(const std::vector<DirectedEdgeStruct>& edges, Side side);
    void AddGraphEdges(const std::vector<DirectedEdgeStruct>& edges, Side side);

    void reset();

    std::vector<NodeId> reachableNodes;
private:
    std::vector<NodeStruct> graphNodes;
    NodeId dfs_root_node{};
    std::vector<DirectedEdgeStruct> possibleEdges;
    std::vector<NodeId> neighborsToVisit;
    std::vector<NodeId> dfs_stack;
    std::vector<DirectedEdgeStruct> left_edges;
    std::vector<DirectedEdgeStruct> right_edges;
    //Used to take neighbor uniformly at random
    std::vector<int> neighborIds;
    std::vector<int> swaps;
    bool print = false;
    bool dfs_tree = false;
    PUNGraph _subgraph;

    void UpdateNodeParameters(NodeStruct &aStruct);
};


#endif //CLOSURES_OUTERPLANARSUBGRAPHDFS_H
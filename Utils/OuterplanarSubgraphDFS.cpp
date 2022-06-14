//
// Created by anonymous on 10.06.2021.
//

#include <iostream>
#include "OuterplanarSubgraphDFS.h"

OuterplanarSubgraphDFS::OuterplanarSubgraphDFS(const PUNGraph graph) : OuterplanarSubgraph(
        graph) {
    int maxDegree = 0;
    graphNodes = std::vector<NodeStruct>(graph->GetNodes(), NodeStruct());
    for (auto node = graph->BegNI(); node != graph->EndNI(); node++){
        graphNodes[node.GetId()] = NodeStruct(node.GetId());
        maxDegree = std::max(maxDegree, node.GetDeg());
    }
    neighborIds = std::vector<int>(maxDegree);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
}

void OuterplanarSubgraphDFS::generate(GraphData& subgraph, std::mt19937_64& gen, bool p) {
    subgraph.graphType = GraphType::OUTERPLANAR;
    this->_subgraph = subgraph.graph();
    this->print = p;
    this->reset();
    dfs_root_node = std::uniform_int_distribution<int>(0, this->_graph->GetNodes() - 1)(gen);
    dfs_stack.emplace_back(dfs_root_node);
    NodeStruct* path_root = &this->graphNodes[dfs_root_node];
    path_root->last_spanned_root = path_root;
    path_root->dfs_parent = path_root;
    path_root->dfs_depth = 0;

    if (p) {
        std::cout << std::endl;
    }
    while (!this->dfs_stack.empty()) {
        get_next_node();
    }
    if (p) {
        std::cout << std::endl;
    }
}

void OuterplanarSubgraphDFS::get_next_node() {
    //Current _node in the dfs
    NodeId currentNodeId = this->dfs_stack.back();
    dfs_stack.pop_back();
    NodeStruct &currentNode = this->graphNodes[currentNodeId];
    //If current get_node is not visited go further down in the dfs and update parameters
    if (!currentNode.visited) {
        //Check if get_node is not the root node
        if (currentNodeId != dfs_root_node) {
            //Add tree edge to the graph
            _subgraph->AddEdge(currentNodeId, currentNode.dfs_parent->node_id);
            UpdateNodeParameters(currentNode);
        }

        //Set get_node to visited
        currentNode.visited = true;
        reachableNodes.emplace_back(currentNodeId);

        //Traverse all the neighbors of currentNode and get cotree edges
        GetCotreeEdges(currentNode);
        if (!neighborsToVisit.empty()) {
            GetValidEdges(currentNode, left_edges, Side::LEFT);
            GetValidEdges(currentNode, right_edges, Side::RIGHT);
            if (left_edges.size() >= right_edges.size()) {
                AddEdges(left_edges, Side::LEFT);
                AddGraphEdges(left_edges, Side::LEFT);
            } else {
                AddEdges(right_edges, Side::RIGHT);
                AddGraphEdges(right_edges, Side::RIGHT);
            }
        }
    }
}


void OuterplanarSubgraphDFS::GetCotreeEdges(NodeStruct& currentNode) {
    //Traverse all neighbors of get_node
    auto node = this->_graph->GetNI(currentNode.node_id);
    int degree = node.GetDeg();
    this->neighborsToVisit.clear();
    int randIdx;
    NodeId neighbor;
    NodeStruct* neighborNode;
    for (int i = 0; i < degree; ++i) {
        randIdx = std::uniform_int_distribution<int>(i, degree - 1)(_gen);
        neighbor = node.GetNbrNId(neighborIds[randIdx]);
        swaps.emplace_back(randIdx);
        std::swap(neighborIds[i], neighborIds[randIdx]);
        if (!currentNode.dfs_parent || neighbor != currentNode.dfs_parent->node_id) {
            if (this->graphNodes[neighbor].visited) {
                //Edge currentNodeId, visitedNeighbor can be some diagonal
                this->neighborsToVisit.emplace_back(neighbor);
            } else {
                //Add node to dfs_stack and update parent get_node
                neighborNode = &this->graphNodes[neighbor];
                neighborNode->dfs_parent = &currentNode;
                this->dfs_stack.emplace_back(neighbor);
            }
        }
    }
    for (int i = (int) swaps.size() - 1; i >= 0; --i) {
        std::swap(neighborIds[i], neighborIds[swaps[i]]);
        swaps.pop_back();
    }
}

bool OuterplanarSubgraphDFS::CheckLeft(const NodeStruct &src,const NodeStruct &dst) {
    return !enclosingPointByLeftEdge(src, dst) && !crossesLeftDiagonal(dst);// || (dst.dfs_depth == src.last_spanned_root->dfs_depth - 1 && !current_left));
}

bool OuterplanarSubgraphDFS::CheckRight(const NodeStruct &src, const NodeStruct &dst) {
    return !enclosingPointByRightEdge(src, dst) && !crossesRightDiagonal(dst);// || (dst.dfs_depth == src.last_spanned_root->dfs_depth - 1 && !current_right));
}

bool OuterplanarSubgraphDFS::crossesLeftDiagonal(const NodeStruct& dst) const {
    return !dst.nodeReachability.left();
}
bool OuterplanarSubgraphDFS::crossesRightDiagonal(const NodeStruct& dst) const {
    return !dst.nodeReachability.right();
}

bool OuterplanarSubgraphDFS::enclosingPointByLeftEdge(const NodeStruct &src, const NodeStruct &dst) const {
    return dst < src.last_L;
}

bool OuterplanarSubgraphDFS::enclosingPointByRightEdge(const NodeStruct &src, const NodeStruct &dst) const {
    return dst < src.last_R;
}

void OuterplanarSubgraphDFS::GetValidEdges(NodeStruct& currentNode, std::vector<DirectedEdgeStruct>& edges, Side side) {
    edges.clear();
    for (int neighborId : this->neighborsToVisit) {
        NodeStruct &neighbor = this->graphNodes[neighborId];
        if (side == Side::LEFT) {
            if (CheckLeft(currentNode, neighbor)) {
                edges.emplace_back(DirectedEdgeStruct(currentNode, neighbor));
            }
        } else {
            if (CheckRight(currentNode, neighbor)) {
                edges.emplace_back(DirectedEdgeStruct(currentNode, neighbor));
            }
        }
    }
}

void OuterplanarSubgraphDFS::AddEdges(const std::vector<DirectedEdgeStruct>& edges, Side side) {
    if (!edges.empty()) {
//        if (side == Side::LEFT){
//            left = true;
//        }
//        else{
//            right = true;
//        }

        for (auto &edge: edges) {
            NodeStruct& src = this->graphNodes[edge.src];
            NodeStruct& dst = this->graphNodes[edge.dst];
            //Update low
//            if (dst < *src.last_spanned_root){
//                left = true;
//                right = true;
//            }
            int i = (int) this->reachableNodes.size() - 2;
            NodeStruct *currentNode = &this->graphNodes[this->reachableNodes.back()];
            if (side == Side::LEFT) {
                currentNode->last_R = currentNode->dfs_parent->dfs_depth;
            } else {
                currentNode->last_L = currentNode->dfs_parent->dfs_depth;
            }
            while (!this->reachableNodes.empty() && this->graphNodes[this->reachableNodes[i]] > dst) {
                currentNode = &this->graphNodes[this->reachableNodes[i]];
                //Mark as side
                if (side == Side::LEFT) {
                    currentNode->nodeReachability.close_left();
                    currentNode->last_R = currentNode->dfs_depth;
                    currentNode->left_spanning = true;
                }
                else{
                    currentNode->nodeReachability.close_right();
                    currentNode->last_L = currentNode->dfs_depth;
                    currentNode->right_spanning = true;
                }
                this->reachableNodes.erase(this->reachableNodes.end() - 2);
                --i;
            }
        }
    }
}




void OuterplanarSubgraphDFS::reset() {
    for (int i = 0; i < this->_graph->GetNodes(); ++i) {
        _subgraph->AddNode();
        this->graphNodes[i].reset();
    }
    dfs_stack.clear();
    reachableNodes.clear();
}

OuterplanarSubgraphDFS::OuterplanarSubgraphDFS(const OuterplanarSubgraphDFS& other) :  OuterplanarSubgraph(other) {
    reachableNodes = other.reachableNodes;
    graphNodes = other.graphNodes;
    dfs_root_node = other.dfs_root_node;
    possibleEdges = other.possibleEdges;
    neighborsToVisit = other.neighborsToVisit;
    dfs_stack = other.dfs_stack;
    left_edges = other.left_edges;
    right_edges = other.right_edges;
    //Used to take neighbor uniformly at random
    neighborIds = other.neighborIds;
    swaps = other.swaps;
    print = other.print;
    dfs_tree = other.dfs_tree;
}

void OuterplanarSubgraphDFS::UpdateNodeParameters(NodeStruct &currentNode) {
    if (print) {
        std::cout << " " << currentNode.dfs_parent->node_id << "-" << currentNode.node_id << " ";
    }
    //get parent node
    NodeStruct &currentParent = *currentNode.dfs_parent;

    //check if current node is in a new path
    if (this->graphNodes[this->reachableNodes.back()] > currentParent) {
        //Update the stack
        while (this->graphNodes[this->reachableNodes.back()] > currentParent) {
            this->reachableNodes.pop_back();
        }
        NodeStruct& path_root = currentParent;
        NodeStruct& path_root_parent = *path_root.dfs_parent;
        //if current branching node is not reachable from left (i.e. there is a left edge over the branching node) then there can be no higher left edge than to the branching get_node
        if (path_root.left_spanning) {
            path_root.last_L = path_root.dfs_depth;
        }
        else{
            path_root.last_L = path_root_parent.last_L;
        }
        //if current branching node is not reacheable from right (i.e there is a right edge over the branching node) then there can be no higher right edge than to the branching get_node
        if (path_root.right_spanning) {
            path_root.last_R = path_root.dfs_depth;
        }
        else{
            path_root.last_R = path_root_parent.last_R;
        }
        //Set reachability for branching get_node to true in branch
        path_root.nodeReachability.set(true, true);
        reachableNodes.emplace_back(path_root.node_id);
    }
    currentNode.last_L = currentParent.last_L;
    currentNode.last_R = currentParent.last_R;

    //Set get_node depth in dfs tree and id of the branch in which it was found
    currentNode.dfs_depth = currentParent.dfs_depth + 1;
}

void OuterplanarSubgraphDFS::AddGraphEdges(const std::vector<DirectedEdgeStruct> &edges, Side side) {
    for (auto &edge: edges) {
        NodeStruct src = this->graphNodes[edge.src];
        NodeStruct dst = this->graphNodes[edge.dst];
        if (print){
            if (side == Side::LEFT) {
                std::cout << " " << src.node_id << "--L" << dst.node_id << " ";
            }
            else{
                std::cout << " " << src.node_id << "--R" << dst.node_id << " ";
            }
        }
        this->_subgraph->AddEdge(src.node_id, dst.node_id);
    }
}

void OuterplanarSubgraphDFS::generate(OuterplanarGraphData &subgraph, std::mt19937_64 &gen, bool p) {
    generate(dynamic_cast<GraphData &>(subgraph), gen, p);
}


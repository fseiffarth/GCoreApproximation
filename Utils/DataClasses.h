//
// Created by anonymous on 15.04.2021.
//

#ifndef HOPS_DATACLASSES_H
#define HOPS_DATACLASSES_H
#include <vector>
#include "Enums.h"
#include "typedefs.h"
#include <map>
#include <set>

namespace GRAPHS {
    class GraphClass {
    public:
        explicit GraphClass(int seed = 0){
            this->_gen = std::mt19937_64(seed);
        };

        explicit GraphClass(const std::string &graph_path){
            if (std::filesystem::path(graph_path).extension() == ".edges") {
                try {
                    NodeId src;
                    NodeId dest;
                    std::string a, b;
                    std::string line;
                    std::ifstream infile(graph_path);
                    while (std::getline(infile, line)) {
                        std::istringstream iss(line);
                        iss >> a;
                        if (a == "#") {
                            continue;
                        } else if (iss >> b) {
                            src = std::stoi(a);
                            dest = std::stoi(b);
                            int max = std::max(src, dest);
                            if (max >= _nodes.size()) {
                                _nodes.resize(max, std::vector<NodeId>());
                            }
                        } else {

                        }
                    }
                    _name = std::filesystem::path(graph_path).stem();
                }
                catch (...) {
                }
            }
        };

        [[nodiscard]] int GetNodes() const { return node_num; };

        [[nodiscard]] int GetEdges() const { return edge_num; };

        void NewSeed(int seed){
            this->_gen = std::mt19937_64(seed);
        };

        void DeleteAllEdges() {
            for (auto &node: _nodes) {
                node.clear();
            }
            edge_num = 0;
        };

        int AddNode() {
            _nodes.emplace_back(std::vector<NodeId>());
            ++node_num;
            return node_num - 1;
        };

        int AddEdge(NodeId src, NodeId dst) {
            _nodes[src].emplace_back(dst);
            ++edge_num;
            return static_cast<int>(_nodes[src].size());
        };

        [[nodiscard]] int GetDegree(NodeId src) const { return static_cast<int>(_nodes[src].size()); };

        int GetRandomNeighbor(NodeId src) {
            int randIdx = std::uniform_int_distribution<int>(0, GetDegree(src) - 1)(_gen);
            return _nodes[src][randIdx];
        };

        int GetRandomNeighborInRange(NodeId src, int minIdx) {
            int randIdx = std::uniform_int_distribution<int>(minIdx, GetDegree(src) - 1)(_gen);
            std::swap(_nodes[src][randIdx], _nodes[src][0]);
            return _nodes[src][0];
        };

        [[nodiscard]] int GetNeighbor(NodeId src, int idx) const {
            return _nodes[src][idx];
        }

        void Save(std::string &path) {
            std::string edgePath = path + this->_name + ".edges";
            std::ofstream file;
            file.open(edgePath);
            for (int i = 0; i < _nodes.size(); ++i) {
                for (NodeId node: _nodes[i]) {
                    file << i << " " << node << std::endl;
                }
            }
            file.close();
        };

    struct NeighborIterator{
        NeighborIterator(const GraphClass& graph, NodeId nodeId, int idx) : _graph(graph), _nodeId(nodeId), _idx(idx){};
        // Prefix increment
        void operator++() { ++_idx;};
        NodeId operator*() const { return _graph.GetNeighbor(_nodeId, _idx); }
        friend bool operator== (const NeighborIterator& a, int b) { return a._idx == b;};
        friend bool operator!= (const NeighborIterator& a, int b) { return a._idx != b;};
    private:
        const GraphClass& _graph;
        NodeId _nodeId;
        int _idx;
    };
    NeighborIterator begin(NodeId nodeId) const { return NeighborIterator(*this, nodeId, 0); }
    int end(NodeId nodeId) const   { return this->GetDegree(nodeId); }

    struct NeighborIteratorRand{
        NeighborIteratorRand(GraphClass& graph, NodeId nodeId, int idx) : _graph(graph), _nodeId(nodeId), _idx(idx){};
        // Prefix increment
        void operator++() { ++_idx;};
        NodeId operator*() const { return _graph.GetRandomNeighborInRange(_nodeId, _idx); }
        friend bool operator== (const NeighborIteratorRand& a, const NeighborIteratorRand& b) { return a._idx == b._idx; };
        friend bool operator!= (const NeighborIteratorRand& a, const NeighborIteratorRand& b) { return a._idx != b._idx; };
    private:
        GraphClass& _graph;
        NodeId _nodeId;
        int _idx;
    };
    NeighborIteratorRand beginRand(NodeId nodeId) { return NeighborIteratorRand(*this, nodeId, 0); }
    NeighborIteratorRand endRand(NodeId nodeId)   { return NeighborIteratorRand(*this, nodeId, this->GetDegree(nodeId)); }

    protected:
        int node_num = 0;
        int edge_num = 0;
        std::vector<std::vector<NodeId>> _nodes;
        std::string _name;
        std::mt19937_64 _gen;
    };
}
#endif //HOPS_DATACLASSES_H

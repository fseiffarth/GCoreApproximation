//
// Created by anonymous on 20.09.21.
//

#include "SimpleGraph.h"
#include "../Utils/GraphFunctions.h"

SimpleGraph::SimpleGraph(const std::string &graphPath) {
    if (std::filesystem::path(graphPath).extension() == ".edges" || std::filesystem::path(graphPath).extension() == ".txt") {
        try {
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(graphPath);
            while(std::getline(infile, line)){
                std::istringstream iss(line);
                iss >> a;
                if (a == "#"){
                    continue;
                }
                else if (iss >> b){
                    src = std::stoi(a);
                    dest = std::stoi(b);
                    AddEdge(src, dest);
                }
                else{

                }
            }
            _name = std::filesystem::path(graphPath).stem();
        }
        catch (...) {
        }
    }
    classNumber = 0;
    for (int i = 0; i < _labels.size(); ++i) {
        if (_labelMap.find(_labels[i]) == _labelMap.end()){
            _labelMap[_labels[i]] = Nodes{i};
            ++classNumber;
        }
        else{
            _labelMap[_labels[i]].emplace_back(i);
        }
    }
    int size = static_cast<int>(this->size());
    DistanceList = std::vector<int>(size, 0);
    ContainmentList = std::vector<bool>(size, false);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

SimpleGraph::SimpleGraph(int size) {
    for (int i = 0; i < size; ++i) {
        this->AddNode();
    }
    DistanceList = std::vector<int>(size, 0);
    ContainmentList = std::vector<bool>(size, false);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

SimpleGraph::SimpleGraph(const PUNGraph &Graph, int size, const std::string &name) : Data<SIMPLE_GRAPH>(SIMPLE_GRAPH(), name) {
    for (auto EdgeIt = Graph->BegEI(); EdgeIt != Graph->EndEI(); EdgeIt++) {
        this->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
    }
    DistanceList = std::vector<int>(size, 0);
    ContainmentList = std::vector<bool>(size, false);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

SimpleGraph::SimpleGraph(const PUNGraph &Graph, const std::string &name) : Data<std::vector<std::vector<int>>>(SIMPLE_GRAPH(), name) {
    for (auto EdgeIt = Graph->BegEI(); EdgeIt != Graph->EndEI(); EdgeIt++) {
        this->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
    }
    DistanceList = std::vector<int>(this->size(), 0);
    ContainmentList = std::vector<bool>(this->size(), false);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

SimpleGraph::SimpleGraph(const PUNGraph &Graph, const std::string &name, Labels &labels): Data<std::vector<std::vector<int>>>(SIMPLE_GRAPH(), name, labels) {
    for (auto EdgeIt = Graph->BegEI(); EdgeIt != Graph->EndEI(); EdgeIt++) {
        this->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
    }
    DistanceList = std::vector<int>(this->size(), 0);
    ContainmentList = std::vector<bool>(this->size(), false);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

SimpleGraph::SimpleGraph(const std::string &graphPath, const std::string &labelPath) {
    if (std::filesystem::path(graphPath).extension() == ".edges" || std::filesystem::path(graphPath).extension() == ".txt") {
        try {
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(graphPath);
            while(std::getline(infile, line)){
                std::istringstream iss(line);
                iss >> a;
                if (a == "#"){
                    continue;
                }
                else if (iss >> b){
                    src = std::stoi(a);
                    dest = std::stoi(b);
                    AddEdge(src, dest);
                }
                else{

                }
            }
            _name = std::filesystem::path(graphPath).stem();
        }
        catch (...) {
        }
    }
    if (std::filesystem::path(labelPath).extension() == ".labels") {
        try {
            std::ifstream infile(labelPath);
            NodeId id;
            Label label;
            while (infile >> id >> label){
                _labels.push_back(label);
            }
        }
        catch (...) {
        }
    }
    classNumber = 0;
    for (int i = 0; i < _labels.size(); ++i) {
        if (_labelMap.find(_labels[i]) == _labelMap.end()){
            _labelMap[_labels[i]] = Nodes{i};
            ++classNumber;
        }
        else{
            _labelMap[_labels[i]].emplace_back(i);
        }
    }
    DistanceList = std::vector<int>(this->size(), 0);
    ContainmentList = std::vector<bool>(this->size(), false);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

int SimpleGraph::degree(NodeId Id) const {
    return static_cast<int>(this->get_graph()[Id].size());
}

NodeId SimpleGraph::neighbor(int NodeId, int NeighborIdx) const {
    return this->get_graph()[NodeId][NeighborIdx];
}

NodeId SimpleGraph::random_neighbor(int NodeId, int minIdx, std::mt19937_64& gen) {
    int randIdx = std::uniform_int_distribution<int>(minIdx, this->degree(NodeId) - 1)(gen);
    std::swap(this->graph()[NodeId][randIdx], this->graph()[NodeId][minIdx]);
    return this->get_graph()[NodeId][minIdx];
}

void SimpleGraph::AddEdge(NodeId src, NodeId dest, bool check) {
    if (this->get_graph().size() <= src){
        this->graph().resize(src);
    }
    if (this->get_graph().size() <= dest){
        this->graph().resize(src);
    }

    if(!check) {
        this->graph()[src].emplace_back(dest);
        this->graph()[dest].emplace_back(src);
    }
    else {
        if (this->degree(src) < this->degree(dest)) {
            if (std::find(this->get_graph()[src].begin(), this->get_graph()[src].end(), dest) != this->get_graph()[src].end()){
                this->graph()[src].emplace_back(dest);
                this->graph()[dest].emplace_back(src);
            }
        } else {
            if (std::find(this->get_graph()[dest].begin(), this->get_graph()[dest].end(), src) != this->get_graph()[dest].end()){
                this->graph()[src].emplace_back(dest);
                this->graph()[dest].emplace_back(src);
            }
        }
    }
}

bool SimpleGraph::FindEdge(NodeId src, NodeId dest, bool directed) {
    if (directed){
        return std::find(this->get_graph()[src].begin(), this->get_graph()[src].end(), dest) != this->get_graph()[src].end();
    }
    else {
        if (this->degree(src) < this->degree(dest)) {
            return std::find(this->get_graph()[src].begin(), this->get_graph()[src].end(), dest) !=
                this->get_graph()[src].end();
        } else {
            return std::find(this->get_graph()[dest].begin(), this->get_graph()[dest].end(), src) !=
                this->get_graph()[dest].end();
        }
    }
}

const SIMPLE_GRAPH & SimpleGraph::get_graph() const {
    return this->get_data();
}

SIMPLE_GRAPH &SimpleGraph::graph() {
    return this->data();
}

const Nodes& SimpleGraph::get_node(NodeId Id) const {
    return this->get_graph()[Id];
}

Nodes &SimpleGraph::node(NodeId Id) {
    return this->graph()[Id];
}

NodeId SimpleGraph::elem(NodeId Id) {
    return Id;
}

size_t SimpleGraph::size() const {
    return this->get_graph().size();
}

void SimpleGraph::save(const std::string &path) const {
    std::string labelPath = path + this->_name + ".labels";
    std::ofstream file;
    file.open(labelPath);
    int id = 0;
    for (int label : this->get_labels()) {
        file << id << " " << label << std::endl;
        ++id;
    }
    file.close();
    save_edges(path);
}

void SimpleGraph::save_edges(const std::string &path) const {
    std::string edgePath = path + this->_name + ".edges";
    std::ofstream file;
    file.open(edgePath);
    for (auto NodeIt = this->BegNI(); NodeIt != this->EndNI(); ++NodeIt) {
        auto const Node = *NodeIt;
        NodeId src = NodeIt._idx;
        for (NodeId dest : Node) {
            if (src < dest){
                file << src << " " << dest << std::endl;
            }
        }
    }
    file.close();
}

void SimpleGraph::clear_edges() {
    for (auto Nodes : this->graph()) {
        Nodes.clear();
    }
}

void SimpleGraph::AddNode() {
    this->graph().emplace_back(std::vector<int>());
}










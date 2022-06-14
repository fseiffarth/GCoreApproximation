//
// Created by anonymous on 12.04.2021.
//

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "GraphData.h"
#include "../Utils/GraphFunctions.h"

//TODO new format first line is node number of the graph (for old format use Update)
GraphData::GraphData(const std::string &graphPath) : Data<PUNGraph>(new TUNGraph()) {
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
                iss >> b;
                if (b.empty()){
                    for (int i = 0; i < std::stoi(a); ++i) {
                        this->graph()->AddNode(i);
                    }
                }
                else{
                    src = std::stoi(a);
                    dest = std::stoi(b);
                    if (!this->graph()->IsNode(src)){
                        this->graph()->AddNode(src);
                    }
                    if (!this->graph()->IsNode(dest)){
                        this->graph()->AddNode(dest);
                    }
                    this->get_graph()->AddEdge(src, dest);
                }
            }
            _name = std::filesystem::path(graphPath).stem().string();
        }
        catch (...) {
        }
    }
    if (std::filesystem::path(graphPath).extension() == ".bin"){
        ReadBinary(graphPath);
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
    DistanceList = std::vector<int>(size, -2);
    ContainmentList = std::vector<int>(size, -1);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

GraphData::GraphData(const std::string& graphPath, const std::string& labelPath) : GraphData(graphPath){
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
    DistanceList = std::vector<int>(this->size(), -2);
    ContainmentList = std::vector<int>(this->size(), -1);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

GraphData::GraphData(const PUNGraph & Graph, const std::string& name) : Data<PUNGraph>(Graph, name) {
    DistanceList = std::vector<int>(this->size(), -2);
    ContainmentList = std::vector<int>(this->size(), -1);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

GraphData::GraphData(const PUNGraph& Graph, const std::string& name, Labels& labels) : Data<PUNGraph>(Graph, name, labels) {
    DistanceList = std::vector<int>(this->size(), -2);
    ContainmentList = std::vector<int>(this->size(), -1);
    PredecessorsList = std::vector<std::vector<int>>(this->size(), std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

GraphData::GraphData(const PUNGraph& graph, int size, const std::string& name) : Data<PUNGraph>(graph, name) {
    DistanceList = std::vector<int>(size, -2);
    ContainmentList = std::vector<int>(size, -1);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

size_t GraphData::size() const {
    return this->get_graph()->GetNodes();
}

TUNGraph::TNodeI GraphData::node(NodeId Id) const {
    return get_graph()->GetNI(Id);
}

int GraphData::degree(NodeId Id) const {
    return node(Id).GetDeg();
}

NodeId GraphData::neighbor(int NodeId, int NeighborIdx) const {
    return node(NodeId).GetNbrNId(NeighborIdx);
}

void GraphData::save(const std::string& path) const {
    std::string labelPath = path + this->_name + ".labels";
    std::string edgePath = path + this->_name + ".edges";
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

void GraphData::save_edges(const std::string &path) const {
    std::string edgePath = path + this->_name + ".edges";
    std::ofstream file;
    int id = 0;
    file.open(edgePath);
    file << nodes() << std::endl;
    for (TUNGraph::TEdgeI edgeI = this->get_graph()->BegEI(); edgeI != this->get_graph()->EndEI(); edgeI++) {
        file << edgeI.GetSrcNId() << " " << edgeI.GetDstNId() << std::endl;
        ++id;
    }
    file.close();
}

NodeId GraphData::elem(NodeId Id) {
    return Id;
}

GraphData GraphData::mergeGraphs(const GraphData& other) const {
    PUNGraph punGraph = new TUNGraph();
    for(TUNGraph::TNodeI node = this->get_graph()->BegNI(); node != this->get_graph()->EndNI(); node++){
        punGraph->AddNode(node);
    }
    for(TUNGraph::TEdgeI edge = this->get_graph()->BegEI(); edge != this->get_graph()->EndEI(); edge++){
        punGraph->AddEdge(edge);
    }
    for(TUNGraph::TNodeI node = other.get_graph()->BegNI(); node != other.get_graph()->EndNI(); node++){
        punGraph->AddNode();
    }
    for(TUNGraph::TEdgeI edge = other.get_graph()->BegEI(); edge != other.get_graph()->EndEI(); edge++){
        punGraph->AddEdge(this->get_graph()->GetNodes() + edge.GetSrcNId(),
                          this->get_graph()->GetNodes() + edge.GetDstNId());
    }
    GraphData graph = GraphData(punGraph, "");
    return graph;
}

void GraphData::mergeGraphs(const GraphData &other) {
    int initialSize = this->size();
    for(TUNGraph::TNodeI node = other.get_graph()->BegNI(); node != other.get_graph()->EndNI(); node++){
        this->get_graph()->AddNode();
    }
    for(TUNGraph::TEdgeI edge = other.get_graph()->BegEI(); edge != other.get_graph()->EndEI(); edge++){
        this->get_graph()->AddEdge(initialSize + edge.GetSrcNId(), initialSize + edge.GetDstNId());
    }
}

PUNGraph GraphData::get_graph() const {
    return this->get_data();
}

PUNGraph GraphData::graph() {
    return this->data();
}

void GraphData::set_graph(PUNGraph graph) {
    set_data(graph);
}

GraphData::GraphData(const GraphData &other) : Data<PUNGraph>(other){
    this->DistanceList = other.DistanceList;
    this->ContainmentList = other.ContainmentList;
    this->PredecessorsList = other.PredecessorsList;
    this->neighborList = other.neighborList;
    this->graphType = other.graphType;
}

void GraphData::print() {
    std::cout << "Edges: " << std::endl;
    for(auto EdgeIt = get_graph()->BegEI(); EdgeIt != get_graph()->EndEI(); EdgeIt++){
        std::cout << EdgeIt.GetSrcNId() << "--" << EdgeIt.GetDstNId() << " ";
    }
    std::cout << std::endl;
}

void GraphData::save_dot(const std::string &path) const {

    std::string outPath = path + this->_name + ".dot";
    std::ofstream file;
    file.open(outPath);
    file << "strict graph {" << std::endl;
    for (auto EdgeIt = get_graph()->BegEI(); EdgeIt != get_graph()->EndEI(); EdgeIt++) {
        file << std::to_string(EdgeIt.GetSrcNId()) << " -- " << std::to_string(EdgeIt.GetDstNId()) << std::endl;
    }
    file << "}" << std::endl;
    file.close();
}

void GraphData::getDegreeDistribution(std::map<int, int> &degree_distribution) const {
    for (auto NodeIt = get_graph()->BegNI(); NodeIt != get_graph()->EndNI(); NodeIt++) {
        int degree = NodeIt.GetDeg();
        if (degree_distribution.find(degree) == degree_distribution.end()) {
            degree_distribution[degree] = 1;
        } else {
            degree_distribution[degree] += 1;
        }
    }
}

void GraphData::getDegreeDistribution(std::map<int, int> &degree_distribution, const std::set<NodeId>& NodeSet) const {
    degree_distribution.clear();
        for (auto const Node : NodeSet) {
            auto NodeIt = this->get_graph()->GetNI(Node);
            int degree = NodeIt.GetDeg();
            if (degree_distribution.find(degree) == degree_distribution.end()) {
                degree_distribution[degree] = 1;
            } else {
                degree_distribution[degree] += 1;
            }
        }
}

void GraphData::getDegreeDistribution(std::map<int, int> &degree_distribution, const TIntV& NodeSet) const {
    degree_distribution.clear();
        for (int i = 0; i < NodeSet.Len(); ++i) {
            auto NodeIt = this->get_graph()->GetNI(NodeSet[i]);
            int degree = NodeIt.GetDeg();
            if (degree_distribution.find(degree) == degree_distribution.end()) {
                degree_distribution[degree] = 1;
            } else {
                degree_distribution[degree] += 1;
            }
        }
}


GraphData::GraphData() {

}

int GraphData::nodes() const {
    return (int) size();
}

int GraphData::edges() const {
    return this->get_graph()->GetEdges();
}

double GraphData::density() const{
    return edges()/ ((double) nodes() * nodes());
}

void GraphData::UpdateSingletons(const std::string &graphPath) {
    GraphData graphData = GraphData(new TUNGraph());
    std::string name = std::filesystem::path(graphPath).stem().string();
    std::cout << "Updating " << name << std::filesystem::path(graphPath).extension() << std::endl;
    std::string extension = std::filesystem::path(graphPath).extension().string();
    if (extension == ".edges" || extension == ".txt") {
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
                iss >> b;
                if (b.empty()){
                    for (int i = 0; i < std::stoi(a); ++i) {
                        graphData.graph()->AddNode(i);
                    }
                }
                else{
                    src = std::stoi(a);
                    dest = std::stoi(b);
                    if (!graphData.graph()->IsNode(src)){
                        graphData.graph()->AddNode(src);
                    }
                    if (!graphData.graph()->IsNode(dest)){
                        graphData.graph()->AddNode(dest);
                    }
                    graphData.graph()->AddEdge(src, dest);
                }
            }

            graphData.setName(name);
        }
        catch (...) {
        }
        for (int i = 0; i < graphData.nodes(); ++i) {
            if (!graphData.get_graph()->IsNode(i)){
                graphData.graph()->AddNode(i);
            }
        }
        if (graphData.nodes() != graphData.get_graph()->GetMxNId()){
            throw std::range_error("Nodes " + std::to_string(graphData.nodes()) + " and Ids (" + std::to_string(graphData.get_graph()->GetMxNId()) +") do not fit for the graph! " + graphData.getName() + "\n");
        }
    }
    std::string parent_path = std::filesystem::path(graphPath).parent_path().string();
    graphData.save_edges(parent_path + "/");

}

//Update txt graphs to new binary format!!!
void GraphData::Update(const std::string &graphPath) {
    GraphData graph = GraphData(new TUNGraph());
    std::string name = std::filesystem::path(graphPath).stem().string();
    std::cout << "Updating " << name << std::filesystem::path(graphPath).extension() << std::endl;
    std::string extension = std::filesystem::path(graphPath).extension().string();
    if (extension == ".edges" || extension == ".txt") {
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
                iss >> b;
                if (b.empty()){
                    for (int i = 0; i < std::stoi(a); ++i) {
                        graph.graph()->AddNode(i);
                    }
                }
                else{
                    src = std::stoi(a);
                    dest = std::stoi(b);
                    if (!graph.graph()->IsNode(src)){
                        graph.graph()->AddNode(src);
                    }
                    if (!graph.graph()->IsNode(dest)){
                        graph.graph()->AddNode(dest);
                    }
                    graph.graph()->AddEdge(src, dest);
                }
            }


        }
        catch (...) {
        }
        GraphData graphData = GraphData(GraphFunctions::ResetGraphIds(graph.get_graph()));
        graphData.setName(name);
        for (int i = 0; i < graphData.nodes(); ++i) {
            if (!graphData.get_graph()->IsNode(i)){
                graphData.graph()->AddNode(i);
            }
        }
        if (graphData.nodes() != graphData.get_graph()->GetMxNId()){
            throw std::range_error("Nodes " + std::to_string(graphData.nodes()) + " and Ids (" + std::to_string(graphData.get_graph()->GetMxNId()) +") do not fit for the graph! " + graphData.getName() + "\n");
        }
        std::string parent_path = std::filesystem::path(graphPath).parent_path().string();
        //graphData.save_edges(parent_path + "/");
        //graphData.save_bin(parent_path + "/");
    }
}

void GraphData::init(int size) {
    DistanceList = std::vector<int>(size, -2);
    ContainmentList = std::vector<int>(size, -1);
    PredecessorsList = std::vector<std::vector<int>>(size, std::vector<int>());
    GraphFunctions::generateNeighborVector(this->get_graph(), neighborList);
    maxDegree = (int) neighborList.size();
}

void GraphData::init() {
    init((int) this->size());
}

void GraphData::WriteBinary(const std::string& graphPath, bool Labeled, bool OnlyGraph) const {
    std::ofstream Out(graphPath + this->getName() + ".bin", std::ios::out | std::ios::binary);
    const std::string Name = this->getName();
    unsigned int stringLength = Name.length();
    GraphType Type = this->graphType;
    Out.write((char*) (&stringLength), sizeof(stringLength));
    Out.write(Name.c_str(), stringLength);
    Out.write((char*) (&Type), sizeof(GraphType));
    int Size = (int) this->size();
    int Edges = this->edges();
    Out.write((char*) (&Size), sizeof(int));
    Out.write((char*) (&Edges), sizeof(int));
    for (auto EdgeIt = this->get_graph()->BegEI(); EdgeIt != this->get_graph()->EndEI(); EdgeIt++) {
        int Src = EdgeIt.GetSrcNId();
        int Dst = EdgeIt.GetDstNId();
        Out.write((char*) (&Src), sizeof(int));
        Out.write((char*) (&Dst), sizeof(int));
    }
    Out.write((char*) (&Labeled), sizeof(bool));
    if (Labeled){
        for (auto const & L : get_labels()) {
            Out.write((char*) (&L), sizeof(Label));
        }
    }
    Out.write((char*) (&OnlyGraph), sizeof(bool));
    if (!OnlyGraph){
    }
    Out.close();
}

void GraphData::ReadBinary(const std::string &graphPath) {
    auto size = std::filesystem::file_size(graphPath);
    std::ifstream In(graphPath, std::ios::in | std::ios::binary);
    // Check if file is good
    if (!In) {
        std::cout << "Error occurred" << std::endl;
        exit(-1);
    }
    std::string Name;
    unsigned int stringLength;
    In.read( (char*)( &stringLength ), sizeof( stringLength ) );
    Name.resize( stringLength );
    In.read( (char*)Name.c_str(), stringLength );
    this->setName(Name);
    In.read((char*) (&this->graphType), sizeof(GraphType));
    int Size = 0;
    int Edges = 0;
    In.read((char*) (&Size), sizeof(int));
    In.read((char*) (&Edges), sizeof(int));
    for (int i = 0; i < Size; ++i) {
        this->graph()->AddNode(i);
    }
    for (int i = 0; i < Edges; ++i) {
        int Src = 0;
        int Dst = 0;
        In.read((char*) (&Src), sizeof(int));
        In.read((char*) (&Dst), sizeof(int));
        this->graph()->AddEdge(Src, Dst);
    }
    bool Labeled;
    In.read((char*) (&Labeled), sizeof(bool));
    if (Labeled){
        Label L;
        for (int i = 0; i < Size; ++i) {
            In.read((char*) (&L), sizeof(Label));
            this->labels().emplace_back(L);
        }
    }
    bool OnlyGraph;
    In.read((char*) (&OnlyGraph), sizeof(bool));
    if (!OnlyGraph){
    }
    In.close();
}

void GraphData::save_bin(const std::string &path) const {
    WriteBinary(path);
}





//
// Created by anonymous on 12.08.21.
//

#include "OuterplanarSubgraph.h"
#include "GraphFunctions.h"

OuterplanarSubgraph::OuterplanarSubgraph(const PUNGraph graph) : _graph(graph){}

OuterplanarSubgraph::OuterplanarSubgraph(const OuterplanarSubgraph &other) : _graph(other._graph), _gen(other._gen), print(other.print) {
}

void OuterplanarGraphStatistics::getStatistics(const PUNGraph &outerplanarGraph) {
    std::vector<PUNGraph> components;
    GraphFunctions::GetBiconnectedComponents(outerplanarGraph, components);
    this->ComponentSizes.clear();
    this->ComponentFaces.clear();
    this->ComponentSizes.emplace_back(std::vector<int>(0));
    this->ComponentFaces.emplace_back(std::vector<int>(0));
    for (const auto& component : components) {
        this->ComponentSizes.back().emplace_back(component->GetNodes());
        this->ComponentFaces.back().emplace_back(0.0);
        GraphFunctions::GetBiconnectedOuterplanarFaceNum(component, this->ComponentFaces.back().back());
        component->Clr();
    }
    this->ComponentNumber.emplace_back(static_cast<int>(this->ComponentSizes.back().size()));
    this->BiggestComponent.emplace_back(*std::max_element(this->ComponentSizes.back().begin(), this->ComponentSizes.back().end()));
    this->FaceNumbers.emplace_back(*std::max_element(this->ComponentFaces.back().begin(), this->ComponentFaces.back().end()));
}

OuterplanarGraphStatistics::OuterplanarGraphStatistics(const PUNGraph outerplanarGraph) {
    getStatistics(outerplanarGraph);
}

void OuterplanarGraphStatistics::evaluate(std::vector<std::string> &headers, std::vector<std::string> &values) const {
    headers.clear();
    values.clear();
    headers = {"Avg. Component Number", "Std. Component Number", "Avg. Biggest Component", "Std. Biggest Component", "Avg. Face Number", "Std. Face Number"};
    values = {std::to_string(StaticFunctions::mean(ComponentNumber)),
              std::to_string(StaticFunctions::standard_deviation(ComponentNumber)),
              std::to_string(StaticFunctions::mean(BiggestComponent)),
              std::to_string(StaticFunctions::standard_deviation(BiggestComponent)),
              std::to_string(StaticFunctions::mean(FaceNumbers)),
              std::to_string(StaticFunctions::standard_deviation(FaceNumbers))};
}
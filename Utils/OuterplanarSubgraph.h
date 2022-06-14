//
// Created by anonymous on 12.08.21.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPH_H
#define CLOSURES_OUTERPLANARSUBGRAPH_H

#include "../Utils/typedefs.h"
#include "../Data/GraphData.h"
#include "StaticFunctions.h"

class OuterplanarSubgraph {
public:
    explicit OuterplanarSubgraph(const PUNGraph graph);

    //copying
    OuterplanarSubgraph(const OuterplanarSubgraph& other);

    virtual void generate(GraphData& subgraph, std::mt19937_64& gen, bool print) = 0;

protected:
    std::mt19937_64 _gen;
    const PUNGraph _graph;
    bool print = false;
};

struct OuterplanarGraphStatistics{
    OuterplanarGraphStatistics() = default;
    explicit OuterplanarGraphStatistics(const PUNGraph outerplanarGraph);
    std::vector<std::vector<int>> ComponentSizes;
    std::vector<std::vector<int>> ComponentFaces;
    std::vector<double> ComponentNumber;
    std::vector<double> BiggestComponent;
    std::vector<double> FaceNumbers;
    void operator +=(const OuterplanarGraphStatistics& other){
        ComponentSizes.insert(ComponentSizes.end(), other.ComponentSizes.begin(),  other.ComponentSizes.end());
        ComponentNumber.insert(ComponentNumber.end(), other.ComponentNumber.begin(),  other.ComponentNumber.end());
        BiggestComponent.insert(BiggestComponent.end(), other.BiggestComponent.begin(),  other.BiggestComponent.end());
        FaceNumbers.insert(FaceNumbers.end(), other.FaceNumbers.begin(), other.FaceNumbers.end());
    };
    void getStatistics(const PUNGraph& outerplanarGraph);
    void evaluate(std::vector<std::string>& headers, std::vector<std::string>& values) const;
};


#endif //CLOSURES_OUTERPLANARSUBGRAPH_H

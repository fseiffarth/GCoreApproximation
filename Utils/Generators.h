//
// Created by anonymous on 08.05.2021.
//

#ifndef DS2021_GENERATORS_H
#define DS2021_GENERATORS_H


#include <string>
#include <utility>
#include <random>
#include "../Data/GraphData.h"

struct Properties{
public:
    Properties(int sizeA, double densityA, int sizeB, double densityB, int connectivity, const std::string& name)
            : sizeA(sizeA), densityA(densityA), sizeB(sizeB), densityB(densityB), connectivity(connectivity) {
        this->_size = sizeA + sizeB;
        this->name = name + "_A_" + std::to_string(sizeA)  + "_dA_" + std::to_string((int) (densityA*10)) + "_B_" + std::to_string(sizeB)  + "_dB_" + std::to_string((int) (densityB*10)) + "_C_" + std::to_string(connectivity) + "_";
    }
    Properties(int sizeA, int sizeB, double densityA, const std::string& name)
            : sizeA(sizeA), densityA(densityA), sizeB(sizeB){
        this->_size = sizeA + sizeB;
        this->name = name + "_A_" + std::to_string(sizeA)  + "_B_" + std::to_string(sizeB)  + "_d_" + std::to_string((int) (densityA*10)) + "_";
    }
    Properties(int size, double density, const std::string& name)
            : _size(size), _density(density){
        this->name = name + "_Size_" + std::to_string(_size)  + "_Density_" + std::to_string((int) (_density*10)) + "_";
    }
    int _size;
    double _density;

    int sizeA;
    double densityA;
    int sizeB;
    double densityB;
    int connectivity;
    std::string name{};
};

class Generators {
public:
    static void generateConnectedRandomGraphs(const std::string& path, std::vector<GraphData>& graphs, const Properties& properties, int number, std::mt19937_64& gen, bool save = true);
    static void generateTwoComponentGraphs(const std::string& path, std::vector<GraphData>& graphs, const Properties& properties, int number, std::mt19937_64& gen, bool random_connections = false, bool save = true);
    static void generateOneComponentGraphs(const std::string& path, std::vector<GraphData>& graphs, const Properties& properties, int number, std::mt19937_64& gen, bool save = true);

    static void generateConnectedRandomOuterplanarGraphs(const std::string &path, std::vector<GraphData> &graphs,
                                                         const Properties &properties, int number, std::mt19937_64 &gen,
                                                         bool save = true);
};


#endif //DS2021_GENERATORS_H

//
// Created by anonymous on 08.05.2021.
//

#include "Generators.h"
#include "graph_generators.h"
#include "drawing.h"
#include "OuterplanarSubgraphDFS.h"
#include <iostream>

void Generators::generateTwoComponentGraphs(const std::string& path, std::vector<GraphData>& graphs, const Properties& properties, int number, std::mt19937_64& gen, bool random_connections, bool save) {
    int counter = 0;
    while(counter < number) {
        bool valid;
        GraphClosureSP graphClosureSp("ShortestPath");
        GraphData graph = GraphGeneration::generateSeparableComponents(graphClosureSp, properties.sizeA,
                                                                       properties.densityA, properties.sizeB,
                                                                       properties.densityB, properties.connectivity,
                                                                       valid, random_connections, gen);
        std::cout << omp_get_thread_num() << "Is _valid? " << valid << std::endl;
        if (valid) {
            std::cout << "Graph Size: " << properties.sizeA + properties.sizeB << ", Density: " << properties.densityA << ", Connections: " << properties.connectivity << std::endl;
            std::cout << counter + 1 << "/" << number << std::endl;
            if (save) {
                graph.save(path + properties.name + std::to_string(counter));
            }
            graph.setName(properties.name + std::to_string(counter));
            graphs.emplace_back(graph);
            ++counter;
        }
    }
}

void Generators::generateOneComponentGraphs(const std::string& path, std::vector<GraphData>& graphs, const Properties& properties, int number, std::mt19937_64& gen, bool save) {
    int counter = 0;
    while(counter < number) {
        bool valid;
        GraphClosureSP graphClosureSp("ShortestPath");
        GraphData graph = GraphGeneration::generateOneComponent(graphClosureSp, properties.sizeA, properties.sizeB,
                                                                       properties.densityA,
                                                                       valid, gen);
        std::cout << omp_get_thread_num() << "Is _valid? " << valid << std::endl;
        if (valid) {
            std::cout << "Graph Size: " << properties.sizeA + properties.sizeB << ", Density: " << properties.densityA << ", Connections: " << properties.connectivity << std::endl;
            std::cout << counter + 1 << "/" << number << std::endl;
            if (save) {
                graph.save(path + properties.name + std::to_string(counter));
            }
            graph.setName(properties.name + std::to_string(counter));
            graphs.emplace_back(graph);
            ++counter;
        }
    }
}

void Generators::generateConnectedRandomGraphs(const std::string &path, std::vector<GraphData> &graphs,
                                               const Properties &properties, int number, std::mt19937_64 &gen,
                                               bool save) {
    int counter = 0;
    while(counter < number) {
        GraphData graph = GraphGeneration::generateGraph(properties._size, properties._density, gen);
        std::cout << "Graph Size: " << properties._size << ", Density: " << properties._density << std::endl;
        std::cout << counter + 1 << "/" << number << std::endl;
        if (save) {
            graph.save(path + "Graph_" + properties.name + std::to_string(counter));
        }
        graph.setName(properties.name + std::to_string(counter));
        graphs.emplace_back(graph);
        ++counter;
    }
}

void Generators::generateConnectedRandomOuterplanarGraphs(const std::string &path, std::vector<GraphData> &graphs,
                                                          const Properties &properties, int number, std::mt19937_64 &gen,
                                                          bool save) {
    int counter = 0;

    while(counter < number) {
        GraphData subgraph = GraphData(new TUNGraph());
        GraphData graph = GraphGeneration::generateGraph(properties._size, properties._density, gen);
        OuterplanarSubgraphDFS subgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
        subgraphDfs.generate(subgraph, gen, false);
        GraphData subgraphData = GraphData(subgraph);
        std::cout << "Graph Size: " << properties._size << ", Density: " << properties._density << std::endl;
        std::cout << counter + 1 << "/" << number << std::endl;
        if (save) {
            subgraphData.save(path + "Graph_" + properties.name + std::to_string(counter));
        }
        subgraphData.setName(properties.name + std::to_string(counter));
        graphs.emplace_back(subgraphData);
        ++counter;
    }
}

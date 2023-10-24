//
// Created by florian on 15.08.23.
//

#include <random>
#include <unordered_set>
#include <iostream>
#include "CoreGrowAlgorithm.h"
#include "typedefs.h"
#include "../ClosureOperators/GraphClosures.h"


void CoreGrowAlgorithm::Run(CoreGrowAlgorithmOutputParameters& outputParameters, const CoreGrowAlgorithmInputParameters& inputParameters)
{
    // set start time
    auto start = std::chrono::high_resolution_clock::now();
    // set the generator using the seed
    std::mt19937 generator(inputParameters.seed);
    // neighbors of the core nodes (in fact the out_edges thus duplicates are possible)
    std::vector<NodeId> core_neighbors = std::vector<NodeId>();
    // vector measuring the size evolution of the core
    std::vector<int> core_size_evolution = std::vector<int>(inputParameters.grow_steps + 1, 0);
    core_size_evolution[0] = 1;
    std::vector<int> core_nodes = std::vector<int>(graph.nodes(), 0);
    // iterate over the number of runs
    for (int i = 0; i < inputParameters.num_runs; ++i)
    {
        // get the start vertex
        NodeId start_vertex = std::uniform_int_distribution<NodeId>(0, graph.nodes() - 1)(generator);
        core_neighbors = graph.get_neighbors(start_vertex);

        // create the graph closure
        GraphClosureSP graphClosureSP = GraphClosureSP();
        // get the closure of the start vertex
        std::set<NodeId> start_set = std::set<NodeId>({start_vertex});
        ClosureParameters closureParameters = ClosureParameters(start_set);
        graphClosureSP.closure(graph, closureParameters);


        // iterate over the grow steps
        for (int j = 0; j < inputParameters.grow_steps; ++j)
        {
            // print the grow step
            if (inputParameters.print)
            {
                std::cout << "Run: " << i << " Grow step: " << j << std::endl;
            }
            // delete all core neighbors that are already in the closureParameters.closed_set
            std::vector<NodeId> core_neighbors_copy;
            for (auto node : core_neighbors) {
                if (closureParameters.closed_set.find(node) == closureParameters.closed_set.end())
                {
                    core_neighbors_copy.emplace_back(node);
                }

            }
            core_neighbors = core_neighbors_copy;
            // add the neighbors of all added elements to the core neighbors if the neighbor is not already in the core nodes
            for (auto added_element : closureParameters.added_elements)
            {
                // get the neighbors of the added element
                std::vector<NodeId> added_element_neighbors = graph.get_neighbors(added_element);
                // iterate over the neighbors of the added element
                for (auto added_element_neighbor : added_element_neighbors)
                {
                    // add the neighbor to the core neighbors if it is not already in the core nodes
                    if (closureParameters.closed_set.find(added_element_neighbor) == closureParameters.closed_set.end())
                    {
                        core_neighbors.emplace_back(added_element_neighbor);
                    }
                }
            }

            // get a random element from the core neighbors
            NodeId random_element = core_neighbors[std::uniform_int_distribution<NodeId>(0, (NodeId) core_neighbors.size() - 1)(generator)];
            // calculate the closure of the core + the random element
            closureParameters.input_set = closureParameters.closed_set;
            closureParameters.element_to_add = random_element;
            graphClosureSP.closure(graph, closureParameters);

            // track the core growth
            core_size_evolution[i+1] = (int) closureParameters.added_elements.size();
        }
        // print the closure size
        if (inputParameters.print)
        {
            std::cout << "Closure size: " << closureParameters.closed_set.size() << std::endl;
        }
        // add the core nodes to the output parameters
        for (auto node : closureParameters.closed_set)
        {
            core_nodes[node] += 1;
        }
    }
    // print the core nodes vector
    if (inputParameters.print)
    {
        std::cout << "Core nodes: " << std::endl;
        for (int i = 0; i < core_nodes.size(); ++i)
        {
            std::cout << i << ": " << core_nodes[i] << std::endl;
        }
    }
    // iterate over core_nodes and add to output parameters if the node is in core node percentage of the runs
    for (int i = 0; i < core_nodes.size(); ++i)
    {
        if (core_nodes[i] >= int(inputParameters.core_percentage * inputParameters.num_runs + 0.99))
        {
            outputParameters.core_nodes.emplace_back(i);
        }
    }
    // set the end time
    outputParameters.runtime = (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000.0;
}
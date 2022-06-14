//
// Created by anonymous on 07.11.21.
//

#include "Experiments.h"
#include "../Utils/GraphFunctions.h"
#include "../Utils/OuterplanarSubgraphDFS.h"
#include "../Utils/StaticFunctions.h"
#include "../Utils/Enums.h"
#include "../Data/OuterplanarGraphData.h"

void Experiments::GetOuterplanarSamples(const std::pair<int, int>& size, int number, std::vector<GraphData>& outerplanar, std::vector<GraphData>& graphs, std::vector<GraphData>& trees) {
    std::vector<int> neighbors = std::vector<int>(size.first);
    std::iota(neighbors.begin(), neighbors.end(), 0);
    int nonOuterPlanarNumber = 0;
    std::vector<int> nonOuterPlanarSeeds;
    int nonMaximalNumber = 0;
    std::vector<int> nonMaximalSeeds;
    TRnd rnd = TInt::Rnd;
    std::vector<int> maximalEdges;
    double nonMaximalEdgeNumPerCent = 0.0;
    int outerPlanarEdgeNum = 0;
    int maximalEdgesNum = 0;
    std::vector<int> algorithmMissingEdges;
    int missingEdgesNum = 0;
    size_t outerplanar_runtime = 0;
    size_t tree_runtime = 0;
    OuterplanarGraphStatistics privateStatistics = OuterplanarGraphStatistics();
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    rnd.PutSeed(0);
    for (int i = 0; i < number; ++i) {
        auto graph = TSnap::GenRndGnm<PUNGraph>(size.first, size.second, false, rnd);
        if (TSnap::IsConnected(graph)) {
            //std::cout << "\tSeed: " << i << std::endl;
            std::mt19937_64 gen(i);
            GraphData subgraph = GraphData(new TUNGraph(), graph->GetNodes());
            auto start = std::chrono::high_resolution_clock::now();
            OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graph);
            subgraphGeneration.generate(subgraph, gen, false);
            outerplanar_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start).count();
            outerPlanarEdgeNum += subgraph.get_graph()->GetEdges();
            outerplanar.emplace_back(subgraph);
            int i = subgraph.edges();
            GraphData subtree = GraphData(new TUNGraph, subgraph.nodes());
            GraphFunctions::bfsSubtree(subgraph.get_graph(), subtree, neighbors, gen);
            i = subtree.edges();
            trees.emplace_back(subtree);
            GraphData rand_graph = GraphData(new TUNGraph, graph->GetNodes());
            for(auto NodeIt = subtree.get_graph()->BegNI(); NodeIt != subtree.get_graph()->EndNI(); NodeIt++){
                rand_graph.graph()->AddNode(NodeIt.GetId());
            }
            for(auto EdgeIt = subtree.get_graph()->BegEI(); EdgeIt != subtree.get_graph()->EndEI(); EdgeIt++){
                rand_graph.graph()->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
            }
            GraphFunctions::addRandomEdges(rand_graph, subgraph.edges() - subtree.edges(), i);
            i = rand_graph.edges();
            graphs.emplace_back(rand_graph);
        } else {
            std::cout << "Not connected!" << std::endl;
            --i;
        }
        graph.Clr();
    }
}

void Experiments::OuterplanarSampling(const std::string& out_path, Experiments::TestParameter& testParameter) {
    std::vector<int> neighbors = std::vector<int>();
    for (auto const &graph_info: testParameter.graph_sizes) {
        FileEvaluation fileEvaluation = FileEvaluation(out_path, "");
        int nonOuterPlanarNumber = 0;
        std::vector<int> nonOuterPlanarSeeds;
        int nonMaximalNumber = 0;
        std::vector<int> nonMaximalSeeds;
        TRnd rnd = TInt::Rnd;

        std::vector<double> maximalEdges;
        std::vector<double> nonMaximalEdgeNumPerCent;
        std::vector<double> outerPlanarEdgeNum;
        std::vector<double> maximalEdgesNum;
        std::vector<double> algorithmMissingEdges;
        std::vector<double> missingEdgesNum;
        size_t outerplanar_runtime = 0;
        size_t outerplanar_structure_runtime = 0;
        size_t bfs_tree_runtime = 0;
        size_t dfs_tree_runtime = 0;
        size_t snap_bfs_runtime = 0;
        OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
        int iterations = testParameter.seedMax + 1 - testParameter.seedMin;
        rnd.PutSeed(0);
        for (int seed = testParameter.seedMin; seed < testParameter.seedMax + 1; ++seed) {
            auto graph = TSnap::GenRndGnm<PUNGraph>(graph_info.first, graph_info.second, false, rnd);
            if (TSnap::IsConnected(graph)) {
                std::cout << "Seed: " << seed << std::endl;
                if (seed >= testParameter.seedMin) {
                    std::mt19937_64 gen(seed);
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        GraphData subgraph = GraphData(new TUNGraph(), graph->GetNodes());
                        GraphFunctions::generateNeighborVector(graph, neighbors);
                        GraphFunctions::bfsSubtree(graph, subgraph, neighbors, gen);
                        bfs_tree_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - start).count();
                    }
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        GraphData subgraph = GraphData(new TUNGraph(), graph->GetNodes());
                        GraphFunctions::generateNeighborVector(graph, neighbors);
                        GraphFunctions::dfsSubtree(graph, subgraph, neighbors, gen);
                        dfs_tree_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - start).count();
                    }
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        TSnap::GetBfsTree(graph, std::uniform_int_distribution<int>(0, graph->GetNodes() - 1)(gen),
                                          true,
                                          true);
                        snap_bfs_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - start).count();
                    }
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        GraphData subgraph = GraphData(new TUNGraph(), graph->GetNodes());
                        OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graph);
                        subgraphGeneration.generate(subgraph, gen, false);

                        outerplanar_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - start).count();
                        start = std::chrono::high_resolution_clock::now();
                        OuterplanarGraphData outerplanarSubgraph = OuterplanarGraphData(subgraph);
                        outerplanar_structure_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - start).count();

                        //Check outerplanarity + quality
                        if (testParameter.outerplanarity_check) {
                            GraphFunctions::checkingOuterplanarity(graph, subgraph.get_graph(), nonOuterPlanarNumber,
                                                                   nonMaximalNumber,
                                                                   nonOuterPlanarSeeds,
                                                                   nonMaximalSeeds, algorithmMissingEdges, maximalEdges,
                                                                   seed);
                            maximalEdgesNum.emplace_back(subgraph.get_graph()->GetEdges() + algorithmMissingEdges.back());
                            nonMaximalEdgeNumPerCent.emplace_back(
                                    (double) (subgraph.get_graph()->GetEdges()) /
                                            ((double) subgraph.get_graph()->GetEdges() + algorithmMissingEdges.back())*100.0);
                            missingEdgesNum.emplace_back(algorithmMissingEdges.back());
                        }
                        statistics += OuterplanarGraphStatistics(subgraph.get_graph());
                        outerPlanarEdgeNum.emplace_back(subgraph.get_graph()->GetEdges());
                    }
                }
            }
            else {
                std::cout << "Not connected!" << std::endl;
                --seed;
            }
        }
        std::vector<std::string> stat_headers;
        std::vector<std::string> stat_values;
        statistics.evaluate(stat_headers, stat_values);

        std::vector<std::string> headers = {"Size",
                                            "Edges",
                                            "Density",
                                            "O1 Time",
                                            "O2 Time",
                                            "Outerplanar Struct Time",
                                            "BFS Time",
                                            "Snap BFS",
                                            "DFS Time",
                                            "Factor",
                                            "Samples",
                                            "Not Outerplanar Samplings",
                                            "Not Maximal Samplings",
                                            "Avg. Maximal Output Edges",
                                            "Std. Maximal Output Edges",
                                            "Avg. Output Edges",
                                            "Std. Output Edges",
                                            "Avg. Missing Edges per Graph",
                                            "Std. Missing Edges per Graph",
                                            "Avg. % of Maximal outerplanar graph",
                                            "Std. % of Maximal outerplanar graph"};

        std::vector<std::string> values = {std::to_string(graph_info.first),
                                           std::to_string(graph_info.second),
                                           std::to_string(graph_info.second / ((double) graph_info.first/2 * (graph_info.first - 1))),
                                           std::to_string((double) outerplanar_runtime / 1000000),
                                           std::to_string((double) outerplanar_runtime / 1000000 + (double) outerplanar_structure_runtime / 1000000),
                                           std::to_string((double) outerplanar_structure_runtime / 1000000),
                                           std::to_string((double) bfs_tree_runtime / 1000000),
                                           std::to_string((double) snap_bfs_runtime / 1000000),
                                           std::to_string((double) dfs_tree_runtime / 1000000),
                                           std::to_string((double) outerplanar_runtime / (double) bfs_tree_runtime),
                                           std::to_string(iterations),
                                           std::to_string(nonOuterPlanarNumber),
                                           std::to_string(nonMaximalNumber),
                                           std::to_string(StaticFunctions::mean(maximalEdgesNum)),
                                           std::to_string(StaticFunctions::standard_deviation(maximalEdgesNum)),
                                           std::to_string(StaticFunctions::mean(outerPlanarEdgeNum)),
                                           std::to_string(StaticFunctions::standard_deviation(outerPlanarEdgeNum)),
                                           std::to_string(StaticFunctions::mean(missingEdgesNum)),
                                           std::to_string(StaticFunctions::standard_deviation(missingEdgesNum)),
                                           std::to_string(StaticFunctions::mean(nonMaximalEdgeNumPerCent)),
                                           std::to_string(StaticFunctions::standard_deviation(nonMaximalEdgeNumPerCent))};

        fileEvaluation.headerValueInsert(headers, values);
        fileEvaluation.headerValueInsert(stat_headers, stat_values);
#pragma omp critical
        {
            fileEvaluation.save();

            std::cout << std::endl << "******************************* FileEvaluation "
                      << " *********************************" << std::endl;
            std::cout << "Runtime: " << (double) outerplanar_runtime << " microseconds" << ", Factor: "
                      << (double) outerplanar_runtime / (double) bfs_tree_runtime << std::endl;
            std::cout << "BFSTime: " << (double) bfs_tree_runtime << " microseconds" << std::endl;
            std::cout << "Nodes: " << graph_info.first << std::endl;
            std::cout << "Edges: " << graph_info.second << std::endl;
            std::cout << "Density: " << graph_info.second / ((double) graph_info.first / 2 * (graph_info.first - 1))
                      << std::endl;
            std::cout << "Samples: " << testParameter.seedMax - testParameter.seedMin + 1 << std::endl;
            std::cout << "\tNot Outerplanar: " << nonOuterPlanarNumber << " "
                      << StaticFunctions::print<std::vector<int>, int>(nonOuterPlanarSeeds) << std::endl;
            std::cout << "\tNot Maximal: " << nonMaximalNumber << " "
                      << StaticFunctions::print<std::vector<int>, int>(nonMaximalSeeds) << std::endl;
            std::cout << "Average Edges: " << StaticFunctions::mean(outerPlanarEdgeNum) << std::endl;
            std::cout << "Average Maximal Edges: " << std::accumulate(maximalEdges.begin(), maximalEdges.end(), 0.0) /
                       (testParameter.seedMax - testParameter.seedMin + 1) << std::endl;
            std::cout << "Missing Edges per Graph: " << std::accumulate(algorithmMissingEdges.begin(), algorithmMissingEdges.end(), 0.0) /
                       (double) algorithmMissingEdges.size() << std::endl;
            std::cout << "Edges per Graph in %: "
                      << StaticFunctions::mean(nonMaximalEdgeNumPerCent) << std::endl;
        }
    }
}

void Experiments::GetOuterplanarSamplesBiconnected(const std::pair<int, int> &size, int number,
                                                   std::vector<GraphData> &outerplanar, std::vector<GraphData> &graphs,
                                                   std::vector<GraphData> &trees) {
    std::vector<int> neighbors = std::vector<int>(size.first);
    std::iota(neighbors.begin(), neighbors.end(), 0);
    int nonOuterPlanarNumber = 0;
    std::vector<int> nonOuterPlanarSeeds;
    int nonMaximalNumber = 0;
    std::vector<int> nonMaximalSeeds;
    TRnd rnd = TInt::Rnd;
    std::vector<int> maximalEdges;
    double nonMaximalEdgeNumPerCent = 0.0;
    int outerPlanarEdgeNum = 0;
    int maximalEdgesNum = 0;
    std::vector<int> algorithmMissingEdges;
    int missingEdgesNum = 0;
    size_t outerplanar_runtime = 0;
    size_t tree_runtime = 0;
    OuterplanarGraphStatistics privateStatistics = OuterplanarGraphStatistics();
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    rnd.PutSeed(0);
    for (int i = 0; i < number; ++i) {
        auto graph = TSnap::GenRndGnm<PUNGraph>(size.first, size.second, false, rnd);
        if (TSnap::IsConnected(graph)) {
            std::cout << "Seed: " << i << std::endl;
            std::mt19937_64 gen(i);
            OuterplanarGraphData subgraph = OuterplanarGraphData(new TUNGraph(), graph->GetNodes());
            auto start = std::chrono::high_resolution_clock::now();
            OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graph);
            subgraphGeneration.generate(subgraph, gen, false);
            outerplanar_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start).count();
            outerPlanarEdgeNum += subgraph.get_graph()->GetEdges();
            outerplanar.emplace_back(subgraph);
            int i = subgraph.edges();
            GraphData subtree = GraphData(new TUNGraph, subgraph.nodes());
            GraphFunctions::bfsSubtree(subgraph.get_graph(), subtree, neighbors, gen);
            i = subtree.edges();
            trees.emplace_back(subtree);
            GraphData rand_graph = GraphData(new TUNGraph, graph->GetNodes());
            for(auto NodeIt = subtree.get_graph()->BegNI(); NodeIt != subtree.get_graph()->EndNI(); NodeIt++){
                rand_graph.graph()->AddNode(NodeIt.GetId());
            }
            for(auto EdgeIt = subtree.get_graph()->BegEI(); EdgeIt != subtree.get_graph()->EndEI(); EdgeIt++){
                rand_graph.graph()->AddEdge(EdgeIt.GetSrcNId(), EdgeIt.GetDstNId());
            }
            GraphFunctions::addRandomEdges(rand_graph, subgraph.edges() - subtree.edges(), i);
            i = rand_graph.edges();
            graphs.emplace_back(rand_graph);
        } else {
            std::cout << "Not connected!" << std::endl;
            --i;
        }
        graph.Clr();
    }
}




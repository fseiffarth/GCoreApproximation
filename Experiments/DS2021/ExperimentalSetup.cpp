//
// Created by anonymous on 10.05.2021.
//

#include <string>
#include <iostream>
#include "ExperimentalSetup.h"
#include "../../Data/GraphData.h"
#include "../../Utils/graphio.h"
#include "../../Utils/pattern_generator.h"
#include "../../Operators/GraphClosures.h"
#include "../../Evaluations/evaluations.h"
#include "../../Algorithms/MCSS.h"
#include "../../Utils/Generators.h"
#include <omp.h>


ExperimentalSetup::ExperimentalSetup(std::string resultsName, std::string resultPath, std::string graphPath, std::string treePath,
                                     std::string outerplanarPath, int graphNumber,
                                     const std::vector<int> &blockSizes, std::vector<double> blockDensities,
                                     std::vector<int> blockConnections, std::vector<int> trainingSizes,
                                     int patterns, int pattern_out_stepsize, int pattern_step_ int trainingConfigurations, bool randomConnections, int threads)
        : resultsName(resultsName), blockSizes(blockSizes), blockDensities(blockDensities), blockConnections(blockConnections),
          trainingSizes(trainingSizes), patterns(patterns), trainingConfigurations(trainingConfigurations), resultPath(resultPath),
          graphPath(graphPath), graphNumber(graphNumber), treePath(treePath), outerplanarPath(outerplanarPath), _out_stepsize(pattern_out_stepsize), _random_connections(randomConnections){
    omp_set_num_threads(threads);
}


void ExperimentalSetup::generateGraphs(int seed) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int blockSize : blockSizes) {
        std::vector<std::pair<double, int>> pairs;
        for (double blockDensity : blockDensities) {
            if (blockConnections.empty()){
                pairs.emplace_back(std::pair<double, int>{blockDensity, blockSize/2});
            }
            else {
                for (int blockConnection : blockConnections) {
                    pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
                }
            }
        }
        auto gen = std::mt19937_64(seed);
#pragma omp parallel for default(none) shared(pairs, blockSize) firstprivate(gen)
        for (auto const &[blockDensity, blockConnection] : pairs) {
            std::vector<GraphData> graphs;
            Generators::generateGraphs(graphPath, graphs, Properties(blockSize, blockDensity, blockSize, blockDensity,
                                                  blockConnection, "Graph"), graphNumber, gen);
            for (GraphData &graph : graphs) {
                std::vector<GraphData> trees = PatternGenerator::generatePatterns(graph, PatternType::BFS_TREE,
                                                                                  patterns, gen);
                for (GraphData &pattern : trees) {
                    pattern.save(treePath + graph.name + "_");
                }
            }
            for (GraphData &graph : graphs) {
                std::vector<GraphData> outerplanar = PatternGenerator::generatePatterns(graph, PatternType::OUTERPLANAR,
                                                                                        patterns, gen);
                for (GraphData &pattern : outerplanar) {
                    pattern.save(outerplanarPath + graph.name + "_");
                }
            }
        }
    }
    runtimes.generation = std::chrono::duration_cast<std::chrono::seconds>(
            (std::chrono::high_resolution_clock::now() - start));
}


void ExperimentalSetup::runExperimentalSetup(int seed) {

    for (int blockSize : blockSizes) {
        std::vector<std::pair<double, int>> pairs;
        for (double blockDensity : blockDensities) {
            for (int blockConnection : blockConnections) {
                pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
            }
        }
        auto gen = std::mt19937_64(seed);
#pragma omp parallel for default(none) shared(pairs, blockSize) firstprivate(gen)
        for (auto const &[blockDensity, blockConnection] : pairs) {
            for (int trainingSize : trainingSizes) {
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<GraphData> graphs;
                std::vector<std::vector<GraphData>> treePatterns;
                std::vector<std::vector<GraphData>> outerplanarPatterns;
                LoadGraphData(graphPath, graphs, LoadProperties({{"_A_",  blockSize},
                                                                 {"_dA_", (int) (blockDensity * 10)},
                                                                 {"_C_",  blockConnection}}, graphNumber));
                for (GraphData &graph : graphs) {
                    LoadGraphData(treePath, treePatterns, LoadProperties({{graph.name, -1}}, patterns));
                    LoadGraphData(outerplanarPath, outerplanarPatterns, LoadProperties({{graph.name, -1}}, patterns));
                }
#pragma omp critical
                runtimes.Loading += std::chrono::duration_cast<std::chrono::seconds>(
                        (std::chrono::high_resolution_clock::now() - start));
                start = std::chrono::high_resolution_clock::now();
                singleRun(graphs, treePatterns, outerplanarPatterns, trainingSize, gen);
#pragma omp critical
                runtimes.all += std::chrono::duration_cast<std::chrono::seconds>(
                        (std::chrono::high_resolution_clock::now() - start));
            }
        }
        std::cout << "Runtime: " << runtimes.all.count() << std::endl;
        std::cout << "\t" << "Trees: " << (double) runtimes.trees.count() / 1000000 << std::endl;
        std::cout << "\t" << "Outerplanar: " << (double) runtimes.outerplanar.count() / 1000000 << std::endl;
        std::cout << "\t" << "Graphs: " << (double) runtimes.graphs.count() / 1000000 << std::endl;
        std::cout << "Load Time: " << runtimes.Loading.count() << std::endl;
        std::cout << "Single Run Time: " << runtimes.singleRun.count() << std::endl;
    }
}

void ExperimentalSetup::singleRun(std::vector<GraphData>& graphs, std::vector<std::vector<GraphData>>& treePatterns, std::vector<std::vector<GraphData>>& outerplanarPatterns, int trainingSize, std::mt19937_64& gen) {
    GraphClosureSP gc = GraphClosureSP();
    Evaluations treeEval = Evaluations(graphs[0], graphs.size(), trainingConfigurations, "SpanningTrees");
    Evaluations outerEval = Evaluations(graphs[0], graphs.size(), trainingConfigurations, "Outerplanar Graphs");
    Evaluations graphEval = Evaluations(graphs[0], graphs.size(), trainingConfigurations, "Graphs");
    for (int n = 0; n < graphs.size(); ++n) {
        auto& trees = treePatterns[n];
        auto& outerplanars = outerplanarPatterns[n];
        auto& graph = graphs[n];
        int runs = trees.size();
        Evaluations treeEvalPerGraph = Evaluations(graph, trees.size(), trainingConfigurations, "SpanningTrees");
        Evaluations outerEvalPerGraph = Evaluations(graph, outerplanars.size(), trainingConfigurations, "Outerplanar Graphs");
        Evaluations graphEvalPerGraph = Evaluations(graph, runs, trainingConfigurations, "Graphs");
        for (int i = 0; i < trainingConfigurations; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            graph.setTrainingNodes(trainingSize, gen);
            std::vector<Labels> predictions;
            for (GraphData &tree : trees) {
                if (trainingSize > 1) {
                    for (std::set<NodeId>& nodes : graph.trainingSet) {
                        for (NodeId node : nodes) {
                            tree.trainingSet.emplace_back(std::set<NodeId>{node});
                        }
                    }
                    MCSS algorithm = MCSS(tree, gc);
                    Labels labels = algorithm.predict(tree.trainingSet, gen);
                    for (Label& label : labels) {
                        label = (int) label/trainingSize;
                    }
                    predictions.emplace_back(labels);
                }
                else{
                    tree.trainingSet = graph.trainingSet;
                    MCSS algorithm = MCSS(tree, gc);
                    Labels labels = algorithm.predict(tree.trainingSet, gen);
                    predictions.emplace_back(labels);
                }
            }
            treeEvalPerGraph.evaluate(graph, predictions);
            predictions.clear();
            //runtimes.trees += std::chrono::duration_cast<std::chrono::microseconds>((std::chrono::high_resolution_clock::now() - start));
            start = std::chrono::high_resolution_clock::now();
            for (GraphData &outerplanar : outerplanars) {
                if (trainingSize > 1) {
                    for (std::set<NodeId>& nodes : graph.trainingSet) {
                        for (NodeId node : nodes) {
                            outerplanar.trainingSet.emplace_back(std::set<NodeId>{node});
                        }
                    }
                    MCSS algorithm = MCSS(outerplanar, gc);
                    Labels labels = algorithm.predict(outerplanar.trainingSet, gen);
                    for (Label& label : labels) {
                        label = (int) label/trainingSize;
                    }
                    predictions.emplace_back(labels);
                }
                else {
                    outerplanar.trainingSet = graph.trainingSet;
                    MCSS algorithm = MCSS(outerplanar, gc);
                    Labels labels = algorithm.predict(outerplanar.trainingSet, gen);
                    predictions.emplace_back(labels);
                }
            }
            outerEvalPerGraph.evaluate(graph, predictions);
            predictions.clear();
            //runtimes.outerplanar += std::chrono::duration_cast<std::chrono::microseconds>((std::chrono::high_resolution_clock::now() - start));
//            start = std::chrono::high_resolution_clock::now();
//            for (int j = 0; j < runs; ++j) {
//                MCSS algorithm = MCSS(graph, gc);
//                Labels labels = algorithm.predict(graph.trainingSet, gen);
//                predictions.emplace_back(labels);
//            }
//            graphEvalPerGraph.evaluate(graph, predictions);
            //runtimes.graphs += std::chrono::duration_cast<std::chrono::microseconds>((std::chrono::high_resolution_clock::now() - start));
        }
#pragma omp critical
        treeEvalPerGraph.save(resultPath + resultsName + ".csv");
#pragma omp critical
        outerEvalPerGraph.save(resultPath + resultsName + ".csv");
#pragma omp critical
        graphEvalPerGraph.save(resultPath + resultsName + ".csv");
        treeEval.accuracies.insert(treeEval.accuracies.end(), treeEvalPerGraph.accuracies.begin(), treeEvalPerGraph.accuracies.end());
        outerEval.accuracies.insert(outerEval.accuracies.end(), outerEvalPerGraph.accuracies.begin(), outerEvalPerGraph.accuracies.end());
        graphEval.accuracies.insert(graphEval.accuracies.end(), graphEvalPerGraph.accuracies.begin(), graphEvalPerGraph.accuracies.end());
    }
#pragma omp critical
    treeEval.save(resultPath + resultsName + "_average.csv");
#pragma omp critical
    outerEval.save(resultPath + resultsName + "_average.csv");
#pragma omp critical
    graphEval.save(resultPath + resultsName + "_average.csv");
}


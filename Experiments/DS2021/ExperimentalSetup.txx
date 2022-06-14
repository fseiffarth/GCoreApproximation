//
// Created by anonymous on 10.05.2021.
//

#include <string>
#include <iostream>
#include "ExperimentalSetup.h"
#include "../../Data/GraphData.h"
#include "../../Utils/graphio.h"
#include "../../Utils/pattern_generator.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Evaluations/evaluations.h"
#include "../../Algorithms/MCSS.h"
#include "../../Utils/Generators.h"
#include "../../InitUpdate/InitUpdateBase.h"
#include <omp.h>

class ParallelStream{
    std::ostringstream stdStream;
public:
    ParallelStream(){}
    template <class T>
            ParallelStream& operator<<(const T& inData){
        stdStream << inData;
        return *this;
    }
    std::string toString() const{
        return stdStream.str();
    }
};

ExperimentalSetup::ExperimentalSetup(std::string resultsName, std::string resultPath, std::string graphPath, std::string treePath,
                                     std::string outerplanarPath, int graphNumber,
                                     const std::vector<int> &blockSizes, std::vector<double> blockDensities,
                                     std::vector<int> blockConnections, std::vector<int> trainingSizes,
                                     int patterns, int pattern_out_stepsize, int trainingConfigurations, bool randomConnections, int threads, bool testMode)
        : resultsName(resultsName), blockSizes(blockSizes), blockDensities(blockDensities), blockConnections(blockConnections),
          trainingSizes(trainingSizes), patternNumber(patterns), trainingConfigurations(trainingConfigurations), resultPath(resultPath),
          graphPath(graphPath), graphNumber(graphNumber), treePath(treePath), outerplanarPath(outerplanarPath), _out_stepsize(pattern_out_stepsize), _random_connections(randomConnections), _threads(threads), _testMode(testMode){
    omp_set_num_threads(threads);
}

void ExperimentalSetup::generateGraphs(int seed, GenerationType generationType) {
    auto start = std::chrono::high_resolution_clock::now();
    switch (generationType) {
        case GenerationType::TWO_COMPONENTS:
            for (int blockSize : blockSizes) {
                std::vector<std::pair<double, int>> pairs;
                for (double blockDensity : blockDensities) {
                    if (blockConnections.empty()){
                        pairs.emplace_back(std::pair<double, int>{blockDensity, blockSize * blockDensity / 10});
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
                    Generators::generateTwoComponentGraphs(graphPath, graphs,
                                                           Properties(blockSize, blockDensity, blockSize, blockDensity,
                                                                      blockConnection, "TwoComponents"), graphNumber, gen, _random_connections);
                }
            }
            break;
        case GenerationType::ONE_COMPONENT:
            for (int blockSize : blockSizes) {
                std::vector<std::pair<double, int>> pairs;
                auto gen = std::mt19937_64(seed);
#pragma omp parallel for default(none) shared(pairs, blockSize) firstprivate(gen)
                for (double blockDensity : blockDensities) {
                    std::vector<GraphData> graphs;
                    Generators::generateOneComponentGraphs(graphPath, graphs,
                                                           Properties(blockSize, blockDensity, blockSize, blockDensity,
                                                                      0, "OneComponent"), graphNumber,
                                                           gen);
                }
            }
            break;
        default:
            break;
    }
    runtimes.generation = std::chrono::duration_cast<std::chrono::seconds>(
            (std::chrono::high_resolution_clock::now() - start));
}

template <typename T>
void ExperimentalSetup::generatePatterns(T& closure, int seed, int numPatterns, GenerationType generationType, PatternType patternType) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int blockSize : blockSizes) {
        std::vector<std::pair<double, int>> pairs;
        for (double blockDensity : blockDensities) {
            switch (generationType) {
                case GenerationType::TWO_COMPONENTS:
                    if (blockConnections.empty()){
                        pairs.emplace_back(std::pair<double, int>{blockDensity, blockSize * blockDensity / 10});
                    }
                    else {
                        for (int blockConnection : blockConnections) {
                            pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
                        }
                    }
                    break;
                case GenerationType::ONE_COMPONENT:
                    pairs.emplace_back(std::pair<double, int>{blockDensity, 0});
                    break;
            }

        }
        auto gen = std::mt19937_64(seed);
        for (auto const &[blockDensity, blockConnection] : pairs) {

            //Load the graph data to generate the patternNumber
            std::vector<GraphData> graphs;
            switch (generationType) {
                case GenerationType::TWO_COMPONENTS:
                    LoadGraphData(graphPath, graphs, LoadProperties({{"TwoComponents", -1},{"_A_",  blockSize},
                                                                     {"_dA_", (int) (blockDensity * 10)},
                                                                     {"_C_",  blockConnection}}, graphNumber));
                    break;
                case GenerationType::ONE_COMPONENT:
                    LoadGraphData(graphPath, graphs, LoadProperties({{"OneComponent", -1}, {"_A_",  blockSize},
                                                                     {"_dA_", (int) (blockDensity * 10)},
                                                                     {"_C_",  0}}, graphNumber));
                    break;
            }
            //Iterate over all loaded graphs and generate patternNumber many patterns of the proposed type
            float overallGenerationTime = 0;
            for (GraphData &graph : graphs) {
                float genTime = 0;
                float perPatternTime = 0;
                auto startTime = std::chrono::high_resolution_clock::now();
                //Generate patterns of pattern type
                std::vector<GraphData> pat = PatternGenerator::generatePatterns<T>(graph, closure, patternType,
                                                                                   numPatterns, gen, false);
                genTime = (float) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000000;
                overallGenerationTime += genTime / (float) graphs.size();
                perPatternTime =  genTime/(float)pat.size();
                switch (patternType) {
                    case PatternType::BFS_TREE:
                        std::cout << "Tree generation time: " << genTime << "s" << " Per Pattern: " << perPatternTime << "s" << std::endl;
                        //Save patternNumber
                        for (GraphData &pattern : pat) {
//                            if (_testMode){
//                                std::cout << "Generated (Tree): " << "Nodes: " << pattern.graph()->GetNodes() << "Edges: " << pattern.graph()->GetEdges() << std::endl;
//                            }
                            pattern.save(treePath + graph.getName() + "_");
                        }
                        break;
                    case PatternType::OUTERPLANAR:
                        std::cout << "Outerplanar generation time: " << genTime << "s" << " Per Pattern: " << perPatternTime << "s" << std::endl;
                        for (GraphData &pattern : pat) {
//                            if (_testMode){
//                                std::cout << "Generated (Outerplanar): " << "Nodes: " << pattern.graph()->GetNodes() << "Edges: " << pattern.graph()->GetEdges() << std::endl;
//                            }
                            pattern.save(outerplanarPath + graph.getName() + "_");
                        }
                        break;
                }
            }
            savePatternGenerationRuntimes(graphs[0], blockDensity, blockConnection, overallGenerationTime, numPatterns, (int) graphs.size(), generationType, patternType);
        }
    }
    runtimes.generation = std::chrono::duration_cast<std::chrono::seconds>(
            (std::chrono::high_resolution_clock::now() - start));
}

template <typename T1, typename T2>
void ExperimentalSetup::run(T1& closure, T2& initUpdate, int seed, int max_patterns, int numGraphs, std::string& folderPath, GenerationType generationType, PatternType patternType) {
    if (folderPath.empty()){
        folderPath = createNewFolder(resultPath);
        std::ifstream  src(resultPath + "/GenerationRuntimes/"+ resultsName + "_GenerationRuntimes.csv", std::ios::binary);
        std::ofstream  dst(resultPath + "/"+ resultsName + "_GenerationRuntimes.csv",   std::ios::binary);
        dst << src.rdbuf();
    }
    for (int blockSize : blockSizes) {
        std::vector<std::pair<double, int>> pairs;
        switch (generationType) {
            case GenerationType::TWO_COMPONENTS:
                for (double blockDensity : blockDensities) {
                    if (blockConnections.empty()){
                        pairs.emplace_back(std::pair<double, int>{blockDensity, blockSize * blockDensity / 10});
                    }
                    else {
                        for (int blockConnection : blockConnections) {
                            pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
                        }
                    }
                }
                break;
            case GenerationType::ONE_COMPONENT:
                for (double blockDensity : blockDensities) {
                    pairs.emplace_back(std::pair<double, int>{blockDensity, 0});
                }
                break;
        }
        auto gen = std::mt19937_64(seed);
        size_t loading = 0;
        size_t runtime = 0;
        size_t time = 0;
        std::vector<Evaluations> allAverageEvaluations;


//#pragma omp parallel for default(none) shared(allAverageEvaluations, pairs, blockSize) firstprivate(closure, initUpdate, gen, generationType, patternType, max_patterns, numGraphs) reduction(+: runtime, loading, time)
        //Iterate over different parameters of the blocks
        for (auto & [blockDensity, blockConnection] : pairs) {
            std::cout << "Size: " << blockSize*2 << " Density: " << blockDensity << " Connections: " << blockConnection << std::endl;
            std::vector<FileEvaluation> averageEvaluations;
            auto start = std::chrono::high_resolution_clock::now();
            //Load all graphs and patterns
            std::vector<GraphData> graphs;
            std::vector<std::vector<GraphData>> patterns;
            LoadGraphsAndPatternsForExperiment(graphs, patterns,
                                               generationType, blockSize,
                                               blockDensity, blockConnection, numGraphs,
                                               patternType, max_patterns);
            loading += std::chrono::duration_cast<std::chrono::seconds>(
                    (std::chrono::high_resolution_clock::now() - start)).count();

            for (int trainingSize : trainingSizes) {
                std::cout << "\t" << "Training Size: " << trainingSize << std::endl;
                start = std::chrono::high_resolution_clock::now();
                trainingConfigurationRun(closure, initUpdate, averageEvaluations, graphs, patterns, max_patterns, time, trainingSize, blockDensity, blockConnection, patternType, gen);
                runtime += std::chrono::duration_cast<std::chrono::seconds>(
                        (std::chrono::high_resolution_clock::now() - start)).count();
            }
//#pragma omp critical
            {
                for (auto const& evaluation : averageEvaluations) {
                    allAverageEvaluations.template emplace_back(evaluation);
                }
            }
        }
        runtimes.loading = loading;
        runtimes.all = runtime;
        runtimes.algorithm = time;
        for (auto& evaluation : allAverageEvaluations) {
            evaluation.generationType = generationType;
            evaluation.save(folderPath + resultsName + ".csv");
        }
        std::cout << "Runtime: " << runtimes.all << std::endl;
        std::cout << "\t" << "Algorithm: " << (double) runtimes.algorithm / 1000000 << std::endl;
        std::cout << "Load Time: " << runtimes.loading << std::endl;
        std::cout << "Single Run Time: " << runtimes.singleRun.count() << std::endl;
    }
}

template<typename T1, typename T2>
void ExperimentalSetup::run(T1 &closure, T2 &initUpdate, int seed, int numPatterns, int numGraphs, bool withClosureOfTraining) {
    if (!withClosureOfTraining){
        run(closure, initUpdate, seed);
    }
    else{
        for (int blockSize : blockSizes) {
            std::vector<std::pair<double, int>> pairs;
            for (double blockDensity : blockDensities) {
                if (blockConnections.empty()){
                    pairs.emplace_back(std::pair<double, int>{blockDensity, blockSize / 5});
                }
                else {
                    for (int blockConnection : blockConnections) {
                        pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
                    }
                }
            }
            auto gen = std::mt19937_64(seed);
            std::vector<Evaluations> allSingleEvaluations;
            std::vector<Evaluations> allAverageEvaluations;
            size_t loading = 0;
            size_t runtime = 0;
            size_t time = 0;
            std::vector<Evaluations> singleEvaluations;
            std::vector<Evaluations> averageEvaluations;
#pragma omp parallel for default(none) shared(allSingleEvaluations, allAverageEvaluations, pairs, blockSize) private(singleEvaluations, averageEvaluations) firstprivate(closure, initUpdate, gen) reduction(+: runtime, loading, time)
                for (auto & [blockDensity, blockConnection] : pairs) {
                    auto start = std::chrono::high_resolution_clock::now();
                    std::vector<GraphData> graphs;
                    std::vector<std::vector<GraphData>> treePatterns;
                    std::vector<std::vector<GraphData>> outerplanarPatterns;
                    LoadGraphData(graphPath, graphs, LoadProperties({{"_A_",  blockSize},
                                                                     {"_dA_", (int) (blockDensity * 10)},
                                                                     {"_C_",  blockConnection}}, numGraphs));
                    loading += std::chrono::duration_cast<std::chrono::seconds>(
                            (std::chrono::high_resolution_clock::now() - start)).count();
                    for (int trainingSize : trainingSizes) {
                        int graphSize = static_cast<int>(graphs.size());
                        Evaluations treeEval = Evaluations(graphs[0], closure.getName(), initUpdate.getName(),
                                                         graphSize, trainingConfigurations,
                                                         "SpanningTrees");
                        treeEval.trainingSize = trainingSize;
                        Evaluations outerEval = Evaluations(graphs[0], closure.getName(), initUpdate.getName(),
                                                          graphSize, trainingConfigurations,
                                                          "OuterPlanarSubgraph");
                        outerEval.trainingSize = trainingSize;
                        Evaluations graphEval = Evaluations(graphs[0], closure.getName(), initUpdate.getName(),
                                                          graphSize, trainingConfigurations,
                                                          "Graphs");
                        graphEval.trainingSize = trainingSize;
                        for (auto &graph: graphs) {
                            Evaluations treeEvalPerGraph = Evaluations(graph, closure.getName(), initUpdate.getName(),
                                                                     numPatterns,
                                                                     trainingConfigurations,
                                                                     "SpanningTrees");
                            Evaluations outerEvalPerGraph = Evaluations(graph, closure.getName(), initUpdate.getName(),
                                                                      numPatterns,
                                                                      trainingConfigurations,
                                                                      "OuterPlanarSubgraph");
                            Evaluations graphEvalPerGraph = Evaluations(graph, closure.getName(), initUpdate.getName(),
                                                                      numPatterns, trainingConfigurations,
                                                                      "Graphs");
                            for (int i = 0; i < trainingConfigurations; ++i) {
                                std::vector<GraphData> patterns;

                                start = std::chrono::high_resolution_clock::now();
                                graph.setTrainingNodes(trainingSize, gen);
                                patterns = PatternGenerator::generatePatterns(graph, closure, PatternType::BFS_TREE, numPatterns, gen,
                                                                           true);
                                patterns = PatternGenerator::generatePatterns(graph, closure, PatternType::OUTERPLANAR,
                                                                                  numPatterns, gen, true);
                                loading += std::chrono::duration_cast<std::chrono::seconds>(
                                        (std::chrono::high_resolution_clock::now() - start)).count();
                                start = std::chrono::high_resolution_clock::now();
                                singleRun(graph, closure, initUpdate, patterns, trainingSize, time, treeEvalPerGraph, outerEvalPerGraph, graphEvalPerGraph, gen);
                                runtime += std::chrono::duration_cast<std::chrono::seconds>(
                                        (std::chrono::high_resolution_clock::now() - start)).count();
                            }
                            singleEvaluations.template emplace_back(treeEvalPerGraph);
                            singleEvaluations.template emplace_back(outerEvalPerGraph);
                            singleEvaluations.template emplace_back(graphEvalPerGraph);
                            treeEval.accuracies.insert(treeEval.accuracies.end(), treeEvalPerGraph.accuracies.begin(),
                                                       treeEvalPerGraph.accuracies.end());
                            outerEval.accuracies.insert(outerEval.accuracies.end(),
                                                        outerEvalPerGraph.accuracies.begin(),
                                                        outerEvalPerGraph.accuracies.end());
                            graphEval.accuracies.insert(graphEval.accuracies.end(),
                                                        graphEvalPerGraph.accuracies.begin(),
                                                        graphEvalPerGraph.accuracies.end());
                        }
                        averageEvaluations.template emplace_back(treeEval);
                        averageEvaluations.template emplace_back(outerEval);
                        averageEvaluations.template emplace_back(graphEval);
                    }
#pragma omp critical
                    {
                        for (auto const &evaluation : singleEvaluations) {
                            allSingleEvaluations.template emplace_back(evaluation);
                        }
                        for (auto const &evaluation : averageEvaluations) {
                            allAverageEvaluations.template emplace_back(evaluation);
                        }
                    }
                }
            runtimes.loading = loading;
            runtimes.all = runtime;
            runtimes.algorithm = time;
            for (auto& evaluation : allSingleEvaluations) {
                evaluation.save(resultPath + resultsName + ".csv");
            }
            for (auto& evaluation : allAverageEvaluations) {
                evaluation.save(resultPath + resultsName + "_average.csv");
            }
            std::cout << "Runtime: " << runtimes.all << std::endl;
            std::cout << "\t" << "Algorithm: " << (double) runtimes.algorithm / 1000000 << std::endl;
            std::cout << "Load Time: " << runtimes.loading << std::endl;
            std::cout << "Single Run Time: " << runtimes.singleRun.count() << std::endl;
        }
    }
}

template<typename T1, typename T2>
void ExperimentalSetup::trainingConfigurationRun(T1& closure, T2& initUpdate, std::vector<FileEvaluation>& averageEvaluations, std::vector<GraphData>& graphs, std::vector<std::vector<GraphData>>& patterns, int maxPatterns, size_t& time, int trainingSize, double density, int connections, PatternType type, std::mt19937_64& gen) {
    int numGraphs = static_cast<int>(graphs.size());

    std::string generationText;
    if (type == PatternType::BFS_TREE){
        generationText = "SpanningTrees";
    }
    else if (type == PatternType::OUTERPLANAR){
        generationText = "OuterPlanarSubgraph";
    }

    FileEvaluation eval = FileEvaluation(graphs[0], closure.getName(), initUpdate.getName(), density, connections, numGraphs,
                                 trainingConfigurations,
                                 generationText);

    for (int n = 0; n < numGraphs; ++n) {
        std::vector<GraphData>* pat = nullptr;
        if (patterns.size() > n) {
            pat = &patterns[n];
        }

        //Set the training sets, because graphs are randomly drawn we can choose fixed training sets to make experiments comparable
        for (int i = 0; i < trainingConfigurations; ++i) {
            std::set<NodeId> trainingA;
            std::set<NodeId> trainingB;
            for (int j = i*trainingSize; j < (i+1)*trainingSize; ++j) {
                trainingA.insert(j);
                trainingB.insert(j + (int) graphs[n].size()/2);
            }
            graphs[n].trainingSet = {trainingA, trainingB};
            //graph.setTrainingNodes(trainingSize, _gen);
            if (pat != nullptr){
                singleRun(graphs[n], closure, initUpdate, pat, trainingSize, time, eval, gen);
            }
        }
        averageEvaluations.template emplace_back(eval);
    }
}

template<typename T1, typename T2>
void ExperimentalSetup::singleRun(GraphData &graph, T1 &closure, T2 &initUpdate, std::vector<GraphData> *patterns,int trainingSize, size_t &patternTime, Evaluations &patternEvaluation,std::mt19937_64 &gen) const {

    std::vector<Labels> predictions;
    std::vector<double> algorithm_runtimes;
    auto start = std::chrono::high_resolution_clock::now();
    for (GraphData &pattern : *patterns) {
        if (trainingSize > 1) {
            pattern.trainingSet.clear();
            for (std::set<NodeId> &nodes : graph.trainingSet) {
                for (NodeId node : nodes) {
                    pattern.trainingSet.emplace_back(std::set<NodeId>{node});
                }
            }
            MCSS algorithm = MCSS<GraphData, T1, T2>(pattern, closure, initUpdate);
            Labels labels = algorithm.predict(pattern.trainingSet, gen);
            if (_testMode){
                //std::cout << "Nodes: " << pattern.graph()->GetNodes() << " Edges: " << pattern.graph()->GetEdges() << std::endl;
                std::vector<std::set<NodeId>> closed_sets;
                GraphClosureSP gc = GraphClosureSP("Test");
                std::vector<std::set<NodeId>> nodes;
                int labelSize = 2*trainingSize;
                for (int i = 0; i < labelSize; ++i) {
                    nodes.template emplace_back(StaticFunctions::getNodesWithLabel(labels, i));
                }
                for (int i = 0; i < labelSize; ++i) {
                    std::vector<std::set<NodeId>*> forbidden_elements;
                    for (int j = 0; j < labelSize; ++j) {
                        if (j!=i){
                            forbidden_elements.template emplace_back(&nodes[j]);
                        }
                    }
                    ClosureParameters closureParameters = ClosureParameters(nodes[i], forbidden_elements);
                    gc.closure(pattern, closureParameters);
                    if (closureParameters.input_set != closureParameters.closed_set || closureParameters.output_intersects_forbidden){
                        if (pattern.graphType == GraphType::TREE){
                            std::cout << "Tree: ";
                        }
                        else{
                            std::cout << "Outerplanar Pattern";
                        }
                        std::cout << "Closure Test fails!" << std::endl;
                    }
                }
            }
            for (Label &label : labels) {
                label = (int) (label / trainingSize);
            }
            predictions.emplace_back(labels);
        } else {
            pattern.trainingSet = graph.trainingSet;
            MCSS algorithm = MCSS<GraphData, T1, T2>(pattern, closure, initUpdate);
            Labels labels = algorithm.predict(pattern.trainingSet, gen);
            predictions.emplace_back(labels);
            if (_testMode){
                std::vector<std::set<NodeId>> closed_sets;
                GraphClosureSP gc = GraphClosureSP("Test");
                std::vector<std::set<NodeId>> nodes;
                int labelSize = 2;
                for (int i = 0; i < labelSize; ++i) {
                    nodes.template emplace_back(StaticFunctions::getNodesWithLabel(labels, i));
                }
                for (int i = 0; i < labelSize; ++i) {
                    std::vector<std::set<NodeId>*> forbidden_elements;
                    for (int j = 0; j < labelSize; ++j) {
                        if (j!=i){
                            forbidden_elements.template emplace_back(&nodes[j]);
                        }
                    }
                    ClosureParameters closureParameters = ClosureParameters(nodes[i], forbidden_elements);
                    gc.closure(pattern, closureParameters);
                    if (closureParameters.input_set != closureParameters.closed_set || closureParameters.output_intersects_forbidden){
                        if (pattern.graphType == GraphType::TREE){
                            std::cout << "Tree: ";
                        }
                        else{
                            std::cout << "Outerplanar Pattern";
                        }
                        std::cout << "Closure Test fails!" << std::endl;
                    }
                }
            }
        }
        algorithm_runtimes.template emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(
                (std::chrono::high_resolution_clock::now() - start)).count());

        if (_testMode && std::find(predictions.back().begin(), predictions.back().end(), -1) != predictions.back().end()){
            std::cout << "Algorithm is invalid unclassified nodes" << std::endl;
            if (pattern.graphType == GraphType::TREE){
                std::cout << "\tTree" << std::endl;;
            }
            else{
                std::cout << "\tOuterplanar" << std::endl;
            }
            std::cout << "\t" << pattern.getName() << std::endl;
            std::string t_set;
            int c = 0;
            for (const auto & set : pattern.trainingSet) {
                t_set.append("\t");
                for (const auto nodeId : set) {
                    t_set.append(std::to_string(nodeId));
                }
            }
            std::cout << std::endl;
            std::cout << "\tTraining Set: " << t_set << std::endl;
            std::cout << "\tPrediction: " << StaticFunctions::print<std::vector<int>, int>(predictions.back()) << std::endl;
        }
    }


    patternEvaluation.evaluate(graph, predictions, _out_stepsize, algorithm_runtimes);
    predictions.clear();
    patternTime += std::chrono::duration_cast<std::chrono::microseconds>(
            (std::chrono::high_resolution_clock::now() - start)).count();


    //    if (outerplanars != nullptr) {
    //        //Start Outerplanar run
    //        algorithm_runtimes.clear();
    //        predictions.clear();
    //        start = std::chrono::high_resolution_clock::now();
    //        for (GraphData &outerplanar : *outerplanars) {
    //            if (trainingSize > 1) {
    //                outerplanar.trainingSet.clear();
    //                for (std::set<NodeId> &nodes : graph.trainingSet) {
    //                    for (NodeId get_node : nodes) {
    //                        outerplanar.trainingSet.emplace_back(std::set<NodeId>{get_node});
    //                    }
    //                }
    //                MCSS algorithm = MCSS<GraphData, T1, T2>(outerplanar, closure, initUpdate);
    //                Labels labels = algorithm.predict(outerplanar.trainingSet, gen);
    //                for (Label &label : labels) {
    //                    label = (int) label / trainingSize;
    //                }
    //                predictions.emplace_back(labels);
    //            } else {
    //                outerplanar.trainingSet = graph.trainingSet;
    //                MCSS algorithm = MCSS<GraphData, T1, T2>(outerplanar, closure, initUpdate);
    //                Labels labels = algorithm.predict(outerplanar.trainingSet, gen);
    //                predictions.emplace_back(labels);
    //            }
    //            algorithm_runtimes.template emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(
    //                    (std::chrono::high_resolution_clock::now() - start)).count());
    //        }
    //        outerplanarEval.evaluate(graph, predictions, _out_stepsize, algorithm_runtimes);
    //        predictions.clear();
    //        outerplanarTime += std::chrono::duration_cast<std::chrono::microseconds>(
    //                (std::chrono::high_resolution_clock::now() - start)).count();
    //    }
    //    start = std::chrono::high_resolution_clock::now();
    ////            for (int j = 0; j < runs; ++j) {
    ////                MCSS algorithm = MCSS(graph, gc);
    ////                Labels labels = algorithm.predict(graph.trainingSet, _gen);
    ////                predictions.emplace_back(labels);
    ////            }
    ////            graphEval.calculate_accuracies(graph, predictions, (double) std::chrono::duration_cast<std::chrono::microseconds>((std::chrono::high_resolution_clock::now() - start)).count()/1000000);
    //    graphsTime += std::chrono::duration_cast<std::chrono::microseconds>((std::chrono::high_resolution_clock::now() - start)).count();
}

void ExperimentalSetup::analyzeFaces(GenerationType generationType) {
    for (int blockSize : blockSizes) {
        std::vector<std::pair<double, int>> pairs;
        switch (generationType) {
            case GenerationType::TWO_COMPONENTS:
                for (double blockDensity : blockDensities) {
                    for (int blockConnection : blockConnections) {
                        pairs.emplace_back(std::pair<double, int>{blockDensity, blockConnection});
                    }
                }
                break;
            case GenerationType::ONE_COMPONENT:
                for (double blockDensity : blockDensities) {
                    pairs.emplace_back(std::pair<double, int>{blockDensity, 0});
                }
                break;
        }
        for (auto &[blockDensity, blockConnection] : pairs) {
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<GraphData> graphs;
            std::vector<std::vector<GraphData>> outerplanarPatterns;
            switch (generationType) {
                case GenerationType::TWO_COMPONENTS:
                    LoadGraphData(graphPath, graphs, LoadProperties({{"TwoComponents", -1},
                                                                     {"_A_",           blockSize},
                                                                     {"_dA_",          (int) (blockDensity * 10)},
                                                                     {"_C_",           blockConnection}}, graphNumber));
                    break;
                case GenerationType::ONE_COMPONENT:
                    LoadGraphData(graphPath, graphs, LoadProperties({{"OneComponent", -1},
                                                                     {"_A_",          blockSize},
                                                                     {"_dA_",         (int) (blockDensity * 10)},
                                                                     {"_C_",          0}}, graphNumber));
                    break;
            }
            for (GraphData &graph : graphs) {
                LoadGraphData(outerplanarPath, outerplanarPatterns,
                              LoadProperties({{graph.getName(), -1}}, patternNumber));
                int face_sum = 0;
                int face_size_sum = 0;
                for (auto &pattern : outerplanarPatterns.back()) {
                    if (GraphFunctions::IsOuterPlanar(pattern.get_data())){
                        std::cout << "Outerplanar! " << pattern.getName() << std::endl;
                    }
                    else{
                        std::cout << "Not Outerplanar! " << pattern.getName() << std::endl;
                    }
                    pattern.graphType = GraphType::OUTERPLANAR;
                    int size = (int) pattern.size();
                    int overall_faces = 0;
                    std::vector<std::pair<int, int>> size_face_num_pair;
                    std::vector<PUNGraph> components;
                    GraphFunctions::GetBiconnectedComponents(pattern.get_data(), components);
                    for (auto& component : components) {
                        int compSize = component->GetNodes();
                        if (!GraphFunctions::IsTree(component)) {
                            int face_num;
                            GraphFunctions::GetBiconnectedOuterplanarFaceNum(component, face_num);
                            overall_faces += face_num;
                            face_sum += face_num;
                            face_size_sum += compSize*face_num;
                        }
                        else{
                            face_size_sum += compSize*1;
                        }
                    }
                    std::cout << "Size: " << size << std::endl;
                    StaticFunctions::print<std::vector<std::pair<int, int>>, std::pair<int, int>>(size_face_num_pair, true);
                    std::cout << "FaceNums: " << overall_faces << std::endl;
                }
                std::cout << "Size: " << graph.size() << std::endl;
                std::cout << "Average FaceNums: " << (double) face_sum/outerplanarPatterns.back().size() << std::endl;
                std::cout << "Average Sum Face Size: " << (double) face_size_sum/outerplanarPatterns.back().size() << std::endl;
                std::cout << "Average Sum Face Size/Size: " << (double) (face_size_sum/outerplanarPatterns.back().size())/graph.size() << std::endl;
                std::cout << "Average Sum Face Size/Size^2: " << (double) (face_size_sum/outerplanarPatterns.back().size())/(graph.size()*graph.size()) << std::endl;
            }
        }
    }
}

void ExperimentalSetup::savePatternGenerationRuntimes(const GraphData& graphData, double density, int connection, float perPatternTime,
                                                      int patternNum, int graphNum, GenerationType generationType, PatternType patternType) {
    std::string path = resultPath + "/GenerationRuntimes/"+ resultsName + "_GenerationRuntimes.csv";
    bool newFile = std::filesystem::exists(path);
    std::ofstream fs;
    fs.open(path);
    // write the file headers
    if (!newFile) {
        fs << ",Graphs" << "," << "GraphSize" << "," << "GraphEdgeNum" << "," << "BlockDensity" << "," << "BlockConnections" << "," << "PatternType" << "," << "MeanRuntimePerPattern" << "," << "NumGraphs" << ","
           << "Iterations" << "," << "Generation" << std::endl;
    }
    std::string EstimationString;
    std::string generationTypeString;
    switch (generationType) {
        case GenerationType::ONE_COMPONENT:
            generationTypeString = "OneComponent";
            break;
        case GenerationType::TWO_COMPONENTS:
            generationTypeString = "TwoComponents";
            break;
    }
    std::string patternTypeString;
    switch (patternType) {
        case PatternType::BFS_TREE:
            patternTypeString = "SpanningTrees";
            break;
        case PatternType::OUTERPLANAR:
            patternTypeString = "OuterPlanarSubgraph";
            break;
    }

    fs << std::fixed << "," << graphData.getName() << "," << graphData.size() << "," << graphData.get_data()->GetEdges() << "," << density << "," << connection << "," << patternTypeString << "," << perPatternTime / (float) patternNum << "," << graphNum << ","
       << patternNum << "," << generationTypeString << std::endl;

    fs << std::scientific;
    fs.close();

}

std::string ExperimentalSetup::createNewFolder(const std::string& pathString) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y_%m_%d_%H_%M_%S");
    auto str = oss.str();
    std::filesystem::create_directories(pathString + "/" + str);
    return pathString + "/" + str + "/";
}

void
ExperimentalSetup::LoadGraphsAndPatternsForExperiment(std::vector<GraphData> &graphs,
                                                      std::vector<std::vector<GraphData>> &patterns,
                                                      GenerationType generationType, int blockSize,
                                                      double blockDensity, int blockConnection, int numGraphs,
                                                      PatternType patternType, int max_patterns) {

    switch (generationType) {
        case GenerationType::TWO_COMPONENTS:
            LoadGraphData(graphPath, graphs, LoadProperties({{"TwoComponents", -1}, {"_A_",  blockSize},
                                                             {"_dA_", (int) (blockDensity * 10)},
                                                             {"_C_",  blockConnection}}, numGraphs));
            break;
            case GenerationType::ONE_COMPONENT:
                LoadGraphData(graphPath, graphs, LoadProperties({{"OneComponent", -1}, {"_A_",  blockSize},
                                                                 {"_dA_", (int) (blockDensity * 10)},
                                                                 {"_C_",  0}}, numGraphs));
                break;
    }


    for (GraphData &graph : graphs) {
        if (patternType == PatternType::BFS_TREE) {
            //Load trees
            LoadGraphData(treePath, patterns, LoadProperties({{graph.getName(), -1}}, max_patterns));
            for (auto &pattern : patterns.back()) {
//                if (_testMode){
//                    std::cout << "Tree: " "Nodes: " << pattern.graph()->GetNodes() << " Edges: "  << pattern.graph()->GetEdges() << std::endl;
//                }
                pattern.graphType = GraphType::TREE;
            }
        }
        if (patternType == PatternType::OUTERPLANAR) {
            //Load outerplanar graphs
            LoadGraphData(outerplanarPath, patterns,
                          LoadProperties({{graph.getName(), -1}}, max_patterns));
            for (auto &pattern : patterns.back()) {
//                if (_testMode){
//                    std::cout << "Outerplanar: " "Nodes: " << pattern.graph()->GetNodes() << " Edges: "  << pattern.graph()->GetEdges() << std::endl;
//                }
                pattern.graphType = GraphType::OUTERPLANAR;
            }
        }
    }
}




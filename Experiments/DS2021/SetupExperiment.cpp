//
// Created by anonymous on 25.08.21.
//

#include <set>
#include "SetupExperiment.h"
#include "../../Utils/StaticFunctions.h"

#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

void SetupExperiment::samplePreprocessing(GraphData& graph, std::vector<GraphData>& treeSamples, std::pair<double, double>& best_threshold_tree,
                                          std::vector<GraphData>& outerplanarSamples, std::pair<double, double>& best_threshold_outerplanar,
                                          int& simple_approx_iteration_number, std::vector<int>& variation_check,
                                          const std::vector<int> &generator_sizes, int num_training_sets, double training_elements, int seed,
                                          int thresholds_max_number, double lower_bound, int max_iterations, void *generator_elements_choosing,
                                          bool completeSampling) {
    this->t_max_number = thresholds_max_number;
    GraphClosureSP gc = GraphClosureSP("Closure");
    ClosureParameters closureParameters = ClosureParameters();
    OuterplanarSubgraphDFS outerplanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
    std::vector<NodeId> neighborIds;
    GraphFunctions::generateNeighborVector(graph.get_graph(), neighborIds);
    int size = graph.get_graph()->GetNodes();
    std::vector<int> approximationVector = std::vector<int>(size, 0);
    std::vector<bool> realClosure = std::vector<bool>(size, false);
    variation_check = std::vector<int>(size, 0);
    //CorrectClosure
    std::vector<std::string> header = {"Name", "Size", "Edges", "Density", "InputSize", "ClosedSetSize", "Type",
                                       "Jaccard", "TP",
                                       "TN", "FP", "FN", "Time", "Threads", "Closure Time", "Sampling Time",
                                       "Iterations", "Threshold", "TP weights (max/mean/min)",
                                       "TN weights (max/mean/min)", "FP weights (max/mean/min)",
                                       "FN weights (max/mean/min)", "Approximation Frequencies"};

    std::vector<std::pair<double, double>> best_thresholds_trees;
    std::vector<std::pair<double, double>> best_thresholds_outerplanar;

    for (int i = 0; i < num_training_sets; ++i) {
        //Get generator set
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_sizes[0], seed);
        if (!completeSampling) {
            CorrectClosure(gc, graph, closureParameters, 0, header, realClosure, true);
            for (auto elem: closureParameters.closed_set) {
                variation_check[elem] = 1;
            }
            int closureSize = (int) closureParameters.closed_set.size();
            //Set Iterations
            SetApproximationIterations(gc, graph, best_threshold_outerplanar, outerplanarSamples, neighborIds,
                                       &outerplanarSubgraphDfs, "Outerplanar", PatternType::OUTERPLANAR,
                                       closureParameters,
                                       realClosure, approximationVector, max_iterations, lower_bound);
            SetApproximationIterations(gc, graph, best_threshold_tree, treeSamples, neighborIds,
                                       nullptr, "Tree", PatternType::BFS_TREE, closureParameters,
                                       realClosure, approximationVector, max_iterations, lower_bound);
            SetSimpleApproximationIterations(gc, graph, simple_approx_iteration_number, closureParameters, closureSize,
                                             lower_bound);
        }
        else {
            OuterplanarSubgraphDFS outerPlanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
            double runtime = 0;
            GraphFunctions::GetSamples(graph, PatternType::BFS_TREE, treeSamples, nullptr, neighborIds, max_iterations, _seed, runtime);
            GraphFunctions::GetSamples(graph, PatternType::OUTERPLANAR, outerplanarSamples, &outerplanarSubgraphDfs, neighborIds, max_iterations, _seed, runtime);
        }
    }
}
void SetupExperiment::runApproximation(GraphData& graph, std::vector<GraphData>& treeSamples, const std::pair<double, double>& best_threshold_tree, std::vector<GraphData>& outerplanarSamples, const std::pair<double, double>& best_threshold_outerplanar, int simple_approx_iterations, std::vector<int>& variation_check,  const std::vector<int>& generator_sizes, int seed) {
    GraphClosureSP gc = GraphClosureSP("Closure");
    ClosureParameters closureParameters = ClosureParameters();
    OuterplanarSubgraphDFS outerplanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
    std::vector<NodeId> neighborIds;
    GraphFunctions::generateNeighborVector(graph.get_graph(), neighborIds);
    int size = graph.get_graph()->GetNodes();
    std::vector<int> approximationVector = std::vector<int>(size, 0);
    std::vector<bool> realClosure = std::vector<bool>(size, false);
    //Get generator set

    StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_sizes[0], seed);


    //CorrectClosure
    std::vector<std::string> header = {"Name", "Size", "Edges", "Density", "InputSize", "ClosedSetSize", "Type", "Jaccard", "TP",
                                       "TN", "FP", "FN", "Time", "Threads", "Closure Time", "Sampling Time", "Iterations", "Threshold", "TP weights (max/mean/min)", "TN weights (max/mean/min)", "FP weights (max/mean/min)", "FN weights (max/mean/min)",  "Approximation Frequencies"};

    CorrectClosure(gc, graph, closureParameters, 0, header, realClosure, true);
    //Sampling Approximations
    double runtime;
    Triples triples;
    SamplingApproximation(gc, graph, outerplanarSamples, neighborIds, &outerplanarSubgraphDfs, "Outerplanar Approximation", PatternType::OUTERPLANAR,
                          closureParameters, 0, {}, realClosure, approximationVector, runtime, triples, {best_threshold_outerplanar.first}, true);
    SamplingApproximation(gc, graph, treeSamples, neighborIds, nullptr, "Tree Approximation", PatternType::BFS_TREE,
                          closureParameters, 0, {}, realClosure, approximationVector, runtime, triples, {best_threshold_tree.first}, true);
    SimpleApproximation(gc, graph, closureParameters, 0, {}, realClosure, approximationVector, triples, false, simple_approx_iterations, -1, true);

    //Check variations
    VariationCheck(graph, closureParameters, realClosure, variation_check, triples);
}


void SetupExperiment::runApproximation(const std::vector<PatternType>& patternTypes,const std::vector<double>& thresholds, const std::vector<std::string>& graph_paths, int seed) {
    std::vector<GraphData> graphs;
    std::set<NodeId> inputSet;
    GraphClosureSP gc = GraphClosureSP("Closure");
    if (graph_paths.empty()) {
        for (int size: _graphSizes) {
            for (double density: _densities) {
                graphs.clear();
                LoadGraphData(_path + "/Graphs/", graphs, LoadProperties({{"Graph",     -1},
                                                                          {"_Size_",    size},
                                                                          {"_Density_", (int) (density * 10)}},
                                                                         _graphNumber));
                run(gc, graphs, inputSet, density, thresholds, patternTypes);
            }
        }
    } else {
        for (auto const &graph_path: graph_paths) {
            graphs.clear();
            graphs.emplace_back(GraphData(graph_path));

            GraphFunctions::GetLargestComponent(graphs.back());
            run(gc, graphs, inputSet, 0, thresholds, patternTypes);
        }
    }
}



void SetupExperiment::generateGraphs(std::vector<GraphData> &graphs)  {
    std::string name;
    for (int size : _graphSizes) {
        for (double density : _densities) {
            //TSnap::GenRndGnm<PUNGraph>(size, size*)
            Generators::generateConnectedRandomGraphs(_path + "/Graphs/", graphs, Properties(size, density, name), _graphNumber, _gen, false);
            graphs.back().save_dot(_path + "/Graphs/");
            //Generators::generateConnectedRandomOuterplanarGraphs(_path + "/Graphs/", graphs, Properties(size, density, name), _graphNumber, _gen);
        }
    }
}

void SetupExperiment::run(GraphClosureSP& gc, std::vector<GraphData> &graphs, std::set<NodeId>& inputSet, double density,const std::vector<double>& thresholds, const std::vector<PatternType>& patternTypes) {
    std::vector<std::string> header = {"Name", "Size", "Edges", "Density", "InputSize", "ClosedSetSize", "Type", "Jaccard", "TP",
                                       "TN", "FP", "FN", "Time", "Threads", "Closure Time", "Sampling Time", "Iterations", "Threshold", "TP weights (max/mean/min)", "TN weights (max/mean/min)", "FP weights (max/mean/min)", "FN weights (max/mean/min)",  "Approximation Frequencies"};

    for (auto & graph : graphs) {
        std::vector<int> approximation(graph.size(), 0);
        std::vector<int> private_approximations(graph.size(), 0);
        Triples triples;
        ClosureParameters closureParameters = ClosureParameters(inputSet);
        std::vector<bool> realClosure(graph.size(), false);
        OuterplanarSubgraphDFS outerPlanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
        std::vector<NodeId> neighborIds;
        GraphFunctions::generateNeighborVector(graph.get_graph(), neighborIds);
        std::vector<GraphData> subgraphData = {GraphData(new TUNGraph(), (int) graph.size())};

        int inputSeed = 36543843;
        for (int inputSize: _generatorElements) {
            double runtime = 0;
            OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
            //Inner Loop
            StaticFunctions::generateInputSet(closureParameters.input_set, graph, inputSize, inputSeed);

            //CorrectClosure
            CorrectClosure(gc, graph, closureParameters, density, header, realClosure);
            //Approximation by PreClosure
            SimpleApproximation(gc, graph, closureParameters, density, {}, realClosure, approximation, triples, true);
            //Sampling Approximations
            SamplingApproximation(gc, graph, subgraphData, neighborIds, nullptr, "Tree Approximation", PatternType::BFS_TREE,
                                  closureParameters, density, {}, realClosure, approximation, runtime, triples, thresholds);
            SamplingApproximation(gc, graph, subgraphData,neighborIds, &outerPlanarSubgraphDfs, "Outerplanar Approximation", PatternType::OUTERPLANAR,
                                  closureParameters, density, {}, realClosure, approximation, runtime, triples, thresholds);
            //Approximation by Incremental BFS
            SimpleApproximation(gc, graph, closureParameters, density, {}, realClosure, approximation, triples, false, -1, runtime);
            ++inputSeed;
        }
    }
}


void SetupExperiment::evaluateApproximation(std::vector<int>& approximation, const std::vector<bool>& realClosure, int iterations, double threshold, double& similarity, int& tp, int& tn, int& fp, int& fn, Triples& triples) {
    tp = 0; fp = 0; tn = 0; fn = 0;
    triples.tp_triple = ValTriple{0, 0, std::numeric_limits<int>::max()};
    triples.fp_triple = ValTriple{0, 0, std::numeric_limits<int>::max()};
    triples.tn_triple = ValTriple{0, 0, std::numeric_limits<int>::max()};
    triples.fn_triple = ValTriple{0, 0, std::numeric_limits<int>::max()};
    for (int i = 0; i < approximation.size(); ++i) {
        bool closed_approximated = ((((double) approximation[i] / (double) iterations) / threshold) >= 1);
        bool closed = realClosure[i];
        if (closed_approximated && closed) {
            triples.tp_triple.mean += approximation[i];
            triples.tp_triple.max = std::max(triples.tp_triple.max, approximation[i]);
            triples.tp_triple.min = std::min(triples.tp_triple.min, approximation[i]);
            tp += 1;
        } else if (closed_approximated) {
            triples.fp_triple.mean += approximation[i];
            triples.fp_triple.max = std::max(triples.fp_triple.max, approximation[i]);
            triples.fp_triple.min = std::min(triples.fp_triple.min, approximation[i]);
            fp += 1;
        } else if (closed) {
            triples.fn_triple.mean += approximation[i];
            triples.fn_triple.max = std::max(triples.fn_triple.max, approximation[i]);
            triples.fn_triple.min = std::min(triples.fn_triple.min, approximation[i]);
            fn += 1;
        } else {
            triples.tn_triple.mean += approximation[i];
            triples.tn_triple.max = std::max(triples.tn_triple.max, approximation[i]);
            triples.tn_triple.min = std::min(triples.tn_triple.min, approximation[i]);
            tn += 1;
        }
    }
    triples.tp_triple.mean/=std::max(tp, 1);
    triples.tn_triple.mean/=std::max(tn, 1);
    triples.fp_triple.mean/=std::max(fp, 1);
    triples.fn_triple.mean/=std::max(fn, 1);

    triples.tp_triple/=iterations;
    triples.tn_triple/=iterations;
    triples.fp_triple/=iterations;
    triples.fn_triple/=iterations;
    similarity = (double) tp / (double) (tp + fp + fn);
}

void SetupExperiment::analyzeApproximations(std::map<int, int>& approximationFrequency, std::vector<int>& approximation){
    for (int val : approximation) {
        if (approximationFrequency.find(val) != approximationFrequency.end()){
            ++approximationFrequency[val];
        }
        else{
            approximationFrequency[val] = 1;
        }
    }
}

void SetupExperiment::CorrectClosure(GraphClosureSP& gc, GraphData& graph, ClosureParameters& closureParameters, double density, std::vector<std::string>& header, std::vector<bool>& realClosure, bool save) {
    closureParameters.clear();
    auto start = std::chrono::high_resolution_clock::now();
    gc.naive_closure(graph, closureParameters);
    double time = (double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
    std::cout << std::endl;
    for(auto elem: closureParameters.closed_set){
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::vector<std::string> values = {graph.getName(), std::to_string(graph.size()), std::to_string(graph.get_graph()->GetEdges()),
                                       std::to_string(density), std::to_string(closureParameters.input_set.size()),
                                       std::to_string(closureParameters.closed_set.size()), "Correct Closure",
                                       std::to_string(1),
                                       std::to_string(closureParameters.closed_set.size()),
                                       std::to_string(graph.get_graph()->GetNodes() - closureParameters.closed_set.size()),
                                       std::to_string(0),
                                       std::to_string(0),
                                       std::to_string(time), "1", "",
                                       "","", "", "", "", "", "", ""};
    std::vector<std::string> stat_headers;
    std::vector<std::string> stat_values;
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    statistics.evaluate(stat_headers, stat_values);
    header.insert(header.end(), stat_headers.begin(), stat_headers.end());
    values.insert(values.end(), stat_values.begin(), stat_values.end());

    if(save) {
        StaticFunctions::saveValuesToFile(_out_path, header, values);
    }
    std::fill(realClosure.begin(), realClosure.end(), false);
    for (NodeId id: closureParameters.closed_set) {
        realClosure[id] = true;
    }
}

void
SetupExperiment::SimpleApproximation(GraphClosureSP &gc, GraphData &graph, ClosureParameters closureParameters,
                                               double density, const std::vector<std::string>& header,
                                               const std::vector<bool>& realClosure, std::vector<int>& approximation, Triples& triples, bool preClosure,
                                               int incrementalApproximations, double timeConstraint, bool save) {
    closureParameters.clear();
    closureParameters.onlyPreClosure = preClosure;
    closureParameters.incrementalCount = incrementalApproximations;
    closureParameters.timeConstraint = timeConstraint;
    auto start = std::chrono::high_resolution_clock::now();
    gc.naive_closure(graph, closureParameters);
    double time = (double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
    double similarity = 0;
    int tp = 0, tn = 0, fp = 0, fn = 0;
    std::fill(approximation.begin(), approximation.end(), 0);
    for (auto elem : closureParameters.closed_set) {
        approximation[elem] = 1;
    }
    evaluateApproximation(approximation, realClosure, 1, 1, similarity, tp, tn, fp, fn, triples);
    std::vector<std::string> values = {graph.getName(),
                                       std::to_string(graph.size()),
                                       std::to_string(graph.get_graph()->GetEdges()),
                                       std::to_string(density),
                                       std::to_string(closureParameters.input_set.size()),
                                       std::to_string(closureParameters.closed_set.size()),
                                       "Incremental Count",
                                       std::to_string(similarity),
                                       std::to_string(tp),
                                       std::to_string(tn),
                                       std::to_string(fp),
                                       std::to_string(fn),
                                       std::to_string(time), "1", "",
                                       "", std::to_string(incrementalApproximations), "", "", "", "", "", ""};
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    std::vector<std::string> stat_headers;
    std::vector<std::string> stat_values;
    statistics.evaluate(stat_headers, stat_values);
    values.insert(values.end(), stat_values.begin(), stat_values.end());
    if (save) {
        StaticFunctions::saveValuesToFile(_out_path, {}, values);
    }
}

void SetupExperiment::SamplingApproximation(GraphClosureSP &gc, GraphData &graph, std::vector<GraphData>& subgraphData, std::vector<NodeId>& neighborIds, OuterplanarSubgraphDFS* outerPlanarSubgraphDfs,
                                                      const std::string& type, PatternType patternType,
                                                      ClosureParameters closureParameters, double density,
                                                      const std::vector<std::string>& header, const std::vector<bool>& realClosure,
                                                      std::vector<int>& approximation, double& runtime,
                                                      SetupExperiment::Triples triples,const std::vector<double>& thresholds, bool save) {
    closureParameters.clear();
    //Approximation vector
    double closure_time = 0;
    double generation_time = 0;
    std::fill(approximation.begin(), approximation.end(), 0);
    double edges = 0;
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    if (subgraphData.size() == 1) {
        GraphData& sD = subgraphData[0];
        for (int i = 0; i < _approximation_iterations; ++i) {
            PUNGraph subgraph = new TUNGraph();
            sD.set_graph(subgraph);
            std::mt19937_64 generator(i + _approximation_iterations * _seed);
            std::chrono::time_point<std::chrono::system_clock> start_generation;
            switch (patternType) {
                case PatternType::BFS_TREE:
                    start_generation = std::chrono::high_resolution_clock::now();
                    GraphFunctions::bfsSubtree(graph.get_graph(), sD, neighborIds, generator);
                    sD.graphType = GraphType::TREE;
                    generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - start_generation).count() /
                                        1000000.0);
                    break;
                case PatternType::OUTERPLANAR:
                    start_generation = std::chrono::high_resolution_clock::now();
                    outerPlanarSubgraphDfs->generate(sD, generator, false);
                    sD.graphType = GraphType::OUTERPLANAR;
                    generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - start_generation).count() /
                                        1000000.0);
                    statistics += OuterplanarGraphStatistics(sD.get_graph());
                    break;
            }
            edges += sD.get_graph()->GetEdges();
            auto start_closure = std::chrono::high_resolution_clock::now();
            closureParameters.closed_set.clear();
            gc.naive_closure(sD, closureParameters);
            for (NodeId id: closureParameters.closed_set) {
                approximation[id] += 1;
            }
            closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
        }
    }
    else{
        _approximation_iterations = (int) subgraphData.size();
        int sample_number;
        gc.approx_closure(graph, subgraphData, sample_number, closureParameters, approximation, closure_time);
        for (auto& sample : subgraphData) {
            edges += sample.get_graph()->GetEdges();
        }
    }
    generation_time /= this->_threads;
    closure_time /= this->_threads;
    edges /= this->_threads;
    runtime = closure_time + generation_time;
    double similarity = 0;
    int tp = 0, tn = 0, fp = 0, fn = 0;
    std::map<int, int> approximationFreqs;
    analyzeApproximations(approximationFreqs, approximation);
    for (double threshold: thresholds) {
        evaluateApproximation(approximation, realClosure, _approximation_iterations, threshold, similarity, tp, tn, fp, fn, triples);
        std::vector<std::string> stat_headers;
        std::vector<std::string> stat_values;
        statistics.evaluate(stat_headers, stat_values);
        std::vector<std::string> values = {graph.getName(), std::to_string(graph.size()),
                  std::to_string(
                          edges /
                          (double) _approximation_iterations),
                  std::to_string(density),
                  std::to_string(closureParameters.input_set.size()),
                  std::to_string(tp + fp),
                  type,
                  std::to_string(similarity),
                  std::to_string(tp),
                  std::to_string(tn),
                  std::to_string(fp),
                  std::to_string(fn),
                  std::to_string(runtime),
                  std::to_string(this->_threads),
                  std::to_string(closure_time),
                  std::to_string(generation_time),
                  std::to_string(_approximation_iterations),
                  std::to_string(threshold),
                  triples.tp_triple.to_string(),
                  triples.tn_triple.to_string(),
                  triples.fp_triple.to_string(),
                  triples.fn_triple.to_string(),
                  StaticFunctions::printMap(approximationFreqs)};
        values.insert(values.end(), stat_values.begin(), stat_values.end());
        if(save) {
            StaticFunctions::saveValuesToFile(_out_path, {}, values);
        }
    }
}

void SetupExperiment::SetApproximationIterations(GraphClosureSP &gc, GraphData &graph, std::pair<double, double>& max_threshold_similarity,
                                                 std::vector<GraphData>& subgraphData, std::vector<NodeId>& neighborIds,
                                                 OuterplanarSubgraphDFS* outerPlanarSubgraphDfs,
                                                 const std::string& type, PatternType patternType,
                                                 ClosureParameters& closureParameters, const std::vector<bool>& realClosure,
                                                 std::vector<int>& approximationVector, int max_iterations, double lower_bound) {
    closureParameters.clear();
    //Approximation vector
    double closure_time = 0;
    double generation_time = 0;
    double runtime = 0;
    SetupExperiment::Triples triples = Triples();
    max_threshold_similarity = {0, 0};
    int best_sample_num = 0;
    std::fill(approximationVector.begin(), approximationVector.end(), 0);
    double edges = 0;
    for (int i = 0; i < max_iterations; ++i) {
        subgraphData.emplace_back(GraphData(new TUNGraph, (int) graph.size()));
        std::mt19937_64 generator(i + max_iterations * _seed);
        std::chrono::time_point<std::chrono::system_clock> start_generation;
        switch (patternType) {
            case PatternType::BFS_TREE:
                start_generation = std::chrono::high_resolution_clock::now();
                GraphFunctions::bfsSubtree(graph.get_graph(), subgraphData.back(), neighborIds, generator);
                subgraphData.back().graphType = GraphType::TREE;
                generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start_generation).count() /
                                    1000000.0);
                break;
            case PatternType::OUTERPLANAR:
                start_generation = std::chrono::high_resolution_clock::now();
                outerPlanarSubgraphDfs->generate(subgraphData.back(), generator, false);
                subgraphData.back().save_dot(_path + "/Graphs/Graph_Outerplanar_" + "_Size_" + std::to_string(subgraphData.back().size()) + "_GraphEdges_" + std::to_string(graph.get_graph()->GetEdges()) + "_Sample_" + std::to_string(i));
                subgraphData.back().graphType = GraphType::OUTERPLANAR;
                generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start_generation).count() /
                                    1000000.0);
                break;
        }
        edges += subgraphData.back().get_graph()->GetEdges();
        auto start_closure = std::chrono::high_resolution_clock::now();
        gc.naive_closure(subgraphData.back(), closureParameters);
        for (NodeId id: closureParameters.closed_set) {
            approximationVector[id] += 1;
        }
        closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);

        generation_time /= this->_threads;
        closure_time /= this->_threads;
        edges /= this->_threads;
        runtime += (closure_time + generation_time);
        double similarity = 0;
        int tp = 0, tn = 0, fp = 0, fn = 0;
        std::vector<double> thresholds;
        thresholds.reserve(this->t_max_number);
        for (int j = 0; j < this->t_max_number; ++j) {
            thresholds.emplace_back(((double) (j+1)/ (double) t_max_number));
        }
        for (double threshold: thresholds) {
            evaluateApproximation(approximationVector, realClosure, i + 1, threshold, similarity, tp, tn, fp, fn, triples);
            if (similarity > max_threshold_similarity.second){
                max_threshold_similarity.first = threshold;
                max_threshold_similarity.second = similarity;
                best_sample_num = i+1;
            }
        }
        if (max_threshold_similarity.second > lower_bound){
            break;
        }
    }
    while (subgraphData.size() > best_sample_num) {
        subgraphData.pop_back();
    }
}

void SetupExperiment::SamplingApproximation_parallel(GraphClosureSP &gc, GraphData &graph, GraphData &subgraphData, std::vector<NodeId>& neighborIds, OuterplanarSubgraphDFS* outerPlanarSubgraphDfs,
                                            const std::string& type, PatternType patternType,
                                            ClosureParameters closureParameters, double density,
                                            const std::vector<std::string>& header, const std::vector<bool>& realClosure,
                                            std::vector<int>& approximation, double& runtime,
                                            SetupExperiment::Triples triples,const std::vector<double>& thresholds, bool save) {
    closureParameters.clear();
    //Approximation vector
    double closure_time = 0;
    double generation_time = 0;
    std::fill(approximation.begin(), approximation.end(), 0);
    double edges = 0;
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
#pragma omp parallel shared(patternType, graph, gc, statistics, approximation) default(none) firstprivate(subgraphData, neighborIds, outerPlanarSubgraphDfs, closureParameters) reduction(+: generation_time, edges, closure_time)
    {
        for (int i = 0; i < _approximation_iterations; ++i) {
            PUNGraph subgraph = new TUNGraph();
            subgraphData.set_graph(subgraph);
            std::mt19937_64 generator(i + _approximation_iterations * _seed);
            auto start_generation = std::chrono::high_resolution_clock::now();
            switch (patternType) {
                case PatternType::BFS_TREE:
                    start_generation = std::chrono::high_resolution_clock::now();
                    GraphFunctions::bfsSubtree(graph.get_graph(), subgraphData, neighborIds, generator);
                    subgraphData.graphType = GraphType::TREE;
                    generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - start_generation).count() /
                                        1000000.0);
                    break;
                case PatternType::OUTERPLANAR:
                    start_generation = std::chrono::high_resolution_clock::now();
                    outerPlanarSubgraphDfs->generate(subgraphData, generator, false);
                    subgraphData.graphType = GraphType::OUTERPLANAR;
                    generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - start_generation).count() /
                                        1000000.0);
#pragma omp critical
                    {
                        statistics += OuterplanarGraphStatistics(subgraphData.get_graph());
                    }
                    break;
            }
            edges += subgraphData.get_graph()->GetEdges();
            auto start_closure = std::chrono::high_resolution_clock::now();
            closureParameters.closed_set.clear();
            gc.naive_closure(subgraphData, closureParameters);
#pragma omp critical
            {
                for (NodeId id: closureParameters.closed_set) {
                    approximation[id] += 1;
                }
            }
            closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
        }
    }
    generation_time /= this->_threads;
    closure_time /= this->_threads;
    edges /= this->_threads;
    runtime = closure_time + generation_time;
    double similarity = 0;
    int tp = 0, tn = 0, fp = 0, fn = 0;
    std::map<int, int> approximationFreqs;
    analyzeApproximations(approximationFreqs, approximation);
    for (double threshold: thresholds) {
        evaluateApproximation(approximation, realClosure, _approximation_iterations, threshold, similarity, tp, tn, fp, fn, triples);
        std::vector<std::string> stat_headers;
        std::vector<std::string> stat_values;
        statistics.evaluate(stat_headers, stat_values);
        std::vector<std::string> values = {graph.getName(), std::to_string(graph.size()),
                                           std::to_string(
                                                   edges /
                                                   (double) _approximation_iterations),
                                           std::to_string(density),
                                           std::to_string(closureParameters.input_set.size()),
                                           std::to_string(tp + fp),
                                           type,
                                           std::to_string(similarity),
                                           std::to_string(tp),
                                           std::to_string(tn),
                                           std::to_string(fp),
                                           std::to_string(fn),
                                           std::to_string(runtime),
                                           std::to_string(this->_threads),
                                           std::to_string(closure_time),
                                           std::to_string(generation_time),
                                           std::to_string(_approximation_iterations),
                                           std::to_string(threshold),
                                           triples.tp_triple.to_string(),
                                           triples.tn_triple.to_string(),
                                           triples.fp_triple.to_string(),
                                           triples.fn_triple.to_string(),
                                           StaticFunctions::printMap(approximationFreqs)};
        values.insert(values.end(), stat_values.begin(), stat_values.end());
        if (save) {
            StaticFunctions::saveValuesToFile(_out_path, {}, values);
        }
    }
}

void SetupExperiment::SetSimpleApproximationIterations(GraphClosureSP& gc, GraphData& graph, int &optimal_iteration_number,
                                                       ClosureParameters& closureParameters,
                                                       int correctClosureSize, double lower_bound) {
    closureParameters.clear();
    closureParameters.lower_bound = lower_bound;
    closureParameters.target_set_size = correctClosureSize;
    closureParameters.approximate = true;
    gc.naive_closure(graph, closureParameters);
    optimal_iteration_number = closureParameters.iteration_number;
}

void SetupExperiment::VariationCheck(const GraphData& graph, const ClosureParameters& closureParameters, const std::vector<bool>& realClosure, std::vector<int> &variation_check,
                                     Triples& triples) {
    double similarity = 0;
    int tp = 0, tn = 0, fp = 0, fn = 0;
    evaluateApproximation(variation_check, realClosure, 1, 1, similarity, tp, tn, fp, fn, triples);



    int size = std::count(variation_check.begin(), variation_check.end(), 1);
    std::vector<std::string> values = {graph.getName(),
                                       std::to_string(graph.size()),
                                       std::to_string(graph.get_graph()->GetEdges()),
                                       std::to_string(0),
                                       std::to_string(closureParameters.input_set.size()),
                                       std::to_string(size),
                                       "Variation Check",
                                       std::to_string(similarity),
                                       std::to_string(tp),
                                       std::to_string(tn),
                                       std::to_string(fp),
                                       std::to_string(fn),
                                       std::to_string(0), "1", "",
                                       "", std::to_string(0), "", "", "", "", "", ""};
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    std::vector<std::string> stat_headers;
    std::vector<std::string> stat_values;
    statistics.evaluate(stat_headers, stat_values);
    values.insert(values.end(), stat_values.begin(), stat_values.end());
    StaticFunctions::saveValuesToFile(_out_path, {}, values);
}







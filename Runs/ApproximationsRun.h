//
// Created by anonymous on 25.08.21.
//

#ifndef CLOSURES_APPROXIMATIONSRUN_H
#define CLOSURES_APPROXIMATIONSRUN_H

#include "../Experiments/DS2021/SetupExperiment.h"

static void run_approximations_test(const std::string& input_path, const std::string& out_path, int seed, int iterations = 10, int threads = 1, bool synthetic = false){
    std::mt19937_64 gen(seed);
    SetupExperiment experiment = SetupExperiment("../out/Experiments", out_path, {0}, {0}, {4, 2},
                                                 gen, 1, iterations, 1, threads);
    SetupExperiment synthetic_experiment = SetupExperiment("../out/Experiments", out_path, {10000},
                                              {1.0, 1.1, 1.2, 1.3, 1.5, 2, 3}, {3, 4, 5, 6, 8}, gen, 1,
                                              100, 1, threads);

    std::vector<GraphData> outerplanarSamples;
    std::vector<GraphData> treeSamples;
    std::pair<double, double> best_threshold_tree;
    std::pair<double, double> best_threshold_outerplanar;
    int simple_approx_iterations;
    std::vector<int> variation_check;
    int input_seed = 654968413;
    int generator_size = 200;
    int threshold_steps = 100;
    int max_iterations = 100;
    double lower_bound = 0.95;
    std::vector<std::string> graph_paths;
    std::vector<GraphData> graphs;


    if (synthetic){
        synthetic_experiment.generateGraphs(graphs);
        for(auto& graph : graphs) {
            //Get BiggestComponent
            GraphFunctions::GetLargestComponent(graph);
            synthetic_experiment.samplePreprocessing(graph, treeSamples, best_threshold_tree,
                                           outerplanarSamples, best_threshold_outerplanar, simple_approx_iterations,
                                           variation_check, {generator_size}, 1, 0.1, input_seed, threshold_steps, lower_bound, max_iterations, nullptr);
            std::cout << "Generated " << std::to_string(treeSamples.size()) << " Tree Samples!" << std::endl;
            std::cout << "Generated " << std::to_string(outerplanarSamples.size()) << " Outerplanar Samples!" << std::endl;
            for (int i = 0; i < iterations; ++i) {
                synthetic_experiment.runApproximation(graph, treeSamples, best_threshold_tree,
                                            outerplanarSamples, best_threshold_outerplanar, simple_approx_iterations,
                                            variation_check, {generator_size}, i);
            }
        }
    }
    else {
        for (const auto &entry: std::filesystem::directory_iterator(input_path)) {
            if (entry.path().extension() == ".edges" &&
                GraphData(entry.path().string()).get_data()->GetNodes() < 100000) {
                graph_paths.emplace_back(entry.path().string());
            }
        }
        graph_paths = {"../../GraphData/CA-CondMat.edges"};
        for (const auto &graph_path: graph_paths) {
            GraphData graph = GraphData(graph_path);
            //Get BiggestComponent
            GraphFunctions::GetLargestComponent(graph);
            experiment.samplePreprocessing(graph, treeSamples, best_threshold_tree,
                                           outerplanarSamples, best_threshold_outerplanar, simple_approx_iterations,
                                           variation_check, {generator_size}, 1, 0.1, input_seed, threshold_steps,
                                           lower_bound, max_iterations, nullptr);
            std::cout << "Generated " << std::to_string(treeSamples.size()) << " Tree Samples!" << std::endl;
            std::cout << "Generated " << std::to_string(outerplanarSamples.size()) << " Outerplanar Samples!"
                      << std::endl;
            for (int i = 0; i < iterations; ++i) {
                experiment.runApproximation(graph, treeSamples, best_threshold_tree,
                                            outerplanarSamples, best_threshold_outerplanar, simple_approx_iterations,
                                            variation_check, {generator_size}, i);
            }
        }
    }
    //experiment.runApproximation({PatternType::BFS_TREE, PatternType::OUTERPLANAR}, {0.02}, {"../../GraphData/CA-AstroPh.edges"}, generator_seed);
}

//Experiments with real graphs smaller than 100000 nodes
static void run_approximations_real_graphs_small(const std::string& out_path, int threads = 1){
    std::mt19937_64 gen(0);
    SetupExperiment sCA = SetupExperiment("../out/Experiments", out_path, {0}, {0}, {2, 4, 6},
                                                              gen, 1, 100, 1, threads);

    std::vector<std::string> graph_paths;
    for (const auto & entry : std::filesystem::directory_iterator("../../GraphData/")){
        if (entry.path().extension() == ".edges" && GraphData(entry.path().string()).get_data()->GetNodes() < 100000) {
            graph_paths.emplace_back(entry.path().string());
        }
    }
    sCA.runApproximation({PatternType::BFS_TREE, PatternType::OUTERPLANAR}, {0.01, 0.02, 0.03, 0.04, 0.05, 0.06}, graph_paths, 0);
}

static void run_approximations_road_subgraphs(const std::string& out_path, int threads = 1){
    std::mt19937_64 gen(0);
    int approximation_iterations = 100;
    SetupExperiment sCA = SetupExperiment("../out/Experiments", out_path, {0}, {0}, {2, 4, 6},
                                                              gen, 1, approximation_iterations, 1, threads);

    std::vector<std::string> graph_paths = {"../../GraphData/roadNet-CA_seed_0_size_200000.edges", "../../GraphData/roadNet-CA_seed_0_size_200000.edges"};
    std::vector<double> thresholds(approximation_iterations, 0.0);
    double val = 1.0/approximation_iterations;
    for (double & threshold : thresholds) {
        threshold = val;
        val+=1.0/approximation_iterations;
    }
    sCA.runApproximation({PatternType::BFS_TREE, PatternType::OUTERPLANAR}, thresholds, graph_paths, 0);
}

static void run_approximations_synthetic_graphs(const std::string& out_path, int threads = 1){
    std::mt19937_64 gen(0);
    SetupExperiment sCA = SetupExperiment("../out/Experiments", out_path, {10000},
                                                              {1, 1.1, 1.2, 1.3, 1.5, 2, 3}, {3, 4, 5, 6, 8}, gen, 1,
                                                              100, 1, threads);
    std::vector<GraphData> graphs;
    sCA.generateGraphs(graphs);
    sCA.runApproximation({PatternType::BFS_TREE, PatternType::OUTERPLANAR}, {0.02});
}
#endif //CLOSURES_APPROXIMATIONSRUN_H

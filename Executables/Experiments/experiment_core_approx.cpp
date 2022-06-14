//
// Created by anonymous on 30.09.21.
//

#include <iostream>
#include "../../Data/GraphData.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Utils/StaticFunctions.h"
#include "../../Utils/GraphFunctions.h"

// Comparison function to sort the vector elements
// by second element of tuples
bool sortbysecond(const std::tuple<std::string, int>& a,
               const std::tuple<std::string, int>& b)
{
    return (std::get<1>(a) < std::get<1>(b));
}

struct OverlapApproxParams{
    //Paths
    std::string in_path = "../../GraphData/";
    std::string out_path = "../out/Approximation/";
    std::string eval_name = "kdd_2022";

    int thread_num = 4;

    int max_nodes = std::numeric_limits<int>::max();
    int max_edges = std::numeric_limits<int>::max();
    int min_nodes = 0;
    int min_edges = 0;

    std::vector<int> generator_size = {10};
    int generator_seed = 0;
    int coreIterations = 3;
    int sample_number = 100;
    std::vector<double> threshold = {1.0 / (double) this->sample_number};

    int outer_loop = 0;
    double overall_threshold = 1.0 / (double) this->sample_number;
    bool overall = false;
    int sample_seed = 32487643;
    bool outerplanar_new = true;
    
    //Eval
    bool periphery = false;
    bool core_iteration = false;
    bool exact = true;
    bool approx_core = true;

    int runtime = 0;
    int simple_approx_iterations = 0;

    bool core_stats = false;
    bool outerplanar_statistics = false;
    bool threshold_percentage = false;
    bool tree_eval = true;
    bool outer_fixed_point = false;
    bool save_load_samples = false;
    bool small_graphs = true;
    bool large_graphs = false;
    bool save_outerplanar_approx = true;
    bool delete_samples = true;

    //inline variables
    std::string stripped_path;
};

void simple_approx(GraphData &graph, int generator_size, int core_iterations, int seed, std::map<int, int>& degree_distribution,
                   TIntV& coreNodes, double& runtime, double time_constraint = -1, int iteration_constraint = 0) {
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<NodeId> v_intersection;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < core_iterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        closureParameters.iteration_number = iteration_constraint;
        closureParameters.timeConstraint = time_constraint/core_iterations;
        gc.naive_closure(graph, closureParameters);
        std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;
        if (i == 0) {
            overlap = closureParameters.closed_set;
        } else {
            v_intersection.clear();
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closureParameters.closed_set.begin(), closureParameters.closed_set.end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            if (i == core_iterations -1) {
                closureParameters.closed_set.clear();
                closureParameters.closed_set.insert(v_intersection.begin(), v_intersection.end());
            }
            else {
                for (auto elem: v_intersection) {
                    overlap.insert(elem);
                }
            }
        }
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: closureParameters.closed_set) {
        coreNodes.Add(elem);
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
}

void approximate_core(GraphData& graph, std::vector<GraphData>& samples, int generator_size, OverlapApproxParams& params, std::map<int, int>& degree_distribution, std::vector<TIntV>& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<int> overall_approximation = std::vector<int>(graph.size(), 0);
    double closure_time = 0;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    std::vector<int> approximation = std::vector<int>(graph.size(), 0);

    std::vector<std::set<NodeId>> coreSets = std::vector<std::set<NodeId>>(params.threshold.size(), std::set<NodeId>());
    for (int i = 0; i < params.coreIterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + params.coreIterations * params.generator_seed);
        int sample_number = (int) samples.size();
        gc.approx_closure(graph, samples, sample_number, closureParameters, approximation, closure_time, params.save_load_samples);
        if (!params.overall) {
            for (int j = 0; j < params.threshold.size(); ++j) {
                StaticFunctions::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                            (int) sample_number, params.threshold[j]);

                for (int k = 0; k < params.outer_loop; ++k) {
                    closureParameters.input_set = closureParameters.closed_set;
                    gc.approx_closure(graph, samples, sample_number, closureParameters, approximation,
                                      closure_time,
                                      params.save_load_samples);
                    StaticFunctions::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                                (int) sample_number,params.threshold[j]);
                }

                std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;

                if (i == 0) {
                    coreSets[j] = closureParameters.closed_set;
                } else {
                    //Intersect with previous found
                    auto it1 = coreSets[j].begin();
                    auto it2 = closureParameters.closed_set.begin();
                    while ((it1 != coreSets[j].end()) && (it2 != closureParameters.closed_set.end())) {
                        if (*it1 < *it2) {
                            coreSets[j].erase(it1++);
                        } else if (*it2 < *it1) {
                            ++it2;
                        } else { // *it1 == *it2
                            ++it1;
                            ++it2;
                        }
                    }
                    coreSets[j].erase(it1, coreSets[j].end());
                }
                std::cout << "\tCore Size: " << coreSets[j].size() << std::endl;
            }
        }
        else{
            for(int j = 0; j < approximation.size(); ++j){
                overall_approximation[j] += approximation[j];
            }
        }
    }
    if (params.overall){
        StaticFunctions::getClosedFromApproximation(overall_approximation, closureParameters.closed_set, (int) samples.size(), params.overall_threshold);
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.clear();
    for(int j = 0; j < coreSets.size(); ++j) {
        coreNodes.emplace_back(TIntV());
        for (auto elem: coreSets[j]) {
            coreNodes[j].Add(elem);
        }
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
}

void approximate_core(GraphData& graph, std::vector<OuterplanarGraphData>& samples, int generator_size, OverlapApproxParams& params, std::map<int, int>& degree_distribution, std::vector<TIntV>& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<int> overall_approximation = std::vector<int>(graph.size(), 0);
    double closure_time = 0;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    std::vector<int> approximation = std::vector<int>(graph.size(), 0);

    std::vector<std::set<NodeId>> coreSets = std::vector<std::set<NodeId>>(params.threshold.size(), std::set<NodeId>());
    for (int i = 0; i < params.coreIterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + params.coreIterations * params.generator_seed);
        int sample_number = (int) samples.size();
        gc.approx_closure(graph, samples, sample_number, closureParameters, approximation, closure_time, params.save_load_samples);
        if (!params.overall) {
            for (int j = 0; j < params.threshold.size(); ++j) {
                StaticFunctions::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                            (int) sample_number, params.threshold[j]);
                for (int k = 0; k < params.outer_loop; ++k) {
                    closureParameters.input_set = closureParameters.closed_set;
                    gc.approx_closure(graph, samples, sample_number, closureParameters, approximation,
                                      closure_time,
                                      params.save_load_samples);
                    StaticFunctions::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                                (int) sample_number,params.threshold[j]);
                }

                std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;

                if (i == 0) {
                    coreSets[j] = closureParameters.closed_set;
                } else {
                    //Intersect with previous found
                    auto it1 = coreSets[j].begin();
                    auto it2 = closureParameters.closed_set.begin();
                    while ((it1 != coreSets[j].end()) && (it2 != closureParameters.closed_set.end())) {
                        if (*it1 < *it2) {
                            coreSets[j].erase(it1++);
                        } else if (*it2 < *it1) {
                            ++it2;
                        } else { // *it1 == *it2
                            ++it1;
                            ++it2;
                        }
                    }
                    coreSets[j].erase(it1, coreSets[j].end());
                }
                std::cout << "\tCore Size: " << coreSets[j].size() << std::endl;
            }
        }
        else{
            for(int j = 0; j < approximation.size(); ++j){
                overall_approximation[j] += approximation[j];
            }
        }
    }
    if (params.overall){
        StaticFunctions::getClosedFromApproximation(overall_approximation, closureParameters.closed_set, (int) samples.size(), params.overall_threshold);
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.clear();
    for(int j = 0; j < coreSets.size(); ++j) {
        coreNodes.emplace_back(TIntV());
        for (auto elem: coreSets[j]) {
            coreNodes[j].Add(elem);
        }
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
    if (params.save_outerplanar_approx) {
        for(int j = 0; j < coreNodes.size(); ++j) {
            StaticFunctions::save(params.stripped_path, coreNodes[j],
                                  "approx_core_" + std::to_string((int) (params.threshold[j] * 100)));
        }
    }
}



void overlap_eval(OverlapApproxParams& params) {
    omp_set_num_threads(params.thread_num);
    std::vector<std::tuple<std::string, int>> paths_and_sizes;
    std::vector<std::string> stripped_paths;

    for (const auto &entry: std::filesystem::recursive_directory_iterator(params.in_path)) {
        std::string stripped_path =
                entry.path().parent_path().string() + "/" + entry.path().stem().string();
        if(!params.exact){
            if(entry.path().extension() == ".bin" && !std::filesystem::is_regular_file(stripped_path + ".core"))
            {
                stripped_paths.emplace_back(stripped_path);
            }
        }
        else if (entry.path().extension() == ".bin" && std::filesystem::is_regular_file(stripped_path + ".core")) {
            stripped_paths.emplace_back(stripped_path);
        }
    }

#pragma omp parallel for default(none) shared(stripped_paths, paths_and_sizes, params) schedule(dynamic,1)
    for (auto const &stripped_path: stripped_paths) {
        const GraphData &data = GraphData(stripped_path + ".bin");
        int nodes = data.nodes();
        int edges = data.edges();

        if (params.small_graphs && GraphFunctions::SmallGraph(data)) {
            if (nodes < params.max_nodes && edges < params.max_edges && nodes > params.min_nodes &&
                edges > params.min_edges) {
#pragma omp critical
                paths_and_sizes.emplace_back(std::tuple<std::string, int>(stripped_path + ".bin", edges));
            }
        }
        else if (params.large_graphs && !GraphFunctions::SmallGraph(data)) {
            if (nodes < params.max_nodes && edges < params.max_edges && nodes > params.min_nodes &&
                edges > params.min_edges) {
#pragma omp critical
                paths_and_sizes.emplace_back(std::tuple<std::string, int>(stripped_path + ".bin", edges));
            }
        }
    }
    std::sort(paths_and_sizes.begin(), paths_and_sizes.end(), sortbysecond);

#pragma omp parallel for default(none) shared(paths_and_sizes) firstprivate(params) schedule(dynamic,1)
    for (const auto &[path, size]: paths_and_sizes) {
        params.stripped_path =
                std::filesystem::path(path).parent_path().string() + "/" + std::filesystem::path(path).stem().string();
        TIntV exact_core;
        StaticFunctions::load(params.stripped_path + ".core", exact_core);
        bool exact_computation;
        if (exact_core.Len() > 0) {
            exact_computation = true;
        } else {
            exact_computation = false;
        }
        if (!std::filesystem::is_directory(params.out_path)){
            std::filesystem::create_directory(params.out_path);
        }
        FileEvaluation eval = FileEvaluation(params.out_path, params.eval_name, ".csv");

        double exact_runtime = -1;
        std::vector<std::vector<std::string>> info_array;
        if (exact_computation) {
            StaticFunctions::load_csv(params.stripped_path + ".core_info", info_array);
            exact_runtime = std::stoi(info_array[1][9]);
        }
        GraphData graph = GraphData(path);
        //std::cout << "Graph Nodes: " << graph.nodes() << std::endl;
        GraphFunctions::GetLargestComponent(graph);
        //GraphFunctions::analyse_graph(graph.get_graph(), graph.getName(), false, &eval);
        eval.headerValueInsert({"Graph", "Nodes", "Edges", "Density", "Core Size", "Exact Runtime"},
                               {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                                std::to_string(graph.density()),
                                std::to_string(exact_core.Len()),
                                std::to_string(exact_runtime)}, -2, true, true);

        for (int i=0; i < params.generator_size.size(); ++i) {
            int generator_size = params.generator_size[i];
            for (int j = 0; j < params.threshold.size(); ++j) {
                double threshold = params.threshold[j];
                eval.headerValueInsert({"Threshold", "Generator Size", "OuterLoop"},
                                       {std::to_string(threshold),
                                        std::to_string(generator_size), std::to_string(params.outer_loop)}, i*(int) params.threshold.size()+ j,
                                       true, true);
            }
        }
        std::map<int, int> core_degree_distribution;
        if (exact_computation) {
            GraphFunctions::GetCoreGraphStats(graph, exact_core, exact_runtime, "Exact Core", &eval);
        }

        double tree_sampling_runtime = 0, outerplanar_sampling_runtime = 0;
        double tree_core_runtime = 0, outerplanar_core_runtime = 0;
        double conversion_runtime = 0;
        std::vector<std::string> out_stat_headers;
        std::vector<std::string> out_stat_values;
        if (params.approx_core) {
            std::vector<int> variation_check;
            std::vector<NodeId> neighborIds;
            GraphFunctions::generateNeighborVector(graph.get_graph(), neighborIds);
            OuterplanarGraphStatistics statistic;

            //Approximate cores
            double overall_tree = 0, overall_outerplanar = 0;
            std::vector<TIntV> approx_cores;
            bool save_load_samples = params.save_load_samples && !(graph.nodes() < 350000 || graph.edges() < 2000000);
            for (auto type: {PatternType::BFS_TREE, PatternType::OUTERPLANAR}) {
                if (type == PatternType::BFS_TREE && params.tree_eval) {
                    //Tree samples
                    std::vector<GraphData> treeSamples;
                    GraphFunctions::GetSamples(graph, PatternType::BFS_TREE, treeSamples, nullptr, neighborIds,
                                               params.sample_number, params.sample_seed, tree_sampling_runtime,
                                               save_load_samples);
                    //std::cout << "Tree Sampling Runtime: " << runtime << "s" << std::endl;

                    //std::cout << std::endl << "Approximate Tree Core:" << std::endl;
                    for (int i=0; i < params.generator_size.size(); ++i) {
                        int generator_size = params.generator_size[i];
                        approximate_core(graph, treeSamples, generator_size, params,
                                         core_degree_distribution,
                                         approx_cores, tree_core_runtime);
                        for (int j=0; j<approx_cores.size(); ++j) {
                            TIntV& approx_core = approx_cores[i];
                            if (params.core_stats) {
                                GraphFunctions::GetCoreGraphStats(graph, approx_core, tree_core_runtime,
                                                                  "Tree Core",
                                                                  &eval);
                            } else {
                                eval.headerValueInsert(
                                        {"Tree Core Nodes", "Tree Core Relative Nodes", "Tree Core Edges",
                                         "Tree Core Relative Edges", "Tree Core Out Edges"},
                                        {std::to_string(approx_core.Len()),
                                         std::to_string((double) graph.nodes() / approx_core.Len()),
                                         std::to_string(0),
                                         std::to_string(0),
                                         std::to_string(0)}, i*(int) params.threshold.size()+ j);
                            }
                        }
                        overall_tree = tree_sampling_runtime + tree_core_runtime;
                        //std::cout << "Tree Core Runtime: " << tree_core_runtime << "s" << std::endl;
                        //std::cout << "Tree Runtime: " << overall_tree << "s" << std::endl;
                        //GraphFunctions::analyse_graph(core_graph.get_graph(), "Approximation Core Tree");

                        for (int j = 0; j < approx_cores.size(); ++j) {
                            eval.headerValueInsert({"Tree Core Runtime", "Tree Runtime"},
                                                   {std::to_string(tree_core_runtime),
                                                    std::to_string(overall_tree)}, i*(int) params.threshold.size()+ j, true,
                                                   true);
                            if (exact_computation) {
                                int intersection_length = approx_cores[j].IntrsLen(exact_core);
                                int union_length = approx_cores[j].UnionLen(exact_core);
                                double recall = (double) intersection_length / exact_core.Len();
                                double similarity = (double) intersection_length / union_length;
                                //std::cout << "Tree Recall: " << recall << std::endl;
                                //std::cout << "Tree Similarity: " << similarity << std::endl;
                                //Approx values
                                eval.headerValueInsert(
                                        {"Tree Recall", "Tree Similarity"},
                                        {
                                                std::to_string(recall),
                                                std::to_string(similarity)}, i*(int) params.threshold.size()+ j, true, true);
                            } else {
                                eval.headerValueInsert(
                                        {"Tree Recall", "Tree Similarity"},
                                        {
                                                std::to_string(0),
                                                std::to_string(0)}, i*(int) params.threshold.size()+ j, true, true);
                            }
                        }
                    }
                    for (const auto &entry: std::filesystem::recursive_directory_iterator(
                            "../out/Samples/" + graph.getName() + "/")) {
                        std::filesystem::remove(entry);
                    }
                } else if (type == PatternType::OUTERPLANAR) {
                    std::vector<GraphData> outerplanarSamples;
                    std::vector<OuterplanarGraphData> outerplanarSamples_new;
                    //Outerplanar samples
                    OuterplanarSubgraphDFS outerPlanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
                    if (!params.outerplanar_new) {
                        GraphFunctions::GetSamples(graph, PatternType::OUTERPLANAR, outerplanarSamples,
                                                   &outerPlanarSubgraphDfs, neighborIds, params.sample_number,
                                                   params.sample_seed, outerplanar_sampling_runtime,
                                                   save_load_samples);
                        if (params.outerplanar_statistics) {
                            for (auto const &out_graph: outerplanarSamples) {
                                statistic += OuterplanarGraphStatistics(out_graph.get_graph());
                            }
                        }
                    } else {
                        GraphFunctions::GetOuterplanarSamples(graph, PatternType::OUTERPLANAR,
                                                              outerplanarSamples_new,
                                                              &outerPlanarSubgraphDfs,
                                                              params.sample_number,
                                                              params.sample_seed, outerplanar_sampling_runtime,
                                                              conversion_runtime, true, params.save_load_samples);
                        if (params.outerplanar_statistics) {
                            for (auto const &out_graph: outerplanarSamples_new) {
                                statistic += OuterplanarGraphStatistics(out_graph.get_graph());
                            }
                        }
                    }

                    //std::cout << std::endl << "Approximate Outerplanar Core:" << std::endl;
                    for (int i = 0; i < params.generator_size.size(); ++i) {
                        if (!params.outerplanar_new) {
                            approximate_core(graph, outerplanarSamples, params.generator_size[i], params,
                                             core_degree_distribution,
                                             approx_cores, outerplanar_core_runtime);
                        } else {
                            approximate_core(graph, outerplanarSamples_new, params.generator_size[i], params,
                                             core_degree_distribution,
                                             approx_cores, outerplanar_core_runtime);
                        }
                        for (int j = 0; j<approx_cores.size();++j) {
                            if (params.core_stats) {
                                GraphFunctions::GetCoreGraphStats(graph, approx_cores[j],
                                                                  outerplanar_core_runtime,
                                                                  "Outerplanar Core", &eval);
                            } else {
                                eval.headerValueInsert({"Outerplanar Core Nodes", "Outerplanar Core Relative Nodes",
                                                        "Outerplanar Core Edges", "Outerplanar Core Relative Edges",
                                                        "Outerplanar Core Out Edges"},
                                                       {std::to_string(approx_cores[j].Len()),
                                                        std::to_string((double) graph.nodes() / approx_cores[j].Len()),
                                                        std::to_string(0),
                                                        std::to_string(0),
                                                        std::to_string(0)}, i*(int) params.threshold.size()+ j);
                            }
                        }
                        overall_outerplanar = outerplanar_sampling_runtime + conversion_runtime + outerplanar_core_runtime;
                        //std::cout << "Outerplanar Core Runtime: " << outerplanar_core_runtime << "s" << std::endl;
                        //std::cout << "Outerplanar Runtime: " << overall_outerplanar << "s" << std::endl;
                        //GraphFunctions::analyse_graph(core_graph.get_graph(), "Approximation Core Outerplanar");


                        for (int j = 0; j < approx_cores.size(); ++j) {
                            eval.headerValueInsert({"Outerplanar Core Runtime", "Outerplanar Runtime"},
                                                   {std::to_string(outerplanar_core_runtime),
                                                    std::to_string(overall_outerplanar)}, i*(int) params.threshold.size()+ j, true, true);

                            if (exact_computation) {
                                if (exact_computation) {
                                    int intersection_length = approx_cores[j].IntrsLen(exact_core);
                                    int union_length = approx_cores[j].UnionLen(exact_core);
                                    double recall = (double) intersection_length / exact_core.Len();
                                    double similarity = (double) intersection_length / union_length;
                                    //std::cout << "Outerplanar Recall: " << recall << std::endl;
                                    //std::cout << "Outerplanar Similarity: " << similarity << std::endl;
                                    //Approx values
                                    eval.headerValueInsert(
                                            {"Outerplanar Recall",
                                             "Outerplanar Similarity"},
                                            {
                                                    std::to_string(recall),
                                                    std::to_string(similarity)}, i*(int) params.threshold.size()+ j,
                                            true, true);
                                } else {
                                    eval.headerValueInsert(
                                            {"Outerplanar Recall",
                                             "Outerplanar Similarity"},
                                            {
                                                    std::to_string(0),
                                                    std::to_string(0)}, i*(int) params.threshold.size()+ j,
                                            true, true);
                                }
                            }
                        }
                    }
                    if (params.delete_samples) {
                        if (std::filesystem::exists("../out/Samples/" + graph.getName() + "/")) {
                            for (const auto &entry: std::filesystem::recursive_directory_iterator(
                                    "../out/Samples/" + graph.getName() + "/")) {
                                std::filesystem::remove(entry);
                            }
                        }
                    }
                }
                if (params.outerplanar_statistics) {
                    statistic.evaluate(out_stat_headers, out_stat_values);
                }
            }
            if (params.outerplanar_statistics) {
                eval.headerValueInsert(out_stat_headers, out_stat_values, -2, true, true);
            }
            graph.get_graph().Clr();
#pragma omp critical
            eval.save(true, true);
        }

    }
}


int main(int argc, char *argv[]) {
    //std::string out_path = "../out/ICDE_2022/";

    OverlapApproxParams params;


    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0){
            params.in_path = std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "-o") == 0){
            params.out_path = std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--file_name") == 0){
            params.eval_name = std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators") == 0){
            int j = i+1;
            params.generator_size.clear();
            while (j < argc && std::string(argv[j]).find('-') == std::string::npos){
                params.generator_size.emplace_back(atoi(argv[j]));
                ++j;
            }
        }
        if (std::strcmp(argv[i], "--threshold") == 0){
            int j = i+1;
            params.threshold.clear();
            while (j < argc && std::string(argv[j]).find('-') == std::string::npos){
                params.threshold.emplace_back(atof(argv[j]));
                ++j;
            }
        }
        if (std::strcmp(argv[i], "--generators_seed") == 0){
            params.generator_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--core_iterations") == 0){
            params.coreIterations = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--samples") == 0){
            params.sample_number = atoi(argv[i+1]);
        }

        if (std::strcmp(argv[i], "--overall_threshold") == 0){
            params.overall_threshold = atof(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_nodes") == 0){
            params.max_nodes = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_edges") == 0){
            params.max_edges = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--min_nodes") == 0){
            params.min_nodes = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--min_edges") == 0){
            params.min_edges = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--overall") == 0){
            params.overall = true;
        }
        if (std::strcmp(argv[i], "--no_tree") == 0){
            params.tree_eval = false;
        }
        if (std::strcmp(argv[i], "--sample_seed") == 0){
            params.sample_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--naive_outerplanar") == 0){
            params.outerplanar_new = false;
        }
        if (std::strcmp(argv[i], "--outerplanar_stats") == 0){
            params.outerplanar_statistics = true;
        }
        if (std::strcmp(argv[i], "--core_stats") == 0){
            params.core_stats = true;
        }
        if (std::strcmp(argv[i], "--threads") == 0){
            params.thread_num = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--percentage") == 0){
            params.threshold_percentage = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--outer_loop") == 0){
            params.outer_loop = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--outer_fixed_point") == 0){
            params.outer_fixed_point = true;
        }
        if (std::strcmp(argv[i], "--no_core") == 0){
            params.exact = false;
        }
        if (std::strcmp(argv[i], "--save_load_samples") == 0){
            params.save_load_samples = true;
        }
        if (std::strcmp(argv[i], "--small_graphs") == 0){
            params.small_graphs = true;
            params.large_graphs = false;
        }
        if (std::strcmp(argv[i], "--large_graphs") == 0){
            params.small_graphs = false;
            params.large_graphs = true;
            params.save_load_samples = true;
        }
        if (std::strcmp(argv[i], "--save_outerplanar_approx") == 0){
            params.save_outerplanar_approx = true;
        }
        if (std::strcmp(argv[i], "--keep_samples_saved") == 0){
            params.delete_samples = false;
        }

    }
        overlap_eval(params);
}
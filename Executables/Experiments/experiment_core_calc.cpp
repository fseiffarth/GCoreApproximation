//
// Created by anonymous on 30.09.21.
//

#include <iostream>
#include "../../Data/GraphData.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Utils/StaticFunctions.h"
#include "../../Experiments/DS2021/SetupExperiment.h"



struct OverlapApproxParams{
    //Paths
    std::string in_path = "../../GraphData/";
    std::string out_path = "../../GraphData/";

    int thread_num = 4;
    int max_nodes = std::numeric_limits<int>::max();
    int max_edges = std::numeric_limits<int>::max();
    int min_nodes = 0;
    int min_edges = 0;

    int generator_size = 10;
    int generator_seed = 0;
    int coreIterations = 3;
    int sample_number = 100;
    double threshold = 1.0 / (double) this->sample_number;
    double overall_threshold = 1.0 / (double) this->sample_number;
    bool overall = false;
    int sample_seed = 32487643;
    bool recalculate = false;

    //Eval
    bool periphery = false;
    bool core_iteration = false;
    bool exact = true;
    bool approx_core = false;
};



void get_core(GraphData& graph, int generator_size, int core_iterations, int seed, std::vector<int>& intersection_loss, TIntV& coreNodes, double& runtime, double& avg_preclosure_steps){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<std::set<NodeId>> closures;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    avg_preclosure_steps = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < core_iterations; ++i) {
        std::cout << "Iteration " << std::to_string(i) << " of graph " << graph.getName() << " started" << std::endl;
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        gc.naive_closure(graph, closureParameters);
        avg_preclosure_steps += closureParameters.preclosure_steps;
        closures.emplace_back(closureParameters.closed_set);
        if (i == 0) {
            overlap = closureParameters.closed_set;
        } else {
            int overlap_size = (int) overlap.size();
            std::vector<NodeId> v_intersection;
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closures.back().begin(), closures.back().end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            overlap.insert(v_intersection.begin(), v_intersection.end());
            intersection_loss.emplace_back(overlap_size - overlap.size());
        }
        std::cout << "Iteration " << std::to_string(i) << " of graph " << graph.getName() << " finished after " << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count() /1000000.0)) << "s"  << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: overlap) {
        coreNodes.Add(elem);
    }
    avg_preclosure_steps /= core_iterations;
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);

}

void core_calc(OverlapApproxParams& params) {
    std::vector<std::tuple<std::string, int>> paths_and_sizes;
    for (const auto &entry: std::filesystem::recursive_directory_iterator(params.in_path)) {
        if (entry.path().extension() == ".bin") {
            std::string stripped_path = entry.path().parent_path().string() + "/" + entry.path().stem().string();
            if (params.recalculate || !std::filesystem::exists(stripped_path + ".core")) {
                const GraphData &data = GraphData(entry.path().string());
                int nodes = data.nodes();
                int edges = data.edges();
                if (nodes < params.max_nodes && edges < params.max_edges && nodes > params.min_nodes &&
                    edges > params.min_edges) {
                    paths_and_sizes.emplace_back(std::tuple<std::string, int>(entry.path().string(), edges));
                    std::cout << "Computing core of graph: " << data.getName() << std::endl;
                }
            }
        }
    }
    std::sort(paths_and_sizes.begin(), paths_and_sizes.end(), StaticFunctions::sortbysecond);

#pragma omp parallel for default(none) shared(paths_and_sizes, params)
    for (const auto &[path, size]: paths_and_sizes) {
        std::string stripped_path = std::filesystem::path(path).parent_path().string() + "/" + std::filesystem::path(path).stem().string();

        FileEvaluation eval = FileEvaluation(stripped_path, "", ".core_info");
        GraphData graph = GraphData(path);
        GraphFunctions::GetLargestComponent(graph);
        eval.headerValueInsert({"Graph", "Nodes", "Edges", "Density"}, {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()), std::to_string(graph.density())});

        std::vector<int> intersection_loss;
        TIntV exact_core;
        double exact_runtime = 0;
        double avg_preclosure_depth = 0;
        get_core(graph, params.generator_size, params.coreIterations, params.generator_seed, intersection_loss, exact_core, exact_runtime, avg_preclosure_depth);
        StaticFunctions::save(stripped_path, exact_core);

        eval.headerValueInsert( {"Core Size", "Generators", "CoreIterations", "Intersection Loss", "Exact Runtime", "Avg. Preclosure Depth"},
                                {std::to_string(exact_core.Len()),
                                 std::to_string(params.generator_size),
                                 std::to_string(params.coreIterations),
                                 StaticFunctions::print<std::vector<int>, int>(intersection_loss),
                                 std::to_string(exact_runtime),
                                 std::to_string(avg_preclosure_depth)});
        //GraphFunctions::analyse_graph(core_graph.get_graph(), "Core");

        std::stringstream sstream;
        sstream << "Graph " << graph.getName() << " Exact Runtime: " << exact_runtime << "s" << std::endl;
        StaticFunctions::PrintStream(sstream);
        eval.save();
    }
}



int main(int argc, char *argv[]) {
    omp_set_num_threads(omp_get_max_threads());

    OverlapApproxParams params;

    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0){
            params.in_path = std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "-o") == 0){
            params.out_path = std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--threads") == 0){
            omp_set_num_threads(atoi(argv[i+1]));
        }
        if (std::strcmp(argv[i], "--generators") == 0){
            params.generator_size = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_seed") == 0){
            params.generator_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--core_iterations") == 0){
            params.coreIterations = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--min_nodes") == 0){
            params.min_nodes = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--min_edges") == 0){
            params.min_edges = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_nodes") == 0){
            params.max_nodes = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_edges") == 0){
            params.max_edges = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--overall") == 0){
            params.overall = true;
        }
        if (std::strcmp(argv[i], "--recalculate") == 0){
            params.recalculate = true;
        }
        if (std::strcmp(argv[i], "--sample_seed") == 0){
            params.sample_seed = atoi(argv[i+1]);
        }
    }

    core_calc(params);
}
//
// Created by florian on 15.08.23.
//

#include <string>
#include "../Utils/CoreGrowAlgorithm.h"
#include "../Utils/StaticFunctions.h"
#include "../ClosureOperators/GraphClosures.h"
#include "../Utils/FileEvaluation.h"

void get_core(GraphData& graph, int generator_size, int core_iterations, int seed, std::vector<int>& intersection_loss, TIntV & coreNodes, double& runtime, double& avg_preclosure_steps){
    GraphClosureSP gc = GraphClosureSP();
    ClosureParameters closureParameters;
    std::vector<std::set<NodeId>> closures;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    avg_preclosure_steps = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto c_iterations = core_iterations;
    if (c_iterations == -1) {
        c_iterations = graph.nodes();
    }
    for (int i = 0; i < c_iterations; ++i) {
        std::cout << "\tIteration " << std::to_string(i) << " of graph " << graph.getName() << " started"
                  << std::endl;

        // generate input set
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        // compute the closure
        gc.closure(graph, closureParameters);
        closures.emplace_back(closureParameters.closed_set);
        if (i == 0) {
            overlap = closureParameters.closed_set;

            std::cout << "\tIteration " << std::to_string(i) << " of graph " << graph.getName() << " finished after "
                      << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                              std::chrono::high_resolution_clock::now() - start).count() / 1000000.0)) << "s"
                      << std::endl << std::endl;
        } else {
            int overlap_size = (int) overlap.size();
            std::vector<NodeId> v_intersection;
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closures.back().begin(), closures.back().end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            overlap.insert(v_intersection.begin(), v_intersection.end());
            intersection_loss.emplace_back(overlap_size - overlap.size());
            std::cout << "\tIteration " << std::to_string(i) << " of graph " << graph.getName() << " finished after "
                      << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                              std::chrono::high_resolution_clock::now() - start).count() / 1000000.0)) << "s"
                      << std::endl << std::endl;

            if (core_iterations == -1 && intersection_loss.back() == 0) {
                break;
            }
        }

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



void grow_vs_intersection_core()
{
    // num runs
    int num_runs =  10;
    int num_steps = 20;
    double threshold = 0.9;
    // iterate over several graphs
    std::string directory = "Collaboration";
    // paths
    std::string directory_path = "../../../GraphData/RealWorld/" + directory + "/";
    // get all files with extension .edges in the directory using std::filesystem
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == ".edges") {
            // check if entry path contains "CA-GrQc"
            if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                continue;
            }
            graph_paths.emplace_back(entry.path().string());
        }
    }

    for (const auto &path: graph_paths) {
        auto inputParameters = CoreGrowAlgorithmInputParameters{num_runs, num_steps, threshold, 0,
                                                                true, true};
        CoreGrowAlgorithmOutputParameters outputParameters;
        GraphData graph = GraphData(path);
        // print the graph information (nodes, edges, density)
        std::cout << "Nodes: " << graph.nodes() << std::endl;
        std::cout << "Edges: " << graph.edges() << std::endl;
        std::cout << "Density: " << graph.density()  << std::endl;
        // run the core grow algorithm
        // measure the runtime
        auto start = std::chrono::high_resolution_clock::now();
        auto CGA = CoreGrowAlgorithm(graph);
        CGA.Run(outputParameters, inputParameters);
        double grow_runtime = outputParameters.runtime;
        // print the core size
        std::cout << "Core size: " << outputParameters.core_nodes.size() << std::endl;

        // Calculate the exact core
        std::vector<int> intersection_loss;
        TIntV coreNodes;
        double runtime;
        double avg_preclosure_steps;
        get_core(graph, num_steps, 3, 0, intersection_loss, coreNodes, runtime, avg_preclosure_steps);
        // exact core nodes as std::vector
        std::vector<NodeId> exact_core_nodes;
        for (int i = 0; i < coreNodes.Len(); ++i) {
            exact_core_nodes.emplace_back(coreNodes[i]);
        }
        // compare outputParameters.core_nodes with coreNodes
        std::cout << "Exact Core size: " << coreNodes.Len() << std::endl;
        // compute Jaccard similarity between outputParameters.core_nodes and coreNodes i.e. the size of the intersection divided by the size of the union
        // intersection of outputParameters.core_nodes and coreNodes
        std::vector<NodeId> _intersection;
        // union of outputParameters.core_nodes and coreNodes
        std::vector<NodeId> _union;
        // compute the intersection and the union
        std::set_intersection(outputParameters.core_nodes.begin(), outputParameters.core_nodes.end(),
                              exact_core_nodes.begin(), exact_core_nodes.end(), std::back_inserter(_intersection));
        std::set_union(outputParameters.core_nodes.begin(), outputParameters.core_nodes.end(), exact_core_nodes.begin(),
                       exact_core_nodes.end(), std::back_inserter(_union));
        // compute the Jaccard similarity
        double jaccard_similarity = (double) _intersection.size() / (double) _union.size();
        // print the Jaccard similarity
        std::cout << "Jaccard similarity: " << jaccard_similarity << std::endl;

        // save results to file using FileEvaluation
        FileEvaluation fileEvaluation = FileEvaluation("../results/", directory, ".csv");
        fileEvaluation.headerValueInsert(
                {"Graph", "Nodes", "Edges", "Density", "Grow Core Size", "Grow Runtime", "Exact Core Size", "Exact Runtime", "Jaccard similarity",
                 "Avg. preclosure steps"},
                {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                 std::to_string(graph.density()), std::to_string(outputParameters.core_nodes.size()), std::to_string(grow_runtime),
                 std::to_string(coreNodes.Len()), std::to_string(runtime), std::to_string(jaccard_similarity),
                 std::to_string(avg_preclosure_steps)});
        fileEvaluation.save();
    }

}

void stability_grow_core() {
    // num runs
    auto num_runs_list = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    auto num_steps_list = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

    for (auto num_runs: num_runs_list) {
        for (auto num_steps: num_steps_list) {

            // iterate over several graphs
            std::string directory = "Collaboration";
            // paths
            std::string directory_path = "../../../GraphData/RealWorld/" + directory + "/";
            // get all files with extension .edges in the directory using std::filesystem
            std::vector<std::string> graph_paths;
            for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
                if (entry.path().extension() == ".edges") {
                    // check if entry path contains "CA-GrQc"
                    if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                        continue;
                    }
                    graph_paths.emplace_back(entry.path().string());
                }
            }

            for (const auto &path: graph_paths) {
                auto inputParameters = CoreGrowAlgorithmInputParameters{num_runs, num_steps, 0.9, 0,
                                                                        true, true};
                CoreGrowAlgorithmOutputParameters outputParameters;
                GraphData graph = GraphData(path);
                // print the graph information (nodes, edges, density)
                std::cout << "Nodes: " << graph.nodes() << std::endl;
                std::cout << "Edges: " << graph.edges() << std::endl;
                std::cout << "Density: " << graph.density() << std::endl;
                // run the core grow algorithm
                // measure the runtime
                auto start = std::chrono::high_resolution_clock::now();
                auto CGA = CoreGrowAlgorithm(graph);
                //CGA.Run(outputParameters, inputParameters);
                double grow_runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count() /
                                       1000000.0);
                // print the core size
                std::cout << "Core size: " << outputParameters.core_nodes.size() << std::endl;

                // save results to file using FileEvaluation
                FileEvaluation fileEvaluation = FileEvaluation("../results/", directory + "_stability_grow", ".csv");
                fileEvaluation.headerValueInsert(
                        {"Graph", "Nodes", "Edges", "Density", "Grow Core Size", "Grow Runtime", "Num Runs",
                         "Num Steps"},
                        {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                         std::to_string(graph.density()), std::to_string(outputParameters.core_nodes.size()),
                         std::to_string(grow_runtime), std::to_string(num_runs), std::to_string(num_steps)});
                fileEvaluation.save();
            }
        }
    }
}

void stabilility_intersection_core(){
    // num runs
    auto num_steps_list = {5}; //, 10, 15, 20, 25, 30, 35, 40, 45, 50};
        for (auto num_steps: num_steps_list) {

            // iterate over several graphs
            std::string directory = "Collaboration";
            // paths
            std::string directory_path = "../../../GraphData/RealWorld/" + directory + "/";
            // get all files with extension .edges in the directory using std::filesystem
            std::vector<std::string> graph_paths;
            for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
                if (entry.path().extension() == ".edges") {
                    // check if entry path contains "CA-GrQc"
                    if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                        continue;
                    }
                    graph_paths.emplace_back(entry.path().string());
                }
            }

            for (const auto &path: graph_paths) {
                GraphData graph = GraphData(path);
                // print the graph information (nodes, edges, density)
                std::cout << "Nodes: " << graph.nodes() << std::endl;
                std::cout << "Edges: " << graph.edges() << std::endl;
                std::cout << "Density: " << graph.density() << std::endl;
                // Calculate the exact core
                std::vector<int> intersection_loss;
                TIntV coreNodes;
                double runtime;
                double avg_preclosure_steps;
                get_core(graph, num_steps, -1, 0, intersection_loss, coreNodes, runtime, avg_preclosure_steps);
                // exact core nodes as std::vector
                std::vector<NodeId> exact_core_nodes;
                for (int i = 0; i < coreNodes.Len(); ++i) {
                    exact_core_nodes.emplace_back(coreNodes[i]);
                }
                // compare outputParameters.core_nodes with coreNodes
                std::cout << "Exact Core size: " << coreNodes.Len() << std::endl;


                // save results to file using FileEvaluation
                FileEvaluation fileEvaluation = FileEvaluation("../results/", directory + "_stability_intersection", ".csv");
                fileEvaluation.headerValueInsert(
                        {"Graph", "Nodes", "Edges", "Density", "Intersection Core Size", "Exact Runtime", "Num Steps"},
                        {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                         std::to_string(graph.density()), std::to_string(coreNodes.Len()),
                         std::to_string(runtime), std::to_string(num_steps)});
                fileEvaluation.save();
            }
        }
}


int main(int argc, char *argv[]) {
    //grow_vs_intersection_core();
    //stability_grow_core();
    stabilility_intersection_core();
}
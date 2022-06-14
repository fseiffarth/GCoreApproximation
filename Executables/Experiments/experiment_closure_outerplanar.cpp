//
// Created by anonymous on 31.08.21.
//


#include "../../Data/GraphData.h"
#include "../../Experiments/Experiments.h"
#include "../../Utils/FileEvaluation.h"
#include "../../Experiments/DS2021/SetupExperiment.h"
#include "../../Utils/StaticFunctions.h"

void closure_runtime(int number) {
    int n = 1000;
    std::string out_name = "closure_runtime";

    std::vector<std::string> headers = {"Size",
                                        "Edges",
                                        "Density",
                                        "Input Size",
                                        "C1 Time",
                                        "C2 Time",
                                        "O2 Generation Time",
                                        "O2 Closure Time",
                                        "CTree Time",
                                        "CGraph Time",
                                        "Samples"};

#pragma omp parallel for default(none) shared(n, headers, number, out_name)
    for (int j = 0; j < 8; ++j) {
        double p = 0.006 + 0.002*j;
        std::vector<std::pair<int, int>> graph_sizes;
        for (int i = 1; i <= 10; ++i) {
            graph_sizes.emplace_back(std::pair<int, int>{i * n, (int) (i * n * (i * n - 1) / 2 * p)});
        }
        for (auto const &size: graph_sizes) {
            std::stringstream sstream;
            sstream << "Nodes " << size.first << " Edges " << size.second  << " Edge Probability: " << p << std::endl;
            StaticFunctions::PrintStream(sstream);
            ClosureParameters closureParameters;
            GraphClosureSP cl;
            std::vector<GraphData> outerplanar;
            std::vector<GraphData> graphs;
            std::vector<GraphData> trees;
            size_t o1 = 0;
            size_t o2_generation = 0;
            size_t o2_closure = 0;
            size_t tree_time = 0;
            size_t cgraph_time = 0;
            int input_size = size.first / 100;
            if (!std::filesystem::is_directory("../out/Closure/")){
                std::filesystem::create_directory("../out/Closure/");
            }
            FileEvaluation fileEvaluation = FileEvaluation("../out/Closure/", out_name);
            OuterplanarGraphStatistics statistics;
            std::vector<std::string> values;
            outerplanar.clear();
            graphs.clear();
            trees.clear();
            sstream.clear();
            sstream << "Get Samples" << std::endl;
            StaticFunctions::PrintStream(sstream);
            Experiments::GetOuterplanarSamples(size, number, outerplanar, graphs, trees);
            sstream.clear();
            sstream << "Calculate Closures" << std::endl;
            StaticFunctions::PrintStream(sstream);
            int edges = 0;
            for (int i = 0; i < outerplanar.size(); ++i) {
                statistics += OuterplanarGraphStatistics(outerplanar[i].get_graph());
                std::set<NodeId> input_set;
                StaticFunctions::generateInputSet(closureParameters.input_set, outerplanar[i], input_size, i);
                edges += outerplanar[i].edges();

                auto start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(outerplanar[i], closureParameters);
                o1 += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                OuterplanarGraphData outerplanar_new = OuterplanarGraphData(outerplanar[i]);
                o2_generation += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(outerplanar_new, closureParameters);
                o2_closure += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(graphs[i], closureParameters);
                cgraph_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(trees[i], closureParameters);
                tree_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();
            }
            std::vector<std::string> stat_headers;
            std::vector<std::string> stat_values;
            statistics.evaluate(stat_headers, stat_values);
            edges /= (int) outerplanar.size();
            fileEvaluation.headerValueInsert(headers, {std::to_string(size.first),
                                                       std::to_string(edges),
                                                       std::to_string(size.second / ((double) size.first/2 * (size.first - 1))),
                                                       std::to_string(input_size),
                                                       std::to_string((double) o1 / 1000000.0),
                                                       std::to_string((double) o2_generation / 1000000.0 + (double) o2_closure / 1000000.0),
                                                       std::to_string((double) o2_generation / 1000000.0),
                                                       std::to_string((double) o2_closure / 1000000.0),
                                                       std::to_string((double) tree_time/1000000.0),
                                                       std::to_string((double) cgraph_time / 1000000.0),
                                                       std::to_string(outerplanar.size())});
            fileEvaluation.headerValueInsert(stat_headers, stat_values);
#pragma omp critical
            fileEvaluation.save();
        }
    }
}


int main(){
    int number = 100;
    closure_runtime(number);
}
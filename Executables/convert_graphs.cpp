//
// Created by anonymous on 24.11.2021.
//

#include <string>
#include <sstream>
#include <filesystem>
#include <iostream>
#include "../Utils/FileEvaluation.h"
#include "../Data/GraphData.h"
#include "../Utils/GraphFunctions.h"

int main(int argc, char *argv[]) {
    omp_set_num_threads(omp_get_thread_num()/2);
    std::string in_path = "../../GraphData/";
    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0) {
            in_path = std::string(argv[i + 1]);
        }
    }
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::recursive_directory_iterator(in_path)) {
        std::string extension;
        if (entry.path().extension() == ".edges") {
            extension = ".edges";
        } else if (entry.path().extension() == ".txt") {
            extension = ".txt";
        }
        if (!extension.empty()) {
            graph_paths.emplace_back(entry.path());
            std::cout << "Found " << graph_paths.back() << std::endl;
        }
    }
#pragma omp parallel for shared(graph_paths) default(none)
    for (auto const & path : graph_paths) {
        auto const file_path = std::filesystem::path(path);
        std::string extension;
        if (file_path.extension() == ".edges"){
            extension = ".edges";
        }
        else if(file_path.extension() == ".txt"){
            extension = ".txt";
        }
        if (!extension.empty()) {
            if (!std::filesystem::is_directory("../out/Other/")) {
                std::filesystem::create_directory("../out/Other/");
            }
            FileEvaluation fileEvaluation = FileEvaluation("../out/Other/", "converted_graphs");
            std::string stripped_path =
                    file_path.parent_path().string() + "/" + file_path.stem().string();
            //GraphData::Update(stripped_path + ".txt");
            auto start = std::chrono::system_clock::now();
            GraphData data = GraphData(stripped_path + extension);
            GraphFunctions::GetLargestComponent(data, true);
            data.setName(file_path.stem().string() + "_component");
            //data.save_edges(entry.path().parent_path().string() + "/");
            data.save_bin(file_path.parent_path().string() + "/");
            size_t nodes = data.nodes();
            size_t edges = data.edges();
            double density = (double) edges / ((double) (nodes * (nodes - 1)) / 2);
            std::stringstream sstream;
            double time = (double) std::chrono::duration_cast<std::chrono::milliseconds>((std::chrono::system_clock::now() - start)).count() / 1000.0;
            sstream << "Finished conversion of " << data.getName() << " in " << time << "s" << std::endl;
            sstream << "Written to " << file_path.parent_path().string() + "/" + data.getName() << ".bin" << std::endl;
            StaticFunctions::PrintStream(sstream);
            fileEvaluation.headerValueInsert({"Graph", "Size", "Edges", "Density", "Conversion Time"},
                                             {data.getName(),
                                              std::to_string(nodes),
                                              std::to_string(edges),
                                              std::to_string(density),
                                              std::to_string(time),
                                             });
            fileEvaluation.save();
        }
    }
}
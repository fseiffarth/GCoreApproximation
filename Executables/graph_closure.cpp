//
// Created by florian on 03.05.23.
//

#include <string>
#include "../ClosureOperators/GraphClosures.h"
#include "../Utils/StaticFunctions.h"
#include "Snap.h"
#include "../Utils/GraphFunctions.h"

int main(int argc, char *argv[]) {
    std::string in_path;
    std::string out_path;
    std::vector<std::string> node_id_str;
    std::vector<NodeId> node_ids;
    std::string help = "Parameters: \n -i \t Input path for the graph \n -o \t Output path for the closure \n -ids \t list of vertex ids for calculating the closure, e.g. <0 1 2> or path to file of vertex ids for calculating the closure (one id per line)";

    for (int i = 0; i < argc; ++i) {
        std::string str = argv[i];
        if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "-help") == 0 || std::strcmp(argv[i], "--h") == 0 || std::strcmp(argv[i], "--help") == 0){
            std::cout << help << std::endl;
        }
        if (std::strcmp(argv[i], "-i") == 0) {
            in_path =argv[i + 1];
        }
        if (std::strcmp(argv[i], "-o") == 0) {
            out_path = argv[i + 1];
        }
        //Define vertex ids for closure computation
        if (std::strcmp(argv[i], "-ids") == 0) {
            std::string arg = argv[i+1];
            if (StaticFunctions::load(arg, node_ids)){
            }
            else {
                StaticFunctions::ArgumentList(argv, argc, i, node_id_str);
                for (auto &val: node_id_str) {
                    node_ids.emplace_back(std::stoi(val));
                }
            }
        }
    }

    if(!in_path.empty() && !out_path.empty()) {
        GraphClosureSP gc;
        ClosureParameters closureParameters;
        for (auto & val : node_ids) {
            closureParameters.input_set.emplace(val);
        }
        GraphData g = GraphData(in_path);
        GraphData graph = GraphData(GraphFunctions::ResetGraphIds(g.get_graph()));
        gc.naive_closure(graph, closureParameters);
        std::set<NodeId> closed_set = closureParameters.closed_set;
        StaticFunctions::save(out_path, closed_set);
    }
    else{
        std::cout << help << std::endl;
    }


}
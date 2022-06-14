//
// Created by anonymous on 07.11.21.
//

#ifndef CLOSURES_EXPERIMENTS_H
#define CLOSURES_EXPERIMENTS_H

#include <vector>
#include <iostream>
#include "../Data/GraphData.h"


class Experiments {

public:
    struct TestParameter{
        const std::vector<std::pair<int, int>> graph_sizes = {{10, 10}};
        int seedMin = 0;
        int seedMax = 1000;
        bool outerplanarity_check = true;
        int threads = 1;
    };

    static void GetOuterplanarSamples(const std::pair<int, int>& size, int number, std::vector<GraphData>& outerplanar, std::vector<GraphData>& graphs, std::vector<GraphData>& trees);
    static void GetOuterplanarSamplesBiconnected(const std::pair<int, int>& size, int number, std::vector<GraphData>& outerplanar, std::vector<GraphData>& graphs, std::vector<GraphData>& trees);
    static void OuterplanarSampling(const std::string& out_path,  Experiments::TestParameter& testParameter);

};


#endif //CLOSURES_EXPERIMENTS_H

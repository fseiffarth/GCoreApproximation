//
// Created by anonymous on 31.08.21.
//

#include "../../Experiments/Experiments.h"

void runtime_experiment() {
    int n = 1000;
    std::string out_name = "sampling_runtime";
#pragma omp parallel for default(none) shared(n, out_name)
    for (int j = 0; j < 8; ++j) {
        double p = 0.006 + 0.002*j;
        std::vector<std::pair<int, int>> graph_sizes;
        for (int i = 1; i <= 10; ++i) {
            graph_sizes.emplace_back(std::pair<int, int>{i * n, (int) (i * n * (i * n - 1) / 2 * p)});
        }
        auto params = Experiments::TestParameter{graph_sizes, 0, 99, false, 1};
        if (!std::filesystem::is_directory("../out/Sampling/")){
            std::filesystem::create_directory("../out/Sampling/");
        }
        Experiments::OuterplanarSampling("../out/Sampling/" + out_name, params);
    }
}

int main(){
    runtime_experiment();
}
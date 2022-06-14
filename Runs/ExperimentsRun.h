//
// Created by anonymous on 09.06.2021.
//
#ifndef CLOSURES_EXPERIMENTSRUN_H
#define CLOSURES_EXPERIMENTSRUN_H

#include "../Experiments/DS2021/ExperimentalSetup.h"

static void run_experiments(){
    GraphClosureSP gc("shortestPath");
    InitUpdateGreedy<GraphData> initUpdateGreedy("greedy");
    InitUpdateGreedyRandom<GraphData> initUpdateGreedyRandom("greedyRandom");

    ExperimentalSetup setup("Results", "../out/results/", "../out/graphs/", "../out/patterns/trees/",
                            "../out/patterns/outerplanar/", 1, {5},
                            {2}, {}, {1, 2, 4},
                            100, 1, 1, false, 1, true);

    //Graph and pattern generation
    setup.generateGraphs(10, GenerationType::TWO_COMPONENTS);
    std::cout << "Finished graph two component generation!" << std::endl;
    setup.generatePatterns(gc, 20, 100, GenerationType::TWO_COMPONENTS, PatternType::BFS_TREE);
    std::cout << "Finished pattern two component generation!" << std::endl;
    setup.generatePatterns(gc, 20, 100, GenerationType::TWO_COMPONENTS, PatternType::OUTERPLANAR);
    std::cout << "Finished pattern two component generation!" << std::endl;

//    //Graph and pattern generation
//    setup.generateGraphs(10, GenerationType::ONE_COMPONENT);
//    std::cout << "Finished graph two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 500, GenerationType::ONE_COMPONENT, PatternType::BFS_TREE);
//    std::cout << "Finished pattern two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 100, GenerationType::ONE_COMPONENT, PatternType::OUTERPLANAR);
//    std::cout << "Finished pattern two component generation!" << std::endl;

    //Algorithm run
    std::string path;

    std::cout << "Trees" << std::endl;
    setup.run(gc, initUpdateGreedyRandom, 5, 100, 1, path, GenerationType::TWO_COMPONENTS, PatternType::BFS_TREE);
    std::cout << "Outerplanar" << std::endl;
    setup.run(gc, initUpdateGreedyRandom, 5, 100, 1, path, GenerationType::TWO_COMPONENTS, PatternType::OUTERPLANAR);

}
#endif //CLOSURES_EXPERIMENTSRUN_H

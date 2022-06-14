//
// Created by anonymous on 10.05.2021.
//

#ifndef CLOSURES_EXPERIMENTALSETUP_H
#define CLOSURES_EXPERIMENTALSETUP_H


#include "../../Utils/typedefs.h"
#include "../../Data/GraphData.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Evaluations/evaluations.h"
#include "../../Utils/pattern_generator.h"



class ExperimentalSetup {
public:
    ExperimentalSetup(std::string resultsName, std::string resultPath, std::string graphPath, std::string treePath,
                      std::string outerplanarPath, int graphNumber,
                      const std::vector<int> &blockSizes, std::vector<double> blockDensities,
                      std::vector<int> blockConnections, std::vector<int> trainingSizes = {1, 2, 4},
                      int patterns = 500, int pattern_out_stepsize = 1, int trainingConfigurations = 10, bool randomConnections = false, int threads = 1, bool testMode = false);

    template <typename T1, typename T2>
    void run(T1& closure, T2& initUpdate, int seed, int max_patterns, int numGraphs, std::string& folderPath, GenerationType generationType = GenerationType::TWO_COMPONENTS, PatternType patternType = PatternType::BFS_TREE);
    template <typename T1, typename T2>
    void run(T1& closure, T2& initUpdate, int seed, int numPatterns, int numGraphs, bool withClosureOfTraining);
    void generateGraphs(int seed, GenerationType generationType = GenerationType::TWO_COMPONENTS);
    template <typename T>
    void generatePatterns(T& closure, int seed, int numPatterns, GenerationType generationType = GenerationType::TWO_COMPONENTS, PatternType patternType = PatternType::BFS_TREE);
    void analyzeFaces(GenerationType generationType = GenerationType::TWO_COMPONENTS);
private:
    std::string resultPath;
    std::string resultsName;
    std::string graphPath;
    std::string treePath;
    std::string outerplanarPath;
    int graphNumber;
    std::vector<int> blockSizes;
    std::vector<double> blockDensities;
    std::vector<int> blockConnections;
    std::vector<int> trainingSizes;
    int patternNumber;
    int trainingConfigurations;
    int _out_stepsize;
    bool _random_connections;
    int _threads;
    bool _testMode;
    template <typename T1, typename T2>
    void trainingConfigurationRun(T1& closure, T2& initUpdate, std::vector<FileEvaluation>& averageEvaluations, std::vector<GraphData>& graphs, std::vector<std::vector<GraphData>>& patterns, int maxPatterns, size_t& treeTime, int trainingSize, double density, int connections, PatternType type, std::mt19937_64& gen);
    template <typename T1, typename T2>
    void singleRun(GraphData &graph, T1 &closure, T2 &initUpdate, std::vector<GraphData> *patterns,int trainingSize, size_t &patternTime, Evaluations &patternEvaluation, std::mt19937_64 &gen) const;

    struct Runtimes{
        //Time vars
        size_t all{};
        size_t algorithm{};
        std::chrono::seconds singleRun{};
        size_t loading{};
        std::chrono::seconds generation{};
    };
    Runtimes runtimes;

    void savePatternGenerationRuntimes(const GraphData& graphData, double density, int connection, float perPatternTime, int patternNum, int graphNum, GenerationType generationType, PatternType patternType);

    static std::string createNewFolder(const std::string& pathString);


    void LoadGraphsAndPatternsForExperiment(std::vector<GraphData> &graphs,
                                            std::vector<std::vector<GraphData>> &patterns,
                                            GenerationType generationType, int blockSize,
                                            double blockDensity, int blockConnection, int numGraphs,
                                            PatternType patternType, int max_patterns);


};

#include "ExperimentalSetup.txx"


#endif //CLOSURES_EXPERIMENTALSETUP_H

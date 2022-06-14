//
// Created by anonymous on 25.08.21.
//

#ifndef CLOSURES_SETUPEXPERIMENT_H
#define CLOSURES_SETUPEXPERIMENT_H


#include <vector>
#include "../../Utils/Generators.h"
#include "../../Utils/pattern_generator.h"
#include "../../Utils/graphio.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../ClosureOperators/GraphClosures.h"


struct LoopParameters{
    int size;
    double density;
    int input_size;
};

class SetupExperiment {
public:
    SetupExperiment(const std::string &path, const std::string &out_path, const std::vector<int> &graphSizes,
                              const std::vector<double> &densities, const std::vector<int> &generatorElements,
                              const std::mt19937_64 &gen, int repetitions = 1, int approximation_iterations = 100,
                              int graphNumber = 1, int threads = 1)
            : _path(path), _out_path(out_path), _graphSizes(graphSizes), _densities(densities), _generatorElements(generatorElements), _gen(gen), _repetitions(repetitions), _approximation_iterations(approximation_iterations), _graphNumber(graphNumber), _threads(threads){
        if (threads != -1) {
            omp_set_num_threads(_threads);
        }
        else{
            omp_set_num_threads(omp_get_max_threads());
        }
    }


    void generateGraphs(std::vector<GraphData>& graphs);

    static void generateInputSet(std::set<NodeId>& input_set, const GraphData& graphData, int input_size, int seed);

    void runApproximation(const std::vector<PatternType>& patternTypes,const std::vector<double>& threshold = {0.5}, const std::vector<std::string>& graph_paths = {}, int seed = 0);

    struct ValTriple{
        ValTriple()= default;;
        int max = 0;
        double mean = 0;
        int min = 0;
        [[nodiscard]] std::string to_string() const{
            return "(" + std::to_string(max) + " : " + std::to_string(mean) + " : " + std::to_string(min) + ")";
        }
        void operator/=(int value){this->max/=value, this->mean/=value, this->min/=value;}
    };
    struct Triples{
        ValTriple tp_triple;
        ValTriple tn_triple;
        ValTriple fp_triple;
        ValTriple fn_triple;
    };


    void samplePreprocessing(GraphData& graph, std::vector<GraphData>& treeSamples, std::pair<double, double>& best_threshold_tree,
                             std::vector<GraphData>& outerplanarSamples, std::pair<double, double>& best_threshold_outerplanar,
                             int& simple_approx_iteration_number, std::vector<int>& variation_check,
                             const std::vector<int> &generator_sizes = {2}, int num_training_sets = 1, double training_elements = 0.1,
                             int seed = 0, int thresholds_max_number = 100, double lower_bound = 0.75, int max_iterations = 100,
                             void *generator_elements_choosing = nullptr, bool completeSampling = false);

    void
    runApproximation(GraphData& graph, std::vector<GraphData> &treeSamples, const std::pair<double, double>& best_threshold_tree, std::vector<GraphData>& outerplanarSamples, const std::pair<double, double>& best_threshold_outerplanar, int simple_approx_iteration_number,std::vector<int>& variation_check,
                     const std::vector<int> &generator_sizes, int seed);



    static void GetSamples(GraphData &data, PatternType type, std::vector<GraphData> &subgraphs,
                           OuterplanarSubgraphDFS *outerplanarSubgraphDFS, std::vector<NodeId>& neighborIds, int iterations, int seed, double& runtime);

private:
    std::string _path;
    std::string _out_path;
    std::vector<int> _graphSizes;
    std::vector<double> _densities;
    std::vector<int> _generatorElements;
    int _graphNumber;
    int _repetitions;
    int _approximation_iterations;
    std::mt19937_64 _gen;
    int _seed{};
    int _threads;

    //Sample data
    std::vector<PUNGraph> _outerplanarSamples;
    std::vector<PUNGraph> _treeSamples;
    int simple_approximation_iterations;
    int t_max_number;

    void estimate_iteration_threshold();

    void run(GraphClosureSP& gc, std::vector<GraphData> &graphs, std::set<NodeId>& inputSet, double density,const std::vector<double>& threshold, const std::vector<PatternType>& patternType);

    static void evaluateApproximation(std::vector<int>& approximation, const std::vector<bool>& realClosure, int iterations, double threshold, double& similarity, int& tp, int& tn, int& fp, int& fn, Triples& triples) ;

    static void analyzeApproximations(std::map<int, int> &approximationFrequency, std::vector<int> &approximation);


    void CorrectClosure(GraphClosureSP& gc, GraphData& graph, ClosureParameters& closureParameters, double density, std::vector<std::string>& header, std::vector<bool>& realClosure, bool save = false);

    void SimpleApproximation(GraphClosureSP &gc, GraphData &graph, ClosureParameters closureParameters, double density,
                             const std::vector<std::string>& header, const std::vector<bool>& realClosure, std::vector<int>& approximation, Triples& triples, bool preClosure = false,
                             int incrementalApproximations = -1, double timeConstraint = -1, bool save = false);

    void SamplingApproximation(GraphClosureSP &gc, GraphData &graph, std::vector<GraphData>& subgraphData, std::vector<int>& neighborIds,OuterplanarSubgraphDFS* outerPlanarSubgraphDfs,
                               const std::string& type, PatternType patternType,
                               ClosureParameters closureParameters, double density,
                               const std::vector<std::string>& header, const std::vector<bool>& realClosure,
                               std::vector<int>& approximation, double& runtime,
                               SetupExperiment::Triples triples,const std::vector<double>& thresholds, bool save = false);
    void SamplingApproximation_parallel(GraphClosureSP &gc, GraphData &graph, GraphData &subgraphData, std::vector<int>& neighborIds,OuterplanarSubgraphDFS* outerPlanarSubgraphDfs,
                               const std::string& type, PatternType patternType,
                               ClosureParameters closureParameters, double density,
                               const std::vector<std::string>& header, const std::vector<bool>& realClosure,
                               std::vector<int>& approximation, double& runtime,
                               SetupExperiment::Triples triples,const std::vector<double>& thresholds, bool save = false);

    void SetApproximationIterations(GraphClosureSP &gc, GraphData &graph, std::pair<double, double>& best_threshold, std::vector<GraphData> &subgraphData,
                                    std::vector<NodeId> &neighborIds, OuterplanarSubgraphDFS *outerPlanarSubgraphDfs,
                                    const std::string &type, PatternType patternType,
                                    ClosureParameters& closureParameters,
                                    const std::vector<bool> &realClosure, std::vector<int> &approximationVector, int max_iterations, double lower_bound = 0.75);

    void SetSimpleApproximationIterations(GraphClosureSP& sp, GraphData& data, int &optimal_iteration_number,
                                          ClosureParameters& closureParameters, int correctClosureSize,
                                          double lower_bound = 0.75);

    void VariationCheck(const GraphData& graph, const ClosureParameters& closureParameters ,const std::vector<bool>& realClosure, std::vector<int> &variation_check, Triples& triples);

};




#endif //CLOSURES_SETUPEXPERIMENT_H

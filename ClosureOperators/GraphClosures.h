//
// Created by anonymous on 12.04.2021.
//

#ifndef CLOSURES_GRAPHCLOSURES_H
#define CLOSURES_GRAPHCLOSURES_H


#include <string>
#include <vector>
#include <deque>
#include <set>
#include <unordered_set>
#include "../Data/GraphData.h"
#include "BaseOperator.h"
#include "../Data/OuterplanarGraphData.h"

typedef void (*ClosureFunction)(GraphData&, ClosureParameters&);

class GraphClosureSP : public BaseOperator<GraphData> {
public:
    GraphClosureSP();

    explicit GraphClosureSP(const std::string& name, int threshold = std::numeric_limits<int>::max());
    void closure(GraphData &dataObject, ClosureParameters &closureParameters) override;
    void naive_closure(GraphData &graph, ClosureParameters& closureParameters);
    void fast_closure(GraphData &graph, ClosureParameters& closureParameters);
    void outerplanar_closure(OuterplanarGraphData &outerplanarGraphData, ClosureParameters& closureParameters);
    void approx_closure(const GraphData& graph, std::vector<GraphData> &samples, int& sample_number, ClosureParameters &parameters, std::vector<int> &approximation,
                        double& closure_time, bool load = false, ClosureFunction closure = nullptr);
    void approx_closure(const GraphData& graph, std::vector<OuterplanarGraphData> &samples, int& sample_number, ClosureParameters &parameters, std::vector<int> &approximation,
                        double& closure_time, bool load = false, ClosureFunction closure = nullptr);


private:
    bool preClosure(GraphData& graph, NodeId element, std::set<NodeId>& preClosedSet, const std::vector<std::set<NodeId>*>& forbidden_elements, bool& output_intersects_forbidden, std::set<NodeId>& new_elements);
    void forward_step(GraphData &graph, std::set<NodeId> &closed_interval_set, NodeId element, std::deque<NodeId>& bfsQueue, const std::vector<std::set<NodeId> *> &forbidden_elements,
                      bool &output_intersects_forbidden) const;
    static bool forward_step_break_condition(const GraphData &graph, const std::deque<NodeId> &bfsQueue,
                                      int visitedElementsInInterval,
                                      int closedIntervalSetSize, int possiblePaths);
    bool backward_step(GraphData &DataObject, NodeId element, std::set<NodeId> &closed_interval_set, std::deque<NodeId>& bfsQueue, const std::vector<std::set<NodeId> *> &forbidden_elements, bool &output_intersects_forbidden,
                       std::set<NodeId> &new_elements);
    int _threshold;
    std::chrono::time_point<std::chrono::system_clock> time;

    bool incrementalStep(GraphData &graph, NodeId newElement, ClosureParameters &closureParameters);

    void bfs_forward(GraphData &graph, std::set<NodeId>& target_set, NodeId bfs_start,
                std::deque<NodeId> &bfsQueue) const;
    void bfs_fast_forward(GraphData &graph, std::set<NodeId>& target_set, NodeId bfs_start,
                     std::deque<NodeId> &bfsQueue) const;
    void bfs_backward(GraphData &graph, std::set<NodeId>& target_set, ClosureParameters &closureParameters,
    std::deque<NodeId> &bfsQueue, std::deque<NodeId>& bfsElements, std::unordered_set<NodeId>& generatedElements) const;

    void addToClosed(GraphData &graph, NodeId elem, ClosureParameters &closureParameters, std::deque<NodeId>& bfsQueue, std::unordered_set<NodeId>& generatedElements);

    void set_break_condition(bool &condition, std::deque<NodeId> &newElements, ClosureParameters &closureParameters,
                             int iterations);
    void get_generators(OuterplanarComponent& outerplanarComponent, std::set<NodeId> &input_set,
                        std::set<NodeId> &generator_set);


};


#endif //CLOSURES_GRAPHCLOSURES_H

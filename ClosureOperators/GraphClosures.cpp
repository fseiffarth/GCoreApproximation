//
// Created by anonymous on 12.04.2021.
//

#include "GraphClosures.h"
#include "../Utils/GraphFunctions.h"
#include <vector>
#include <deque>

GraphClosureSP::GraphClosureSP(const std::string& name, int threshold) : BaseOperator<GraphData>(name), _threshold(threshold){}

GraphClosureSP::GraphClosureSP() {}

/// gets the closures in the graphData object respective the closures settings
/// \param dataObject
/// \param closureParameters
//TODO there are some bugs
void GraphClosureSP::closure(GraphData &dataObject, ClosureParameters &closureParameters) {
    this->time = std::chrono::system_clock::now();
    closureParameters.clear();
    if(closureParameters.timeConstraint > 0) {
        closureParameters.timeConstraint = (int) (closureParameters.timeConstraint * 1000000);
    }
    else{
        closureParameters.timeConstraint = -1;
    }
    if (closureParameters.element_to_add != -1){
        closureParameters.added_elements.insert(closureParameters.element_to_add);
        //Add only one element to input_set (the result is only closed if closureParameters.input set is closed
        closureParameters.closed_set.insert(closureParameters.input_set.begin(), closureParameters.input_set.end());
        if (!incrementalStep(dataObject,closureParameters.element_to_add, closureParameters)){
            closureParameters.output_intersects_forbidden = true;
            return;
        }
    }
    else{
        size_t i = 0;
        for (NodeId elem : closureParameters.input_set) {
            if (i == 0){
                closureParameters.closed_set.insert(elem);
            }
            else{
                if (!incrementalStep(dataObject, elem, closureParameters)){
                    closureParameters.output_is_all = true;
                    return;
                }
            }
            ++i;
        }
    }
}




bool GraphClosureSP::preClosure(GraphData& graph, NodeId element, std::set<NodeId> &preClosedSet, const std::vector<std::set<NodeId>*>& forbidden_elements, bool& output_intersects_forbidden, std::set<NodeId> &new_elements) {
    preClosedSet.insert(element);
    new_elements.clear();
    std::fill(graph.DistanceList.begin(), graph.DistanceList.end(), -2);
    std::fill(graph.PredecessorsList.begin(), graph.PredecessorsList.end(), std::vector<NodeId>());
    std::deque<NodeId> bfsQueue;
    forward_step(graph, preClosedSet, element, bfsQueue,forbidden_elements, output_intersects_forbidden);
    return backward_step(graph, element, preClosedSet, bfsQueue, forbidden_elements, output_intersects_forbidden, new_elements);
}

bool GraphClosureSP::incrementalStep(GraphData& graph, NodeId newElement, ClosureParameters& closureParameters) {
    std::set<NodeId> new_closed = std::set<NodeId>{newElement};
    std::set<NodeId> new_elements = std::set<NodeId>();
    bool proceedWhile = !new_closed.empty() && closureParameters.incrementalCount != 0 && (closureParameters.timeConstraint == -1 || closureParameters.timeConstraint > std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count());
    if (closureParameters.approximate && closureParameters.target_set_size != 0){
        if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
            proceedWhile = false;
        }
    }
    while(proceedWhile){
        NodeId next_element = *new_closed.begin();
        new_closed.erase(new_closed.begin());
        if (!preClosure(graph, next_element, closureParameters.closed_set, closureParameters.forbidden_elements, closureParameters.output_intersects_forbidden, new_elements)){
            return false;
        }
        new_closed.insert(new_elements.begin(), new_elements.end());
        closureParameters.added_elements.insert(new_elements.begin(), new_elements.end());
        if (graph.graphType == GraphType::TREE || graph.graphType == GraphType::OUTERPLANAR || closureParameters.onlyPreClosure){
            break;
        }
        --closureParameters.incrementalCount;
        ++closureParameters.iteration_number;
        proceedWhile = !new_closed.empty() && closureParameters.incrementalCount != 0 && (closureParameters.timeConstraint == -1 || closureParameters.timeConstraint > std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count());
        if (closureParameters.approximate && closureParameters.target_set_size != 0){
            if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
                proceedWhile = false;
            }
        }
    }
    return true;
}

void GraphClosureSP::bfs_forward(GraphData &graph, std::set<NodeId>& target_set, NodeId bfs_start, std::deque<NodeId>& bfsQueue) const {

    //Forward Search through graph O(m)
    //all elements in closed set
    ++graph.Id;
    if (graph.Id == 0){
        std::fill(graph.ContainmentList.begin(), graph.ContainmentList.end(), 0);
        ++graph.Id;
    }
    for (NodeId elem : target_set) {
        graph.ContainmentList[elem] = graph.Id;
    }
    for (auto & x : graph.PredecessorsList) {
        x.clear();
    }
    std::fill(graph.DistanceList.begin(), graph.DistanceList.end(), -2);
    graph.DistanceList[bfs_start] = 0;
    bfsQueue.clear();
    bfsQueue.push_back(bfs_start);
    int visitedSize = 1;
    int lastDistance = std::numeric_limits<int>::max();
    while (!bfsQueue.empty()) {
        NodeId currentNodeId = bfsQueue.back();
        bfsQueue.pop_back();
        int currentDistance = graph.DistanceList[currentNodeId];
        if ((visitedSize == target_set.size() && graph.DistanceList[currentNodeId] > lastDistance)) {
            break;
        }
        TUNGraph::TNodeI Node = graph.get_graph()->GetNI(currentNodeId);
        for (int i = 0; i < Node.GetDeg(); ++i) {
            NodeId neighborId = Node.GetNbrNId(i);
            //Neighbor is unvisited by bfs
            if (graph.DistanceList[neighborId] < 0) {
                graph.DistanceList[neighborId] = currentDistance + 1;
                bfsQueue.push_front(neighborId);
                //Neighbor is in closed set
                if (graph.ContainmentList[neighborId] == graph.Id) {
                    ++visitedSize;
                }
            }
            if (graph.DistanceList[neighborId] == currentDistance + 1) {
                graph.PredecessorsList[neighborId].push_back(currentNodeId);
            }
        }
        if (visitedSize == target_set.size()) {
            lastDistance = currentDistance;
        }
    }
}

void GraphClosureSP::bfs_backward(GraphData &graph, std::set<NodeId>& target_set, ClosureParameters &closureParameters,
                                  std::deque<NodeId> &bfsQueue, std::deque<NodeId>& bfsElements, std::unordered_set<NodeId>& generatedElements) const {
    //Backward search
    ++graph.Id;
    if (graph.Id == 0){
        std::fill(graph.ContainmentList.begin(), graph.ContainmentList.end(), 0);
        ++graph.Id;
    }
    bfsQueue.clear();
    if (graph.graphType == GraphType::OUTERPLANAR){
        for (NodeId elem: closureParameters.input_set) {
            graph.ContainmentList[elem] = graph.Id;
            bfsQueue.push_back(elem);
        }
    }
    else {
        for (NodeId elem: target_set) {
            graph.ContainmentList[elem] = graph.Id;
            bfsQueue.push_back(elem);
        }
    }

    while (!bfsQueue.empty()){
        NodeId NextElementId = bfsQueue.back();
         bfsQueue.pop_back();
        for (NodeId NeighborNodeId : graph.PredecessorsList[NextElementId]) {
            if (graph.ContainmentList[NeighborNodeId] != graph.Id) {
                graph.ContainmentList[NeighborNodeId] = graph.Id;
                closureParameters.closed_set.insert(NeighborNodeId);
                bfsQueue.push_back(NeighborNodeId);

                //New elements should be only considered in the case that the graph is of general type
                if (graph.graphType != GraphType::TREE && graph.graphType != GraphType::OUTERPLANAR) {
                    bfsElements.push_back(NeighborNodeId);
                }
            }
            else{
                if (graph.graphType == GraphType::TREE || graph.graphType == GraphType::OUTERPLANAR){
                    generatedElements.insert(NeighborNodeId);
                }
            }
        }
    }

}


void GraphClosureSP::forward_step(GraphData &graph, std::set<NodeId> &closed_interval_set, NodeId element, std::deque<NodeId>& bfsQueue,
                                  const std::vector<std::set<NodeId>*> &forbidden_elements,
                                  bool &output_intersects_forbidden) const {

    //Forward Search through graph O(m)
    //all elements in closed set
    ++graph.Id;
    if (graph.Id == 0){
        std::fill(graph.ContainmentList.begin(), graph.ContainmentList.end(), 0);
        ++graph.Id;
    }
    for (auto & x : graph.PredecessorsList) {
        x.clear();
    }
    std::fill(graph.DistanceList.begin(), graph.DistanceList.end(), -2);
    for (NodeId elem : closed_interval_set) {
        graph.ContainmentList[elem] = graph.Id;
    }
    graph.DistanceList[element] = 0;
    bfsQueue.push_back(element);
    std::vector<NodeId> visitedElementsInInterval = std::vector<NodeId>(1, element);
    int possiblePaths = 1;
    int lastDistance = std::numeric_limits<int>::max();
    bool lastDistanceSet = false;
    while (!forward_step_break_condition(graph, bfsQueue, static_cast<int>(visitedElementsInInterval.size()),
                                        static_cast<int>(closed_interval_set.size()), possiblePaths)) {
        NodeId currentNodeId = bfsQueue.back();
        bfsQueue.pop_back();
        if (this->_threshold < graph.DistanceList[currentNodeId] || lastDistance < graph.DistanceList[currentNodeId]) {
            break;
        }
        bool currentNodeInClosed = graph.ContainmentList[currentNodeId] == graph.Id;
        bool currentNeighborInClosed;
        //Current node is already in closure and is not the first element
        if (currentNodeInClosed && currentNodeId != element){
            if (!lastDistanceSet) {
                TUNGraph::TNodeI Node = graph.get_graph()->GetNI(currentNodeId);
                for (int i = 0; i < Node.GetDeg(); ++i) {
                    NodeId neighborId = Node.GetNbrNId(i);
                    currentNeighborInClosed = graph.ContainmentList[neighborId] == graph.Id;
                    if (currentNeighborInClosed) {
                        if (graph.DistanceList[neighborId] < 0) {
                            visitedElementsInInterval.push_back(neighborId);
                            int next_distance = graph.DistanceList[currentNodeId] + 1;
                            graph.DistanceList[neighborId] = next_distance;
                            bfsQueue.push_front(neighborId);
                            if (!lastDistanceSet && visitedElementsInInterval.size() == closed_interval_set.size()) {
                                lastDistance = next_distance;
                                lastDistanceSet = graph.DistanceList[currentNodeId];
                            }
                        }
                        if (graph.DistanceList[neighborId] == graph.DistanceList[currentNodeId] + 1) {
                            graph.PredecessorsList[neighborId].push_back(currentNodeId);
                        }
                    }
                }
            }
        }
        //Current get_node is not in closure
        else{
            --possiblePaths;
            for (auto it = graph.begin(currentNodeId); it != graph.end(currentNodeId); ++it) {
                NodeId neighborId = *it;
                currentNeighborInClosed = graph.ContainmentList[neighborId] == graph.Id;
                if (graph.DistanceList[neighborId] < 0) {
                    int next_distance = graph.DistanceList[currentNodeId] + 1;
                    graph.DistanceList[neighborId] = next_distance;
                    bfsQueue.push_front(neighborId);
                    if (currentNeighborInClosed) {
                        visitedElementsInInterval.push_back(neighborId);
                        if (!lastDistanceSet && visitedElementsInInterval.size() == closed_interval_set.size()){
                            lastDistance = graph.DistanceList[currentNodeId];
                            lastDistanceSet = true;
                        }
                    }
                    else {
                        ++possiblePaths;
                    }
                }
                if (graph.DistanceList[neighborId] == graph.DistanceList[currentNodeId] + 1) {
                    graph.PredecessorsList[neighborId].push_back(currentNodeId);
                }
            }
        }
    }
}

bool GraphClosureSP::forward_step_break_condition(const GraphData &graph, const std::deque<NodeId> &bfsQueue,
                                                  int visitedElementsInInterval,
                                                  int closedIntervalSetSize, int possiblePaths) {
    switch (graph.graphType) {
        case GraphType::TREE:
            return bfsQueue.empty() || (visitedElementsInInterval >= 2) || possiblePaths == 0;
            break;
        default:
            return bfsQueue.empty() || possiblePaths == 0;
            break;
    }
}

bool GraphClosureSP::backward_step(GraphData &graph, NodeId element, std::set<NodeId> &closed_interval_set, std::deque<NodeId>& bfsQueue,
                                   const std::vector<std::set<NodeId> *> &forbidden_elements, bool &output_intersects_forbidden,
                                   std::set<NodeId> &new_elements) {
    //Backward search
    ++graph.Id;
    bfsQueue.clear();
    for (NodeId elem : closed_interval_set) {
        graph.ContainmentList[elem] = graph.Id;
        bfsQueue.push_back(elem);
    }
    while (!bfsQueue.empty()){
        NodeId NextElementId = bfsQueue.back();
        bfsQueue.pop_back();
        for (NodeId NeighborNodeId : graph.PredecessorsList[NextElementId]) {
            if (graph.ContainmentList[NeighborNodeId] != graph.Id){
                if (IsForbidden(NeighborNodeId, forbidden_elements)){
                    output_intersects_forbidden = true;
                    return false;
                }
                else{
                    graph.ContainmentList[NeighborNodeId] = graph.Id;
                    closed_interval_set.insert(NeighborNodeId);
                    new_elements.insert(NeighborNodeId);
                    bfsQueue.push_front(NeighborNodeId);
                }
            }
        }
    }
    return true;
}

void GraphClosureSP::naive_closure(GraphData &graph, ClosureParameters &closureParameters) {
    if (dynamic_cast<OuterplanarGraphData*>(&graph) != nullptr && graph.ContainmentList.empty()){
        graph.init();
    }
    if (closureParameters.input_set.empty()){
        return;
    }
    if (dynamic_cast<OuterplanarGraphData*>(&graph) != nullptr){
        outerplanar_closure(dynamic_cast<OuterplanarGraphData&>(graph), closureParameters);
        return;
    }
    std::unordered_set<NodeId> generatedElements;
    std::deque<NodeId> bfsQueue;
    //Clear closed set before calculating closure
    closureParameters.closed_set.clear();
    //Add generators to closure
    closureParameters.closed_set.insert(closureParameters.input_set.begin(), closureParameters.input_set.end());
    //Variable initialization
    std::deque<NodeId> newElements = std::deque<NodeId>(closureParameters.input_set.begin(), closureParameters.input_set.end());
    closureParameters.preclosure_steps = 0;
    NodeId LastElemOfPreviousStep = newElements.back();
    //Approximate coreIterations
    if (closureParameters.approximate && closureParameters.target_set_size != 0){
        bool proceedWhile = !newElements.empty();
        if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
            proceedWhile = false;
        }
        while (proceedWhile){
            NodeId bfsStart = newElements.front();
            if (bfsStart == LastElemOfPreviousStep){
                ++closureParameters.preclosure_steps;
                LastElemOfPreviousStep = newElements.back();
            }
            newElements.pop_front();
            //Graph is outerplanar or tree
            if (graph.graphType == GraphType::TREE || graph.graphType == GraphType::OUTERPLANAR){
                auto const it = generatedElements.find(bfsStart);
                if (it == generatedElements.end()){
                    bfs_forward(graph, closureParameters.closed_set, bfsStart, bfsQueue);
                    bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, newElements, generatedElements);
                    proceedWhile = !newElements.empty();
                    if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
                        proceedWhile = false;
                    }
                }
            }
            else {
                bfs_forward(graph, closureParameters.closed_set, bfsStart, bfsQueue);
                bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, newElements, generatedElements);
                proceedWhile = !newElements.empty();
                if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) >
                    closureParameters.lower_bound) {
                    proceedWhile = false;
                }
            }
            if (graph.graphType == GraphType::TREE){
                break;
            }
        }
    }
    else {
        int iterations = 0;
        bool break_condition;
        this->time = std::chrono::system_clock::now();
        set_break_condition(break_condition, newElements, closureParameters, iterations);
        while (!break_condition) {
            NodeId bfsStart = newElements.front();
            newElements.pop_front();
            if (graph.graphType == GraphType::OUTERPLANAR) {
                auto const it = generatedElements.find(bfsStart);
                if (it == generatedElements.end()) {
                    bfs_forward(graph, closureParameters.closed_set, bfsStart, bfsQueue);
                    bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, newElements,
                                 generatedElements);
                    ++iterations;
                }
            } else {
                bfs_forward(graph, closureParameters.closed_set, bfsStart, bfsQueue);
                bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, newElements,
                             generatedElements);
                ++iterations;
            }
            if (graph.graphType == GraphType::TREE) {
                break;
            }
            set_break_condition(break_condition, newElements, closureParameters, iterations);
        }
    }
}

void GraphClosureSP::fast_closure(GraphData &graph, ClosureParameters &closureParameters) {
    this->time = std::chrono::system_clock::now();
    closureParameters.clear();
    std::deque<NodeId> bfsQueue;
    if(closureParameters.timeConstraint > 0) {
        closureParameters.timeConstraint = (int) (closureParameters.timeConstraint * 1000000);
    }
    else{
        closureParameters.timeConstraint = -1;
    }
    size_t i = 0;
    std::unordered_set<NodeId> generated_elements;

    //incremental closure
    for (NodeId elem : closureParameters.input_set) {
        //if tree do only one bfs
        if (graph.graphType == GraphType::TREE){
            closureParameters.closed_set.insert(closureParameters.input_set.begin(), closureParameters.input_set.end());
            std::deque<NodeId> bfsElements;
            bfs_forward(graph, closureParameters.closed_set, elem, bfsQueue);
            bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, bfsElements, generated_elements);
            break;
        }
        else {
            if (i == 0) {
                closureParameters.closed_set.insert(elem);
            } else {
                // if outerplanar do not start bfs from already generated elements
                if (graph.graphType == GraphType::OUTERPLANAR &&
                    generated_elements.find(elem) == generated_elements.end()) {
                    addToClosed(graph, elem, closureParameters, bfsQueue, generated_elements);
                } else {
                    addToClosed(graph, elem, closureParameters, bfsQueue, generated_elements);
                }
            }
            ++i;
        }
    }

}



void GraphClosureSP::bfs_fast_forward(GraphData &graph, std::set<NodeId>& target_set, NodeId bfs_start,
                                      std::deque<NodeId> &bfsQueue) const {
    //Forward Search through graph O(m)
    //all elements in closed set
    target_set.insert(bfs_start);
    ++graph.Id;
    if (graph.Id == 0){
        std::fill(graph.ContainmentList.begin(), graph.ContainmentList.end(), 0);
        ++graph.Id;
    }
    for (NodeId elem : target_set) {
        graph.ContainmentList[elem] = graph.Id;
    }
    for (auto & x : graph.PredecessorsList) {
        x.clear();
    }
    std::fill(graph.DistanceList.begin(), graph.DistanceList.end(), -2);

    graph.DistanceList[bfs_start] = 0;
    graph.ContainmentList[bfs_start] = graph.Id;
    bfsQueue.push_back(bfs_start);
    int visitedFromClosed = 1;
    int possiblePaths = 1;
    int lastDistance = std::numeric_limits<int>::max();
    bool lastDistanceSet = false;
    while (!forward_step_break_condition(graph, bfsQueue, visitedFromClosed, (int) target_set.size(), possiblePaths)) {
        NodeId currentNodeId = bfsQueue.front();
        bfsQueue.pop_front();
        if (this->_threshold < graph.DistanceList[currentNodeId] || lastDistance < graph.DistanceList[currentNodeId]) {
            break;
        }
        bool currentNodeInClosed = graph.ContainmentList[currentNodeId] == graph.Id;
        bool currentNeighborInClosed;
        //Current node is already in closure and is not the first element
        if (currentNodeInClosed && currentNodeId != bfs_start){
            if (!lastDistanceSet) {
                TUNGraph::TNodeI Node = graph.get_graph()->GetNI(currentNodeId);
                for (int i = 0; i < Node.GetDeg(); ++i) {
                    NodeId neighborId = Node.GetNbrNId(i);
                    currentNeighborInClosed = graph.ContainmentList[neighborId] == graph.Id;
                    if (currentNeighborInClosed) {
                        if (graph.DistanceList[neighborId] < 0) {
                            ++visitedFromClosed;
                            int next_distance = graph.DistanceList[currentNodeId] + 1;
                            graph.DistanceList[neighborId] = next_distance;
                            bfsQueue.push_back(neighborId);
                            if (!lastDistanceSet && visitedFromClosed == target_set.size()) {
                                lastDistance = next_distance;
                                lastDistanceSet = graph.DistanceList[currentNodeId];
                            }
                        }
                        if (graph.DistanceList[neighborId] == graph.DistanceList[currentNodeId] + 1) {
                            graph.PredecessorsList[neighborId].push_back(currentNodeId);
                        }
                    }
                }
            }
        }
        //Current get_node is not in closure
        else{
            --possiblePaths;
            int next_distance = graph.DistanceList[currentNodeId] + 1;
            for (auto it = graph.begin(currentNodeId); it != graph.end(currentNodeId); ++it) {
                NodeId neighborId = *it;
                if (graph.DistanceList[neighborId] < 0) {
                    currentNeighborInClosed = graph.ContainmentList[neighborId] == graph.Id;
                    if (lastDistanceSet && !currentNeighborInClosed){
                        continue;
                    }
                    graph.DistanceList[neighborId] = next_distance;
                    bfsQueue.push_back(neighborId);
                    if (currentNeighborInClosed) {
                        ++visitedFromClosed;
                        if (!lastDistanceSet && visitedFromClosed == (int) target_set.size()){
                            lastDistance = graph.DistanceList[currentNodeId];
                            lastDistanceSet = true;
                        }
                    }
                    else {
                        ++possiblePaths;
                    }
                }
                if (graph.DistanceList[neighborId] == graph.DistanceList[currentNodeId] + 1) {
                    graph.PredecessorsList[neighborId].push_back(currentNodeId);
                }
            }
        }
    }
}

void GraphClosureSP::addToClosed(GraphData &graph, NodeId elem, ClosureParameters &closureParameters, std::deque<NodeId>& bfsQueue, std::unordered_set<NodeId>& generatedElements) {
    std::set<NodeId> addedInTheStep;
    std::deque<NodeId> unconsideredForBFS;

    bfs_fast_forward(graph, closureParameters.closed_set, elem, bfsQueue);
    bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, unconsideredForBFS, generatedElements);
    for (auto const & x: unconsideredForBFS) {
        addedInTheStep.insert(x);
    }
    bool proceedWhile = !unconsideredForBFS.empty() && closureParameters.incrementalCount != 0 && (closureParameters.timeConstraint == -1 || closureParameters.timeConstraint > std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count());
    if (closureParameters.approximate && closureParameters.target_set_size != 0){
        if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
            proceedWhile = false;
        }
    }
    while(proceedWhile){
        NodeId next_element = unconsideredForBFS.front();
        unconsideredForBFS.pop_front();

        bfs_forward(graph, closureParameters.closed_set, next_element, bfsQueue);
        bfs_backward(graph, closureParameters.closed_set, closureParameters, bfsQueue, unconsideredForBFS, generatedElements);

        if (graph.graphType == GraphType::TREE || graph.graphType == GraphType::OUTERPLANAR || closureParameters.onlyPreClosure){
            break;
        }
        --closureParameters.incrementalCount;
        ++closureParameters.iteration_number;
        proceedWhile = !unconsideredForBFS.empty() && closureParameters.incrementalCount != 0 && (closureParameters.timeConstraint == -1 || closureParameters.timeConstraint > std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count());
        if (closureParameters.approximate && closureParameters.target_set_size != 0){
            if (((double) closureParameters.closed_set.size() / closureParameters.target_set_size) > closureParameters.lower_bound){
                proceedWhile = false;
            }
        }
    }
}


void GraphClosureSP::approx_closure(const GraphData& graph, std::vector<GraphData> &samples, int& sample_number, ClosureParameters &closureParameters,
                                    std::vector<int> &approximation,  double& closure_time, bool load, ClosureFunction closureFunction) {
    std::fill(approximation.begin(), approximation.end(), 0);
    if (samples.empty() && load){
        auto start_closure = std::chrono::high_resolution_clock::now();
        std::vector<std::string> paths;
#pragma omp critical
        {
            for (const auto &entry: std::filesystem::recursive_directory_iterator("../out/Samples/" + graph.getName() + "/")) {
                paths.emplace_back(entry.path().string());
            }
        }
        sample_number = paths.size();
        for (auto const & path : paths) {
            GraphData sample = GraphData(path);
            closureParameters.closed_set.clear();
            if (closureFunction == nullptr) {
                this->naive_closure(sample, closureParameters);
            } else {
                closureFunction(sample, closureParameters);
            }
            for (NodeId id: closureParameters.closed_set) {
                approximation[id] += 1;
            }
        }
        closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
    }
    else {
        sample_number = samples.size();
        auto start_closure = std::chrono::high_resolution_clock::now();
        for (auto &sample: samples) {
            closureParameters.closed_set.clear();
            if (closureFunction == nullptr) {
                this->naive_closure(sample, closureParameters);
            } else {
                closureFunction(sample, closureParameters);
            }
            for (NodeId id: closureParameters.closed_set) {
                approximation[id] += 1;
            }
        }
        closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
    }
}

void GraphClosureSP::approx_closure(const GraphData& graph, std::vector<OuterplanarGraphData> &samples, int& sample_number, ClosureParameters &closureParameters,
                                    std::vector<int> &approximation,  double& closure_time, bool load, ClosureFunction closureFunction) {
    std::fill(approximation.begin(), approximation.end(), 0);
    if (samples.empty() && load){
        std::vector<std::string> paths;
#pragma omp critical
        {
            for (const auto &entry: std::filesystem::recursive_directory_iterator("../out/Samples/" + graph.getName() + "/")) {
                paths.emplace_back(entry.path().string());
            }
        }
        sample_number = paths.size();
        for (auto const & path : paths) {
            auto start_closure = std::chrono::high_resolution_clock::now();
            OuterplanarGraphData sample = OuterplanarGraphData(path);
            closureParameters.closed_set.clear();
            if (closureFunction == nullptr) {
                this->naive_closure(sample, closureParameters);
            } else {
                closureFunction(sample, closureParameters);
            }
            for (NodeId id: closureParameters.closed_set) {
                approximation[id] += 1;
            }
            closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
        }
    }
    else {
        sample_number = samples.size();
        auto start_closure = std::chrono::high_resolution_clock::now();
        for (auto &sample: samples) {
            closureParameters.closed_set.clear();
            if (closureFunction == nullptr) {
                this->naive_closure(sample, closureParameters);
            } else {
                closureFunction(sample, closureParameters);
            }
            for (NodeId id: closureParameters.closed_set) {
                approximation[id] += 1;
            }
        }
        closure_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start_closure).count() / 1000000.0);
    }
}

void GraphClosureSP::set_break_condition(bool &break_condition, std::deque<NodeId> &newElements,
                                         ClosureParameters &closureParameters, int iterations) {
    break_condition = newElements.empty();
    if (closureParameters.iteration_number > 0){
        break_condition = break_condition || iterations > closureParameters.iteration_number;
    }
    if (closureParameters.timeConstraint != -1){
        double time_count = (double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count();
        double max_time =  closureParameters.timeConstraint * 1000000;
        break_condition = break_condition || time_count > max_time ;
    }
}

void GraphClosureSP::get_generators(OuterplanarComponent& outerplanarComponent, std::set<NodeId> &input_set,
                                    std::set<NodeId> &generator_set) {
    generator_set.clear();
    std::vector<PUNGraph> faces;
    std::vector<std::vector<NodeId>> nodeToFaces;
    generator_set = input_set;
}

void
GraphClosureSP::outerplanar_closure(OuterplanarGraphData &outerplanarGraphData, ClosureParameters &closureParameters) {
    ClosureParameters bbTreeParameters;
    for (const auto & elem : closureParameters.input_set) {
        outerplanarGraphData.get_bb_tree_ids(elem, bbTreeParameters.input_set);
    }
    naive_closure(outerplanarGraphData.get_bbTree(), bbTreeParameters);
    int nodeId = 0;
    int componentId = 0;
    std::vector<std::set<int>> componentInput = std::vector<std::set<int>>(outerplanarGraphData.Components.size(), std::set<int>());
    std::vector<NodeOrComponent> nodeOrComponents;
    for (auto const & nodeId : closureParameters.input_set) {
        outerplanarGraphData.get_components(nodeId, nodeOrComponents);
        for (auto const &nodeOrComp: nodeOrComponents) {
            if (nodeOrComp.is_component(componentId)) {
                OuterplanarComponent &currentComponent = outerplanarGraphData.Components[componentId];
                componentInput[componentId].insert(currentComponent.NodeIdToComponentNodeId[nodeId]);
            }
        }
    }
    for (auto const & elem : bbTreeParameters.closed_set) {
        NodeOrComponent nodeComponent = outerplanarGraphData.get_bbNodeOrComponent(elem);
        if (nodeComponent.is_node(nodeId)){
            closureParameters.closed_set.insert(nodeId);
            outerplanarGraphData.get_components(nodeId, nodeOrComponents);
            for (auto const & nodeOrComp : nodeOrComponents) {
                if (nodeOrComp.is_component(componentId)) {
                    OuterplanarComponent &currentComponent = outerplanarGraphData.Components[componentId];
                    componentInput[componentId].insert(currentComponent.componentId(nodeId));
                }
            }
        }
    }

    for (int i = 0; i < outerplanarGraphData.Components.size(); ++i) {
        OuterplanarComponent& component = outerplanarGraphData.Components[i];
        ClosureParameters componentClosure;
        if (componentInput[i].size() > 1) {
            //TODO get_generators
            GraphClosureSP::get_generators(component, componentInput[i],
                                           componentClosure.input_set);
            naive_closure(component.component, componentClosure);
            for (auto const & elem : componentClosure.closed_set) {
                closureParameters.closed_set.insert(component.nodeId(elem));
            }
        }
    }
}










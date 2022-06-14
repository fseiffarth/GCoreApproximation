//
// Created by anonymous on 08.05.2021.
//

#include "drawing.h"
#include "../Data/GraphData.h"
#include <iostream>

void Draw::printSnap(GraphData& graph) {
    std::cout << "Edge List: " << std::endl;
    for (TUNGraph::TEdgeI Edge = graph.get_graph()->BegEI(); Edge != graph.get_graph()->EndEI(); Edge++) {
        std::cout << Edge.GetSrcNId() << " " << Edge.GetDstNId() << std::endl;
    }
}

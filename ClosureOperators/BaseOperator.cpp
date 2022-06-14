//
// Created by anonymous on 12.04.2021.
//

#include <string>
#include <iostream>
#include "BaseOperator.h"
#include "../Utils/StaticFunctions.h"

ClosureParameters::ClosureParameters(std::set<NodeId> &input_set, NodeId element_to_add) : input_set(input_set), closed_set(closed_set), element_to_add(element_to_add) {

}

ClosureParameters::ClosureParameters(std::set<NodeId> &input_set, std::vector<std::set<NodeId>*>& forbidden_elements, NodeId element_to_add) : input_set(input_set), closed_set(closed_set), forbidden_elements(forbidden_elements), element_to_add(element_to_add) {

}

void ClosureParameters::print() {
    std::cout << "Size: " << closed_set.size() << std::endl;
    std::cout << "Input Set: " << StaticFunctions::print<std::set<NodeId>, NodeId>(this->input_set) << std::endl;
    std::cout << "Closed Set: " << StaticFunctions::print<std::set<NodeId>, NodeId>(this->closed_set) << std::endl;
    std::cout << "Added Elements: " << StaticFunctions::print<std::set<NodeId>, NodeId>(this->added_elements) << std::endl;
}

void ClosureParameters::reset() {
    this->output_intersects_forbidden = false;
    this->output_is_all = false;
    lower_bound = 0.75;
    approximate = false;
    target_set_size = 0;
    preclosure_steps = 0;
}

void ClosureParameters::clear() {
    this->output_intersects_forbidden = false;
    this->output_is_all = false;
    this->onlyPreClosure = false;
    this->incrementalCount = -1;
    this->timeConstraint = -1;
    this->closed_set.clear();
    this->added_elements.clear();
    lower_bound = 0.75;
    approximate = false;
    target_set_size = 0;
    iteration_number = 0;
}

ClosureParameters::ClosureParameters() {

}

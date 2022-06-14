//
// Created by anonymous on 12.04.2021.
//

#ifndef CLOSURES_BASEOPERATOR_H
#define CLOSURES_BASEOPERATOR_H

#include <set>
#include <utility>
#include "../Data/Data.h"

class ClosureParameters{
public:
    ClosureParameters();

//!
    //! \param input_set reference to the input_set for which the closures should be computed
    //! \param element_to_add
    explicit ClosureParameters(std::set<NodeId>& input_set, NodeId element_to_add = -1);
    ClosureParameters(std::set<NodeId>& input_set, std::vector<std::set<NodeId>*>& forbidden_elements, NodeId element_to_add = -1);
    std::set<NodeId> input_set;
    std::set<NodeId> closed_set = std::set<NodeId>();
    std::vector<std::set<NodeId>*> forbidden_elements;
    std::set<NodeId> added_elements = std::set<NodeId>();
    bool output_is_all = false;
    bool output_intersects_forbidden = false;
    NodeId element_to_add = -1;
    bool onlyPreClosure = false;
    int incrementalCount = -1;
    double timeConstraint = -1;
    double lower_bound = 0.75;
    bool approximate = false;
    int iteration_number = 0;
    int target_set_size = 0;

    int preclosure_steps = 0;
    void print();
    void reset();


    void clear();


};

template <class T> class BaseOperator {
public:
    explicit BaseOperator(std::string name) : _name(std::move(name)){}

    BaseOperator();;
    std::string getName(){return _name;};

    virtual void closure(T& DataObject, ClosureParameters& ClosureOutput) = 0;


protected:
    [[nodiscard]] bool IsForbidden(NodeId Id, const std::vector<std::set<NodeId>*>& forbidden_elements) const;

private:
    std::string _name;

};

template<class T>
bool BaseOperator<T>::IsForbidden(NodeId Id, const std::vector<std::set<NodeId>*>& forbidden_elements) const {
    for (std::set<NodeId>* nodeSet : forbidden_elements) {
        if (nodeSet->find(Id) != nodeSet->end()){
            return true;
        }
    }
    return false;
}

template<class T>
BaseOperator<T>::BaseOperator() = default;


#endif //CLOSURES_BASEOPERATOR_H

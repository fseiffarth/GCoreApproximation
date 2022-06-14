//
// Created by anonymous on 12.04.2021.
//

#include <utility>

#include "Data.h"

template<typename T>
Data<T>::Data(const T& dataObject, std::string name) : dataObject(dataObject), _name(std::move(name)) {}

template<typename T>
Data<T>::Data(const T& dataObject, std::string name, Labels& labels) : dataObject(dataObject), _name(std::move(name)), _labels(labels){
    classNumber = 0;
    for (int i = 0; i < _labels.size(); ++i) {
        if (_labelMap.find(_labels[i]) == _labelMap.end()){
            _labelMap[_labels[i]] = Nodes{i};
            ++classNumber;
        }
        else{
            _labelMap[_labels[i]].emplace_back(i);
        }
    }
}


template<typename T>
const T& Data<T>::get_data() const {
    return dataObject;
}

template<typename T>
T & Data<T>::data() {
    return dataObject;
}

template<typename T>
std::string Data<T>::getName() const {
    return _name;
}

template<typename T>
void Data<T>::setTrainingNodes(int trainingSizePerClass, std::mt19937_64& gen) {
    this->trainingSet = std::vector<std::set<NodeId>>(this->classNumber, std::set<NodeId>());
    for (int i = 0; i < this->classNumber; ++i) {
        for (int j = 0; j < trainingSizePerClass; ++j) {
            int rand_num = std::uniform_int_distribution<int>(j, this->elementsOfLabel(i).size() - 1)(gen);
            trainingSet[i].insert(this->elementsOfLabel(i)[rand_num]);
            std::swap(this->elementsOfLabel(i)[j], this->elementsOfLabel(i)[rand_num]);
        }
    }
}

template<typename T>
void Data<T>::setName(const std::string& name) {this->_name = name;

}

template<typename T>
Data<T>::Data(const Data<T> &other) {
    this->dataObject = other.dataObject;
    this->_name = other._name;
    this->_labelMap = other._labelMap;
    this->_labels = other._labels;
    this->classNumber = other.classNumber;
}

template<typename T>
void Data<T>::set_data(T &data) {
dataObject = data;
}

template<typename T>
Data<T>::~Data() {
    _labelMap.clear();
    labels().clear();
}



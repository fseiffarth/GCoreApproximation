//
// Created by anonymous on 12.04.2021.
//

#ifndef CLOSURES_DATA_H
#define CLOSURES_DATA_H

#include "../Utils/typedefs.h"

template <typename T>
class Data {
public:
    Data()= default;
    ~Data();
    explicit Data(const T& dataObject, std::string name = "");
    Data(const T& dataObject, std::string name, Labels& labels);
    //Copy data
    explicit Data(const Data<T>& other);

    const T& get_data() const;
    T & data();
    void set_data(T& data);
    [[nodiscard]] virtual size_t size() const = 0;
    size_t class_num(){return classNumber;};
    virtual NodeId elem(NodeId Id) = 0;
    Nodes& elementsOfLabel(Label label) {return _labelMap[label];};
    [[nodiscard]] std::string getName() const;
    void setName(const std::string& name);
    Labels& labels() {return _labels;};
    [[nodiscard]] const Labels& get_labels() const {return _labels;};
    void setTrainingNodes(int trainingSizePerClass, std::mt19937_64& gen);
    std::vector<std::set<NodeId>> trainingSet;
protected:
    T dataObject;
    std::string _name;
    std::map<Label, Nodes> _labelMap;
    Labels _labels;
    size_t classNumber = 0;
};

#include "Data.txx"

#endif //CLOSURES_DATA_H


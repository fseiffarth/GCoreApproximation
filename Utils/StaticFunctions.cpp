//
// Created by anonymous on 13.04.2021.
//

#include <iostream>
#include "StaticFunctions.h"

std::set<NodeId> StaticFunctions::getNodesWithLabel(Labels& labels, Label targetLabel){
    std::set<NodeId> nodes = std::set<NodeId>();
    for (int i = 0; i < labels.size(); ++i) {
        if (labels[i] == targetLabel){
            nodes.insert(i);
        }
    }
    return nodes;
}

void StaticFunctions::saveValuesToFile(const std::string& path, const std::vector<std::string>& header, const std::vector<std::string>& values, std::_Ios_Openmode mode) {
    bool newFile = std::filesystem::exists(path);
    std::ofstream fs;
    fs.open(path, mode);
    // write the file headers

    if (!newFile) {
        for (const auto& string : header) {
            fs << "," << string;
        }
        fs << std::endl;
    }
    fs << std::fixed;
    for (const auto& string : values) {
        fs << "," << string;
    }
    fs << std::endl;
    fs << std::scientific;
    fs.close();
}

std::string StaticFunctions::printMap(const std::map<int, int> &map) {
    std::string out = "{";
    for (auto const & [key, value] : map) {
        out += "(";
        out += std::to_string(key);
        out += ":";
        out += std::to_string(value);
        out += ")";
    }
    out += "}";
    return out;
}

void StaticFunctions::headerValueInsert(std::vector<std::string> &header, std::vector<std::string> &values,
                                        const std::vector<std::string> &new_header,
                                        const std::vector<std::string> &new_values) {
    header.insert(header.end(), new_header.begin(), new_header.end());
    values.insert(values.end(), new_values.begin(), new_values.end());
}

void
StaticFunctions::generateInputSet(std::set<NodeId> &input_set, const GraphData &graphData, int input_size, int seed) {
    std::mt19937_64 generator(seed);
    input_set.clear();
    std::vector<NodeId> nodes(graphData.size());
    std::iota(nodes.begin(),  nodes.end(), 0);
    //std::cout << std::endl;
    //std::cout << "\tInput Set: ";
    for (int i = 0; i < input_size; ++i) {
        int rand_idx = std::uniform_int_distribution<int>(i, ((int) nodes.size()) - 1)(generator);
        input_set.insert(nodes[rand_idx]);
        std::swap(nodes[rand_idx], nodes[i]);
        //std::cout << nodes[i] << " ";
    }
    //std::cout << std::endl;
}

void StaticFunctions::printComponents(TCnComV &components) {
    std::cout << "Components: " << std::endl;
    for (int i= 0; i < components.Len(); ++i) {
        for (int j = 0; j < components[i].Len(); ++j) {
            std::cout << components[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void StaticFunctions::save(const std::string &path, const TIntV &NodeIds, const std::string& extension) {
    std::string edgePath = path + "." + extension;
    std::ofstream file;
    file.open(edgePath);
    for (int i = 0; i < NodeIds.Len(); ++i) {
        file << NodeIds[i] << std::endl;
    }
    file.close();
}

void StaticFunctions::save(const std::string &path,const std::vector<double> &values, const std::string& extension) {
    std::string complete_path = path + extension;

    if (extension.find("bin") != std::string::npos){
        std::ofstream file(complete_path, std::ios::out | std::ios::binary);
        size_t size = (values.size());
        file.write((char*) (&size), sizeof(size_t));
        for (auto const val: values) {
            file.write((char*) (&val), sizeof(double));
        }
        file.close();
    }
    else {
        std::ofstream file(complete_path, std::ios::out);
        for (auto const val: values) {
            file << val << std::endl;
        }
        file.close();
    }
}
void StaticFunctions::load(const std::string &path, std::vector<double> &values) {
    std::string extension = std::filesystem::path(path).extension();
    if (extension.find("bin") != std::string::npos){
        std::ifstream file(path, std::ios::in | std::ios::binary);
        size_t size = 0;
        file.read((char*) (&size), sizeof(size_t));
        values = std::vector<double>(size, 0);
        for (auto & val: values) {
            file.read((char*) (&val), sizeof(double));
        }
        file.close();
    }
    else {
        std::ifstream file(path, std::ios::in);
        std::string row, item;
        values.clear();
        while(std::getline(file, row)) {
            std::istringstream ss(row);
            while (std::getline(ss, item)) {
                values.emplace_back(std::stod(item));
            }
        }
        file.close();
    }
}

void StaticFunctions::load(const std::string& path, TIntV &NodeIds){
    if (std::filesystem::exists(path)) {
        int a;
        std::string line;
        std::ifstream infile(path);
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            iss >> a;
            NodeIds.Add(a);
        }
        infile.close();
    }
}

void StaticFunctions::load_csv(const std::string &path, std::vector<std::vector<std::string>>& out, const char& delimiter) {
    std::string row, item;
    std::ifstream in(path);
    std::vector<std::string> R;
    while(std::getline(in, row))
    {
        R.clear();
        std::istringstream ss(row);
        while (std::getline(ss, item, delimiter)) {
            R.push_back(item);
        }
        out.push_back( R );
    }
    in.close();
}

void StaticFunctions::PrintStream(std::stringstream& stringstream) {
    std::cout << stringstream.str();
}

void StaticFunctions::getClosedFromApproximation(std::vector<NodeId>& approximation, std::set<NodeId>& closedSet, int iterations, double threshold){
    closedSet.clear();
    for (int i = 0; i < approximation.size(); ++i) {
        bool closed_approximated = ((((double) approximation[i] / (double) iterations) / threshold) >= 1);
        if (closed_approximated){
            closedSet.emplace(i);
        }
    }
}

//
// Created by anonymous on 07.05.2021.
//

#ifndef DS2021_GRAPHIO_H
#define DS2021_GRAPHIO_H

#include "typedefs.h"
#include "../Data/GraphData.h"


static void removeCharsFromString(std::string &str,const char* charsToRemove ) {
    for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
        str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
    }
}

static void labels_from_gml(std::string input_path, std::string output_path){
    std::ifstream infile(input_path);
    std::string line;
    size_t counter = 0;
    std::ofstream out(output_path);
    std::string key;
    std::string value;
    size_t label;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> key >> value;
        if (key == "label"){
            removeCharsFromString(value, "\"");
            std::stringstream sstream(value);
            sstream >> label;
            out << counter << " " << label << std::endl;
            ++counter;
        }
    }
    out.close();
}

static bool LoadGraph(const std::string& Path, PUNGraph &Graph) {
    if (std::filesystem::path(Path).extension() == ".graph") {
        TFIn fin(Path.c_str());
        try {

            Graph = TUNGraph::Load(fin);
        }
        catch (...) {
            fin.Reset();
            Graph.Clr();
            return false;
        }
        return true;
    }
    return false;
}

static void LoadGraphsFromPath(const std::string& Path, std::vector<PUNGraph> &graphs,
                                  std::vector<std::string> &names, std::set<int>* graphSizes, int patternNum) {
    PUNGraph UndirectedGraph;
    std::map<int, int> sizesNumMap;
    if (graphSizes != nullptr) {
        for (int size : *graphSizes) {
            sizesNumMap.insert({size, 0});
        }
    }

    for (const auto &entry : std::filesystem::directory_iterator(Path)) {
        if (LoadGraph(entry.path().string(), UndirectedGraph)) {
            if (graphSizes == nullptr || graphSizes->find(UndirectedGraph->GetNodes()) != graphSizes->end()) {
                if (patternNum == -1 || graphSizes == nullptr || sizesNumMap[UndirectedGraph->GetNodes()] < patternNum) {
                    graphs.push_back(UndirectedGraph);
                    names.push_back(entry.path().stem().string());
                    if (graphSizes != nullptr) {
                        ++sizesNumMap[UndirectedGraph->GetNodes()];
                    }
                }
            }
        }
    }
}

static bool LoadLabels(const std::string& Path, Labels &LabelVector) {
    LabelVector.clear();
    if (std::filesystem::path(Path).extension() == ".labels") {
        try {
            std::ifstream infile(Path);
            NodeId id;
            Label label;
            while (infile >> id >> label){
                LabelVector.push_back(label);
            }
        }
        catch (...) {
            return false;
        }
        return true;
    }
    return false;
}

static void LoadLabelsFromPath(const std::string& Path, std::vector<Labels> &Labels, std::map<size_t, std::string> &LabelNames, std::set<int>* graphSizes, int patternNum) {
    std::vector<Label> LabelVector;
    std::map<int, int> sizesNumMap;
    if (graphSizes != nullptr) {
        for (int size : *graphSizes) {
            sizesNumMap.insert({size, 0});
        }
    }

    for (const auto &entry : std::filesystem::directory_iterator(Path)) {
        if (LoadLabels(entry.path().string(), LabelVector)) {
            LabelNames.insert(std::pair<size_t, std::string>(Labels.size(), entry.path().stem().string()));
            if (graphSizes == nullptr || graphSizes->find(LabelVector.size()) != graphSizes->end()) {
                if (patternNum == -1 || graphSizes == nullptr || sizesNumMap[LabelVector.size()] < patternNum) {
                    Labels.push_back(LabelVector);
                    if (graphSizes != nullptr) {
                        ++sizesNumMap[LabelVector.size()];
                    }
                }
            }

        }
    }
}

static void LoadGraphData(const std::string& path, std::vector<GraphData>& data){
    std::vector<std::pair<std::string, std::string>> graphLabelPaths;
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (std::filesystem::path(entry).extension() == ".edges") {
            std::pair<std::string, std::string> pair = {std::filesystem::path(entry).string(),std::filesystem::path(entry).parent_path().string() + "/" + std::filesystem::path(entry).stem().string() + ".labels"};
            graphLabelPaths.emplace_back(pair);
        }
    }
    for (const auto& [graphPath, labelPath] : graphLabelPaths) {
        data.emplace_back(GraphData(graphPath, labelPath));
    }
}

static void LoadGraphData(const std::string& path, std::vector<std::vector<GraphData>>& data){
    data.emplace_back();
    LoadGraphData(path, data.back());
}



struct LoadProperties{
public:
    explicit LoadProperties(const std::map<std::string, int> &propertyMap, int number);
    bool checkProperty(std::string path) const;
    std::map<std::string, int> propertyMap;
    int number;
};

static void LoadGraphData(const std::string& path, std::vector<GraphData>& data, const LoadProperties& properties){
    std::vector<std::pair<std::string, std::string>> graphLabelPaths;
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (std::filesystem::path(entry).extension() == ".edges" && properties.checkProperty(std::filesystem::path(entry))) {
            std::pair<std::string, std::string> pair = {std::filesystem::path(entry).string(),std::filesystem::path(entry).parent_path().string() + "/" + std::filesystem::path(entry).stem().string() + ".labels"};
            graphLabelPaths.emplace_back(pair);
        }
    }
    int counter = 0;
    for (const auto& [graphPath, labelPath] : graphLabelPaths) {
        if (properties.number == -1 || properties.number > counter) {
            data.emplace_back(GraphData(graphPath, labelPath));
        }
        ++counter;
    }
}

static void LoadGraphData(const std::string& path, std::vector<std::vector<GraphData>>& data, const LoadProperties& properties){
    data.emplace_back();
    LoadGraphData(path, data.back(), properties);
}



#endif //DS2021_GRAPHIO_H

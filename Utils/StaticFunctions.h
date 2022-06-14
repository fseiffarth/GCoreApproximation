//
// Created by anonymous on 13.04.2021.
//

#ifndef CLOSURES_STATICFUNCTIONS_H
#define CLOSURES_STATICFUNCTIONS_H


#include <string>
#include <numeric>
#include "typedefs.h"
#include "../Data/GraphData.h"
#include "iostream"

class StaticFunctions {
public:
    template<class T>
    static T mean(const std::vector<T>& vector);

    template<class T>
    static T standard_deviation(const std::vector<T>& vector);

    template<class T1, class T2>
    static std::string print(T1 &Object);

    template<class T1, class T2>
    static std::string print(T1 &Object, bool pair);

    template<class T2>
    static std::string pairToString(const T2 &object);

    static std::string printMap(const std::map<int, int> &map);

    static std::set<NodeId> getNodesWithLabel(Labels& labels, Label targetLabel);

    static void saveValuesToFile(const std::string& path, const std::vector<std::string>& header, const std::vector<std::string>& values, std::_Ios_Openmode mode = std::ios_base::app);
    static void headerValueInsert(std::vector<std::string>& header, std::vector<std::string>&values, const std::vector<std::string>& new_header, const std::vector<std::string>& new_values);
    static void generateInputSet(std::set<NodeId>& input_set, const GraphData& graphData, int input_size, int seed);
    static void printComponents(TCnComV& components);
    static void save(const std::string& path,const TIntV& NodeIds, const std::string& extension = "core");

    static void load(const std::string &path, TIntV &NodeIds);
    static void load_csv(const std::string &path, std::vector<std::vector<std::string>>& out, const char& delimiter = ',');

    // Comparison function to sort the vector elements
    // by second element of tuples
    static bool sortbysecond(const std::tuple<std::string, int>& a,
                      const std::tuple<std::string, int>& b)
    {
        return (std::get<1>(a) < std::get<1>(b));
    }

    static void save(const std::string &path, const std::vector<double> &values, const std::string &extension);

    static void load(const std::string &path, std::vector<double> &values);

    static void PrintStream(std::stringstream& stringstream);

    static void getClosedFromApproximation(std::vector<NodeId> &approximation, std::set<NodeId> &closedSet, int iterations,
                                           double threshold);
};

template <class T1, class T2>
std::string StaticFunctions::print(T1& Object) {
    return std::accumulate(std::begin(Object),
                           std::end(Object),
                           std::string{},
                           [](const std::string &a, const T2 &b) {
                               return a.empty() ? '"' + std::to_string(b)
                                                : a + "|" + std::to_string(b) + '"';
                           });
}

template<class T1, class T2>
std::string StaticFunctions::print(T1 &Object, bool pair) {
    if (pair) {
        return std::accumulate(std::begin(Object),
                               std::end(Object),
                               std::string{},
                               [](const std::string &a, const T2 &b) {
                                   return a.empty() ? '"' + pairToString(b) + '"'
                                                    : a + ", " + '"' + pairToString(b) + '"';
                               });
    }
    return "";
}

template<class T2>
std::string StaticFunctions::pairToString(const T2 &object) {
    return "(" + std::to_string(object.first) + "," + std::to_string(object.second) + ")";
}

template<class T>
T StaticFunctions::mean(const std::vector<T> &vector) {
    return std::accumulate(vector.begin(),  vector.end(), 0.0)/(double) vector.size();
}

template<class T>
T StaticFunctions::standard_deviation(const std::vector<T> &vector) {
    T m = StaticFunctions::mean(vector);
    double accum = 0.0;
    std::for_each (vector.begin(),  vector.end(), [&](const double d) {
        accum += (d - m) * (d - m);
    });
    return std::sqrt(accum / (vector.size()-1));
}


#endif //CLOSURES_STATICFUNCTIONS_H

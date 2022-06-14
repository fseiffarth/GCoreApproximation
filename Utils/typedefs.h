//
// Created by anonymous on 27.04.2021.
//

#ifndef HOPS_TYPEDEFS_H
#define HOPS_TYPEDEFS_H

#include <omp.h>
#include <Snap.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <set>
#include <map>
#include <random>
#include <filesystem>
#include <chrono>
#include <vector>

typedef int Label;
typedef std::vector<Label> Labels;
typedef int NodeId;
typedef std::vector<NodeId> Nodes;
typedef std::vector<std::pair<NodeId, NodeId>> EdgesSrcDst;
typedef unsigned long long int UInt64;
typedef std::vector<std::vector<int>> SIMPLE_GRAPH;

#endif //HOPS_TYPEDEFS_H

cmake_minimum_required(VERSION 3.16)
project(Closures)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(../../libs/Snap-6.0/snap-core)
link_directories(../../libs/Snap-6.0/snap-core)
include_directories(../../libs/Snap-6.0/glib-core)
link_directories(../../libs/Snap-6.0/glib-core)


#Tests
add_executable(Test test.cpp ../ClosureOperators/BaseOperator.cpp ../ClosureOperators/BaseOperator.h ../Data/Data.txx ../Data/Data.h ../ClosureOperators/GraphClosures.cpp ../ClosureOperators/GraphClosures.h ../Data/GraphData.cpp ../Data/GraphData.h ../../libs/Snap-6.0/snap-core/Snap.cpp ../Utils/StaticFunctions.cpp ../Utils/StaticFunctions.h ../Utils/graphio.h ../Utils/graphio.cpp ../Utils/drawing.h ../Utils/drawing.cpp ../Utils/graph_generators.h ../Utils/graph_generators.cpp ../Utils/subgraphs.h ../Utils/subgraphs.cpp ../Utils/typedefs.h ../Utils/Generators.h ../Utils/Generators.cpp ../Utils/pattern_generator.txx ../Utils/pattern_generator.h ../Utils/GraphFunctions.cpp ../Utils/GraphFunctions.h ../Utils/GraphStructs.cpp ../Utils/GraphStructs.h ../Utils/Graphs.cpp ../Utils/Graphs.h ../Utils/subgraphs.txx ../Utils/GraphFunctions.txx  ../Utils/directed_graph.h ../Utils/directed_graph.cpp  ../Utils/OuterplanarSubgraphDFS.cpp ../Utils/OuterplanarSubgraphDFS.h   ../Utils/OuterplanarSubgraph.cpp ../Utils/OuterplanarSubgraph.h ../Utils/Enums.h  ../Utils/DataClasses.h ../Utils/DataClasses.cpp   ../Utils/FileEvaluation.h ../Utils/FileEvaluation.cpp ../Data/OuterplanarGraphData.h ../Data/OuterplanarGraphData.cpp ../Utils/CoreGrowAlgorithm.cpp ../Utils/CoreGrowAlgorithm.h)

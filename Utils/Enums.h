//
// Created by anonymous on 12.08.21.
//

#ifndef CLOSURES_ENUMS_H
#define CLOSURES_ENUMS_H

enum RootNodeCondition {
    TREE_GIVEN,
    NODE_GIVEN,
    MAX_DEGREE,
    MIN_DEGREE,
    MIN_LABEL_BIG_GRAPH,
};

enum LabelType {
    UNLABELED,
    SPARSE_LABELED,
    LABELED,
};

enum class PatternType{
    BFS_TREE,
    OUTERPLANAR,
};

#endif //CLOSURES_ENUMS_H

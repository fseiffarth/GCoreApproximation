//
// Created by anonymous on 10.05.2021.
//

#ifndef CLOSURES_GRAPHSTRUCTS_H
#define CLOSURES_GRAPHSTRUCTS_H


#include <algorithm>
#include <iostream>


struct NodePair{
    NodePair() = default;
    NodePair(int src, int dst){
        source = std::min(src, dst);
        destination = std::max(src, dst);
    };
    bool operator == (const NodePair& nodePair) const{
        return source == nodePair.source && destination == nodePair.destination;
    }
    bool operator < (const NodePair& nodePair) const{
        return source < nodePair.source || source == nodePair.source && destination < nodePair.destination;
    }
    [[nodiscard]] int first() const{return source;};
    [[nodiscard]] int second() const{return destination;};

    void print() const{
        std::cout << " " << source << "-" << destination << " ";
    }

private:
    int source;
    int destination;
};

// The specialized hash function for `unordered_map` keys
struct hashNodePair
{
    std::size_t operator() (const NodePair &nodePair) const
    {
        std::size_t h1 = std::hash<int>()(nodePair.first());
        std::size_t h2 = std::hash<int>()(nodePair.second());
        return h1 ^ h2;
    }
};

#endif //CLOSURES_GRAPHSTRUCTS_H

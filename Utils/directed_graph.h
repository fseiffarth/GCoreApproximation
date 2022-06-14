/*
Copyright (c) 2020 Paul Stahr

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <iostream>
#include <vector>

namespace DIRGRAPH
{

template <typename T, typename Meaning>
struct Explicit
{
    Explicit(){ }
    Explicit(T value) : value(value){ }
    inline operator T () const { return value; }

    T value;
};

typedef Explicit<std::size_t, struct node_id_t> NodeId;
typedef Explicit<std::size_t, struct edge_id_t> EdgeId;

class Graph;
class Node
{
    friend Graph;
public:
    Node();
    std::vector<EdgeId> const & in_edges() const;
    std::vector<EdgeId> const & out_edges() const;
    int out_degree() const;
    int in_degree() const;
    bool is_valid() const;
private:
    void add_incoming_edge(EdgeId);
    void add_outgoing_edge(EdgeId);
    std::vector<EdgeId> _incoming;
    std::vector<EdgeId> _outgoing;
    bool _is_valid;
};

class Edge
{
    friend Graph;
public:
    Edge();
    Edge(NodeId tail, NodeId head);

    NodeId get_tail() const;
    NodeId get_head() const;
    bool is_valid() const;
private:
    NodeId _tail;
    NodeId _head;
    bool _is_valid;
};

class Graph
{
public:
    explicit Graph(size_t num_nodes);
    Graph();

    EdgeId add_edge(NodeId tail, NodeId head);
    NodeId add_node();
    NodeId add_nodes(size_t num_nodes);

    void remove_node(NodeId node_id);
    void remove_edge(EdgeId edge_id);

    void reasign_ids(std::vector<NodeId> & node_id_old_to_new, std::vector<EdgeId> & edge_id_old_to_new);

    size_t num_nodes() const;
    size_t num_edges() const;

    std::vector<Node>::const_iterator node_cbegin() const;
    std::vector<Node>::const_iterator node_cend() const;

    std::vector<Edge>::const_iterator edge_cbegin() const;
    std::vector<Edge>::const_iterator edge_cend() const;

    std::vector<Node>::iterator node_begin();
    std::vector<Node>::iterator node_end();

    std::vector<Edge>::iterator edge_begin();
    std::vector<Edge>::iterator edge_end();


    Node & get_node(NodeId);
    Node const & get_node(NodeId) const;
    int get_node_id(std::vector<Node>::const_iterator node_iterator) const ;

    Edge & get_edge(EdgeId);
    Edge const & get_edge(EdgeId) const;

    bool neighbor(NodeId nodeId, int Idx, NodeId& neighbor) const;
    bool neighbor(const Node& node, int Idx, NodeId& neighbor) const;

    static const EdgeId invalid_edge;
    static const NodeId invalid_node;

    void clean();
    void reassign_clean(std::vector<NodeId> & node_ref, std::vector<EdgeId> & edge_ref);

    void clear();

    template <typename function>void for_all_nodes(function f){std::for_each(node_begin(), node_end(), f);}
    template <typename function>void for_all_edges(function f){std::for_each(edge_begin(), edge_end(), f);}

private:
      std::vector<Node> _nodes;
      std::vector<Edge> _edges;
};

void print_graph(Graph);
Graph read_graph(char const * filename);
}
#endif

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

#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include "directed_graph.h"

namespace DIRGRAPH
{
EdgeId const Graph::invalid_edge = std::numeric_limits<EdgeId>::max();
NodeId const Graph::invalid_node = std::numeric_limits<NodeId>::max();

Graph::Graph(){}

Graph::Graph(size_t num_nodes)
 : _nodes(num_nodes){}

NodeId Graph::add_node()
{
      _nodes.push_back(Node());
      return _nodes.size() - 1;
}

NodeId Graph::add_nodes(size_t num_nodes)
{
    _nodes.reserve(_nodes.size() + num_nodes);
    std::fill_n(std::back_inserter(_nodes), num_nodes, Node());
    return _nodes.size() - 1;
}

EdgeId Graph::add_edge(NodeId tail, NodeId head)
{
    if ((tail >= num_nodes()) or (head >= num_nodes()))
    {
       throw std::runtime_error("Edge cannot be added due to undefined endpoint.");
       return invalid_edge;
    }
    else
    {
        _edges.push_back(Edge(tail, head));
        EdgeId edge_id = _edges.size() - 1;
        _nodes[tail].add_outgoing_edge(edge_id);
        _nodes[head].add_incoming_edge(edge_id);
        return edge_id;
    }
}

void Node::add_incoming_edge(EdgeId edge){_incoming.push_back(edge);}
void Node::add_outgoing_edge(EdgeId edge){_outgoing.push_back(edge);}

std::vector<Node>::const_iterator Graph::node_cbegin()  const{return _nodes.cbegin();}
std::vector<Node>::const_iterator Graph::node_cend()    const{return _nodes.cend();}
std::vector<Edge>::const_iterator Graph::edge_cbegin()  const{return _edges.cbegin();}
std::vector<Edge>::const_iterator Graph::edge_cend()    const{return _edges.cend();}

std::vector<Node>::iterator Graph::node_begin()              {return _nodes.begin();}
std::vector<Node>::iterator Graph::node_end()                {return _nodes.end();}
std::vector<Edge>::iterator Graph::edge_begin()              {return _edges.begin();}
std::vector<Edge>::iterator Graph::edge_end()                {return _edges.end();}

void Graph::clean()
{
    for_all_nodes([this](Node & item)
    {
        if (!item.is_valid())
        {
            item._incoming.clear();
            item._outgoing.clear();
        }
        auto newend = std::remove_if(item._incoming.begin(), item._incoming.end(),  [this](EdgeId eId){return !this->_edges[eId].is_valid();});
        item._incoming.resize(std::distance(item._incoming.begin(), newend));
        newend = std::remove_if(item._outgoing.begin(), item._outgoing.end(),  [this](EdgeId eId){return !this->_edges[eId].is_valid();});
        item._outgoing.resize(std::distance(item._outgoing.begin(), newend));
    });
    for_all_edges([this](Edge & item)
    {
        if (!item.is_valid())
        {
            item._head = invalid_node;
            item._tail = invalid_node;
        }
    });
}

void Graph::reassign_clean(std::vector<NodeId> & node_ref, std::vector<EdgeId> & edge_ref)
{
    clean();

    auto ebegin = edge_begin();
    auto nbegin = node_begin();
    size_t edge_write_index = 0;
    for (size_t read_index = 0; ebegin + read_index != edge_end(); ++read_index)
    {
        if (ebegin[read_index].is_valid())
        {
            ebegin[edge_write_index] = ebegin[read_index];
            edge_ref.push_back(edge_write_index++);
        }
        else
        {
            edge_ref.push_back(invalid_edge);
        }
    }
    size_t node_write_index = 0;
    for (size_t read_index = 0; nbegin + read_index != node_end(); ++read_index)
    {
        if (nbegin[read_index].is_valid())
        {
            nbegin[node_write_index] = nbegin[read_index];
            node_ref.push_back(node_write_index++);
        }
        else
        {
            node_ref.push_back(invalid_node);
        }
    }
}


Node::Node() : _is_valid(true){}

Edge::Edge()
 : _tail(Graph::invalid_node),
   _head(Graph::invalid_node),
   _is_valid(true)
{}

Edge::Edge(NodeId tail, NodeId head)
 : _tail(tail),
   _head(head),
   _is_valid(true)
{}

bool Node::is_valid() const{return _is_valid;}
bool Edge::is_valid() const{return _is_valid;}

NodeId Edge::get_tail() const {return _tail;}
NodeId Edge::get_head() const {return _head;}

std::vector<EdgeId> const & Node::in_edges()  const{return _incoming;}
std::vector<EdgeId> const & Node::out_edges() const{return _outgoing;}

int Node::out_degree() const {
        return (int) out_edges().size();
}
int Node::in_degree() const {
    return (int) out_edges().size();
}



    size_t Graph::num_nodes() const{return _nodes.size();}
size_t Graph::num_edges() const{return _edges.size();}

Node        & Graph::get_node(NodeId node)      {return _nodes[node];}
Node const  & Graph::get_node(NodeId node) const{return _nodes[node];}
Edge        & Graph::get_edge(EdgeId edge)      {return _edges[edge];}
Edge const  & Graph::get_edge(EdgeId edge) const{return _edges[edge];}

void Graph::clear()
{
    _nodes.clear();
    _edges.clear();
}

void print_graph(Graph g)
{
      for (std::size_t v = 0; v < g.num_nodes(); ++v)
      {
            std::cout << "The following edge(s) leave(s) vertex " << v << ":\n";
            for (std::size_t e = 0; e < g.get_node(v).out_edges().size(); ++e)
            {
                  std::cout << g.get_edge(g.get_node(v).out_edges()[e]).get_tail()
                            << " -> "
                            << g.get_edge(g.get_node(v).out_edges()[e]).get_head()
                            << "\n";
            }
      }
}

void Graph::remove_node(NodeId node){_nodes[node]._is_valid = false;}
void Graph::remove_edge(EdgeId edge){_edges[edge]._is_valid = false;}

    int Graph::get_node_id(std::vector<Node>::const_iterator node) const {
        return std::distance(_nodes.begin(), node);
    }

    bool Graph::neighbor(NodeId nodeId, int Idx, NodeId& neighbor) const {
        neighbor = -1;
        int counter = -1;
        for (auto edgeId : this->get_node(nodeId).out_edges()) {
            auto edge = this->get_edge(edgeId);
            if (edge.is_valid()){
                ++counter;
            }
            if (counter == Idx){
                neighbor = edge.get_head();
                return true;
            }
        }
        for (auto edgeId : this->get_node(nodeId).in_edges()) {
            auto edge = this->get_edge(edgeId);
            if (edge.is_valid()){
                ++counter;
            }
            if (counter == Idx){
                neighbor = edge.get_tail();
                return true;
            }
        }
        return false;
    }

    bool Graph::neighbor(const Node& node, int Idx, NodeId& neighbor) const {
        if (node.out_edges().size() <= Idx || Idx < 0)
        {
            return false;
        }
        auto edge = this->get_edge(node.out_edges()[Idx]);
        if (!edge.is_valid()){
            return false;
        }
        neighbor = edge._head;
        return true;
    }

    Graph read_graph(char const * filename)
{
      std::fstream file(filename);
      if (not file)
      {
            throw std::runtime_error("Cannot open file.");
      }

      size_t num_nodes = 0;
      std::string line;
      std::getline(file, line);
      std::stringstream ss(line);
      ss >> num_nodes;
      if (not ss)
      {
            throw std::runtime_error("Invalid file format.");
            return Graph();
      }
      Graph graph(num_nodes);

      while (std::getline(file, line))
      {
            size_t head = Graph::invalid_node, tail = Graph::invalid_node;
            std::stringstream ss(line);
            ss >> tail;
            if (not ss)
            {
                  throw std::runtime_error("Invalid file format.");
                  return Graph();
            }
            ss >> head;
            if (not ss)
            {
                  throw std::runtime_error("Invalid file format.");
                  return Graph();
            }
            if (tail != head)
            {
                  graph.add_edge(NodeId(tail), NodeId(head));
            }
            //else throw std::runtime_error("Invalid file format: loops not allowed.");
      }
      return graph;
}
}

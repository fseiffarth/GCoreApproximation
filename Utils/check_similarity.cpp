#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include "directed_graph.h"

/*void calculate_reachable_vertices(DIRGRAPH::Graph const & g, std::vector<std::vector<size_t> > const & reachable_vertices)
{
    std::vector<size_t
}*/

struct breadth_search
{
    std::vector<DIRGRAPH::NodeId> _reached_vertices;
    std::vector<bool> _node_reached;
    size_t _read;
    std::vector<DIRGRAPH::EdgeId > _reached_edges;
    DIRGRAPH::Graph const & _g;
    breadth_search(DIRGRAPH::Graph const & g_) :  _node_reached(g_.num_nodes(), false), _read(0), _g(g_)
    {
    }

    void add_reached_vertex(DIRGRAPH::NodeId next_node_id)
    {
        _reached_vertices.push_back(next_node_id);
        _node_reached[next_node_id] = true;
    }

    void operator()()
    {
        size_t oldsize = _reached_vertices.size();
        for (; _read != oldsize ; ++_read)
        {
            DIRGRAPH::NodeId node_id = _reached_vertices[_read];
            for (DIRGRAPH::EdgeId edge_id : _g.get_node(node_id).out_edges())
            {
                DIRGRAPH::NodeId next_node_id = _g.get_edge(edge_id).get_head();
                if (!_node_reached[next_node_id])
                {
                    _reached_vertices.push_back(next_node_id);
                    _node_reached[next_node_id] = true;
                }
            }
        }

    }
};



uint8_t graphs_similar(std::vector<DIRGRAPH::Graph *> const & graphs)
{
    size_t num_graphs = graphs.size();
    if (num_graphs < 2)
    {
        return 1;
    }
    size_t num_nodes = graphs[0]->num_nodes();
    size_t num_edges = graphs[0]->num_edges();
    for (size_t i = 0; i < graphs.size(); ++i)
    {
        if (graphs[i]->num_nodes() != num_nodes || graphs[i]->num_edges() != num_edges)
        {
            return - 1;
        }
    }
    std::vector<std::vector<size_t > > node_to_visited_count(num_graphs, std::vector<size_t>(num_nodes));

    std::vector<std::vector<size_t> > node_to_id(num_graphs, std::vector<size_t>(num_nodes, 0));
    std::vector<std::vector<breadth_search> > b;
    std::vector<std::vector<size_t> > permutation_to_node(num_graphs, std::vector<size_t>(num_nodes, 0));
    for (size_t i = 0; i < num_graphs; ++i)
    {
        b.emplace_back(std::vector<breadth_search>(num_nodes, *graphs[i]));
        std::iota(permutation_to_node[i].begin(), permutation_to_node[i].end(), 0);
        for (size_t j = 0; j < num_nodes; ++j)
        {
            b.back()[j].add_reached_vertex(j);
            std::cout <<"  " <<b.back()[j]._reached_vertices.size() << std::endl;
        }

    }

    bool terminate = false;
    do
    {
        terminate = true;
        for (size_t gi= 0; gi < num_graphs; ++gi)
        {
            std::vector<size_t> & current_node_to_visited_count = node_to_visited_count[gi];
            std::cout << "graph " << gi << std::endl;
            for (size_t i = 0; i < num_nodes; ++i)
            {
                size_t old_visited_count = b[gi][i]._reached_vertices.size();
                b[gi][i]();
                current_node_to_visited_count[i] = b[gi][i]._reached_vertices.size();
                std::cout << ' ' << old_visited_count << "-> " << current_node_to_visited_count[i] << std::endl;
                if (current_node_to_visited_count[i] != old_visited_count)
                {
                    terminate = false;
                }
            }

            size_t begin = 0;
            size_t current_id = 0;
            std::vector<size_t> & current_node_to_id = node_to_id[gi];
            std::vector<size_t> & current_permutation_to_node = permutation_to_node[gi];
            for (size_t end = 0; end < num_nodes; ++end)
            {
                while (end < num_nodes && current_node_to_id[begin] == current_node_to_id[end])
                {
                    ++end;
                }
                std::stable_sort(current_permutation_to_node.begin() + begin, current_permutation_to_node.begin() + end, [&current_node_to_visited_count](size_t lhs, size_t rhs){return current_node_to_visited_count[lhs] < current_node_to_visited_count[rhs];});
                size_t current_node_count = std::numeric_limits<size_t>::max();
                while (begin < end)
                {
                    size_t current_node = current_permutation_to_node[begin];
                    current_node_to_id[current_node] = current_id;
                    if (current_node_to_visited_count[current_node] != current_node_count)
                    {
                        current_node_count = current_node_to_visited_count[current_node];
                        ++current_id;
                    }
                    ++begin;
                }
            }
        }
        for (size_t gi = 1; gi < num_graphs; ++gi)
        {
            std::cout << node_to_id[0][0] << std::endl;
            std::vector<size_t> & last_node_to_visited_count = node_to_visited_count[gi - 1];
            std::vector<size_t> & last_permutation_to_node = permutation_to_node[gi - 1];
            std::vector<size_t> & current_node_to_visited_count = node_to_visited_count[gi];
            std::vector<size_t> & current_permutation_to_node = permutation_to_node[gi];
            if (!std::equal(node_to_id[0].begin(), node_to_id[0].end(), node_to_id[gi].begin()))
            {
                return -1;
            }

            for (size_t i = 0; i < num_nodes; ++i)
            {
                if (current_node_to_visited_count[current_permutation_to_node[i]] != last_node_to_visited_count[last_permutation_to_node[i]])
                {
                    return -1;
                }
            }

        }
    }while(!terminate);
    return 0;
}

std::istream & read_graph(std::istream & stream, DIRGRAPH::Graph & g)
{
    g.clear();
    if (! stream)
    {
        throw std::runtime_error("Cannot open file");
    }

    size_t num_nodes = std::numeric_limits<size_t>::max();
    std::string line;
    std::getline(stream, line);
    std::istringstream ss(line);
    if (! (ss >> num_nodes))
    {
        throw std::runtime_error("Invalid file format");
    }
    if (num_nodes == std::numeric_limits<size_t>::max())
    {
        throw std::runtime_error("num nodes invalid");
    }

    g.add_nodes(num_nodes);
    while (std::getline(stream, line))
    {
        ss.str(line);
        ss.clear();
        try
        {
            size_t head, tail;
            ss >> head;
            ss >> tail;
            g.add_edge(head, tail);
        }
        catch (std::runtime_error const & e)
        {
            throw std::runtime_error("line \"" + line + "\" not readable, error:" + std::string(e.what()));
        }
    }
    return stream;
}


int main(int argc, char *argv[])
{
    std::vector<DIRGRAPH::Graph> graphs(argc - 1);
    std::vector<DIRGRAPH::Graph *> graphsp(argc - 1);
    for (int i = 1; i < argc; ++i)
    {
        std::fstream input(argv[i]);
        read_graph(input, graphs[i-1]);
        input.close();
        graphsp[i-1] = &graphs[i- 1];
    }
    std::cout << size_t(graphs_similar(graphsp)) << std::endl;
    return 0;
}

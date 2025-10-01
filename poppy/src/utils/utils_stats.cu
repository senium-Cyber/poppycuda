#include "include/utils_stats.cuh"
#include "include/utils.cuh"
/**
 * @brief Get a adjacency matrix of the number of edges between each partition
 *
 * @param g
 * @param partition
 * @param directed
 * @return std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> get_edges_stats(Graph g, std::map<int, int> partition, bool directed) {
    int number_partition =
        std::max_element(
            partition.begin(),
            partition.end(),
            [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })
            .operator*()
            .second;
    std::vector<std::vector<int>> stats(number_partition + 1, std::vector<int>(number_partition + 1, 0));
    if (directed) {
        for (auto arc : g.arcs) {
            stats[partition[arc.first]][partition[arc.second]]++;
        }
    } else {
        for (auto arc : g.arcs) {
            stats[partition[arc.first]][partition[arc.second]]++;
            stats[partition[arc.second]][partition[arc.first]]++;
        }
    }
    return stats;
}

/**
 * @brief Get a adjacency matrix of the number of edges between each partition
 *
 * @param graphs
 * @param partition
 * @param directed
 * @return std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> get_edges_stats(std::vector<Graph> graphs, std::map<int, int> partition, bool directed) {
    int number_partition =
        std::max_element(
            partition.begin(),
            partition.end(),
            [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })
            .operator*()
            .second;
    std::vector<std::vector<int>> stats(number_partition + 1, std::vector<int>(number_partition + 1, 0));
    if (directed) {
        for (int i = 0; i < graphs.size(); i++) {
            for (auto arc : graphs[i].arcs) {
                if (partition[arc.first] == i) {
                    stats[partition[arc.first]][partition[arc.second]]++;
                }
            }
        }
    } else {
        for (int i = 0; i < graphs.size(); i++) {
            for (auto arc : graphs[i].arcs) {
                if (partition[arc.first] == i) {
                    stats[partition[arc.first]][partition[arc.second]]++;
                    stats[partition[arc.second]][partition[arc.first]]++;
                }
            }
        }
    }
    return stats;
}

/**
 * @brief Average number of inter DC communication by synchro
 * 1 communication = 1 cipher text sent to another DC
 *
 * @param stats
 * @param batchSize
 * @return double
 */
double average_number_inter_DC_communication(std::vector<Graph> graphs, std::map<Node, int> partition, int batchSize) {
    std::vector<double> average;
    for (int i = 0; i < graphs.size(); i++) {
        double sum = 0;
        double count = 0;
        std::map<Node, Node> map_out;
        for (auto arc : graphs[i].arcs) {
            map_out.insert(arc);
        }
        for (int j = 0; j < graphs.size(); j++) {
            if (i != j) {
                int stats = 0;
                // Count the number of out nodes to j from j
                for (auto node : graphs[i].nodes) {
                    auto range = map_out.equal_range(node);
                    for (auto it = range.first; it != range.second; ++it) {
                        if (partition[it->second] == j) {
                            stats++;
                            break;  // Count once
                        }
                    }
                }
                if (stats > 0) {
                    sum += (stats / batchSize) + 1;
                    count++;
                }
            }
        }
        average.push_back((double)sum / (double)count);
    }
    return std::accumulate(average.begin(), average.end(), 0.0) / average.size();
}

/**
 * @brief Get the percentage of nodes in each partition
 *
 * @param partition
 * @return std::vector<double>
 */
std::vector<double> distribution_nodes(std::map<int, int> partition) {
    int number_partition =
        std::max_element(
            partition.begin(),
            partition.end(),
            [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })
            .operator*()
            .second;
    std::vector<double> stats(number_partition + 1, 0);
    for (auto node : partition)
        stats[node.second]++;
    for (auto& stat : stats)
        stat /= partition.size();

    return stats;
}

/**
 * @brief Get the number of dangling nodes in each graph
 *
 * @param graphs
 * @param partition
 * @return std::vector<int>
 */
std::vector<int> number_dangling_nodes(std::vector<Graph> graphs, std::map<int, int> partition) {
    std::vector<int> stats(graphs.size(), 0);
    for (int i = 0; i < graphs.size(); i++) {
        std::multimap<Node, Node> map_out;
        for (auto arc : graphs[i].arcs) {
            map_out.insert(arc);
        }
        for (auto node : graphs[i].nodes) {
            if (map_out.count(node) == 0 && partition[node] == i) {
                stats[i]++;
            }
        }
    }
    return stats;
}

/**
 * @brief Get the number of constant nodes in each graph
 *
 * @param graphs
 * @param partition
 * @return std::vector<int>
 */
std::vector<int> number_constant_nodes(std::vector<Graph> graphs, std::map<int, int> partition) {
    std::vector<int> stats(graphs.size(), 0);
    for (int i = 0; i < graphs.size(); i++) {
        std::multimap<Node, Node> map_in;
        for (auto arc : graphs[i].arcs) {
            map_in.insert({arc.second, arc.first});
        }
        for (auto node : graphs[i].nodes) {
            if (map_in.count(node) == 0 && partition[node] == i) {
                stats[i]++;
            }
        }
    }
    return stats;
}

/**
 * @brief Get the number of active nodes in each graph
 *
 * @param graphs
 * @param partition
 * @return std::vector<int>
 */
std::vector<int> number_active_nodes(std::vector<Graph> graphs, std::map<int, int> partition) {
    std::vector<int> stats(graphs.size(), 0);
    for (int i = 0; i < graphs.size(); i++) {
        std::multimap<Node, Node> map_in;
        std::multimap<Node, Node> map_out;
        for (auto arc : graphs[i].arcs) {
            map_in.insert({arc.second, arc.first});
            map_out.insert(arc);
        }
        for (auto node : graphs[i].nodes) {
            if (map_in.count(node) != 0 && map_out.count(node) != 0 && partition[node] == i) {
                stats[i]++;
            }
        }
    }
    return stats;
}

/**
 * @brief Get the number of active nodes and the neightbor's one useful in each graph
 *
 * @param graphs
 * @param partition
 * @return std::vector<int>
 */
std::vector<int> number_active_and_distant_nodes(std::vector<Graph> graphs, std::map<int, int> partition) {
    std::vector<int> stats(graphs.size(), 0);
    for (int i = 0; i < graphs.size(); i++) {
        std::multimap<Node, Node> map_in;
        std::multimap<Node, Node> map_out;
        for (auto arc : graphs[i].arcs) {
            map_in.insert({arc.second, arc.first});
            map_out.insert(arc);
        }
        std::set<Node> distant_nodes;
        for (auto node : graphs[i].nodes) {
            if (map_in.count(node) != 0 && map_out.count(node) != 0 && partition[node] == i) {
                stats[i]++;
            }
            auto range = map_in.equal_range(node);
            // Count the distant nodes needed to compute PR
            for (auto it = range.first; it != range.second; ++it) {
                if (partition[it->second] != i) {
                    if (distant_nodes.count(it->second) == 0) {
                        distant_nodes.insert(it->second);
                        stats[i]++;
                    }
                }
            }
        }
    }
    return stats;
}

/**
 * @brief Output the stats in a readable format
 *
 * @param stats
 */
void show_stats(std::vector<std::vector<int>> stats) {
    std::cout << "[" << std::endl;
    for (auto sub : stats) {
        std::cout << "[ ";
        for (auto elt : sub) {
            std::cout << elt << ", ";
        }
        std::cout << "]," << std::endl;
    }
    std::cout << "]" << std::endl;
}

/**
 * @brief Display the stats of a graph
 *
 * @param filename
 */
void multi_dc_graph_stats(std::string filename) {
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC(filename, partition);
    std::vector<std::vector<int>> stats = get_edges_stats(graphs, partition, true);
    show_stats(stats);
    std::vector<double> dist = distribution_nodes(partition);
    std::cout << "Distribution of nodes in each partition" << std::endl;
    std::cout << dist << std::endl;
    std::vector<int> constant = number_constant_nodes(graphs, partition);
    std::cout << "Number of constant nodes in each graph" << std::endl;
    std::cout << constant << std::endl;
    std::vector<int> dangling = number_dangling_nodes(graphs, partition);
    std::cout << "Number of dangling nodes in each graph" << std::endl;
    std::cout << dangling << std::endl;
    std::vector<int> active = number_active_nodes(graphs, partition);
    std::cout << "Number of local active nodes in each graph" << std::endl;
    std::cout << active << std::endl;
    std::vector<int> active_distant = number_active_and_distant_nodes(graphs, partition);
    std::cout << "Number of active nodes and distant nodes needed in each graph (size adjacency matrix in POPPYm)" << std::endl;
    std::cout << active_distant << std::endl;
    std::cout << "Average number of inter DC communication" << std::endl;
    std::cout << average_number_inter_DC_communication(
                     graphs, partition, nextPowerOfTwo(*std::max_element(active.begin(), active.end())))
              << std::endl;
}
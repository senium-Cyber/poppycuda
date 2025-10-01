#include "include/utils_graph.cuh"

std::map<int, int> split_uniformaly_ascendant(Graph g, int number_partition) {
    std::map<int, int> partition;
    for (auto node : g.nodes) {
        partition[node] = node % number_partition;
    }
    return partition;
}

std::map<int, int> split_uniformaly_random(Graph g, int number_partition) {
    std::map<int, int> partition;
    for (auto node : g.nodes) {
        int rng = rand() % number_partition;
        partition[node] = rng;
    }
    return partition;
}

/**
 * @brief Creation of partition map.
 * Keys are the nodes and the values are the partition
 *
 * @param g
 * @param stats
 * @return std::map<int, int>
 */
std::map<int, int> split_following_stats_uniformaly_random(Graph g, std::vector<double> stats) {
    double threshold = 0.01;
    double acc = std::accumulate(stats.begin(), stats.end(), 0.0);
    // Checks if the sum of the stats is around 0.99 and 1.1
    assert(acc > 1 - threshold && acc < 1 + threshold);

    std::map<int, int> partition;

    std::vector<double> cumulative_stats;  // Cumulative bound of each iteration from 0 to 100
    cumulative_stats.push_back(stats[0] * 100);
    for (int i = 1; i < stats.size(); i++) {
        cumulative_stats.push_back(cumulative_stats[i - 1] + stats[i] * 100);
    }

    for (auto node : g.nodes) {
        // Generate number beteween 0 and 99
        int rng = rand() % 100;

        // Checks for the good partition
        for (int i = 0; i < cumulative_stats.size(); i++) {
            if (rng < cumulative_stats[i]) {
                partition[node] = i;
                break;
            }
        }
    }
    return partition;
}

/**
 * @brief Creation of partition map.
 * Keys are the nodes and the values are the partition
 *
 * @param g
 * @param stats
 * @return std::map<int, int>
 */
std::map<int, int> split_following_stats_ascendant(Graph g, std::vector<double> stats) {
    double threshold = 0.01;
    double acc = std::accumulate(stats.begin(), stats.end(), 0.0);
    // Checks if the sum of the stats is around 0.99 and 1.1
    assert(acc > 1 - threshold && acc < 1 + threshold);

    std::map<int, int> partition;

    std::vector<int> node_counts;
    int total_nodes = g.nodes.size();
    for (double stat : stats) {
        node_counts.push_back(static_cast<int>(std::round(stat * total_nodes)));
    }

    // Adjust total number of nodes of the rounding
    int sum_counts = std::accumulate(node_counts.begin(), node_counts.end(), 0);
    if (sum_counts != total_nodes) {
        node_counts.back() += total_nodes - sum_counts;
    }

    int current_partition = 0;
    int nodes_assigned = 0;
    for (auto node : g.nodes) {
        partition[node] = current_partition;
        nodes_assigned++;

        // Change partition when partition is full
        if (nodes_assigned >= node_counts[current_partition]) {
            current_partition++;
            nodes_assigned = 0;
        }
    }
    
    return partition;
}

std::vector<Graph> split_following_partition(Graph g, std::map<int, int> partition) {
    std::vector<Graph> graphs;
    // Get max element of partition values
    int max = std::max_element(
                  partition.begin(),
                  partition.end(),
                  [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })
                  .operator*()
                  .second;
    for (int i = 0; i < max + 1; i++) {
        graphs.emplace_back();
    }
    for (auto node : partition) {
        graphs[node.second].nodes.insert(node.first);
    }

    for (auto arc : g.arcs) {
        if (partition[arc.first] == partition[arc.second]) {
            graphs[partition[arc.first]].arcs.insert(arc);
        } else {
            graphs[partition[arc.first]].arcs.insert(arc);
            graphs[partition[arc.second]].arcs.insert(arc);
            graphs[partition[arc.first]].nodes.insert(arc.second);
            graphs[partition[arc.first]].nodes.insert(arc.first);
            graphs[partition[arc.second]].nodes.insert(arc.first);
            graphs[partition[arc.second]].nodes.insert(arc.second);
        }
    }
    return graphs;
}
#pragma once

#include "utils.cuh"


std::vector<std::vector<int>> get_edges_stats(Graph g, std::map<int, int> partition, bool directed);
std::vector<std::vector<int>> get_edges_stats(std::vector<Graph> graphs, std::map<int, int> partition, bool directed);
double average_number_inter_DC_communication(std::vector<Graph> graphs, std::map<Node, int> partition, int batchSize);
std::vector<int> number_active_and_distant_nodes(std::vector<Graph> graphs, std::map<int, int> partition);
std::vector<int> number_constant_nodes(std::vector<Graph> graphs, std::map<int, int> partition);
std::vector<int> number_dangling_nodes(std::vector<Graph> graphs, std::map<int, int> partition);
void show_stats(std::vector<std::vector<int>> stats);
void multi_dc_graph_stats(std::string filename);
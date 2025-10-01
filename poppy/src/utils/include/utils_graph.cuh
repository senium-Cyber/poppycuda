#pragma once

#include "utils.cuh"

std::map<int, int> split_uniformaly_ascendant(Graph g, int number_partition);

std::map<int, int> split_uniformaly_random(Graph g, int number_partition);

std::map<int, int> split_following_stats_uniformaly_random(Graph g, std::vector<double> stats);
std::map<int, int> split_following_stats_ascendant(Graph g, std::vector<double> stats);

std::vector<Graph> split_following_partition(Graph g, std::map<int, int> partition);
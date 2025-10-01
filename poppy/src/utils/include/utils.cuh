#pragma once

#include <dirent.h>
#include <openfhe.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "../../nauty/nausparse.h"

using std::string;
using namespace lbcrypto;

#define GRAPH_DIR "../Graphs/"
#define DOT_DIR "../dot/"

typedef int Node;
typedef int Orbit;

struct doublet {
    Node node;
    Ciphertext<DCRTPoly> cipher;
};

struct triplet {
    unsigned int indice;
    Orbit vertex;
    int deg_out;

    void print_triplet() { std::cout << indice << ' ' << vertex << ' ' << deg_out << ", "; }
};

struct quadruplet {
    int node_pr; // Node corresponding to a column
    unsigned int indice; // Indice in uj
    Orbit vertex; // vertex with an outgoing edge to node_pr
    int deg_out; // deg_out of vertex

    void print_quadruplet() { std::cout << node_pr << ' ' << indice << ' ' << vertex << ' ' << deg_out << ", "; }
};

struct Graph {
    std::set<int> nodes;
    std::set<std::pair<int, int>> arcs;
};

class DisjointSet {
   private:
    std::vector<int> parent;

   public:
    DisjointSet(int n) {
        parent.resize(n);
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);  // Path compression
        }
        return parent[x];
    }

    void unite(int x, int y) {
        int px = find(x);
        int py = find(y);

        if (px == py)
            return;

        // Toujours choisir le plus petit nœud comme parent
        if (px < py) {
            parent[py] = px;
        } else {
            parent[px] = py;
        }
    }
};

/* Return the next power of two of n. If n = 0, then return 0. */
uint32_t nextPowerOfTwo(uint32_t n);

/** Free momory allocated to a vector of ciphertexts.
*/
void delete_vector_of_ciphers(std::vector<Ciphertext<DCRTPoly>> vector_of_ciphers);

/** Free momory allocated to a vector of Plaintexts.
*/
void delete_vector_of_plains(std::vector<Plaintext> vector_of_plains);

/** Free momory allocated to map of vectors of Plaintexts.
*/
void delete_map_of_plains(std::map<int, std::vector<Plaintext>> map_of_plains);

void print_map_in_out(std::multimap<int, int> map, int n);

void print_map_u_j(std::multimap<int, triplet> map);

void print_map_reconstruction_u_j(std::vector<std::multimap<int, quadruplet>> map);

void print_PR(
    int n,
    std::multimap<int, triplet> map_u_j,
    std::multimap<int, triplet> map_deg0,
    std::vector<double> PR,
    std::vector<double> PR_deg0,
    int* orbits);

void print_cipher(Ciphertext<DCRTPoly> ct, CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, int s);

Graph read_graph(string name, string sep, bool directed);
string choose_graph();
bool is_number(const std::string& str);

std::vector<Graph> parse_graph_multi_DC(string filename, std::map<int /*node*/, int /*datacenter*/>& partition);

std::map<int, int>
compute_sub_partition(std::map<int, int> partition, std::set<std::pair<int, int>> arcs, int datacenter_id);

std::pair<std::multimap<int, int>, std::multimap<int, int>> process_maps(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    std::map<int, int>& orbits,
    int datacenter_id);

std::pair<std::multimap<int, int>, std::multimap<int, int>> sparse_nauty_routine(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    int& nb_local_nodes_after_pruning,
    std::map<int, int>& orbits,
    int datacenter_id,
    int& duration,
    std::set<int>& set_dangling_nodes
);

std::pair<std::multimap<int, int>, std::multimap<int, int>> equivalent_nodes_routine(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    int& nb_local_nodes_after_pruning,
    std::map<int, int>& orbits,
    int datacenter_id,
    int& duration,
    std::set<int>& set_dangling_nodes);

void export_to_dot(std::vector<Graph> graphs, std::string filename, std::map<int, int> partition);
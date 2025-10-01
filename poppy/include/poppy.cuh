#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include "../../nauty/nausparse.h"
#include "../src/utils/include/utils.cuh"
#include "cassert"
#include "openfhe.h"
using namespace lbcrypto;

#define naive_MULT_DEPTH(nb_iterations) (2 * nb_iterations)
#define POPPY_m_MULT_DEPTH(nb_iterations) (4 * nb_iterations)
#define POPPY_uj_MULT_DEPTH(nb_iterations) (4 * nb_iterations)
#define naive_LEVELS_CONSUMPTION 2
#define POPPY_m_LEVELS_CONSUMPTION 4
#define POPPY_uj_LEVELS_CONSUMPTION 4

typedef int Node;
typedef int Orbit;
struct NodeInCKKSVector {
    Ciphertext<DCRTPoly> vect;
    unsigned int index;
};

struct NodesInCKKSVector {
    Ciphertext<DCRTPoly> vect;
    std::set<std::pair<Node /* Node in vect*/, unsigned int /* indexe in vect */>> indexes;
};

template <typename PageRank>
class DataCenter {
   public:
    virtual ~DataCenter() = default;

    // Printers
    void print_data_center() {
        std::cout << "number_of_nodes : " << number_of_nodes << std::endl;
        std::cout << "map_in size : " << map_in.size() << std::endl;
        std::cout << "map_out size : " << map_out.size() << std::endl;
    }

    void print_map_in() {
        std::cout << "map_in : " << std::endl;
        for (auto it = map_in.begin(); it != map_in.end(); it++) {
            std::cout << it->first << " -> " << it->second << std::endl;
        }
    }
    void print_map_out() {
        std::cout << "map_out : " << std::endl;
        for (auto it = map_out.begin(); it != map_out.end(); it++) {
            std::cout << it->first << " -> " << it->second << std::endl;
        }
    }

    virtual void print_PRs() = 0;

    // Getters
    int get_number_of_nodes() const { return number_of_nodes; }
    std::multimap<int, int> get_map_in() const { return map_in; }
    std::multimap<int, int> get_map_out() const { return map_out; }
    auto get_cc() const { return cc; }
    auto get_keys() const { return keys; }
    auto get_public_key() const { return keys.publicKey; }
    int get_number_datacenters() { return nb_datacenter; };

    void set_number_datacenters(int datacenter_num) { nb_datacenter = datacenter_num; };

   protected:
    // Graph
    int number_of_nodes; // Total number of nodes
    std::multimap<int, int> map_in;
    std::multimap<int, int> map_out;
    int nb_datacenter = 0;

    // Crypto material
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;

    // PRs
    PageRank PRs;
};

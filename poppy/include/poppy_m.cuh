#pragma once

#include <condition_variable>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include "../../nauty/nausparse.h"
#include "../src/utils/include/utils_crypto.cuh"
#include "../src/utils/include/utils_poppy_m.cuh"
#include "poppy.cuh"

class POPPYm : public DataCenter<std::vector<Ciphertext<DCRTPoly>>> {
   public:
    void pre_processing(std::string csv_name, std::string separator, bool directed);
    void init(std::string csv_name, std::string separator = ",", bool directed = true);
    void run(int nb_iteration);

    // generation of the crypto material
    void gen_crypto_context(uint32_t multDepth);
    void gen_crypto_context_and_serialize(uint32_t multDepth, std::string folder);
    void read_crypto_context_and_keys(std::string folder);

    // getters
    Ciphertext<DCRTPoly> get_PRs();

    // printers

    void print_PRs() override;
    void print_data_center();

    std::map<Node, double> export_results();

   private:
    void initialize_PRs();
    void gen_rotation_keys();
    void compute_PR();
    void compute_PR_dangling_node(int nb_iteration);
    void nauty_routine(Graph nodes_and_arcs);
    void sparse_nauty_routine(Graph nodes_and_arcs);
    void sort_active_and_dangling_nodes();

    bool is_initialized = false;
    int s;  // size vectors, next power of 2 of the number of active nodes
    int m;  // size full matrix padded, least common multiple of s and n (number of active nodes)

    int number_active_nodes;
    int nb_dangling_nodes;
    int* orbits;
    std::vector<Node> active_nodes;
    std::vector<Node> active_dangling_nodes;
    std::vector<unsigned int> index_nodes_in_active_nodes;

    // std::vector<int> diags_null;
    // std::vector<int> diags_null_dangling_nodes;override
    // std::vector<Plaintext> diags;
    // std::vector<Plaintext> diags_dangling_nodes;

    std::vector<Diags> subDiags;
    std::vector<Diags> subDiags_dangling_nodes;

    std::vector<std::vector<double>> mask_matrix;
    std::vector<std::vector<double>> mask_matrix_dangling_nodes;

    std::vector<Ciphertext<DCRTPoly>> PRs;
    std::vector<Ciphertext<DCRTPoly>> PRs_deg0;
};

class POPPYmMultiDC : public DataCenter<std::vector<Ciphertext<DCRTPoly>>> {
   public:
    POPPYmMultiDC(int id) { this->datacenter_id = id; }

    void pre_processing();
    unsigned int
    init_local_graph(Graph& subgraph, int number_of_nodes, const std::map<Node, int>& partition, int& max_deg_in);
    void sort_active_and_dangling_nodes();
    void set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size);
    void set_max_dg_in(int max) { this->max_deg_in = max; }

    void compute_PR();
    void compute_PR_dangling_node();
    void reset_value_distant_nodes();

    void bootstrap(int multDepth) {
        for (auto& subPR : this->PRs) {
            if (multDepth - subPR->GetLevel() <= POPPY_m_LEVELS_CONSUMPTION) {
                subPR = cc->EvalBootstrap(subPR);
            }
        }
    }

    void compute_masked_distant_node(int iteration);
    void add_constant_to_distant_nodes(int iteration);
    void copy_previous_masked_distant_node() {
        for (auto pr : this->previous_masked_distant_nodes)
            for (auto vect : pr.second) {
                vect.reset();
            }
        this->previous_masked_distant_nodes.clear();
        this->previous_masked_distant_nodes = this->masked_distant_nodes;
    }

    void synchronize_distant_PRs(std::vector<NodesInCKKSVector> PRs_distant_nodes);

    // Getters
    std::vector<Ciphertext<DCRTPoly>> get_PRs() { return PRs; };
    std::vector<Ciphertext<DCRTPoly>> get_PRs_deg0() { return PRs_deg0; };
    std::pair<std::vector<Ciphertext<DCRTPoly>> /*cipher value*/, std::map<Node, unsigned int> /*indexes*/>
    get_PR_distant_nodes(int iteration, int distant_dc);

    int get_vector_size() { return s; };
    int get_matrix_size() { return m; };
    std::map<Node /*node*/, Orbit /*orbit*/> get_orbits() { return orbits; };
    std::vector<int> get_active_nodes() { return active_nodes; };
    std::vector<int> get_active_dangling_nodes() { return active_dangling_nodes; };
    int get_current_iteration() { return current_iteration; };
    int get_nb_dangling_nodes() { return nb_dangling_nodes; };
    int get_nb_cst_nodes() { return nb_cst_nodes; };
    int get_datacenter_id() { return datacenter_id; };

    // Printers
    void print_PRs() override;
    void printNodeCategories();

    std::map<Node, double> export_results();

    void delete_PRs() { delete_vector_of_ciphers(PRs); };
    void delete_PRs_deg0() { delete_vector_of_ciphers(PRs_deg0); };

    void creation_vect_PR_const();

   private:
    void initialize_PRs();

    int datacenter_id;
    int current_iteration = 0;
    int max_deg_in = 0;

    int nb_dangling_nodes;
    int nb_cst_nodes;
    int number_active_nodes;

    int s;  // size vectors
    int m;  // size full matrix padded

    std::map<Node /*node*/, Orbit /*orbit*/> orbits;
    std::map<Orbit /*orbit*/, Node /*repr*/> representative_orbits;
    std::map<Node /*node*/, int /*datacenter*/> partition;

    // Contains orbits of the nodes
    std::vector<Node> active_nodes;
    std::vector<Node> dangling_nodes;
    std::vector<Node> active_dangling_nodes;

    // Nodes connected to other datacenters and used in the computation by others, contains the real node id (not the
    // orbits)
    std::set<std::pair<Node, int /*data center*/>> nodes_connected_to_distant;
    std::set<std::pair<Node, int /*data center*/>> constant_nodes_connected_to_distant;

    // Nodes from other datacenters appearing in the subgraph
    std::vector<Node> distant_nodes;

    std::vector<std::vector<std::vector<double>>> mask_matrices;
    std::vector<std::vector<std::vector<double>>> mask_matrices_dangling_nodes;
    // Used to mask own PRs before adding the distant PRs
    std::vector<Plaintext> mask_distant_nodes;
    // Used to mask the current PRs before sending them to the other datacenters
    std::map<int /*data center*/, std::vector<Plaintext>> mask_effective_distant_nodes;

    std::map<int /*data center*/, std::vector<Plaintext>> border_init_constant_PRs;
    std::map<int /*data center*/, std::vector<Plaintext>> border_constant_PRs;

    std::map<int /*data center*/, std::map<Node, unsigned int>> index_distant_nodes;
    std::map<int /*data center*/, std::vector<Ciphertext<DCRTPoly>>> masked_distant_nodes;
    std::map<int /*data center*/, std::vector<Ciphertext<DCRTPoly>>> previous_masked_distant_nodes;

    std::vector<Plaintext> PRs_const_it0;
    std::vector<Plaintext> PRs_const;
    std::vector<Plaintext> PRs_dangling_const_it0;
    std::vector<Plaintext> PRs_dangling_const;

    // Vector is of size m / s
    std::vector<Ciphertext<DCRTPoly>> PRs;
    std::vector<Ciphertext<DCRTPoly>> PRs_deg0;
};
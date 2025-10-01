#pragma once

#include <condition_variable>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include "../src/utils/include/utils_crypto.cuh"
#include "../src/utils/include/utils_naive.cuh"
#include "poppy.cuh"

class Naive : public DataCenter<std::map<int, Ciphertext<DCRTPoly>>> {
   public:
    void pre_processing(std::string csv_name, std::string separator = " ", bool directed = true);
    void init(std::string csv_name, std::string separator = ",", bool directed = true);
    void run(int nb_iteration);

    // generation of the crypto material

    void gen_crypto_context(uint32_t multDepth);
    void gen_rotation_keys();

    // getters
    std::map<int, Ciphertext<DCRTPoly>> get_PRs();

    // printers
    void print_PRs() override;

   private:
    void initialize_PRs();
    void compute_PR();

    std::map<int, Ciphertext<DCRTPoly>> PRs;
};

class NaiveMultiDC : public DataCenter<std::map<int, Ciphertext<DCRTPoly>>> {
   public:
    NaiveMultiDC(int id) { this->datacenter_id = id; }  // Constructor
    virtual ~NaiveMultiDC() = default;                  // Destructor

    int datacenter_id;

    void pre_processing();
    void init_local_graph(Graph& subgraph, int number_of_nodes, const std::map<Node, int>& partition, int& max_deg_in);
    void set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size);
    void set_max_dg_in(int max) { this->max_deg_in = max; }
    void create_active_and_cst_nodes();

    void compute_PR();
    void compute_PR_dangling_node();

    void bootstrap(int multDepth) {
        for (auto& [node, subPR] : this->PRs) {
            if (multDepth - subPR->GetLevel() <= naive_LEVELS_CONSUMPTION)
                subPR = cc->EvalBootstrap(subPR);
        }
    }

    void synchronize_distant_PRs(std::vector<std::shared_ptr<NaiveMultiDC>> datacenters, int iteration);
    void compute_masked_distant_node();
    void copy_previous_masked_distant_node() {
        for (auto pr : this->previous_border_PRs)
            pr.second.second.reset();
        this->previous_border_PRs.clear();
        this->previous_border_PRs = border_PRs;
    }

    // Getters
    std::map<int, Ciphertext<DCRTPoly>> get_PRs() { return PRs; };
    std::map<int, Ciphertext<DCRTPoly>> get_PRs_deg0() { return PRs_deg0; };
    std::vector<std::pair<Node, Ciphertext<DCRTPoly>>> get_PR_distant_nodes(int datacenter_id, int iteration);
    int get_datacenter_id() { return datacenter_id; };
    int get_nb_dangling_nodes() { return nb_dangling_nodes; };

    // Printers
    void print_PRs() override;
    void printNodeCategories();

    std::map<Node, double> export_results();

    void delete_PRs() {
        for (int i = 0; i < PRs.size(); i++) {
            PRs[i].reset();  // Reset each ciphertext
        }
        PRs.clear();  // Clear the empty vector
    }
    void delete_PRs_deg0() {
        for (int i = 0; i < PRs_deg0.size(); i++) {
            PRs_deg0[i].reset();  // Reset each ciphertext
        }
        PRs_deg0.clear();  // Clear the empty vector
    }

   private:
    void initialize_PRs();

    int current_iteration = 0;

    int nb_dangling_nodes;
    int nb_cst_nodes;
    int max_deg_in;

    int s;  // size vectors

    std::map<Node /*node*/, Orbit /*orbit*/> orbits;
    std::map<Orbit /*orbit*/, Node /*repr*/> representative_orbits;
    std::map<Node /*node*/, int /*datacenter*/> partition;

    // Nodes from the datacenter
    std::vector<Node> active_nodes;
    std::vector<Node> cst_nodes;
    std::vector<Node> dangling_nodes;

    std::set<int> connected_datacenters;

    // Map for pagerank of active nodes
    std::map<int, Ciphertext<DCRTPoly>> PRs;
    // Map for pagerank of dangling nodes
    std::map<int, Ciphertext<DCRTPoly>> PRs_deg0;

    // Distant PRs we add in our computation at each step
    std::map<Node /* distant node */, Ciphertext<DCRTPoly /* corresponding ciphertext */>> distant_PRs;
    // Local PRs we send to distant data centers
    std::multimap<int /*datacenter*/, std::pair<Node, Ciphertext<DCRTPoly>>> border_PRs;
    std::multimap<int /*datacenter*/, std::pair<Node, Ciphertext<DCRTPoly>>> previous_border_PRs;
};
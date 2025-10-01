#ifndef POPPY_UJ_H
#define POPPY_UJ_H

#include "../../nauty/nausparse.h"
#include "../src/utils/include/utils_crypto.cuh"
#include "../src/utils/include/utils_poppy_uj.cuh"
#include "poppy.cuh"

class POPPYuj : public DataCenter<Ciphertext<DCRTPoly>> {
   public:
    void pre_processing(std::string csv_name, std::string separator = " ", bool directed = true);
    void init(std::string csv_name, std::string separator = ",", bool directed = true);
    void run(int nb_iteration);

    // generation of the crypto material
    void gen_crypto_context(uint32_t multDepth);
    void gen_rotation_keys();
    void gen_crypto_context_and_serialize(uint32_t multDepth, std::string folder);

    // serialization and deserialization of the crypto material
    void serialize_ciphers_PRs(std::string folder);
    void read_crypto_context_and_keys(std::string folder);

    // getters
    Ciphertext<DCRTPoly> get_PRs();
    std::map<int, int> get_orbits();
    // printers
    void print_data_center();
    void print_PRs() override;

    // debuggers
    void decrypt_and_print(std::string ct_name, Ciphertext<DCRTPoly> ct);

    std::map<Node, double> export_results();

   private:
    void creation_vect_PR_const();
    void reconstruction_u_j();
    void compute_PR_dangling_node(int nb_iteration);
    void compute_PR();
    void initialize_uj();

    // Graph
    int* orbits;
    std::multimap<Node, triplet> map_uj;
    std::multimap<Node, triplet> map_deg0;

    // PR poppy preprocess
    std::multimap<Node, doublet> PRs_cst;
    int nb_u_deg0;

    // PR poppy process
    int uj_size;
    int nb_of_uj;
    std::vector<Ciphertext<DCRTPoly>> u;
    Ciphertext<DCRTPoly> PRs;
    Ciphertext<DCRTPoly> PRs_deg0;
};

class POPPYujMultiDC : public DataCenter<Ciphertext<DCRTPoly>> {
   public:
    POPPYujMultiDC(int id) { this->datacenter_id = id; }

    void pre_processing();
    void set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size);
    int init_local_graph(Graph& subgraph, int number_of_nodes, const std::map<Node, int>& partition, int& max_deg_in);
    void set_max_iteration(int nb_iteration) { this->max_iteration = nb_iteration; }
    void set_max_dg_in(int max) { this->max_deg_in = max; }

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

    void compute_PR_dangling_node();
    void compute_PR();

    void bootstrap(int multDepth, bool last_it) {
        for (auto& subPR : this->PRs) {
            if (multDepth - subPR->GetLevel() <= POPPY_uj_LEVELS_CONSUMPTION || last_it)
                subPR = cc->EvalBootstrap(subPR);
        }
    }

    // Getters
    int get_vector_size() { return s; };
    std::vector<Ciphertext<DCRTPoly>> get_PRs() { return PRs; };
    std::vector<Ciphertext<DCRTPoly>> get_PRs_deg0() { return PRs_deg0; };
    int get_current_iteration() { return current_iteration; };
    int get_nb_dangling_nodes() { return nb_dangling_nodes; };
    int get_datacenter_id() { return datacenter_id; };

    // Printers
    void print_data_center() {};
    void print_PRs();
    void printNodeCategories();

    void creation_vect_PR_const();
    void reconstruction_u_j();

    struct DistantNodes {
        Ciphertext<DCRTPoly> masked_PRs;
        std::map<Node, unsigned int> index_nodes;
    };
    std::pair<std::vector<Ciphertext<DCRTPoly>> /*cipher value*/, std::map<Node, unsigned int> /*indexes*/>
    get_PR_distant_nodes(int iteration, int distant_dc);
    void synchronize_distant_PRs(
        std::vector<NodesInCKKSVector> PRs_distant_nodes,
        std::multimap<Node, triplet> map_uj,
        bool active_nodes);

    auto get_map_uj() { return map_uj; }
    auto get_map_deg0() { return map_deg0; }
    bool finished = false;

    // Completion step functions
    std::map<Node, double> export_results();

    void delete_PRs() { delete_vector_of_ciphers(PRs); };
    void delete_PRs_deg0() { delete_vector_of_ciphers(PRs_deg0); };

   private:
    void initialize_uj();
    int current_iteration = 0;
    int datacenter_id;
    int max_iteration = 0;

    // Graph
    int nb_dangling_nodes;

    int s;  // size vectors
    int m;  // size active_nodes padded

    std::map<Node /*node*/, Orbit /*orbit*/> orbits;
    std::map<Orbit /*orbit*/, Node /*repr*/> representative_orbits;
    std::map<Node /*node*/, int /*data center*/> partition;
    std::multimap<Node, triplet> map_uj;
    std::multimap<Node, triplet> map_deg0;
    std::vector<std::multimap<int, quadruplet>> map_reconstruction_uj;
    std::vector<std::multimap<int, quadruplet>> map_construction_u_deg0;

    std::set<std::pair<Node, int /*data center*/>> nodes_connected_to_distant;
    std::set<std::pair<Node, int /*data center*/>> constant_nodes_connected_to_distant;

    // PR poppy preprocess
    int nb_u_deg0;
    int max_deg_in;
    std::set<int> connected_datacenters;
    std::vector<std::pair<Node /*node*/, int /*data center*/>> nodes_needed;

    // PR poppy process
    int uj_size;
    int nb_of_uj;
    // vectors of size m / s
    std::vector<Ciphertext<DCRTPoly>> u;
    std::vector<Ciphertext<DCRTPoly>> u_dangling;
    std::vector<Plaintext> u_const;
    std::vector<Plaintext> u_const_init;
    std::vector<Plaintext> u_dangling_const;
    std::vector<Ciphertext<DCRTPoly>> PRs;
    std::vector<Ciphertext<DCRTPoly>> PRs_deg0;
    // Used to mask the current PRs before sending them to the distant nodes
    std::map<int /*data center*/, std::vector<Plaintext>> mask_effective_distant_nodes;
    std::map<int /*data center*/, std::vector<Ciphertext<DCRTPoly>>> masked_distant_nodes;
    std::map<int /*data center*/, std::vector<Ciphertext<DCRTPoly>>> previous_masked_distant_nodes;

    std::map<int /*data center*/, std::vector<Plaintext>> init_constant_PRs;
    std::map<int /*data center*/, std::vector<Plaintext>> border_constant_PRs;

    std::map<int /*data center*/, std::map<Node, unsigned int>> index_distant_nodes;
};

#endif  // POPPY_UJ_H

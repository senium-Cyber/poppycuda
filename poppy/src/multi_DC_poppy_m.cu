#include <thread>
#include "openfhe.h"
#include "poppy_m.cuh"
#include "/home/comp/csxtchen/cbct/poppy/poppy-icde/openfhe/openfhe-development/src/pke/include/openfhe.h"

#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Ciphertext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Context.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/KeySwitchingKey.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Limb.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Plaintext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/openfhe-interface/RawCiphertext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib//include/ConstantsGPU.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/Math.cuh"
// #include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/ParametrizedTest.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/cpuNTT.hpp"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/cpuNTT_nega.hpp"

//#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/hook.h"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/ApproxModEval.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Bootstrap.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/BootstrapPrecomputation.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/CoeffsToSlots.cuh"
#ifdef EXPE
#include "utils/include/utils_export.cuh"
Experiment expe_poppy_m;
#endif

size_t getMemoryUsage(const std::vector<std::vector<std::vector<double>>>& m) {
    size_t total_size = sizeof(m);  // Size of the first vector

    for (const auto& mat : m) {
        total_size += sizeof(mat);  // Size of each std::vector<std::vector<double>>
        for (const auto& row : mat) {
            total_size +=
                sizeof(row) + (row.capacity() * sizeof(double));  // Size of std::vector<double> and each double
        }
    }

    return total_size;
}

/**
 * @brief Initialize the local graph of the datacenter
 * In this step we compute the map containing for each nodes the in and out neighbors
 * We may apply the nauty algorithm to compute the orbits of the graph (prunning rule)
 * We initialize the data structure of active and dangling nodes
 * @param subgraph The local graph of the datacenter
 * @param number_of_nodes The number of nodes in the global graph
 * @param partition The local partition of the datacenter
 * @return unsigned int The number of active nodes
 */
unsigned int POPPYmMultiDC::init_local_graph(
    Graph& subgraph,
    int number_of_nodes,
    const std::map<int, int>& partition,
    int& max_deg_in) {
    this->number_of_nodes = number_of_nodes;
    this->partition = partition;
    int nb_local_nodes = 0;
    int nb_local_nodes_after_pruning = 0;
    int duration = 0;
    std::set<int> set_dangling_nodes;

#ifdef COMPUTE_ORBITS
    std::tie(map_in, map_out) = sparse_nauty_routine(
        &subgraph,
        &partition,
        nb_local_nodes,
        nb_local_nodes_after_pruning,
        this->orbits,
        this->datacenter_id,
        duration,
        set_dangling_nodes);

    this->nb_dangling_nodes = set_dangling_nodes.size();
    this->dangling_nodes = std::vector<int>(set_dangling_nodes.begin(), set_dangling_nodes.end());
#elif COMPUTE_EQUIVALENT_NODES
    std::tie(map_in, map_out) = equivalent_nodes_routine(
        &subgraph,
        &partition,
        nb_local_nodes,
        nb_local_nodes_after_pruning,
        this->orbits,
        this->datacenter_id,
        duration,
        set_dangling_nodes);

    this->nb_dangling_nodes = set_dangling_nodes.size();
    this->dangling_nodes = std::vector<int>(set_dangling_nodes.begin(), set_dangling_nodes.end());
#else
    std::tie(map_in, map_out) = process_maps(&subgraph, &partition, nb_local_nodes, this->orbits, this->datacenter_id);
    nb_local_nodes_after_pruning = nb_local_nodes;
#endif

#ifdef EXPE
    expe_poppy_m.addValue(
        duration,
        "Time to compute equivalent nodes on DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::MICROSECONDS);
    expe_poppy_m.addValue(
        nb_local_nodes,
        "Number of nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_poppy_m.addValue(
        nb_local_nodes_after_pruning,
        "Number of nodes after equivalent nodes routine in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
#endif

    // Compue the deg_in max in the sub graph
    for (auto it = map_in.begin(); it != map_in.end();) {
        auto range = map_in.equal_range(it->first);
        int count = std::distance(range.first, range.second);

        if (count > max_deg_in)
            max_deg_in = count;

        it = range.second;  // move to the next key
    }

    // We compute the actives and dangling nodes
    sort_active_and_dangling_nodes();

#ifdef EXPE
    expe_poppy_m.addValue(
        this->nb_cst_nodes,
        "Number of constant nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_poppy_m.addValue(
        this->nb_dangling_nodes,
        "Number of dangling nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
#endif

    // Count the number of nodes connected to distant datacenters
    std::unordered_map<int, std::set<Node>> all_nodes_connected_to_distant;
    for (auto arc : map_out) {
        if (partition.at(arc.second) != datacenter_id) {
            // The number of nodes connected to a remote node doesn't depend on the pruning rules.
            all_nodes_connected_to_distant[partition.at(arc.second)].insert(
                arc.first);  // We regroupe by distant DC all node connected to a distant node
            if (map_in.count(arc.first) == 0) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                constant_nodes_connected_to_distant.insert({arc.first, partition.at(arc.second)});
#else
                nodes_connected_to_distant.insert({arc.first, partition.at(arc.second)});
#endif
            } else {
                nodes_connected_to_distant.insert({arc.first, partition.at(arc.second)});
            }
        }
    }
    // We compute the max number of nodes connected to a distant DC
    std::vector<int> nb_nodes_connected_to_distant;
    for (const auto& [datacenter, nodes] : all_nodes_connected_to_distant) {
        nb_nodes_connected_to_distant.push_back(nodes.size());
    }
    int max_nodes_connected_to_distant;
    // If there is no connexion between DCs
    if (nb_nodes_connected_to_distant.size() == 0)
        max_nodes_connected_to_distant = 0;
    else
        max_nodes_connected_to_distant =
            *max_element(nb_nodes_connected_to_distant.begin(), nb_nodes_connected_to_distant.end());

#ifdef DEBUG
    printNodeCategories();
#endif
    return std::max(number_active_nodes, max_nodes_connected_to_distant);
}

/**
 * @brief Preprocessing of the datacenter
 * Create the adjacency matrix of the datacenter
 * Compute the masks
 * Initialize the PRs
 *
 */
void POPPYmMultiDC::pre_processing() {
    auto mask_matrix =
        create_adjacency_matrix(active_nodes, orbits, map_in, map_out, representative_orbits, partition, datacenter_id);
#ifdef DEBUG
    if (mask_matrix.size() <= 16) {
        std::cout << "mask matrix datacenter " << this->datacenter_id << std::endl;
        print_matrix(mask_matrix, mask_matrix.size());
    }
#endif

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::vector<std::vector<double>> mask_matrix_dangling_nodes;
    if (nb_dangling_nodes > 0) {
        mask_matrix_dangling_nodes = create_adjacency_matrix(
            this->active_dangling_nodes, orbits, map_in, map_out, representative_orbits, partition, datacenter_id);
#ifdef DEBUG
        if (mask_matrix_dangling_nodes.size() <= 16) {
            std::cout << "mask matrix dangling nodes datacenter " << this->datacenter_id << std::endl;
            print_matrix(mask_matrix_dangling_nodes, mask_matrix_dangling_nodes.size());
        }
#endif
    }
#endif

    int n = mask_matrix.size();
    // find the least common multiple of n (size of matrix) and s (size of
    // vectors)
    this->m = this->s;
    while (m < n) {
        this->m += this->s;
    }

    // Split the full matrix in sub matrices of size s
    this->mask_matrices = prepareMatrix(mask_matrix, this->s, this->m);
#ifdef DEBUG
    std::cout << "taille de mask matrices : " << getMemoryUsage(mask_matrices) << " octets" << std::endl;
#endif
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    this->mask_matrices_dangling_nodes = prepareMatrix(mask_matrix_dangling_nodes, this->s, this->m);
#endif

    // Compute the mask distant nodes, used to set the value of distant nodes to 0
    std::vector<double> ones(this->s, 1);
    std::vector<std::vector<double>> mask(m / s, ones);
    for (int node : distant_nodes) {
        int index = std::find(active_nodes.begin(), active_nodes.end(), node) - active_nodes.begin();
        assert(index != active_nodes.size());
        mask[index / s][index % s] = 0;
    }
#ifdef DEBUG
    std::cout << "Mask distant nodes" << std::endl;
#endif
    for (auto sub_mask : mask) {
#ifdef DEBUG
        std::cout << sub_mask << std::endl;
#endif
        auto pt = this->cc->MakeCKKSPackedPlaintext(sub_mask);
        this->mask_distant_nodes.push_back(pt);
#ifdef EXPE
        number_encoding++;
#endif
    }
    mask.clear();

    // Compute the mask actives_distant_nodes
    std::vector<double> zeros(this->s, 0.0);
    std::vector<std::vector<double>> masks(m / s, zeros);
    std::map<int /*data center*/, std::vector<std::vector<double>>> border_masks;
    for (auto node : nodes_connected_to_distant) {
        // if the mask for a distant data center is not already creat
        if (border_masks.count(node.second) == 0) {
            border_masks.insert({node.second, masks});
        }
        // if the key of a distant dc is not already in the set of indexes
        if (index_distant_nodes.count(node.second) == 0) {
            index_distant_nodes[node.second] = {};  // creation of an empty map
        }

        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbits[node.first]);
        assert(it != active_nodes.end());
        int index = it - active_nodes.begin();

        std::vector<std::vector<double>> current_mask = border_masks.at(node.second);
        current_mask[index / s][index % s] = 1.0 / map_out.count(node.first);
        index_distant_nodes[node.second][node.first] =  index;  // we add the value node.first -> index in the map in node.second

        border_masks.at(node.second) = current_mask;
    }

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    // Fill empty slots with the constant nodes
    std::map<int /*data center*/, std::vector<std::vector<double>>> vect_border_const_PRs;
    std::map<int /*data center*/, std::vector<std::vector<double>>> vect_border_init_const_PRs;
    std::vector<int> index_constant(this->nb_datacenter, -1);
    for (auto node : constant_nodes_connected_to_distant) {
        // if the mask for a distant data center is not already creat
        int dc_id = node.second;
        if (vect_border_const_PRs.count(dc_id) == 0) {
            vect_border_const_PRs.insert({dc_id, masks});
            vect_border_init_const_PRs.insert({dc_id, masks});
        }

        // if the key of a distant dc is not already in the set of indexes
        if (index_distant_nodes.count(node.second) == 0) {
            index_distant_nodes[node.second] = {};  // creation of an empty map
        }

        for (int i = index_constant[dc_id] + 1; i < m; i++) {
            // If the vector has not been created
            if (border_masks.count(dc_id) == 0) {
                index_constant[dc_id] = i;

                // for initialization
                std::vector<std::vector<double>> current_mask_init = vect_border_init_const_PRs.at(node.second);
                current_mask_init[index_constant[dc_id] / s][index_constant[dc_id] % s] =
                    1.0 / (number_of_nodes * map_out.count(node.first));

                index_distant_nodes[node.second][node.first] = index_constant[dc_id];
                border_masks.insert({node.second, masks});
                vect_border_init_const_PRs.at(node.second) = current_mask_init;

                // for other iterations
                std::vector<std::vector<double>> current_mask = vect_border_const_PRs.at(node.second);
                current_mask[index_constant[dc_id] / s][index_constant[dc_id] % s] =
                    0.15 / (double)(this->max_deg_in * map_out.count(node.first));

                vect_border_const_PRs.at(node.second) = current_mask;
                break;
                // Or if the slot is empty
            } else if (border_masks.at(node.second)[i / s][i % s] == 0) {
                index_constant[dc_id] = i;

                // for initialization
                std::vector<std::vector<double>> current_mask_init = vect_border_init_const_PRs.at(node.second);
                current_mask_init[index_constant[dc_id] / s][index_constant[dc_id] % s] =
                    1.0 / (number_of_nodes * map_out.count(node.first));

                index_distant_nodes[node.second][node.first] = index_constant[dc_id];
                border_masks.insert({node.second, masks});
                vect_border_init_const_PRs.at(node.second) = current_mask_init;

                // for other iterations
                std::vector<std::vector<double>> current_mask = vect_border_const_PRs.at(node.second);
                current_mask[index_constant[dc_id] / s][index_constant[dc_id] % s] =
                    0.15 / (double)(this->max_deg_in * map_out.count(node.first));

                vect_border_const_PRs.at(node.second) = current_mask;
                break;
            }
        }
    }
#endif

#ifdef DEBUG
    std::cout << "Masks PRs DC " << datacenter_id << std::endl;
    for (auto masks : border_masks) {
        std::cout << "connected to DC " << masks.first << std::endl;
        std::cout << masks.second << std::endl;
    }
#endif
    // Encode masks
    for (auto masks : border_masks) {
        std::vector<Plaintext> plains;
        for (auto sub_vect : masks.second) {
            plains.push_back(cc->MakeCKKSPackedPlaintext(sub_vect));
#ifdef EXPE
            number_encoding++;
#endif
        }
        this->mask_effective_distant_nodes.insert({masks.first, plains});
    }

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    // Encode border const vect
    for (auto masks : vect_border_const_PRs) {
        std::vector<Plaintext> plain;
        for (auto sub_vect : masks.second) {
            plain.push_back(cc->MakeCKKSPackedPlaintext(sub_vect));
#ifdef EXPE
            number_encoding++;
#endif
        }
        this->border_constant_PRs.insert({masks.first, plain});
    }

    // Encode init const vect
    for (auto masks : vect_border_init_const_PRs) {
        std::vector<Plaintext> plain;
        for (auto sub_vect : masks.second) {
            plain.push_back(cc->MakeCKKSPackedPlaintext(sub_vect));
#ifdef EXPE
            number_encoding++;
#endif
        }
        this->border_init_constant_PRs.insert({masks.first, plain});
    }

    // Creation vect PR const
    creation_vect_PR_const();
#endif

    initialize_PRs();
    compute_masked_distant_node(0);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    add_constant_to_distant_nodes(0);  // Add constant nodes at initialization (iteration 0)
#endif
    copy_previous_masked_distant_node();
}

/**
 * @brief Set the common crypto context and keys
 *
 * @param cc Common crypto context
 * @param keys
 * @param vector_size Common batch size
 */
void POPPYmMultiDC::set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size) {
    this->cc = cc;
    this->keys = keys;
    this->s = vector_size;
}

/**
 * @brief Create the active_nodes vector
 * Sort the node in active and dangling nodes list
 *
 */
void POPPYmMultiDC::sort_active_and_dangling_nodes() {
    std::set<int> active_nodes_set;  // nodes needed for pagerank computation

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::set<int> constant_nodes_set;  // constant nodes
    std::multimap<int /*dangling node*/, int /*active node*/>
        dangling_actives_nodes;  // active nodes connected to dangling nodes
    std::set<int> dangling_actives_nodes_set;

    // Creation active_nodes and actives_dangling_nodes
    for (auto node : this->orbits) {
        if ((map_out.count(node.first) > 0 && map_in.count(node.first) > 0) ||
            (partition[node.first] != datacenter_id &&
             map_out.count(node.first) > 0)) {  // Not a dangling node nor a constant
            active_nodes_set.insert(node.second);
            representative_orbits[node.second] = node.first;
            if (partition[node.first] != datacenter_id)
                distant_nodes.push_back(node.second);

        } else {
            auto it = find(this->dangling_nodes.begin(), this->dangling_nodes.end(), node.second);
            if (it != this->dangling_nodes.end()) {
                representative_orbits[node.second] = node.first;
                dangling_actives_nodes_set.insert(node.second);

                // Loop on all neighbors of the dangling node
                auto neighbors = map_in.equal_range(node.first);
                for (auto it = neighbors.first; it != neighbors.second; ++it) {
                    // We check if its neighbor is not a local constant node
                    if (map_in.count(it->second) > 0 || partition[it->second] != datacenter_id) {
                        dangling_actives_nodes.insert({node.second, orbits[it->second]});
                        dangling_actives_nodes_set.insert(orbits[it->second]);
                        representative_orbits[orbits[it->second]] = it->second;
                    }
                }
            }
        }
        // Constant node in the current data center
        if (map_in.count(node.first) == 0 && partition[node.first] == datacenter_id)
            constant_nodes_set.insert(node.second);
    }
    this->nb_cst_nodes = constant_nodes_set.size();

#else
    for (auto node : this->orbits) {
        if (partition[node.first] == datacenter_id ||
            (partition[node.first] != datacenter_id && map_out.count(node.first) > 0)) {
            active_nodes_set.insert(node.second);

            if (partition[node.first] != datacenter_id)
                distant_nodes.push_back(node.second);
        }
        representative_orbits[node.second] = node.first;
    }
    this->nb_dangling_nodes = 0;
    this->nb_cst_nodes = 0;
#endif

    // We transform sets into vectors
    this->active_nodes = std::vector<int>(active_nodes_set.begin(), active_nodes_set.end());

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)

    // Set nb_active_nodes equal to the max between active and dangling size
    this->number_active_nodes = std::max(this->active_nodes.size(), dangling_actives_nodes_set.size());
#else
    this->number_active_nodes = this->active_nodes.size();
#endif

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    this->active_dangling_nodes = std::vector<int>(this->number_active_nodes);
    std::vector<bool> completed_cells(this->number_active_nodes, false);

    // We order the active nodes (connected to dangling one) in a manner that the nodes take the same place as in
    // actives_nodes And we place between them the dangling nodes
    for (auto actives_dangling_node : dangling_actives_nodes) {
        auto it = std::find(active_nodes.begin(), active_nodes.end(), actives_dangling_node.second);
        assert(it != active_nodes.end());
        if (it != active_nodes.end()) {
            int index_in_actives = it - active_nodes.begin();
            active_dangling_nodes[index_in_actives] = active_nodes[index_in_actives];
            completed_cells[index_in_actives] = true;
        }
    }
    // We place the dangling nodes in the remaining cells
    int index = 0;
    for (auto dangling_node : dangling_nodes) {
        assert(index < this->number_active_nodes);
        for (int i = index; i < this->number_active_nodes; i++) {
            if (!completed_cells[i]) {
                active_dangling_nodes[i] = dangling_node;
                completed_cells[i] = true;
                index = i;
                break;
            }
        }
    }
#endif
}

void POPPYmMultiDC::printNodeCategories() {
    std::cout << "\n--- Datacenter " << this->datacenter_id << " Node Categories ------------------\n";

    // 1) Print active nodes
    std::cout << "Active nodes (" << active_nodes.size() << "): ";
    for (auto node : active_nodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // 2) Print dangling nodes
    std::cout << "Dangling nodes (" << dangling_nodes.size() << "): ";
    for (auto node : dangling_nodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // 3) Collect constant nodes:
    //    We'll iterate through all nodes owned by this DC (based on 'partition'),
    //    and pick those with in-degree 0.
    std::vector<int> constant_nodes;
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    for (auto& kv : partition) {
        int node = kv.first;
        int dc = kv.second;
        if (dc == this->datacenter_id) {
            // Belongs to this DC
            if (map_in.count(node) == 0) {
                constant_nodes.push_back(node);
            }
        }
    }
#endif
    // 4) Print constant nodes
    std::cout << "Constant nodes (" << constant_nodes.size() << "): ";
    for (auto node : constant_nodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "---------------------------------------------------\n\n";
}

/**
 * @brief Initialize the PageRank of the datacenter, all the nodes have the same PR equal to 1/n,
 * where n is the number of nodes in the global graph
 * At the end of the initialization, the PRs masked that will be sent to the distant datacenters are computed
 *
 */
void POPPYmMultiDC::initialize_PRs() {
    std::vector<double> x(this->s, 1.0 / this->number_of_nodes);

    Plaintext iv = this->cc->MakeCKKSPackedPlaintext(x);
    Ciphertext<DCRTPoly> ct_iv = this->cc->Encrypt(keys.publicKey, iv);
#ifdef EXPE
    number_encoding++;
    number_encryption++;
#endif

    this->PRs = std::vector<Ciphertext<DCRTPoly>>(m / s, ct_iv);

    ct_iv.reset();  // Clear ciphertext
}

/**
 * @brief Set to 0 the PRs of the distant nodes, used at each synchronization
 *
 */
void POPPYmMultiDC::reset_value_distant_nodes() {
    for (int i = 0; i < m / s; i++) {
        
        this->PRs[i] = this->cc->EvalMult(this->PRs[i], this->mask_distant_nodes[i]);
#ifdef EXPE
        number_mult_plain++;
#endif
    }
}

/**
 * @brief We compute the ciphertexts that will be sent to the distant datacenters
 *
 */
void POPPYmMultiDC::compute_masked_distant_node(int iteration) {
    for (auto distant_dc : mask_effective_distant_nodes) {
        if (iteration == 0) {  // We need to creat an empty vector at first iteration
            std::vector<Ciphertext<DCRTPoly>> masked_PRs(m / s);
            masked_distant_nodes.insert({distant_dc.first, masked_PRs});
        }
        std::vector<Ciphertext<DCRTPoly>> current_mask = masked_distant_nodes.at(distant_dc.first);
        for (int i = 0; i < m / s; i++) {
            current_mask[i] = this->cc->EvalMult(this->PRs[i], distant_dc.second[i]);
#ifdef EXPE
            number_mult_plain++;
#endif
        }
        masked_distant_nodes.at(distant_dc.first) = current_mask;
    }
}

void POPPYmMultiDC::add_constant_to_distant_nodes(int iteration) {
    if (iteration == 0) {
        for (auto distant_dc : masked_distant_nodes) {
            std::vector<Ciphertext<DCRTPoly>> current_mask = masked_distant_nodes.at(distant_dc.first);
            for (int i = 0; i < m / s; i++) {
                if (border_init_constant_PRs.count(distant_dc.first) != 0) {
                    current_mask[i] =
                        cc->EvalAdd(border_init_constant_PRs.at(distant_dc.first)[i], distant_dc.second[i]);
#ifdef EXPE
                    number_add++;
#endif
                }
            }
            masked_distant_nodes.at(distant_dc.first) = current_mask;
        }
    } else {
        for (auto distant_dc : masked_distant_nodes) {
            std::vector<Ciphertext<DCRTPoly>> current_mask = masked_distant_nodes.at(distant_dc.first);
            for (int i = 0; i < m / s; i++) {
                if (border_constant_PRs.count(distant_dc.first) != 0) {
                    current_mask[i] = cc->EvalAdd(border_constant_PRs.at(distant_dc.first)[i], distant_dc.second[i]);
#ifdef EXPE
                    number_add++;
#endif
                }
            }
            masked_distant_nodes.at(distant_dc.first) = current_mask;
        }
    }
}

/**
 * @brief Synchronize the PRs of the distant nodes with the current PRs
 * To optimize the composed rotations, we sort the indexes to rotate by the rotation to apply
 * @param PRs_distant_nodes
 */
void POPPYmMultiDC::synchronize_distant_PRs(std::vector<NodesInCKKSVector> PRs_distant_nodes) {
    this->reset_value_distant_nodes();
    // Add distant PRs to the current PRs
    for (auto PR_distant_node : PRs_distant_nodes) {
        auto precomp = cc->EvalFastRotationPrecompute(PR_distant_node.vect);
        auto sorted_indexes =
            std::vector<std::pair<Node, unsigned int>>(PR_distant_node.indexes.begin(), PR_distant_node.indexes.end());
        struct {
            std::vector<Node> active_nodes;
            std::map<Node, int> orbits;
            int nb_slots;
            // Custom comparator to sort the indexes, the less between the both is wich has the smallest
            bool operator()(std::pair<Node, unsigned int> a, std::pair<Node, unsigned int> b) const {
                auto it1 = std::find(active_nodes.begin(), active_nodes.end(), orbits.at(a.first));

                int index1 = it1 - active_nodes.begin();
                auto it2 = std::find(active_nodes.begin(), active_nodes.end(), orbits.at(b.first));

                int index2 = it2 - active_nodes.begin();
                return (a.second + nb_slots - index1) % nb_slots < (b.second + nb_slots - index2) % nb_slots;
            }
        } customLess{this->active_nodes, this->orbits, this->s};

        std::sort(sorted_indexes.begin(), sorted_indexes.end(), customLess);
        std::map<int, Ciphertext<DCRTPoly>> rotated_vects;
        for (auto node : sorted_indexes) {
            auto it = std::find(this->active_nodes.begin(), this->active_nodes.end(), this->orbits[node.first]);
            assert(it != this->active_nodes.end());
            int index = it - this->active_nodes.begin();
            // Rotate the PR of the distant node
            int rot = (node.second + this->s - index) % this->s;
#ifdef DEBUG
            std::cout << "Rot value: " << rot << std::endl;
#endif
            Ciphertext<DCRTPoly> rotated_PR = PR_distant_node.vect;
            if (rot != 0) {
                auto it_rot = rotated_vects.find(rot);
                if (it_rot == rotated_vects.end()) {
#ifdef COEUS
                    // Check if rot is a power of 2
                    if ((rot & (rot - 1)) == 0) {
                        rotated_PR =
                            this->cc->EvalFastRotation(PR_distant_node.vect, rot, 2 * cc->GetRingDimension(), precomp);
                        rotated_vects.insert({rot, rotated_PR});
#ifdef EXPE
                        number_rot++;
#endif
                    } else {
                        int rotated = 0;
                        // Find the largest power of 2 that is smaller than rot
                        int k = 1;
                        while (k <= rot) {
                            k = k << 1;
                        }
                        k = k >> 1;
                        rotated = k;
                        // Check if it has already been computed
                        auto it_rotated = rotated_vects.find(k);
                        if (it_rotated == rotated_vects.end()) {
                            rotated_PR = this->cc->EvalFastRotation(
                                PR_distant_node.vect, k, 2 * cc->GetRingDimension(), precomp);
                            rotated_vects.insert({k, rotated_PR});
#ifdef EXPE
                            number_rot++;
#endif
                        } else {
                            rotated_PR = it_rotated->second;
                        }
                        // Compute the remaining rotation
                        while (rotated != rot) {
                            // Find the largest power of 2 that is smaller than remaining rot
                            int k = 1;
                            while (k <= (rot - rotated)) {
                                k = k << 1;
                            }
                            k = k >> 1;
                            rotated += k;
                            // Rotate the PR
                            rotated_PR = this->cc->EvalRotate(rotated_PR, k);
                            rotated_vects.insert({rotated, rotated_PR});
#ifdef EXPE
                            number_rot++;
#endif
                        }
                    }
#else
                    rotated_PR =
                        this->cc->EvalFastRotation(PR_distant_node.vect, rot, 2 * cc->GetRingDimension(), precomp);
                    rotated_vects.insert({rot, rotated_PR});
#ifdef EXPE
                    number_rot++;
#endif
#endif
                } else {
                    rotated_PR = it_rot->second;
                }
            }
            // Re mask to obtain only one PRs in the vector
            std::vector<double> vect(this->s, 0.0);
            vect[index % s] = 1;
            Plaintext mask = this->cc->MakeCKKSPackedPlaintext(vect);
            auto mult = cc->EvalMult(rotated_PR, mask);
            // Add rotated masked PRs to the current PRs
            cc->EvalAddInPlace(this->PRs[index / s], mult);
#ifdef EXPE
            number_encoding++;
            number_mult_plain++;
            number_add++;
#endif
        }
    }
}

/**
 * @brief Compute the PageRank of the datacenter for the current iteration
 * The computation is done by a matrix multiplication
 * Following by the damping factor
 *
 * @param PRs_distant_nodes
 */
void POPPYmMultiDC::compute_PR() {
#ifdef COEUS
    this->PRs = matmul_coeus_splited(cc, keys, this->mask_matrices, this->PRs);
#else
    this->PRs = matmul_splited(cc, keys, this->mask_matrices, this->PRs);
#endif

    for (int i = 0; i < m / s; i++) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        if (this->current_iteration > 0)
            cc->EvalAddInPlace(PRs[i], PRs_const[i]);
        else
            cc->EvalAddInPlace(PRs[i], PRs_const_it0[i]);
#endif
        cc->EvalMultInPlace(PRs[i], 0.85);
        cc->EvalAddInPlace(PRs[i], 0.15 / this->max_deg_in);
#ifdef EXPE
        number_add_scalar++;
        number_add++;
        number_mult_scalar++;
#endif
    }

    this->current_iteration++;
}

/**
 * @brief Compute PRs for the dangling nodes, called once at the end of the PR algorithm
 *
 * @param PRs_distant_nodes
 */
void POPPYmMultiDC::compute_PR_dangling_node() {
#ifdef COEUS
    PRs_deg0 = matmul_coeus_splited(cc, keys, this->mask_matrices_dangling_nodes, PRs);
#else
    PRs_deg0 = matmul_splited(cc, keys, this->mask_matrices_dangling_nodes, PRs);
#endif

    for (int i = 0; i < m / s; i++) {
        if (this->current_iteration > 0)
            cc->EvalAddInPlace(PRs_deg0[i], PRs_dangling_const[i]);
        else
            cc->EvalAddInPlace(PRs_deg0[i], PRs_dangling_const_it0[i]);

        cc->EvalMultInPlace(PRs_deg0[i], 0.85);
        cc->EvalAddInPlace(PRs_deg0[i], 0.15 / this->max_deg_in);
#ifdef EXPE
        number_add_scalar++;
        number_add++;
        number_mult_scalar++;
#endif
    }
}

/**
 * @brief Create the u vector of constant PRs
 *
 */
void POPPYmMultiDC::creation_vect_PR_const() {
    std::map<Node, double> cst_nodes_it0;
    std::map<Node, double> cst_nodes;
    for (auto orbit : orbits) {
        if (map_in.count(orbit.first) == 0 && partition[orbit.first] == datacenter_id) {
            cst_nodes_it0.insert({orbit.first, 1.0 / (this->number_of_nodes * map_out.count(orbit.first))});
            cst_nodes.insert({orbit.first, 0.15 / (double)(this->max_deg_in * map_out.count(orbit.first))});
        }
    }

    std::vector<std::vector<double>> vect_const_it0(m / s, std::vector<double>(number_active_nodes));
    std::vector<std::vector<double>> vect_const_dangling_it0(m / s, std::vector<double>(number_active_nodes));
    std::vector<std::vector<double>> vect_const(m / s, std::vector<double>(number_active_nodes));
    std::vector<std::vector<double>> vect_const_dangling(m / s, std::vector<double>(number_active_nodes));
    for (auto orbit : orbits) {
        if (partition[orbit.first] == datacenter_id) {
            double sum_cst_value_it0 = 0.0;
            double sum_cst_value = 0.0;
            auto it = std::find(active_nodes.begin(), active_nodes.end(), orbit.second);
            if (it != active_nodes.end()) {
                int indice = it - active_nodes.begin();

                auto neighbors = map_in.equal_range(orbit.first);
                for (auto it = neighbors.first; it != neighbors.second; ++it) {
                    if (cst_nodes.count(it->second) != 0) {
                        sum_cst_value_it0 += cst_nodes_it0.at(it->second);
                        sum_cst_value += cst_nodes.at(it->second);
                    }
                }
                vect_const_it0[indice / s][indice % s] = sum_cst_value_it0;
                vect_const[indice / s][indice % s] = sum_cst_value;

            } else if (nb_dangling_nodes > 0) {
                auto it_dangling_node =
                    std::find(active_dangling_nodes.begin(), active_dangling_nodes.end(), orbit.second);
                if (it_dangling_node != active_dangling_nodes.end()) {
                    int indice = it_dangling_node - active_dangling_nodes.begin();

                    auto neighbors = map_in.equal_range(orbit.first);
                    for (auto it = neighbors.first; it != neighbors.second; ++it) {
                        if (cst_nodes.count(it->second) != 0) {
                            sum_cst_value_it0 += cst_nodes_it0.at(it->second);
                            sum_cst_value += cst_nodes.at(it->second);
                        }
                    }

                    assert(number_active_nodes > indice);
                    vect_const_dangling_it0[indice / s][indice % s] = sum_cst_value_it0;
                    vect_const_dangling[indice / s][indice % s] = sum_cst_value;
                }
            }
        }
    }

    for (int i = 0; i < vect_const_it0.size(); i++) {
        // Encode the constant vectors
        PRs_const.push_back(cc->MakeCKKSPackedPlaintext(vect_const[i]));
        PRs_const_it0.push_back(cc->MakeCKKSPackedPlaintext(vect_const_it0[i]));

        // Encode the constant vectors for dangling nodes
        this->PRs_dangling_const.push_back(cc->MakeCKKSPackedPlaintext(vect_const_dangling[i]));
        this->PRs_dangling_const_it0.push_back(cc->MakeCKKSPackedPlaintext(vect_const_dangling_it0[i]));
#ifdef EXPE
        number_encoding += 4;
#endif
    }
}

std::pair<std::vector<Ciphertext<DCRTPoly>> /*cipher value*/, std::map<Node, unsigned int> /*indexes*/>
POPPYmMultiDC::get_PR_distant_nodes(int iteration, int distant_dc) {
    if (current_iteration > iteration)
        return {previous_masked_distant_nodes.at(distant_dc), index_distant_nodes.at(distant_dc)};
    else
        return {masked_distant_nodes.at(distant_dc), index_distant_nodes.at(distant_dc)};
}

/**
 * @brief Decrypt and print the local PageRank of the datacenter
 *
 */
void POPPYmMultiDC::print_PRs() {
    // Decrypt pagerank of active nodes
    std::vector<double> vect_PRs;
    for (int i = 0; i < this->PRs.size(); i++) {
        Plaintext res;
        cc->Decrypt(keys.secretKey, PRs[i], &res);
        res->SetLength(s);
        auto subVect = res->GetRealPackedValue();
        vect_PRs.insert(vect_PRs.end(), subVect.begin(), subVect.end());
    }

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    // Decrypt pagerank of dangling nodes
    std::vector<double> vect_PRs_deg0;
    if (this->nb_dangling_nodes > 0) {
        for (int i = 0; i < this->PRs_deg0.size(); i++) {
            Plaintext res;
            cc->Decrypt(keys.secretKey, PRs_deg0[i], &res);
            res->SetLength(s);
            auto subVect = res->GetRealPackedValue();
            vect_PRs_deg0.insert(vect_PRs_deg0.end(), subVect.begin(), subVect.end());
        }
    }
#endif

#ifdef DEBUG
    std::cout << "vect PRs: " << vect_PRs << std::endl;
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::cout << "vect PRs deg0: " << vect_PRs_deg0 << std::endl;
#endif
#endif

    // Print the PRs in the order of the nodes
    for (auto orbit : orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;

        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbit.second);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        auto it_dangling_node = std::find(dangling_nodes.begin(), dangling_nodes.end(), orbit.second);
#endif
        if (it != active_nodes.end()) {  // Check if the node is an active node node
            int indice = it - active_nodes.begin();
            assert(indice < vect_PRs.size());
            std::cout << orbit.first << " is " << vect_PRs[indice] << " at " << indice << std::endl;
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (it_dangling_node != dangling_nodes.end()) {  // Check if the node is a dangling node
            auto it_indice_dangling_node =
                std::find(active_dangling_nodes.begin(), active_dangling_nodes.end(), orbit.second);
            int indice = it_indice_dangling_node - active_dangling_nodes.begin();
            assert(indice < vect_PRs_deg0.size());
            std::cout << orbit.first << " is " << vect_PRs_deg0[indice] << " (dangling) at " << indice << std::endl;
        } else  // Otherwise it is a constant node
            std::cout << orbit.first << " is " << 0.15 / this->max_deg_in << " (constant)" << std::endl;
#endif
    }

    // Clear vectors of cleartexts
    vect_PRs.clear();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    vect_PRs_deg0.clear();
#endif
}

std::map<int, double> POPPYmMultiDC::export_results() {
    // Decrypt pagerank of standard nodes
    std::vector<double> vect_PRs;
    for (int i = 0; i < this->PRs.size(); i++) {
        Plaintext res;
        cc->Decrypt(keys.secretKey, PRs[i], &res);
        res->SetLength(s);
        auto subVect = res->GetRealPackedValue();
        vect_PRs.insert(vect_PRs.end(), subVect.begin(), subVect.end());
    }

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    // Decrypt pagerank of dangling nodes
    std::vector<double> vect_PRs_deg0;
    if (this->nb_dangling_nodes > 0) {
        for (int i = 0; i < this->PRs_deg0.size(); i++) {
            Plaintext res;
            cc->Decrypt(keys.secretKey, PRs_deg0[i], &res);
            res->SetLength(s);
            auto subVect = res->GetRealPackedValue();
            vect_PRs_deg0.insert(vect_PRs_deg0.end(), subVect.begin(), subVect.end());
        }
    }
#endif

    std::map<int, double> results;
    for (auto orbit : this->orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;
        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbit.second);
        if (it != active_nodes.end()) {
            int indice = it - active_nodes.begin();
            assert(indice < vect_PRs.size());
            results.insert({orbit.first, vect_PRs[indice]});
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (this->nb_dangling_nodes > 0) {  // Check if there is dangling nodes
            auto it_dangling_node = std::find(active_dangling_nodes.begin(), active_dangling_nodes.end(), orbit.second);
            int indice = it_dangling_node - active_dangling_nodes.begin();
            assert(indice < vect_PRs_deg0.size());
            results.insert({orbit.first, vect_PRs_deg0[indice]});
        } else
            results.insert({orbit.first, 0.15 / this->max_deg_in});
#endif
    }

    // Clear vectors of cleartexts
    vect_PRs.clear();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    vect_PRs_deg0.clear();
#endif

    return results;
}

/**
 * @brief Loop over the distant nodes comming in and get their PRs to add them to the PRs of the current datacenter
 *  This phase is done at the beginning of each iteration to receive the PRs of the distant nodes
 *
 * @param datacenter Datacenter to synchronize
 * @param datacenters Other datacenters
 * @param partition
 * @param current_iteration
 * @return std::vector<NodesInCKKSVector>
 */
std::vector<NodesInCKKSVector> synch_data(
    std::shared_ptr<POPPYmMultiDC> datacenter,
    std::vector<std::shared_ptr<POPPYmMultiDC>> datacenters,
    std::map<Node, int> partition,
    int current_iteration) {
    // Loop over the distant nodes comming in and get their PRs to add them to the PRs of the current datacenter
    // First compute the list of the distant nodes and the datacenters they belong to
    std::set<int> connected_datacenters;
    std::vector<std::pair<Node /*node*/, int /*datacenter*/>> nodes_needed;
    for (auto node : datacenter->get_map_out()) {
        // If it is a distant node
        if (partition[node.first] != datacenter->get_datacenter_id()) {
            connected_datacenters.insert(partition[node.first]);
            nodes_needed.emplace_back(node.first, partition[node.first]);
        }
    }
    // For each datacenter wait for the PRs value of the distant nodes,
    // the datacenter send a masked CKKS vector of PRs, and the corresponding index of the nodes in it
    std::map<int /*datacenter*/, std::pair<std::vector<Ciphertext<DCRTPoly>>, std::map<Node, unsigned int>>>
        distant_values;
    for (int datacenter_id : connected_datacenters) {
        distant_values[datacenter_id] =
            datacenters[datacenter_id]->get_PR_distant_nodes(current_iteration, datacenter->get_datacenter_id());
    }
    // For each distant nodes get the corresponding masked CKKS vector with the corresponding indexes of the nodes
    // in it
    std::vector<NodesInCKKSVector> PRs_distant_nodes;
    for (auto dc : distant_values) {
        int datacenter_id = dc.first;
        auto distant_vect = dc.second.first;
        auto nodesIndex = dc.second.second;
        std::multimap<int /*vect number*/, std::pair<Node, unsigned int>> to_add;
        for (auto node : nodesIndex) {
            bool exist = false;
            for (auto node_needed : nodes_needed) {
                if (node_needed.first == node.first && node_needed.second == datacenter_id) {
                    exist = true;
                    break;
                }
            }
            if (exist) {
                to_add.insert({node.second / datacenters.at(datacenter_id)->get_vector_size(), node});
            }
        }
        auto it = to_add.begin();
        auto it2 = to_add.begin();
        while (it != to_add.end()) {
            std::set<std::pair<Node, unsigned int>> nodes_indexes_set;
            while (it2 != to_add.end() && it2->first == it->first) {
                nodes_indexes_set.insert(it2->second);
                it2++;
            }
            PRs_distant_nodes.push_back({distant_vect[it->first], nodes_indexes_set});
            it = it2;
        }
    }
    return PRs_distant_nodes;
}

/**
 * @brief Global run function to compute the PRs of the datacenters
 * At the beginning of each iteration, we synchronize the PRs of the distant nodes
 * At the last iteration, we compute the PRs of the dangling nodes
 *
 * @param nb_iteration
 * @param datacenters
 * @param partition
 */
void run(
    int nb_iteration,
    std::vector<std::shared_ptr<POPPYmMultiDC>> datacenters,
    std::map<int, int> partition,
    int multDepth) {
    for (int iter = 1; iter < nb_iteration; iter++) {
        std::cout << "iteration " << iter << std::endl;
#ifdef EXPE
        expe_poppy_m.setCurrentRun(iter);
#endif
        for (auto& datacenter : datacenters) {
#ifdef EXPE
            int tmp_add_scalar, tmp_add, tmp_rot, tmp_mult, tmp_mult_scalar, tmp_encryption, tmp_encoding,
                tmp_mult_plain;
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_rot = number_rot;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;
            auto startSync = std::chrono::high_resolution_clock::now();
#endif
            // Synchro
            auto PRs_distant_nodes = synch_data(datacenter, datacenters, partition, iter - 1);
            if (PRs_distant_nodes.size() > 0)
                datacenter->synchronize_distant_PRs(PRs_distant_nodes);
#ifdef EXPE
            auto endSync = std::chrono::high_resolution_clock::now();
            auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
            expe_poppy_m.addValue(
                durationSync.count(), "Time to synchronize distant PRs", Experiment::SYNC, Experiment::SECONDS);
            expe_poppy_m.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition during synchronization " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_add - tmp_add,
                "Number of addition during synchronization " + std::to_string(iter),
                Experiment::Operation::ADD_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during synchronization " + std::to_string(iter),
                Experiment::Operation::ENCODE_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during synchronization " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_rot - tmp_rot,
                "Number of rotation during synchronization " + std::to_string(iter),
                Experiment::Operation::ROT_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                PRs_distant_nodes.size(),
                "Number of received ciphertext during synchronization",
                Experiment::Operation::CIPHER_SYNC,
                Experiment::Unit::UNIT);
#endif
            for (auto& ct : PRs_distant_nodes) {  // Free memory of PRs_distant_nodes
                ct.vect.reset();
            }
            PRs_distant_nodes.clear();
            PRs_distant_nodes.shrink_to_fit();
#ifdef EXPE
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_rot = number_rot;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;
            auto startCompute = std::chrono::high_resolution_clock::now();
#endif
            // Computation
            datacenter->compute_PR();
#ifdef EXPE
            auto endCompute = std::chrono::high_resolution_clock::now();
            auto durationCompute = std::chrono::duration_cast<std::chrono::seconds>(endCompute - startCompute);
            expe_poppy_m.addValue(
                durationCompute.count(), "Time to compute PRs", Experiment::EVAL_OP, Experiment::SECONDS);
            expe_poppy_m.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition after computation " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_add - tmp_add,
                "Number of addition after computation " + std::to_string(iter),
                Experiment::Operation::ADD_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding after computation " + std::to_string(iter),
                Experiment::Operation::ENCODE_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption after computation " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult - tmp_mult,
                "Number of multiplication after computation " + std::to_string(iter),
                Experiment::Operation::MUL_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_rot - tmp_rot,
                "Number of rotation after computation " + std::to_string(iter),
                Experiment::Operation::ROT_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of scalar multiplication after computation " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain after computation " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_COMPUTE,
                Experiment::Unit::UNIT);
#endif

#ifdef BOOTSTRAPPING
#ifdef EXPE
            auto startBootstrap = std::chrono::high_resolution_clock::now();
#endif
            if (iter >= multDepth / POPPY_m_LEVELS_CONSUMPTION)
                datacenter->bootstrap(multDepth);
#ifdef EXPE
            auto endBootstrap = std::chrono::high_resolution_clock::now();
            auto durationBootstrap =
                std::chrono::duration_cast<std::chrono::milliseconds>(endBootstrap - startBootstrap);
            expe_poppy_m.addValue(
                durationBootstrap.count(), "Time to bootstrap", Experiment::BOOTSTRAP, Experiment::MILLISECONDS);
#endif
#endif

#ifdef EXPE
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_rot = number_rot;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;
            auto startPrepare = std::chrono::high_resolution_clock::now();
#endif
            // Preparation next iteration
            datacenter->copy_previous_masked_distant_node();
            datacenter->compute_masked_distant_node(iter);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
            datacenter->add_constant_to_distant_nodes(iter);
#endif

#ifdef EXPE
            auto endPrepare = std::chrono::high_resolution_clock::now();
            auto durationPrepare = std::chrono::duration_cast<std::chrono::milliseconds>(endPrepare - startPrepare);
            expe_poppy_m.addValue(
                durationPrepare.count(), "Time to prepare next iteration", Experiment::PREP, Experiment::MILLISECONDS);
            expe_poppy_m.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of addition during prepare " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_add - tmp_add,
                "Number of addition during prepare " + std::to_string(iter),
                Experiment::Operation::ADD_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during prepare " + std::to_string(iter),
                Experiment::Operation::ENCODE_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during prepare " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_rot - tmp_rot,
                "Number of rotation during prepare " + std::to_string(iter),
                Experiment::Operation::ROT_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_m.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_PREP,
                Experiment::Unit::UNIT);
#endif
        }
#ifdef EXPE
        expe_poppy_m.addValue(
            number_add_scalar,
            "Number of scalar addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD_SCALAR,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_add,
            "Number of addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encoding,
            "Number of encoding after iteration " + std::to_string(iter),
            Experiment::Operation::ENCODE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encryption,
            "Number of encryption after iteration " + std::to_string(iter),
            Experiment::Operation::ENCRYPT,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult,
            "Number of multiplication after iteration " + std::to_string(iter),
            Experiment::Operation::MUL,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_rot,
            "Number of rotation after iteration " + std::to_string(iter),
            Experiment::Operation::ROT,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_scalar,
            "Number of multiplication scalar after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_SCALAR,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_plain,
            "Number of multiplication plain after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_PLAIN,
            Experiment::Unit::UNIT);
#endif
    }
    std::cout << "Last iteration " << nb_iteration << std::endl;
#ifdef EXPE
    expe_poppy_m.setCurrentRun(nb_iteration);
#endif
    for (auto& datacenter : datacenters) {
#ifdef EXPE
        int tmp_add_scalar, tmp_add, tmp_rot, tmp_mult, tmp_mult_scalar, tmp_encoding, tmp_encryption, tmp_mult_plain;
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_rot = number_rot;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
        tmp_mult_plain = number_mult_plain;
        auto startSync = std::chrono::high_resolution_clock::now();
#endif

        // Synchro
        auto PRs_distant_nodes = synch_data(datacenter, datacenters, partition, nb_iteration - 1);
        datacenter->synchronize_distant_PRs(PRs_distant_nodes);
#ifdef EXPE
        auto endSync = std::chrono::high_resolution_clock::now();
        auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
        expe_poppy_m.addValue(
            durationSync.count(),
            "Time to synchronize distant PRs at last iteration",
            Experiment::SYNC,
            Experiment::SECONDS);
        expe_poppy_m.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition during synchronization at last iteration",
            Experiment::Operation::ADD_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_add - tmp_add,
            "Number of addition during synchronization at last iteration",
            Experiment::Operation::ADD_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding during synchronization at last iteration",
            Experiment::Operation::ENCODE_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption during synchronization at last iteration",
            Experiment::Operation::ENCRYPT_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult - tmp_mult,
            "Number of multiplication during synchronization at last iteration",
            Experiment::Operation::MUL_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_rot - tmp_rot,
            "Number of rotation during synchronization at last iteration",
            Experiment::Operation::ROT_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar during synchronization at last iteration",
            Experiment::Operation::MUL_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain during synchronization at last iteration",
            Experiment::Operation::MUL_PLAIN_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            PRs_distant_nodes.size(),
            "Number of received ciphertext during synchronization",
            Experiment::Operation::CIPHER_SYNC,
            Experiment::Unit::UNIT);
#endif
        for (auto& ct : PRs_distant_nodes) {  // Free memory of PRs_distant_nodes
            ct.vect.reset();
        }
        PRs_distant_nodes.clear();
        PRs_distant_nodes.shrink_to_fit();
#ifdef EXPE
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_rot = number_rot;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
        tmp_mult_plain = number_mult_plain;
        auto startCompute = std::chrono::high_resolution_clock::now();
#endif
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        if (datacenter->get_nb_dangling_nodes() > 0) {
            datacenter->compute_PR_dangling_node();
        }
#endif
        datacenter->compute_PR();
#ifdef EXPE
        auto endCompute = std::chrono::high_resolution_clock::now();
        auto durationCompute = std::chrono::duration_cast<std::chrono::seconds>(endCompute - startCompute);
        expe_poppy_m.addValue(
            durationCompute.count(), "Time to compute PRs at last iteration", Experiment::EVAL_OP, Experiment::SECONDS);
        expe_poppy_m.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition after last iteration",
            Experiment::Operation::ADD_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_add - tmp_add,
            "Number of addition after last iteration",
            Experiment::Operation::ADD_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding after last iteration",
            Experiment::Operation::ENCODE_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption after last iteration",
            Experiment::Operation::ENCRYPT_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult - tmp_mult,
            "Number of multiplication after last iteration",
            Experiment::Operation::MUL_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_rot - tmp_rot,
            "Number of rotation after last iteration",
            Experiment::Operation::ROT_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar after last iteration",
            Experiment::Operation::MUL_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain after last iteration",
            Experiment::Operation::MUL_PLAIN_COMPUTE,
            Experiment::Unit::UNIT);
#endif
#ifdef EXPE
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_rot = number_rot;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
        tmp_mult_plain = number_mult_plain;
        auto startPrepare = std::chrono::high_resolution_clock::now();
#endif
        datacenter->copy_previous_masked_distant_node();
#ifdef EXPE
        auto endPrepare = std::chrono::high_resolution_clock::now();
        auto durationPrepare = std::chrono::duration_cast<std::chrono::milliseconds>(endPrepare - startPrepare);
        expe_poppy_m.addValue(
            durationPrepare.count(), "Time to prepare at last iteration", Experiment::PREP, Experiment::MILLISECONDS);
        expe_poppy_m.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition prepare at last iteration",
            Experiment::Operation::ADD_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_add - tmp_add,
            "Number of addition prepare at last iteration",
            Experiment::Operation::ADD_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding prepare at last iteration",
            Experiment::Operation::ENCODE_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption prepare at last iteration",
            Experiment::Operation::ENCRYPT_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult - tmp_mult,
            "Number of multiplication prepare at last iteration",
            Experiment::Operation::MUL_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_rot - tmp_rot,
            "Number of rotation prepare at last iteration",
            Experiment::Operation::ROT_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar prepare at last iteration",
            Experiment::Operation::MUL_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_m.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain prepare at last iteration",
            Experiment::Operation::MUL_PLAIN_PREP,
            Experiment::Unit::UNIT);
#endif
    }
}

/**
 * @brief Main function to compute the PRs of the datacenters in the POPPYm multi DC scenario
 * First we parse graph of the datacenters to obtain multiple graphs
 * We let the datacenters initialize their local graph
 * We generate a common crypto context for all the datacenters
 * We run the PR algorithm
 */
void run_multi_dc_poppy_m(int nb_iteration) {
#ifdef EXPE
    number_add_scalar = 0;
    number_add = 0;
    number_mult = 0;
    number_rot = 0;
    number_mult_scalar = 0;
    number_mult_plain = 0;
    number_encoding = 0;
    number_encryption = 0;
#endif
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC("../dot/students-partition1.dot", partition);
#ifdef DEBUG
    for (auto graph : graphs) {
        std::cout << "Graph" << std::endl;
        for (auto node : graph.nodes) {
            std::cout << node << ", ";
        }
        std::cout << std::endl;
        for (auto edge : graph.arcs) {
            std::cout << edge.first << " -> " << edge.second << std::endl;
        }
    }
#endif
    std::cout << "start multi DC coeus initialization" << std::endl;

    // Initialize the datacenters
    std::vector<int> max_deg_in;
    std::vector<std::shared_ptr<POPPYmMultiDC>> datacenters;
    std::vector<int> number_actives_nodes;
    int datacenter_num = 0;
#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<POPPYmMultiDC>(datacenter_num));
        number_actives_nodes.push_back(datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci));
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }

#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time to init all datacenters",
        Experiment::Operation::INIT_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    std::cout << "start gen_crypto_context" << std::endl;
    // Get the max in number_actives_nodes
    int max_number_actives_nodes = *std::max_element(number_actives_nodes.begin(), number_actives_nodes.end());

    // Compute s, the size of the vectors
    uint32_t nbSlots = nextPowerOfTwo(max_number_actives_nodes);
    int multDepth = POPPY_m_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {3, 3};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = POPPY_m_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif

    auto cc_test = gen_crypto_context_test(multDepth);
    int s;
    int ringDim = cc_test->GetRingDimension();
    if (nbSlots < ringDim) {
        s = nbSlots;
        std::cout << "s : " << s << " (the next power of two of " << max_number_actives_nodes
                  << ", the number of active nodes)" << std::endl;

    } else {
        s = ringDim / 2;
        std::cout << "s : " << s << " (ringDim / 2)" << std::endl;
    }

#ifdef DEBUG
    std::cout << "The common batch size is " << s << std::endl;
    std::cout << "The multiplicative depth is " << multDepth << std::endl;
#endif

    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);
#ifdef COEUS
        gen_rotation_keys_coeus(cc, keys, nbSlots);
#else
        gen_rotation_keys(cc, keys, nbSlots);
#endif
#ifdef BOOTSTRAPPING
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, s);
        cc->EvalBootstrapKeyGen(keys.secretKey, s);
#endif
#ifdef SERIALIZE
        serialize_cc(cc, keys, multDepth);
#endif
    }
#ifdef EXPE
    expe_poppy_m.setRingsize(cc->GetRingDimension());
    start = std::chrono::high_resolution_clock::now();
#endif
    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_number_datacenters(datacenter_num);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time for Preprocessing step all datacenters",
        Experiment::Operation::PREPROC_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    // Run
    std::cout << "start multi DC coeus computation" << std::endl;
#ifdef EXPE
    start = std::chrono::high_resolution_clock::now();
#endif
    run(nb_iteration, datacenters, partition, multDepth);
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count(),
        "Time to compute " + std::to_string(nb_iteration) + " iterations",
        Experiment::Operation::EVAL_DC,
        Experiment::Unit::SECONDS);
#endif

    std::cout << "start multi DC coeus PRs extraction" << std::endl;
    std::map<int, double> results;
    for (auto& datacenter : datacenters) {
        std::cout << "Datacenter " << datacenter->get_datacenter_id() << std::endl;
        results.merge(datacenter->export_results());

        // datacenter->print_PRs();
        datacenter->delete_PRs();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        datacenter->delete_PRs_deg0();
#endif
    }

    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalSumKeys();
}

std::map<int, double> run_multi_dc_poppy_m(int nb_iteration, string filename) {
#ifdef EXPE
    number_add_scalar = 0;
    number_add = 0;
    number_mult = 0;
    number_rot = 0;
    number_mult_scalar = 0;
    number_mult_plain = 0;
    number_encoding = 0;
    number_encryption = 0;
#endif
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC(filename, partition);

    std::cout << "start gen_crypto_context" << std::endl;
    CryptoContext<DCRTPoly> ccTest;
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    int multDepth = POPPY_m_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {3, 3};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = POPPY_m_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif
    ccTest = gen_crypto_context_test(multDepth);

    std::cout << "start multi DC coeus initialization" << std::endl;

    // Initialize the datacenters
    std::vector<int> max_deg_in;
    std::vector<std::shared_ptr<POPPYmMultiDC>> datacenters;
    std::vector<int> number_actives_nodes;
    int datacenter_num = 0;
#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<POPPYmMultiDC>(datacenter_num));
        number_actives_nodes.push_back(datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci));
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }

#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time to init all datacenters",
        Experiment::Operation::INIT_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    // Get the max in number_actives_nodes
    int max_number_actives_nodes = *std::max_element(number_actives_nodes.begin(), number_actives_nodes.end());

    // Compute s, the size of the vectors
    uint32_t nbSlots = nextPowerOfTwo(max_number_actives_nodes);

    usint ringDim = ccTest->GetRingDimension();
    int s;
    if (nbSlots < ringDim) {
        s = nbSlots;
        std::cout << "s : " << s << " (the next power of two of " << max_number_actives_nodes
                  << ", the number of active nodes)" << std::endl;

    } else {
        s = ringDim / 2;
        std::cout << "s : " << s << " (ringDim / 2)" << std::endl;
    }

    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);
#ifdef COEUS
        gen_rotation_keys_coeus(cc, keys, nbSlots);
#else
        gen_rotation_keys(cc, keys, nbSlots);
#endif
#ifdef BOOTSTRAPPING
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, s);
        cc->EvalBootstrapKeyGen(keys.secretKey, s);
#endif
#ifdef SERIALIZE
        serialize_cc(cc, keys, multDepth);
#endif
    }
#ifdef EXPE
    expe_poppy_m.setRingsize(cc->GetRingDimension());
    expe_poppy_m.setCurrentRun(0);
    start = std::chrono::high_resolution_clock::now();
#endif

    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_number_datacenters(datacenter_num);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time for Preprocessing step all datacenters",
        Experiment::Operation::PREPROC_DC,
        Experiment::Unit::MILLISECONDS);
    expe_poppy_m.addValue(
        number_add_scalar,
        "Number of scalar additions after the preprocessing",
        Experiment::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_add, "Number of additions after the preprocessing", Experiment::ADD, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult, "Number of multiplications after the preprocessing", Experiment::MUL, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_rot, "Number of rotations after the preprocessing", Experiment::ROT, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult_scalar,
        "Number of scalar multiplications after the preprocessing",
        Experiment::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult_plain,
        "Number of multiplication plain after the preprocessing",
        Experiment::Operation::MUL_PLAIN,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_encoding, "Number of encoding after the preprocessing", Experiment::ENCODE, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_encryption,
        "Number of encryptions after the preprocessing",
        Experiment::ENCRYPT,
        Experiment::Unit::UNIT);
#endif

    // Run
    std::cout << "start multi DC coeus computation" << std::endl;
#ifdef EXPE
    start = std::chrono::high_resolution_clock::now();
#endif
    run(nb_iteration, datacenters, partition, multDepth);
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_m.addValue(
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count(),
        "Time to compute " + std::to_string(nb_iteration) + " iterations",
        Experiment::Operation::EVAL_DC,
        Experiment::Unit::SECONDS);
    expe_poppy_m.addValue(
        number_add_scalar,
        "Number of scalar addition after last iteration",
        Experiment::Operation::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_add, "Number of addition after last iteration", Experiment::Operation::ADD, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_encoding,
        "Number of encoding after last iteration",
        Experiment::Operation::ENCODE,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_encryption,
        "Number of encryption after last iteration",
        Experiment::Operation::ENCRYPT,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult,
        "Number of multiplication after last iteration",
        Experiment::Operation::MUL,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_rot, "Number of rotation after last iteration", Experiment::Operation::ROT, Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult_scalar,
        "Number of multiplication scalar after last iteration",
        Experiment::Operation::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_m.addValue(
        number_mult_plain,
        "Number of multiplication plain after last iteration",
        Experiment::Operation::MUL_PLAIN,
        Experiment::Unit::UNIT);
#endif

    std::map<int, double> results;
    for (auto& datacenter : datacenters) {
        std::cout << "Datacenter " << datacenter->get_datacenter_id() << std::endl;
        results.merge(datacenter->export_results());

        // datacenter->print_PRs();

        // Free memory
        datacenter->delete_PRs();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        datacenter->delete_PRs_deg0();
#endif
    }

#ifdef EXPE
    expe_poppy_m.disable_crypto_context();
#endif
    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalSumKeys();

    return results;
}

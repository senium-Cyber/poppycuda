#include <sys/resource.h>
#include <iostream>
#include <thread>
#include "openfhe.h"
#include "poppy_uj.cuh"

#ifdef EXPE
#include "utils/include/utils_export.cuh"
Experiment expe_poppy_uj;
#endif

void print_memory(std::string message) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << message << std::endl;
    std::cout << "Mémoire utilisée : " << usage.ru_maxrss / 1024 << " MB" << std::endl;  // ru_maxrss is in MB
}

#include <fstream>
void right_file(string file_name, std::vector<double> result) {
    std::ofstream outFile(file_name, std::ios::app);
    if (outFile.is_open()) {
        outFile << result << "\n";
        outFile.close();
    }
}

int POPPYujMultiDC::init_local_graph(
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

    // Initialize nb_dangling_nodes to zero
    this->nb_dangling_nodes = 0;

#ifdef COMPUTE_ORBITS
    std::tie(map_in, map_out) = sparse_nauty_routine(
        &subgraph, &partition, nb_local_nodes, nb_local_nodes_after_pruning, orbits, this->datacenter_id, duration, set_dangling_nodes);

    this->nb_dangling_nodes = set_dangling_nodes.size();
#elif COMPUTE_EQUIVALENT_NODES
    std::tie(map_in, map_out) = equivalent_nodes_routine(
        &subgraph,
        &partition,
        nb_local_nodes,
        nb_local_nodes_after_pruning,
        orbits,
        this->datacenter_id,
        duration,
        set_dangling_nodes);

    this->nb_dangling_nodes = set_dangling_nodes.size();
#else
    std::tie(map_in, map_out) = process_maps(&subgraph, &partition, nb_local_nodes, orbits, this->datacenter_id);
    nb_local_nodes_after_pruning = nb_local_nodes;
#endif

#ifdef EXPE
    expe_poppy_uj.addValue(
        duration,
        "Time to compute equivalent nodes on DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::MICROSECONDS);
    expe_poppy_uj.addValue(
        nb_local_nodes,
        "Number of nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_poppy_uj.addValue(
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

    // Set representative orbits
    for (auto orbit : orbits) {
        if (representative_orbits.find(orbit.second) == representative_orbits.end())
            representative_orbits.insert({orbit.second, orbit.first});
    }

    // Set the uj size
    int nb_constant_nodes = 0;
    int active_uj_size = set_u_j_size(
        this->orbits, map_in, map_out, partition, datacenter_id, nb_constant_nodes, this->nb_dangling_nodes);
    this->uj_size = (active_uj_size > nb_dangling_nodes) ? active_uj_size : nb_dangling_nodes;

#ifdef EXPE
    expe_poppy_uj.addValue(
        nb_constant_nodes,
        "Number of constant nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_poppy_uj.addValue(
        nb_dangling_nodes,
        "Number of dangling nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
#endif

    // Now we create maps for ujs
    this->nb_of_uj = creation_map_u_j(
        this->map_uj,
        this->orbits,
        map_in,
        map_out,
        partition,
        this->datacenter_id,
        representative_orbits,
        set_dangling_nodes);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    // And the map_deg_out_0
    this->nb_u_deg0 = creation_map_deg_out_0(
        this->map_deg0,
        this->orbits,
        map_in,
        map_out,
        partition,
        this->datacenter_id,
        representative_orbits,
        set_dangling_nodes);
#endif

    // Count the number of nodes connected to distant datacenters
    std::unordered_map<int, std::set<Node>> all_nodes_connected_to_distant;
    for (auto arc : map_out) {
        if (partition.at(arc.second) != datacenter_id) {
            all_nodes_connected_to_distant[partition.at(arc.second)].insert(
                arc.first);                      // We regroupe by distant DC all node connected to a distant node
            if (map_in.count(arc.first) == 0) {  // Constants noeuds are always removed from POPPYuj PRs
                constant_nodes_connected_to_distant.insert({arc.first, partition.at(arc.second)});
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
    std::cout << "Map u_j datacenter " << this->datacenter_id << " : " << std::endl;
    print_map_u_j(map_uj);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::cout << "and Map deg0 datacenter " << this->datacenter_id << " : " << std::endl;
    print_map_u_j(map_deg0);
#endif
#endif

#ifdef DEBUG
    printNodeCategories();
#endif
    // In the inter DC communication, the pagerank of the constant nodes are put in the uj's"
    return std::max(uj_size, max_nodes_connected_to_distant);
}

/**
 * @brief Preprocessing of the datacenter
 * Encrypt all the mask
 * Initialize the uj's
 *
 */
void POPPYujMultiDC::pre_processing() {
    this->m = this->s;
    while (this->m < this->uj_size) {
        this->m += s;
    }

    creation_map_reconstruction_u(
        map_uj,
        map_uj,
        map_reconstruction_uj,
        this->s,
        this->m,
        representative_orbits,
        this->partition,
        this->datacenter_id);

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    creation_map_reconstruction_u(
        this->map_uj,
        this->map_deg0,
        map_construction_u_deg0,
        this->s,
        this->m,
        representative_orbits,
        this->partition,
        this->datacenter_id);
#endif
#ifdef DEBUG
    std::cout << "Map reconstruction u_j datacenter " << this->datacenter_id << " : " << std::endl;
    print_map_reconstruction_u_j(map_reconstruction_uj);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::cout << "and map reconstruction deg0" << std::endl;
    print_map_reconstruction_u_j(map_construction_u_deg0);
#endif
#endif

    // Compute the mask of the effective distant nodes
    std::vector<double> zeros(this->s, 0.0);
    std::map<int /*data center*/, std::vector<std::vector<double>>> border_masks;
    std::vector<std::vector<double>> masks(m / s, zeros);
    std::map<int /*data center*/, std::vector<std::vector<double>>> vect_border_const_PRs;
    std::map<int /*data center*/, std::vector<std::vector<double>>> vect_border_init_const_PRs;

    for (auto node : nodes_connected_to_distant) {
        // if the mask for a distant data center is not already creat
        if (border_masks.count(node.second) == 0) {
            border_masks.insert({node.second, masks});
        }
        // if the key of a distant dc is not already in the set of indexes
        if (index_distant_nodes.count(node.second) == 0) {
            index_distant_nodes[node.second] = {};  // creationx of an empty map
        }
        if (map_uj.find(node.first) == map_uj.end()) {
            std::cout << "Error: node.first " << node.first << " not found in map_uj." << std::endl;
        }
        auto elt = map_uj.find(node.first);
        assert(elt != map_uj.end());
        triplet value = elt->second;
        std::vector<std::vector<double>> current_mask = border_masks.at(node.second);
        current_mask[value.indice / s][value.indice % s] = 1.0 / map_out.count(node.first);
        index_distant_nodes[node.second][node.first] =
            value.indice;  // we add the value node.first -> value.indice in the map in node.second

        border_masks.at(node.second) = current_mask;
    }
    // Fill empty slots with the constant nodes (for initialization and other iterations)
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
        }
        this->mask_effective_distant_nodes.insert({masks.first, plains});
    }

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
        this->init_constant_PRs.insert({masks.first, plain});
    }

    // Initialize the PRs and encrypt them
    std::vector<double> PRs_vect(this->s, 1.0 / number_of_nodes);
    Plaintext plain_PRs = cc->MakeCKKSPackedPlaintext(PRs_vect);
    auto ct_iv = cc->Encrypt(keys.publicKey, plain_PRs);
    this->PRs = std::vector<Ciphertext<DCRTPoly>>(m / s, ct_iv);
#ifdef EXPE
    number_encoding++;
    number_encryption++;
#endif

    // Initialize masked_distant_node to compute masked PRs
    compute_masked_distant_node(0);
    add_constant_to_distant_nodes(0);  // Add constant nodes at initialization (iteration 0)

    // Creation vect PR const
    // Constant nodes are always pruned in POPPY uj
    creation_vect_PR_const();

    // Initialize uj
    initialize_uj();
}

/**
 * @brief Set the common crypto context and keys
 *
 * @param cc Common crypto context
 * @param keys
 * @param vector_size Common batch size
 */
void POPPYujMultiDC::set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size) {
    this->cc = cc;
    this->keys = keys;
    this->s = vector_size;
}

/**
 * @brief Initialize the uj's for the first iteration
 *
 */
void POPPYujMultiDC::initialize_uj() {
    std::vector<double> sub_vect(this->s, 0.0);
    std::vector<std::vector<double>> iv(m / s, sub_vect);

    auto it = map_uj.begin();
    auto it2 = map_uj.begin();
    while (it != map_uj.end()) {
        while (it2 != map_uj.end() && it2->first == it->first) {
            triplet value = it2->second;
            assert(value.indice < uj_size);
            // if distant then can't be init (because we don't know their deg out, synced later)
            // let them to 0
            // if constant we initialize it in u_const
            if (partition.at(representative_orbits.at(value.vertex)) == datacenter_id &&
                map_in.count(representative_orbits.at(value.vertex)) > 0) {
                assert(value.deg_out != 0);
                iv[value.indice / s][value.indice % s] += 1.0 / (number_of_nodes * value.deg_out);
            }
            it2++;
        }

        it = it2;
    }

    // instanciate u a vector of m / s Ciphertexts
    this->u = std::vector<Ciphertext<DCRTPoly>>(m / s);
    for (int i = 0; i < m / s; i++) {
        Plaintext plain = cc->MakeCKKSPackedPlaintext(iv[i]);
        this->u[i] = cc->Encrypt(keys.publicKey, plain);
#ifdef EXPE
        number_encoding++;
        number_encryption++;
#endif
    }
}

/**
 * @brief  Loop over the distant nodes comming in and get their PRs
 *  This phase is done at the beginning of each iteration to receive the PRs of the distant nodes
 *
 * @param datacenter Datacenter to synchronize
 * @param datacenters Other datacenters
 * @param partition
 * @param current_iteration
 * @return std::vector<NodesInCKKSVector>
 */
std::vector<NodesInCKKSVector> synch_data(
    std::shared_ptr<POPPYujMultiDC> datacenter,
    std::vector<std::shared_ptr<POPPYujMultiDC>> datacenters,
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

    // For each distant nodes get the corresponding masked CKKS vector with the corresponding indexes of the nodes in it
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
 * @brief Mask the PRs before sending them to the other datacenters
 *
 */
void POPPYujMultiDC::compute_masked_distant_node(int iteration) {
    for (auto distant_dc : mask_effective_distant_nodes) {
        if (iteration == 0) {  // We need to creat an empty vector at first iteration
            std::vector<Ciphertext<DCRTPoly>> masked_PRs(m / s);
            masked_distant_nodes.insert({distant_dc.first, masked_PRs});
        }
        std::vector<Ciphertext<DCRTPoly>> current_mask = masked_distant_nodes.at(distant_dc.first);
        for (int i = 0; i < m / s; i++) {
            current_mask[i] = cc->EvalMult(this->PRs[i], distant_dc.second[i]);
#ifdef EXPE
            number_mult_plain++;
#endif
        }
        masked_distant_nodes.at(distant_dc.first) = current_mask;
        delete_vector_of_ciphers(current_mask);
    }
}

/**
 * @brief Add to the masked_distant_nodes, the PRs of the constant nodes
 * These Ciphertexts are then sent to the distant DCs
 * @param iteration
 */
void POPPYujMultiDC::add_constant_to_distant_nodes(int iteration) {
    if (iteration == 0) {
        for (auto distant_dc : masked_distant_nodes) {
            std::vector<Ciphertext<DCRTPoly>> current_mask = masked_distant_nodes.at(distant_dc.first);
            for (int i = 0; i < m / s; i++) {
                if (init_constant_PRs.count(distant_dc.first) != 0) {
                    current_mask[i] = cc->EvalAdd(init_constant_PRs.at(distant_dc.first)[i], distant_dc.second[i]);
#ifdef EXPE
                    number_add++;
#endif
                }
            }
            masked_distant_nodes.at(distant_dc.first) = current_mask;
            delete_vector_of_ciphers(current_mask);
        }
        delete_map_of_plains(init_constant_PRs);
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
            delete_vector_of_ciphers(current_mask);
        }
    }
}

void POPPYujMultiDC::printNodeCategories() {
    std::cout << "\n--- Datacenter " << this->datacenter_id << " Node Categories ------------------\n";

    // 1) Collect sets of active
    std::set<int> activeNodes;
    for (auto& kv : map_uj) {
        activeNodes.insert(kv.first);
    }

    // 2) Print active nodes
    std::cout << "Active nodes (" << activeNodes.size() << "): ";
    for (auto node : activeNodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // 3) Collect dangling nodes
    //    We'll iterate through all nodes owned by this DC (based on 'partition'),
    //    and pick those with out-degree 0.
    std::set<int> danglingNodes;
    for (auto& kv : map_deg0) {
        danglingNodes.insert(kv.first);
    }

    // 4) Print dangling nodes
    std::cout << "Dangling nodes (" << danglingNodes.size() << "): ";
    for (auto node : danglingNodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // 5) Collect constant nodes:
    //    We'll iterate through all nodes owned by this DC (based on 'partition'),
    //    and pick those with in-degree 0.
    std::vector<int> constantNodes;
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    for (auto& kv : partition) {
        int node = kv.first;
        int dc = kv.second;
        if (dc == this->datacenter_id) {
            // Belongs to this DC
            if (map_in.count(node) == 0) {
                constantNodes.push_back(node);
            }
        }
    }
#endif
    // 6) Print constant nodes
    std::cout << "Constant nodes (" << constantNodes.size() << "): ";
    for (auto node : constantNodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "---------------------------------------------------\n\n";
}

/**
 * @brief Send the PRs to the other datacenters
 * @param iteration
 * @return std::pair<Ciphertext<DCRTPoly>, std::map<Node, unsigned int>>
 */
std::pair<std::vector<Ciphertext<DCRTPoly>> /*cipher value*/, std::map<Node, unsigned int> /*indexes*/>
POPPYujMultiDC::get_PR_distant_nodes(int iteration, int distant_dc) {
    if (current_iteration > iteration) {
        return {previous_masked_distant_nodes.at(distant_dc), index_distant_nodes.at(distant_dc)};
    } else {
        return {masked_distant_nodes.at(distant_dc), index_distant_nodes.at(distant_dc)};
    }
}

/**
 * @brief Print the PRs
 *
 */
void POPPYujMultiDC::print_PRs() {
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
    std::cout << "Vect PRs: " << vect_PRs << std::endl;
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::cout << "Vect PRs deg0: " << vect_PRs_deg0 << std::endl;
#endif
#endif

    // Print the PRs in the order of the nodes
    for (auto orbit : orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;

        if (map_uj.count(representative_orbits.at(orbit.second)) != 0) {  // Check if the node is an active node
            int indice = map_uj.find(representative_orbits.at(orbit.second))->second.indice;
            std::cout << orbit.first << " is " << vect_PRs[indice] << " at " << indice << std::endl;
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (this->map_deg0.size() > 0) {  // Check if the node is a dangling node
            if (map_deg0.count(representative_orbits[orbit.second]) != 0) {
                int indice = map_deg0.find(representative_orbits[orbit.second])->second.indice;
                std::cout << orbit.first << " is " << vect_PRs_deg0[indice] << " (dangling) at " << indice << std::endl;
            } else
                std::cout << orbit.first << " is " << 0.15 / this->max_deg_in << " (constant)" << std::endl;
        }
#endif
        else  // Otherwise it is a constant node
            std::cout << orbit.first << " is " << 0.15 / this->max_deg_in << " (constant)" << std::endl;
    }

    // Clear vectors of cleartexts
    vect_PRs.clear();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    vect_PRs_deg0.clear();
#endif
}

std::map<int, double> POPPYujMultiDC::export_results() {
    std::map<int, double> results;

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

    for (auto orbit : orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;

        if (map_uj.count(representative_orbits[orbit.second]) != 0) {
            int indice = map_uj.find(representative_orbits[orbit.second])->second.indice;
            results.insert({orbit.first, vect_PRs[indice]});
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (map_deg0.count(representative_orbits[orbit.second]) != 0) {
            if (vect_PRs_deg0.size() > 0) {  // Check if there is dangling nodes
                int indice = map_deg0.find(representative_orbits[orbit.second])->second.indice;
                results.insert({orbit.first, vect_PRs_deg0[indice]});
            }
        }
#endif
        else
            results.insert({orbit.first, 0.15 / this->max_deg_in});
    }
    // Clear vectors of cleartexts
    vect_PRs.clear();
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    vect_PRs_deg0.clear();
#endif

    return results;
}

/**
 * @brief Create the u vector of constant PRs
 *
 */
void POPPYujMultiDC::creation_vect_PR_const() {
    std::vector<double> zero_vect(this->s, 0.0);
    std::vector<std::vector<double>> u_const_clear(m / s, zero_vect);
    std::vector<std::vector<double>> u_const_clear_init(m / s, zero_vect);
    // For all constant node appearing in the uj's at indice create a vector that will be added to the uj's
    for (auto& it : map_uj) {
        triplet value = it.second;
        if (map_uj.count(representative_orbits.at(value.vertex)) == 0 &&
            partition.at(representative_orbits.at(value.vertex)) == datacenter_id) {
            // value.vertex is a constant node
            assert(value.indice < uj_size);
            u_const_clear[value.indice / s][value.indice % s] += 0.15 / (double)(this->max_deg_in * value.deg_out);
            u_const_clear_init[value.indice / s][value.indice % s] += 1.0 / (double)(number_of_nodes * value.deg_out);
        }
    }
#ifdef DEBUG
    std::cout << "Datacenter " << datacenter_id << " u_const : " << u_const_clear << std::endl;
    std::cout << "Datacenter " << datacenter_id << " u_const_init : " << u_const_clear_init << std::endl;
#endif
    // Encode the constant vectors
    this->u_const = std::vector<Plaintext>(m / s);
    for (int i = 0; i < m / s; i++) {
        if (u_const_clear[i] != zero_vect) {
            this->u_const[i] = cc->MakeCKKSPackedPlaintext(u_const_clear[i]);
#ifdef EXPE
            number_encoding++;
#endif
        } else {
            // We don't store a Plaintext of zeros
            this->u_const[i] = nullptr;
        }
    }

    this->u_const_init = std::vector<Plaintext>(m / s);
    for (int i = 0; i < m / s; i++) {
        if (u_const_clear_init[i] != zero_vect) {
            this->u_const_init[i] = cc->MakeCKKSPackedPlaintext(u_const_clear_init[i]);
#ifdef EXPE
            number_encoding++;
#endif
        } else {
            // We don't store a Plaintext of zeros
            this->u_const_init[i] = nullptr;
        }
    }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::vector<std::vector<double>> u_dangling_const_clear(m / s, zero_vect);
    // Same for the u_deg0's
    for (auto& it : map_deg0) {
        triplet value = it.second;
        if (map_uj.count(representative_orbits.at(value.vertex)) == 0 &&
            partition.at(representative_orbits.at(value.vertex)) == datacenter_id) {
            // value.vertex is a constant node
            // The distant node are not in map_deg0 and so don't fill the condition
            assert(value.indice < uj_size);
            // dangling nodes are compute only at last iteration: we only need 1 u_dangling_const
            if (this->max_iteration > 1) {
                u_dangling_const_clear[value.indice / s][value.indice % s] +=
                    0.15 / (double)(this->max_deg_in * value.deg_out);
            } else {
                u_dangling_const_clear[value.indice / s][value.indice % s] +=
                    1.0 / (double)(number_of_nodes * value.deg_out);
            }
        }
    }
#ifdef DEBUG
    std::cout << "Datacenter " << datacenter_id << " u_dangling_const : " << u_dangling_const_clear << std::endl;
#endif
    // Encode the constant vector for the dangling nodes
    this->u_dangling_const = std::vector<Plaintext>(m / s);
    for (int i = 0; i < m / s; i++) {
        if (u_dangling_const_clear[i] != zero_vect) {
            this->u_dangling_const[i] = cc->MakeCKKSPackedPlaintext(u_dangling_const_clear[i]);
#ifdef EXPE
            number_encoding++;
#endif
        } else {
            // We don't store a Plaintext of zeros
            this->u_dangling_const[i] = nullptr;
        }
    }
#endif
}

/**
 * @brief Compute the PRs of the uj's
 * We sum all uj's, we then add the constant PRs and the u distant PRs
 * We multiply by 0.85 and add 0.15 (damping factor)
 *
 * @param u_distant
 */
void POPPYujMultiDC::compute_PR() {
    // Mask and rotate the PRs of the distant nodes and add them to the vector

    // Initialize PRs with u_distant
    for (int i = 0; i < m / s; i++) {
        this->PRs[i].reset();
        this->PRs[i] = u[i]->Clone();
    }
    delete_vector_of_ciphers(u);

    for (int i = 0; i < m / s; i++) {
        if (this->current_iteration > 0) {
            // We add vectors of constant nodes
            if (u_const[i] != nullptr) {
                cc->EvalAddInPlace(PRs[i], u_const[i]);
#ifdef EXPE
                number_add++;
#endif
            }
        } else {
            // We add vectors of constant nodes
            if (u_const_init[i] != nullptr) {
                cc->EvalAddInPlace(PRs[i], u_const_init[i]);
#ifdef EXPE
                number_add++;
#endif
            }
        }

        cc->EvalMultInPlace(PRs[i], 0.85);
        cc->EvalAddInPlace(PRs[i], 0.15 / this->max_deg_in);
#ifdef EXPE
        number_mult_scalar++;
        number_add_scalar++;
#endif
    }
    if (current_iteration == 0) {
        delete_vector_of_plains(u_const_init);
    } else {
        delete_vector_of_plains(u_const);
    }

    this->current_iteration++;
}

/**
 * @brief Compute the PRs of the dangling nodes
 *
 * @param u_distant
 */
void POPPYujMultiDC::compute_PR_dangling_node() {
    // Initialize PRs_deg0
    this->PRs_deg0 = std::vector<Ciphertext<DCRTPoly>>(m / s);
    for (int i = 0; i < m / s; i++) {
        this->PRs_deg0[i] = u_dangling[i]->Clone();
    }
    delete_vector_of_ciphers(u_dangling);

    for (int i = 0; i < m / s; i++) {
        if (u_dangling_const[i] != nullptr) {
            // We add vectors of constant nodes
            cc->EvalAddInPlace(PRs_deg0[i], u_dangling_const[i]);
#ifdef EXPE
            number_add++;
#endif
            delete_vector_of_plains(u_dangling_const);
        }

        cc->EvalMultInPlace(PRs_deg0[i], 0.85);
        cc->EvalAddInPlace(PRs_deg0[i], 0.15 / this->max_deg_in);
#ifdef EXPE
        number_mult_scalar++;
        number_add_scalar++;
#endif
    }
}

/**
 * @brief Sum together the PRs to add to the uj's
 *
 * @param PRs_distant_nodes
 * @param map_uj Used to find the index of the uj where the distant value is needed
 */
void POPPYujMultiDC::synchronize_distant_PRs(
    std::vector<NodesInCKKSVector> PRs_distant_nodes,
    std::multimap<int, triplet> map_uj,
    bool active_nodes) {
    if (!active_nodes) {
        u_dangling = construction_uj(
            this->PRs, m, s, map_construction_u_deg0, keys, cc, representative_orbits, partition, datacenter_id);
    }
    // Add to the uj's the masked vector
    // Before adding we need to rotate the vector to the right index
    for (auto PR_distant_node : PRs_distant_nodes) {
        auto precomp = cc->EvalFastRotationPrecompute(PR_distant_node.vect);
        auto node_indexes =
            std::vector<std::tuple<Node /* node in vect*/, unsigned int /* indexe in vect */, unsigned int>>();

        // Loop over the map_uj to find the index of the uj where the distant value is needed
        for (auto node : PR_distant_node.indexes) {
            for (auto elt : map_uj) {
                if (elt.second.vertex == orbits.at(node.first)) {
                    node_indexes.emplace_back(node.first, node.second, elt.second.indice);
                }
            }
        }
        // Map of rotated vects to re use a rotation already done
        std::map<int, std::unique_ptr<Ciphertext<DCRTPoly>>> rotated_vects;

        for (auto node : node_indexes) {
            Node node_id;             // node id in vector
            unsigned int vect_index;  // indexe in vector
            unsigned int index;       // indexe where we need to place the node
            std::tie(node_id, vect_index, index) = node;
            // Calculation of rotation step
            int rot = ((vect_index % s) + s - (index % s)) % s;
#ifdef DEBUG
            std::cout << "Rot value : " << rot << std::endl;
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
                        rotated_vects.insert({rot, std::make_unique<Ciphertext<DCRTPoly>>(rotated_PR)});
#ifdef EXPE
                        number_rot++;
#endif
                    } else {
                        int rotated = 0;
                        // Find the largest power of 2 that is smaller than rot
                        int k = 1;
                        while (!(rot & k)) {
                            k = k << 1;
                        }
                        rotated = k;
                        // Check if it has already been computed
                        auto it_rotated = rotated_vects.find(k);
                        if (it_rotated == rotated_vects.end()) {
                            rotated_PR = this->cc->EvalFastRotation(
                                PR_distant_node.vect, k, 2 * cc->GetRingDimension(), precomp);
                            rotated_vects.insert({k, std::make_unique<Ciphertext<DCRTPoly>>(rotated_PR)});
#ifdef EXPE
                            number_rot++;
#endif
                        } else {
                            rotated_PR = *(it_rotated->second);
                        }
                        // Compute the remaining rotation
                        while (rotated != rot) {
                            // Find the largest power of 2 that is smaller than remaining rot
                            int k = 1;
                            while (!((rot - rotated) & k)) {
                                k = k << 1;
                            }
                            rotated += k;
                            // Rotate the PR
                            rotated_PR = this->cc->EvalRotate(rotated_PR, k);
                            rotated_vects.insert({rotated, std::make_unique<Ciphertext<DCRTPoly>>(rotated_PR)});
#ifdef EXPE
                            number_rot++;
#endif
                        }
                    }
#else
                    rotated_PR =
                        this->cc->EvalFastRotation(PR_distant_node.vect, rot, 2 * cc->GetRingDimension(), precomp);
                    rotated_vects.insert({rot, std::make_unique<Ciphertext<DCRTPoly>>(rotated_PR)});
#ifdef EXPE
                    number_rot++;
#endif
#endif
                } else {
                    rotated_PR = *(it_rot->second);
                }
            }
            // Mask result to obtain only one PRs in the vector
            std::vector<double> vect(this->s, 0.0);
            vect[index % s] = 1;
            Plaintext mask = this->cc->MakeCKKSPackedPlaintext(vect);
#ifdef EXPE
            number_encoding++;
#endif
            // Add rotated masked PRs to the current PRs
            auto mult = cc->EvalMult(std::move(rotated_PR), mask);
            if (active_nodes)
                cc->EvalAddInPlace(u[index / s], mult);
            else
                cc->EvalAddInPlace(u_dangling[index / s], mult);
#ifdef EXPE
            number_mult_plain++;
            number_add++;
#endif
        }
    }
}

/**
 * @brief Reconstruct the u_j, called at the end of each iteration
 *
 */
void POPPYujMultiDC::reconstruction_u_j() {
    std::vector<double> zero(uj_size);
    Plaintext pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encoding++;
    number_encryption++;
#endif

    // To save memory, we fill each sub uj one by one.
    for (int i = 0; i < m / s; i++) {
        this->u[i] = ct_zero->Clone();  // Reset uj to a vector of zeros

        auto map_reconstruction = map_reconstruction_uj[i];
        std::shared_ptr<std::vector<DCRTPoly>> precomp = cc->EvalFastRotationPrecompute(PRs[i]);

        // Map of rotated vects to re use a rotation already done
        std::map<int, Ciphertext<DCRTPoly>> rotated_vects;
        Ciphertext<DCRTPoly> rotated_PR = PRs[i];

        auto it = map_reconstruction.begin();
        auto it2 = map_reconstruction.begin();

        while (it != map_reconstruction.end()) {
            // for each rotation : save the rotated PRs
            int rot = it->first;
            if (rot != 0) {
                auto it_rot = rotated_vects.find(rot);
                if (it_rot == rotated_vects.end()) {
#ifdef COEUS
                    // Check if rot is a power of 2
                    if ((rot & (rot - 1)) == 0) {
                        rotated_PR = this->cc->EvalFastRotation(PRs[i], rot, 2 * cc->GetRingDimension(), precomp);
                        rotated_vects.insert({rot, Ciphertext<DCRTPoly>(rotated_PR)});
#ifdef EXPE
                        number_rot++;
#endif
                    } else {
                        int rotated = 0;
                        // Find the largest power of 2 that is smaller than rot
                        int k = 1;
                        while (!(rot & k)) {
                            k = k << 1;
                        }
                        rotated = k;
                        // Check if it has already been computed
                        auto it_rotated = rotated_vects.find(k);
                        if (it_rotated == rotated_vects.end()) {
                            rotated_PR = this->cc->EvalFastRotation(PRs[i], k, 2 * cc->GetRingDimension(), precomp);
                            rotated_vects.insert({k, Ciphertext<DCRTPoly>(rotated_PR)});
#ifdef EXPE
                            number_rot++;
#endif
                        } else {
                            rotated_PR = rotated_vects.at(k);
                        }

                        // Compute the remaining rotation
                        while (rotated != rot) {
                            // Find the largest power of 2 that is smaller than remaining rot
                            int k = 1;
                            while (!((rot - rotated) & k)) {
                                k = k << 1;
                            }
                            rotated += k;
                            // Rotate the PR
                            rotated_PR = this->cc->EvalRotate(rotated_PR, k);
                            rotated_vects.insert({rotated, Ciphertext<DCRTPoly>(rotated_PR)});
#ifdef EXPE
                            number_rot++;
#endif
                        }
                    }
#else
                    rotated_PR = this->cc->EvalFastRotation(PRs[i], rot, 2 * cc->GetRingDimension(), precomp);
                    rotated_vects.insert({rot, Ciphertext<DCRTPoly>(rotated_PR)});
#ifdef EXPE
                    number_rot++;
#endif
#endif
                } else {
                    rotated_PR = rotated_vects.at(rot);
                }
            }

            // Now we mask the rotated cipher
            // We creat a unique mask per sub_uj
            std::vector<double> mask_of_zeros(this->s, 0.0);
            std::vector<std::vector<double>> masks(m / s, mask_of_zeros);
            while (it2 != map_reconstruction.end() && it2->first == it->first) {
                auto quadruplet = it2->second;

                assert(quadruplet.deg_out != 0);
                masks[quadruplet.indice / s][quadruplet.indice % s] += 1.0 / quadruplet.deg_out;
                it2++;
            }
            for (int sub_uj = 0; sub_uj < m / s; sub_uj++) {
                Plaintext pl_mask = cc->MakeCKKSPackedPlaintext(masks[sub_uj]);
                auto PR_inter = cc->EvalMult(rotated_PR, pl_mask);
                cc->EvalAddInPlace(this->u[sub_uj], PR_inter);

#ifdef EXPE
                number_encoding++;
                number_add++;
                number_mult_plain++;
#endif
            }
            it = it2;
        }
    }
}

/**
 * @brief Run the POPPYujMultiDC algorithm
 *
 * @param nb_iteration
 * @param datacenters
 * @param partition
 * @param multDepth
 */
void run(
    int nb_iteration,
    std::vector<std::shared_ptr<POPPYujMultiDC>> datacenters,
    std::map<int, int> partition,
    int multDepth) {
    // Boolean for active or dangling nodes synchronization
    bool active_nodes = true;
    std::vector<NodesInCKKSVector> PRs_distant_nodes;
    for (int iter = 1; iter < nb_iteration; iter++) {
        std::cout << "iteration " << iter << std::endl;
#ifdef EXPE
        expe_poppy_uj.setCurrentRun(iter);
#endif
        for (auto datacenter : datacenters) {
#ifdef EXPE
            int tmp_add_scalar, tmp_add, tmp_rot, tmp_mult, tmp_mult_scalar, tmp_encryption, tmp_encoding,
                tmp_mult_plain;
            std::chrono::time_point<std::chrono::high_resolution_clock> startSync;
            startSync = std::chrono::high_resolution_clock::now();
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_rot = number_rot;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;
#endif
            // Synchro
            PRs_distant_nodes = synch_data(datacenter, datacenters, partition, iter - 1);
            datacenter->synchronize_distant_PRs(PRs_distant_nodes, datacenter->get_map_uj(), active_nodes);

#ifdef EXPE
            auto endSync = std::chrono::high_resolution_clock::now();
            auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
            expe_poppy_uj.addValue(
                durationSync.count(), "Time to synchronize distant PRs", Experiment::SYNC, Experiment::SECONDS);
            expe_poppy_uj.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition during synchronization",
                Experiment::Operation::ADD_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_add - tmp_add,
                "Number of addition during synchronization",
                Experiment::Operation::ADD_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during synchronization",
                Experiment::Operation::ENCODE_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during synchronization",
                Experiment::Operation::ENCRYPT_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during synchronization",
                Experiment::Operation::MUL_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_rot - tmp_rot,
                "Number of rotation during synchronization",
                Experiment::Operation::ROT_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during synchronization",
                Experiment::Operation::MUL_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during synchronization",
                Experiment::Operation::MUL_PLAIN_SYNC,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
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
            expe_poppy_uj.addValue(
                durationCompute.count(), "Time to compute one iteration", Experiment::EVAL_OP, Experiment::SECONDS);
            expe_poppy_uj.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition during computation" + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_add - tmp_add,
                "Number of addition during computation" + std::to_string(iter),
                Experiment::Operation::ADD_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during computation" + std::to_string(iter),
                Experiment::Operation::ENCODE_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during computation" + std::to_string(iter),
                Experiment::Operation::ENCRYPT_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during computation" + std::to_string(iter),
                Experiment::Operation::MUL_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_rot - tmp_rot,
                "Number of rotation during computation" + std::to_string(iter),
                Experiment::Operation::ROT_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during computation" + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during computation" + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_COMPUTE,
                Experiment::Unit::UNIT);
#endif

#ifdef BOOTSTRAPPING
#ifdef EXPE
            auto startBootstrap = std::chrono::high_resolution_clock::now();
#endif
            if (iter >= multDepth / POPPY_uj_LEVELS_CONSUMPTION) {
                datacenter->bootstrap(multDepth, false);
            }
#ifdef EXPE
            auto endBootstrap = std::chrono::high_resolution_clock::now();
            auto durationBootstrap =
                std::chrono::duration_cast<std::chrono::milliseconds>(endBootstrap - startBootstrap);
            expe_poppy_uj.addValue(
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
            auto startPrep = std::chrono::high_resolution_clock::now();
#endif
            // Next iteration prep
            datacenter->copy_previous_masked_distant_node();
            datacenter->compute_masked_distant_node(iter);
            datacenter->add_constant_to_distant_nodes(iter);
            datacenter->reconstruction_u_j();
#ifdef EXPE
            auto endPrep = std::chrono::high_resolution_clock::now();
            auto durationPrep = std::chrono::duration_cast<std::chrono::milliseconds>(endPrep - startPrep);
            expe_poppy_uj.addValue(
                durationPrep.count(), "Time to prepare next iteration", Experiment::PREP, Experiment::MILLISECONDS);
            expe_poppy_uj.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition during preparation" + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_add - tmp_add,
                "Number of addition during preparation" + std::to_string(iter),
                Experiment::Operation::ADD_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during preparation" + std::to_string(iter),
                Experiment::Operation::ENCODE_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during preparation" + std::to_string(iter),
                Experiment::Operation::ENCRYPT_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during preparation" + std::to_string(iter),
                Experiment::Operation::MUL_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_rot - tmp_rot,
                "Number of rotation during preparation" + std::to_string(iter),
                Experiment::Operation::ROT_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during preparation" + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_poppy_uj.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during preparation" + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_PREP,
                Experiment::Unit::UNIT);
#endif
        }
#ifdef EXPE
        expe_poppy_uj.addValue(
            number_add_scalar,
            "Number of scalar addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD_SCALAR,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_add,
            "Number of addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encoding,
            "Number of encoding after iteration " + std::to_string(iter),
            Experiment::Operation::ENCODE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encryption,
            "Number of encryption after iteration " + std::to_string(iter),
            Experiment::Operation::ENCRYPT,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult,
            "Number of multiplication after iteration " + std::to_string(iter),
            Experiment::Operation::MUL,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_rot,
            "Number of rotation after iteration " + std::to_string(iter),
            Experiment::Operation::ROT,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_scalar,
            "Number of multiplication scalar after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_SCALAR,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_plain,
            "Number of multiplication plain after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_PLAIN,
            Experiment::Unit::UNIT);
#endif
    }
    std::cout << "Last iteration " << nb_iteration << std::endl;
#ifdef EXPE
    expe_poppy_uj.setCurrentRun(nb_iteration);
#endif
    for (auto datacenter : datacenters) {
#ifdef EXPE
        auto startSync = std::chrono::high_resolution_clock::now();
        int tmp_add_scalar, tmp_add, tmp_rot, tmp_mult, tmp_mult_scalar, tmp_encryption, tmp_encoding, tmp_mult_plain;
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_rot = number_rot;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
        tmp_mult_plain = number_mult_plain;
#endif

        // Synchro
        PRs_distant_nodes = synch_data(datacenter, datacenters, partition, nb_iteration - 1);
        datacenter->synchronize_distant_PRs(PRs_distant_nodes, datacenter->get_map_uj(), active_nodes);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        datacenter->synchronize_distant_PRs(PRs_distant_nodes, datacenter->get_map_deg0(), !active_nodes);
#endif

#ifdef EXPE
        auto endSync = std::chrono::high_resolution_clock::now();
        auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
        expe_poppy_uj.addValue(
            durationSync.count(),
            "Time to synchronize distant PRs at last iteration",
            Experiment::SYNC,
            Experiment::SECONDS);
        expe_poppy_uj.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition during synchronization at last iteration",
            Experiment::Operation::ADD_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_add - tmp_add,
            "Number of addition during synchronization at last iteration",
            Experiment::Operation::ADD_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding during synchronization at last iteration",
            Experiment::Operation::ENCODE_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption during synchronization at last iteration",
            Experiment::Operation::ENCRYPT_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult - tmp_mult,
            "Number of multiplication during synchronization at last iteration",
            Experiment::Operation::MUL_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_rot - tmp_rot,
            "Number of rotation during synchronization at last iteration",
            Experiment::Operation::ROT_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar during synchronization at last iteration",
            Experiment::Operation::MUL_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain during synchronization at last iteration",
            Experiment::Operation::MUL_PLAIN_SYNC,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
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
        auto startCompute = std::chrono::high_resolution_clock::now();
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_rot = number_rot;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
#endif
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        if (datacenter->get_nb_dangling_nodes() > 0) {
            datacenter->compute_PR_dangling_node();
        }
#endif
        datacenter->compute_PR();
        datacenter->finished = true;
#ifdef EXPE
        auto endCompute = std::chrono::high_resolution_clock::now();
        auto durationCompute = std::chrono::duration_cast<std::chrono::seconds>(endCompute - startCompute);
        expe_poppy_uj.addValue(
            durationCompute.count(),
            "Time to compute one iteration at last iteration",
            Experiment::EVAL_OP,
            Experiment::SECONDS);
        expe_poppy_uj.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition during computation at last iteration",
            Experiment::Operation::ADD_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_add - tmp_add,
            "Number of addition during computation at last iteration",
            Experiment::Operation::ADD_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding during computation at last iteration",
            Experiment::Operation::ENCODE_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption during computation at last iteration",
            Experiment::Operation::ENCRYPT_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult - tmp_mult,
            "Number of multiplication during computation at last iteration",
            Experiment::Operation::MUL_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_rot - tmp_rot,
            "Number of rotation during computation at last iteration",
            Experiment::Operation::ROT_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar during computation at last iteration",
            Experiment::Operation::MUL_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain during computation at last iteration",
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
        auto startPrep = std::chrono::high_resolution_clock::now();
#endif
        datacenter->copy_previous_masked_distant_node();
#ifdef EXPE
        auto endPrep = std::chrono::high_resolution_clock::now();
        auto durationPrep = std::chrono::duration_cast<std::chrono::milliseconds>(endPrep - startPrep);
        expe_poppy_uj.addValue(
            durationPrep.count(),
            "Time to prepare next iteration at last iteration",
            Experiment::PREP,
            Experiment::MILLISECONDS);
        expe_poppy_uj.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition during prep at last iteration",
            Experiment::Operation::ADD_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_add - tmp_add,
            "Number of addition during prep at last iteration",
            Experiment::Operation::ADD_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding during prep at last iteration",
            Experiment::Operation::ENCODE_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption during prep at last iteration",
            Experiment::Operation::ENCRYPT_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult - tmp_mult,
            "Number of multiplication during prep at last iteration",
            Experiment::Operation::MUL_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_rot - tmp_rot,
            "Number of rotation during prep at last iteration",
            Experiment::Operation::ROT_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar during prep at last iteration",
            Experiment::Operation::MUL_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_poppy_uj.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain during prep at last iteration",
            Experiment::Operation::MUL_PLAIN_PREP,
            Experiment::Unit::UNIT);
#endif
    }
}

void run_multi_dc_poppy_uj(int nb_iteration) {
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

    CryptoContext<DCRTPoly> ccTest;
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    int multDepth = POPPY_uj_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {3, 3};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = POPPY_uj_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif

    ccTest = gen_crypto_context_test(multDepth);

    std::cout << "start multi DC uj initialization" << std::endl;

    // Initialize the datacenters
    int datacenter_num = 0;
    std::vector<std::shared_ptr<POPPYujMultiDC>> datacenters;
    std::vector<int> uj_sizes;
    std::vector<int> max_deg_in;
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<POPPYujMultiDC>(datacenter_num));
        uj_sizes.push_back(datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci));
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }

    int max_uj_size = *std::max_element(uj_sizes.begin(), uj_sizes.end());
    uint32_t nbSlots = nextPowerOfTwo(max_uj_size);

    usint ringDim = ccTest->GetRingDimension();
    int s;
    if (nbSlots < ringDim) {
        s = nbSlots;
        std::cout << "s : " << s << " (the next power of two of " << max_uj_size << ", the number of active nodes)"
                  << std::endl;

    } else {
        s = ringDim / 2;
        std::cout << "s : " << s << " (ringDim / 2)" << std::endl;
    }

#ifdef DEBUG
    std::cout << "The common batch size is " << s << std::endl;
    std::cout << "The multiplicative depth is " << multDepth << std::endl;
#endif
    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);
#ifdef COEUS
        gen_rotation_keys_coeus(cc, keys, s);
#else
        gen_rotation_keys(cc, keys, s);
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
    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_number_datacenters(datacenter_num);
        datacenter->set_max_iteration(nb_iteration);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }

    // Run
    std::cout << "start multi DC uj computation of " << nb_iteration << " iterations" << std::endl;
    run(nb_iteration, datacenters, partition, multDepth);

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

std::map<int, double> run_multi_dc_poppy_uj(int nb_iteration, std::string filename) {
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

    CryptoContext<DCRTPoly> ccTest;
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    int multDepth = POPPY_uj_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {3, 3};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = POPPY_uj_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif
    ccTest = gen_crypto_context_test(multDepth);

    std::cout << "start multi DC uj initialization" << std::endl;

    // Initialize the datacenters
    int datacenter_num = 0;
#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    std::vector<std::shared_ptr<POPPYujMultiDC>> datacenters;
    std::vector<int> uj_sizes;
    std::vector<int> max_deg_in;
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<POPPYujMultiDC>(datacenter_num));
        uj_sizes.push_back(datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci));
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }



#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    expe_poppy_uj.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time to init all datacenters",
        Experiment::Operation::INIT_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    std::cout << "start gen_crypto_context" << std::endl;
    int max_uj_size = *std::max_element(uj_sizes.begin(), uj_sizes.end());
    uint32_t nbSlots = nextPowerOfTwo(max_uj_size);

    usint ringDim = ccTest->GetRingDimension();
    int s;
    if (nbSlots < ringDim) {
        s = nbSlots;
        std::cout << "s : " << s << " (the next power of two of " << max_uj_size << ", the number of active nodes)"
                  << std::endl;
    } else {
        s = ringDim / 2;
        std::cout << "s : " << s << " (ringDim / 2)" << std::endl;
    }

    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);
#ifdef COEUS
        gen_rotation_keys_coeus(cc, keys, s);
#else
        gen_rotation_keys(cc, keys, s);
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
    expe_poppy_uj.setRingsize(cc->GetRingDimension());
    expe_poppy_uj.setCurrentRun(0);
    start = std::chrono::high_resolution_clock::now();
#endif
    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_number_datacenters(datacenter_num);
        datacenter->set_max_iteration(nb_iteration);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_uj.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time for Preprocessing step all datacenters",
        Experiment::Operation::PREPROC_DC,
        Experiment::Unit::MILLISECONDS);
    expe_poppy_uj.addValue(
        number_add_scalar,
        "Number of scalar addition after the preprocessing",
        Experiment::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_add, "Number of addition after the preprocessing", Experiment::ADD, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult, "Number of multiplication after the preprocessing", Experiment::MUL, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_rot, "Number of rotation after the preprocessing", Experiment::ROT, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult_scalar,
        "Number of scalar multiplication after the preprocessing",
        Experiment::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult_plain,
        "Number of multiplication plain after the preprocessing",
        Experiment::MUL_PLAIN,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_encoding, "Number of encoding after the preprocessing", Experiment::ENCODE, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_encryption, "Number of encryption after the preprocessing", Experiment::ENCRYPT, Experiment::Unit::UNIT);
#endif
    // Run
    std::cout << "start multi DC uj computation of " << nb_iteration << " iterations" << std::endl;
#ifdef EXPE
    start = std::chrono::high_resolution_clock::now();
#endif
    run(nb_iteration, datacenters, partition, multDepth);

#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_poppy_uj.addValue(
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count(),
        "Time to compute " + std::to_string(nb_iteration) + " iterations",
        Experiment::Operation::EVAL_DC,
        Experiment::Unit::SECONDS);
    expe_poppy_uj.addValue(
        number_add_scalar,
        "Number of scalar addition after last iteration",
        Experiment::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_add, "Number of addition after last iteration", Experiment::Operation::ADD, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult, "Number of multiplication after last iteration", Experiment::MUL, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_rot, "Number of rotation after last iteration", Experiment::Operation::ROT, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult_scalar,
        "Number of multiplication scalar after last iteration",
        Experiment::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_mult_plain,
        "Number of multiplication plain after last iteration",
        Experiment::MUL_PLAIN,
        Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_encoding, "Number of encoding after last iteration", Experiment::ENCODE, Experiment::Unit::UNIT);
    expe_poppy_uj.addValue(
        number_encryption, "Number of encryption after last iteration", Experiment::ENCRYPT, Experiment::Unit::UNIT);
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
    expe_poppy_uj.disable_crypto_context();
#endif
    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalSumKeys();

    return results;
}
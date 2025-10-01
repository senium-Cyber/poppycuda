#include <thread>
#include "naive.cuh"
#include "openfhe.h"

#ifdef EXPE
#include "utils/include/utils_export.cuh"
Experiment expe_naive;
#endif

/**
 * @brief Initialize the local graph of the datacenter
 * In this step we compute the map containing for each nodes the in and out neighbors
 * We initialize the data structure of active nodes
 * @param subgraph The local graph of the datacenter
 * @param number_of_nodes The number of nodes in the global graph
 * @param partition The local partition of the datacenter
 */
void NaiveMultiDC::init_local_graph(
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
        orbits,
        this->datacenter_id,
        duration,
        set_dangling_nodes);
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
#else
    std::tie(map_in, map_out) = process_maps(&subgraph, &partition, nb_local_nodes, orbits, this->datacenter_id);
    nb_local_nodes_after_pruning = nb_local_nodes;
#endif
    this->nb_dangling_nodes = set_dangling_nodes.size();
    this->dangling_nodes = std::vector<int>(set_dangling_nodes.begin(), set_dangling_nodes.end());

#ifdef EXPE
    expe_naive.addValue(
        duration,
        "Time to compute equivalent nodes on DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::MICROSECONDS);
    expe_naive.addValue(
        nb_local_nodes,
        "Number of nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_naive.addValue(
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

    create_active_and_cst_nodes();

#ifdef EXPE
    expe_naive.addValue(
        this->cst_nodes.size(),
        "Number of constant nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
    expe_naive.addValue(
        this->nb_dangling_nodes,
        "Number of dangling nodes in DC " + std::to_string(this->datacenter_id),
        Experiment::PRUNING,
        Experiment::UNIT);
#endif

    // Set representative orbits
    for (auto orbit : orbits) {
        if (representative_orbits.find(orbit.second) == representative_orbits.end())
            this->representative_orbits.insert({orbit.second, orbit.first});
    }

    std::set<int> connected_datacenters;
    for (auto node : map_out) {
        // If it is a distant node
        if (partition.at(node.first) != datacenter_id) {
            connected_datacenters.insert(partition.at(node.first));
        }
    }
    this->connected_datacenters = connected_datacenters;

#ifdef DEBUG
    printNodeCategories();
#endif
}

/**
 * @brief Creat the vector of active nodes and the vector of dangling nodes
 *
 */
void NaiveMultiDC::create_active_and_cst_nodes() {
    std::set<Node> active_nodes_set;
    std::set<Node> cst_nodes_set;

    for (auto node : orbits) {
        if (partition[node.first] != datacenter_id)
            continue;

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        if (map_in.count(node.first) != 0) {  // If the node is not constant
            if (map_out.count(node.first) != 0)
#endif
                active_nodes_set.insert(node.second);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        } else
            cst_nodes_set.insert(node.second);
#endif
    }

    // We transform sets into vectors
    this->active_nodes = std::vector<int>(active_nodes_set.begin(), active_nodes_set.end());
    this->cst_nodes = std::vector<int>(cst_nodes_set.begin(), cst_nodes_set.end());

    // Clear the sets of active and dangling nodes
    active_nodes_set.clear();
    cst_nodes_set.clear();
}

/**
 * @brief Initialize the PageRank, all the nodes have the same PR equal to 1/n
 * where n is the number of nodes in the global graph
 */
void NaiveMultiDC::initialize_PRs() {
    int n = this->number_of_nodes;
    auto cc = this->cc;
    auto public_key = this->keys.publicKey;
    std::vector<double> x = {1.0 / n};

    // We create the initialisation
    Plaintext iv = cc->MakeCKKSPackedPlaintext(x);
#ifdef EXPE
    number_encoding++;
#endif

    // We fill the map with the initialisation ciphertext
    for (auto node : active_nodes) {
        auto ct_iv = cc->Encrypt(public_key, iv);
        this->PRs.insert(std::pair<int, Ciphertext<DCRTPoly>>(node, ct_iv));
#ifdef EXPE
        number_encryption++;
#endif
        ct_iv.reset();  // Clear ciphertext
    }
    iv.reset();  // Clear the Plaintext
}

void NaiveMultiDC::printNodeCategories() {
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

    // 3) Print constant nodes
    std::cout << "Constant nodes (" << nb_cst_nodes << "): ";
    for (auto node : cst_nodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "---------------------------------------------------\n\n";
}

/**
 * @brief Preprocessing of the datacenter
 * Initialize the PRs ciphertexts
 *
 */
void NaiveMultiDC::pre_processing() {
    initialize_PRs();

    // Mask the initialized PRs to send the value to distant data centers
    compute_masked_distant_node();
    copy_previous_masked_distant_node();
}

/**
 * @brief Set the common crypto context and keys
 *
 * @param cc Common crypto context
 * @param keys
 * @param vector_size Common batch size
 */
void NaiveMultiDC::set_crypto_context(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, unsigned int vector_size) {
    this->cc = cc;
    this->keys = keys;
    this->s = vector_size;
}

/**
 * @brief Compute the mask multiplication for all the border nodes of the data center
 * Creation of the map border_PRs needed for the synchronization
 * Contains the real node id, not the orbit
 *
 */
void NaiveMultiDC::compute_masked_distant_node() {
    std::multimap<int /*datacenter*/, std::pair<Node, Ciphertext<DCRTPoly>>> border_PRs_current_it;
    for (auto node : map_out) {
        // if it is a border node
        if (partition[node.second] != datacenter_id) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
            auto it = find(active_nodes.begin(), active_nodes.end(), orbits.at(node.first));
            if (it != active_nodes.end()) {
#endif
                auto search = PRs.find(orbits.at(node.first));

                std::pair<Node, Ciphertext<DCRTPoly>> border_node_pr = {
                    node.first, cc->EvalMult(search->second, 1.0 / map_out.count(node.first))};
#ifdef EXPE
                number_mult_scalar++;
#endif
                border_PRs_current_it.insert(
                    std::pair<int, std::pair<Node, Ciphertext<DCRTPoly>>>(partition[node.second], border_node_pr));
                border_node_pr.second.reset();

#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
            } else {  // else node i is a constant node
                std::vector<double> cst_vect;

                if (this->current_iteration > 0)
                    cst_vect = {0.15 / (this->max_deg_in * map_out.count(node.first))};
                else {
                    cst_vect = {1.0 / (this->number_of_nodes * map_out.count(node.first))};
                }

                std::pair<Node, Ciphertext<DCRTPoly>> border_node_pr = {
                    node.first, cc->Encrypt(this->keys.publicKey, cc->MakeCKKSPackedPlaintext(cst_vect))};
                border_PRs_current_it.insert(
                    std::pair<int, std::pair<Node, Ciphertext<DCRTPoly>>>(partition[node.second], border_node_pr));
                border_node_pr.second.reset();

#ifdef EXPE
                number_encoding++;
                number_encryption++;
#endif
            }
#endif
        }
    }

    for (auto pr : this->border_PRs)
        pr.second.second.reset();
    this->border_PRs.clear();
    this->border_PRs = border_PRs_current_it;

    // Clear border_PRs_current_it
    for (auto pr : border_PRs_current_it)
        pr.second.second.reset();
    border_PRs_current_it.clear();
}

/**
 * @brief Extract de nodes with the corresponding ciphertexts from border_PRs
 *
 * @param datacenter_id
 * @return std::vector<std::pair<Node, Ciphertext<DCRTPoly>>>
 */
std::vector<std::pair<Node, Ciphertext<DCRTPoly>>> NaiveMultiDC::get_PR_distant_nodes(
    int datacenter_id,
    int iteration) {
    std::vector<std::pair<Node, Ciphertext<DCRTPoly>>> distant_nodes;
    if (this->current_iteration > iteration) {
        auto nodes = previous_border_PRs.equal_range(datacenter_id);
        for (auto it = nodes.first; it != nodes.second; ++it) {
            distant_nodes.push_back(it->second);
        }
    } else {
        auto nodes = border_PRs.equal_range(datacenter_id);
        for (auto it = nodes.first; it != nodes.second; ++it) {
            distant_nodes.push_back(it->second);
        }
    }
    return distant_nodes;
}

/**
 * @brief Synchronize the distant PRs : add them to the map distant_PRs)
 *
 */
void NaiveMultiDC::synchronize_distant_PRs(std::vector<std::shared_ptr<NaiveMultiDC>> datacenters, int iteration) {
    for (auto pr : this->distant_PRs)
        pr.second.reset();
    this->distant_PRs.clear();
    for (int dc : this->connected_datacenters) {
        // For each dc wait for the PRs value of the distant nodes
        auto nodes_from_dc = datacenters[dc]->get_PR_distant_nodes(datacenter_id, iteration);
        for (auto nodes_and_pr : nodes_from_dc) {
            this->distant_PRs.insert(nodes_and_pr);
        }
    }
}

/**
 * @brief Compute the PageRank of the datacenter for the current iteration
 * The computation is done by scalar multiplication, encrypted additions
 * Following by the damping factor
 *
 */
void NaiveMultiDC::compute_PR() {
    auto cc = this->cc;
    auto keys = this->keys;

    // We create plaintext of 0
    std::vector<double> x_0 = {0.0};
    Plaintext zero = cc->MakeCKKSPackedPlaintext(x_0);
#ifdef EXPE
    number_encoding++;
#endif

    std::map<int, Ciphertext<DCRTPoly>> PR_end_iteration;
    for (auto i : active_nodes) {
        auto nodes_in = map_in.equal_range(representative_orbits[i]);  // To get all the nodes getting in i
        auto sum = cc->Encrypt(keys.publicKey, zero);
#ifdef EXPE
        number_encryption++;
#endif

        for (auto it = nodes_in.first; it != nodes_in.second; ++it) {  // For each nodes_in
            int j = it->second;

            // If node j is one of our nodes
            if (partition[j] == datacenter_id) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                auto it_cst = std::find(active_nodes.begin(), active_nodes.end(), orbits.at(j));
                if (it_cst == active_nodes.end()) {  // if j is a constant node
                    if (this->current_iteration > 0)
                        sum = cc->EvalAdd(sum, 0.15 / (this->max_deg_in * map_out.count(j)));
                    else {
                        sum = cc->EvalAdd(sum, 1.0 / (this->number_of_nodes * map_out.count(j)));
                    }
#ifdef EXPE
                    number_add_scalar++;
#endif
                } else {
#endif

                    if (map_out.count(j) > 1) {
                        auto search = PRs.find(orbits.at(j));
                        auto val = cc->EvalMult(search->second, 1.0 / map_out.count(j));  // PR(j)/deg_out(j)
                        sum = cc->EvalAdd(sum, val);
                        val.reset();
#ifdef EXPE
                        number_add++;
                        number_mult_scalar++;
#endif
                    } else {
                        auto search = PRs.find(orbits.at(j));
                        sum = cc->EvalAdd(sum, search->second);
#ifdef EXPE
                        number_add++;
#endif
                    }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                }
#endif
            }
            // Else we find the pagerank in distant_PRs using the real node id (not the orbit)
            else {
                auto search = distant_PRs.find(j);
                sum = cc->EvalAdd(sum, search->second);
#ifdef EXPE
                number_add++;
#endif
            }
        }

        cc->EvalMultInPlace(sum, 0.85);
        cc->EvalAddInPlace(sum, 0.15 / this->max_deg_in);
#ifdef EXPE
        number_add_scalar++;
        number_mult_scalar++;
#endif

        PR_end_iteration.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, sum));
        sum.reset();
    }

    for (auto pr : PRs)
        pr.second.reset();
    this->PRs.clear();

    this->PRs = PR_end_iteration;
    this->current_iteration++;

    // Clear map of Ciphertexts and the plaintext
    for (auto pr : PR_end_iteration)
        pr.second.reset();
    PR_end_iteration.clear();
    zero.reset();
}

/**
 * @brief Compute the PageRank of the dangling nodes for the last iteration
 * The computation is done by scalar multiplication, encrypted additions
 * Following by the damping factor
 *
 */
void NaiveMultiDC::compute_PR_dangling_node() {
    auto cc = this->cc;
    auto keys = this->keys;

    // We create plaintext of 0
    std::vector<double> x_0 = {0.0};
    Plaintext zero = cc->MakeCKKSPackedPlaintext(x_0);
#ifdef EXPE
    number_encoding++;
#endif

    std::map<int, Ciphertext<DCRTPoly>> PR_end_iteration;
    for (auto i : dangling_nodes) {
        auto nodes_in = map_in.equal_range(representative_orbits[i]);  // To get all the nodes getting in i
        auto sum = cc->Encrypt(keys.publicKey, zero);
#ifdef EXPE
        number_encryption++;
#endif

        for (auto it = nodes_in.first; it != nodes_in.second; ++it) {  // For each nodes_in
            int j = it->second;

            // If node j is one of our nodes
            if (partition[j] == datacenter_id) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                auto it_cst = std::find(active_nodes.begin(), active_nodes.end(), orbits.at(j));
                if (it_cst == active_nodes.end()) {  // if j is a constant node
                    if (this->current_iteration > 0)
                        sum = cc->EvalAdd(sum, 0.15 / (this->max_deg_in * map_out.count(j)));
                    else
                        sum = cc->EvalAdd(sum, 1.0 / (this->number_of_nodes * map_out.count(j)));
#ifdef EXPE
                    number_add_scalar++;
#endif
                } else {
#endif

                    if (map_out.count(j) > 1) {
                        auto search = PRs.find(orbits.at(j));
                        auto val = cc->EvalMult(search->second, 1.0 / map_out.count(j));  // PR(j)/deg_out(j)
                        sum = cc->EvalAdd(sum, val);
                        val.reset();
#ifdef EXPE
                        number_add++;
                        number_mult_scalar++;
#endif
                    } else {
                        auto search = PRs.find(orbits.at(j));
                        sum = cc->EvalAdd(sum, search->second);
#ifdef EXPE
                        number_add++;
#endif
                    }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                }
#endif
            }
            // Else we find the pagerank in distant_PRs using the real node id (not the orbit)
            else {
                auto search = distant_PRs.find(j);
                sum = cc->EvalAdd(sum, search->second);
#ifdef EXPE
                number_add++;
#endif
            }
        }

        cc->EvalMultInPlace(sum, 0.85);
        cc->EvalAddInPlace(sum, 0.15 / this->max_deg_in);
#ifdef EXPE
        number_add_scalar++;
        number_mult_scalar++;
#endif

        PR_end_iteration.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, sum));
        sum.reset();
    }

    this->PRs_deg0 = PR_end_iteration;

    // Clear map of Ciphertexts and the plaintext
    for (auto pr : PR_end_iteration)
        pr.second.reset();
    PR_end_iteration.clear();
    zero.reset();
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
    std::vector<std::shared_ptr<NaiveMultiDC>> datacenters,
    std::map<int, int> partition,
    int multDepth) {
    for (int iter = 1; iter < nb_iteration; iter++) {
        std::cout << "iteration " << iter << std::endl;
#ifdef EXPE
        expe_naive.setCurrentRun(iter);
#endif
        for (auto& datacenter : datacenters) {
#ifdef EXPE
            int tmp_add_scalar, tmp_add, tmp_mult, tmp_mult_scalar, tmp_encoding, tmp_encryption, tmp_mult_plain;

            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;

            auto startSync = std::chrono::high_resolution_clock::now();
#endif
            // Synchro
            datacenter->synchronize_distant_PRs(datacenters, iter - 1);
#ifdef EXPE
            auto endSync = std::chrono::high_resolution_clock::now();
            auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
            expe_naive.addValue(
                durationSync.count(), "Time to synchronize distant PRs", Experiment::SYNC, Experiment::SECONDS);
            expe_naive.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition during synchronization " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_add - tmp_add,
                "Number of addition during synchronization " + std::to_string(iter),
                Experiment::Operation::ADD_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during synchronization " + std::to_string(iter),
                Experiment::Operation::ENCODE_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during synchronization " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_SYNC,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during synchronization " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_SYNC,
                Experiment::Unit::UNIT);

#endif
#ifdef EXPE
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
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
            expe_naive.addValue(
                durationCompute.count(), "Time to compute PRs", Experiment::EVAL_OP, Experiment::SECONDS);
            expe_naive.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of scalar addition after computation " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_add - tmp_add,
                "Number of addition after computation " + std::to_string(iter),
                Experiment::Operation::ADD_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding after computation " + std::to_string(iter),
                Experiment::Operation::ENCODE_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption after computation " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult - tmp_mult,
                "Number of multiplication after computation " + std::to_string(iter),
                Experiment::Operation::MUL_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of scalar multiplication after computation " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_COMPUTE,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain after computation " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_COMPUTE,
                Experiment::Unit::UNIT);

#endif

#ifdef BOOTSTRAPPING
#ifdef EXPE
            auto startBootstrap = std::chrono::high_resolution_clock::now();
#endif
            if (iter >= multDepth / naive_LEVELS_CONSUMPTION) {
                datacenter->bootstrap(multDepth);
            }
#ifdef EXPE
            auto endBootstrap = std::chrono::high_resolution_clock::now();
            auto durationBootstrap =
                std::chrono::duration_cast<std::chrono::milliseconds>(endBootstrap - startBootstrap);
            expe_naive.addValue(
                durationBootstrap.count(), "Time to bootstrap", Experiment::BOOTSTRAP, Experiment::MILLISECONDS);
#endif
#endif

#ifdef EXPE
            tmp_add_scalar = number_add_scalar;
            tmp_add = number_add;
            tmp_mult = number_mult;
            tmp_mult_scalar = number_mult_scalar;
            tmp_encoding = number_encoding;
            tmp_encryption = number_encryption;
            tmp_mult_plain = number_mult_plain;
            auto startPrepare = std::chrono::high_resolution_clock::now();
#endif
            // Prepare next iteration
            datacenter->copy_previous_masked_distant_node();
            datacenter->compute_masked_distant_node();
#ifdef EXPE
            auto endPrepare = std::chrono::high_resolution_clock::now();
            auto durationPrepare = std::chrono::duration_cast<std::chrono::seconds>(endPrepare - startPrepare);
            expe_naive.addValue(
                durationPrepare.count(), "Time to prepare next iteration", Experiment::PREP, Experiment::SECONDS);
            expe_naive.addValue(
                number_add_scalar - tmp_add_scalar,
                "Number of addition during prepare " + std::to_string(iter),
                Experiment::Operation::ADD_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_add - tmp_add,
                "Number of addition during prepare " + std::to_string(iter),
                Experiment::Operation::ADD_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encoding - tmp_encoding,
                "Number of encoding during prepare " + std::to_string(iter),
                Experiment::Operation::ENCODE_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_encryption - tmp_encryption,
                "Number of encryption during prepare " + std::to_string(iter),
                Experiment::Operation::ENCRYPT_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult - tmp_mult,
                "Number of multiplication during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_scalar - tmp_mult_scalar,
                "Number of multiplication scalar during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_SCALAR_PREP,
                Experiment::Unit::UNIT);
            expe_naive.addValue(
                number_mult_plain - tmp_mult_plain,
                "Number of multiplication plain during prepare " + std::to_string(iter),
                Experiment::Operation::MUL_PLAIN_PREP,
                Experiment::Unit::UNIT);
#endif
        }
#ifdef EXPE
        expe_naive.addValue(
            number_add_scalar,
            "Number of scalar addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD_SCALAR,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_add,
            "Number of addition after iteration " + std::to_string(iter),
            Experiment::Operation::ADD,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encoding,
            "Number of encoding after iteration " + std::to_string(iter),
            Experiment::Operation::ENCODE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encryption,
            "Number of encryption after iteration " + std::to_string(iter),
            Experiment::Operation::ENCRYPT,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult,
            "Number of multiplication after iteration " + std::to_string(iter),
            Experiment::Operation::MUL,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_scalar,
            "Number of multiplication scalar after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_SCALAR,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_plain,
            "Number of multiplication plain after iteration " + std::to_string(iter),
            Experiment::Operation::MUL_PLAIN,
            Experiment::Unit::UNIT);
#endif
    }
    std::cout << "Last iteration " << nb_iteration << std::endl;
#ifdef EXPE
    expe_naive.setCurrentRun(nb_iteration);
#endif
    for (auto& datacenter : datacenters) {
#ifdef EXPE
        int tmp_add_scalar, tmp_add, tmp_mult, tmp_mult_scalar, tmp_encoding, tmp_encryption, tmp_mult_plain;
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
        tmp_mult = number_mult;
        tmp_mult_scalar = number_mult_scalar;
        tmp_encoding = number_encoding;
        tmp_encryption = number_encryption;
        tmp_mult_plain = number_mult_plain;
        auto startSync = std::chrono::high_resolution_clock::now();
#endif
        // Synchro
        datacenter->synchronize_distant_PRs(datacenters, nb_iteration - 1);

#ifdef EXPE
        auto endSync = std::chrono::high_resolution_clock::now();
        auto durationSync = std::chrono::duration_cast<std::chrono::seconds>(endSync - startSync);
        expe_naive.addValue(
            durationSync.count(),
            "Time to synchronize distant PRs at last iteration",
            Experiment::SYNC,
            Experiment::SECONDS);
        expe_naive.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition during synchronization at last iteration",
            Experiment::Operation::ADD_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_add - tmp_add,
            "Number of addition during synchronization at last iteration",
            Experiment::Operation::ADD_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding during synchronization at last iteration",
            Experiment::Operation::ENCODE_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption during synchronization at last iteration",
            Experiment::Operation::ENCRYPT_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult - tmp_mult,
            "Number of multiplication during synchronization at last iteration",
            Experiment::Operation::MUL_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar during synchronization at last iteration",
            Experiment::Operation::MUL_SCALAR_SYNC,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain during synchronization at last iteration",
            Experiment::Operation::MUL_PLAIN_SYNC,
            Experiment::Unit::UNIT);
#endif

#ifdef EXPE
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
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
        expe_naive.addValue(
            durationCompute.count(), "Time to compute PRs at last iteration", Experiment::EVAL_OP, Experiment::SECONDS);
        expe_naive.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition after last iteration",
            Experiment::Operation::ADD_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_add - tmp_add,
            "Number of addition after last iteration",
            Experiment::Operation::ADD_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding after last iteration",
            Experiment::Operation::ENCODE_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption after last iteration",
            Experiment::Operation::ENCRYPT_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult - tmp_mult,
            "Number of multiplication after last iteration",
            Experiment::Operation::MUL_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar after last iteration",
            Experiment::Operation::MUL_SCALAR_COMPUTE,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain after last iteration",
            Experiment::Operation::MUL_PLAIN_COMPUTE,
            Experiment::Unit::UNIT);
#endif
#ifdef EXPE
        tmp_add_scalar = number_add_scalar;
        tmp_add = number_add;
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
        auto durationPrepare = std::chrono::duration_cast<std::chrono::seconds>(endPrepare - startPrepare);
        expe_naive.addValue(
            durationPrepare.count(), "Time to prepare at last iteration", Experiment::PREP, Experiment::SECONDS);
        expe_naive.addValue(
            number_add_scalar - tmp_add_scalar,
            "Number of scalar addition prepare at last iteration",
            Experiment::Operation::ADD_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_add - tmp_add,
            "Number of addition prepare at last iteration",
            Experiment::Operation::ADD_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encoding - tmp_encoding,
            "Number of encoding prepare at last iteration",
            Experiment::Operation::ENCODE_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_encryption - tmp_encryption,
            "Number of encryption prepare at last iteration",
            Experiment::Operation::ENCRYPT_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult - tmp_mult,
            "Number of multiplication prepare at last iteration",
            Experiment::Operation::MUL_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_scalar - tmp_mult_scalar,
            "Number of multiplication scalar prepare at last iteration",
            Experiment::Operation::MUL_SCALAR_PREP,
            Experiment::Unit::UNIT);
        expe_naive.addValue(
            number_mult_plain - tmp_mult_plain,
            "Number of multiplication plain prepare at last iteration",
            Experiment::Operation::MUL_PLAIN_PREP,
            Experiment::Unit::UNIT);
#endif
    }
}

/**
 * @brief Decrypt and print the local PageRank of the datacenter
 *
 */
void NaiveMultiDC::print_PRs() {
    // We store the decrypted values in the map to decrypt only one per orbit
    std::map<Node, double> decrypted_PRs;
    for (auto orbit : orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;

        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbit.second);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        auto it_dang = std::find(dangling_nodes.begin(), dangling_nodes.end(), orbit.second);
#endif
        if (it != active_nodes.end()) {
            if (decrypted_PRs.find(orbit.second) !=
                decrypted_PRs.end()) {  // if the orbit of the node is already decrypted
                std::cout << "find" << std::endl;
                std::cout << orbit.first << " is " << decrypted_PRs.at(orbit.second) << std::endl;
            } else {
                Plaintext res;
                std::cout.precision(8);
                this->cc->Decrypt(keys.secretKey, PRs.at(orbit.second), &res);
                res->SetLength(s);
                auto vect_res = res->GetRealPackedValue();
                std::cout << orbit.first << " is " << vect_res[0] << std::endl;

                decrypted_PRs.insert({orbit.second, vect_res[0]});
            }
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (it_dang != dangling_nodes.end()) {
            if (decrypted_PRs.find(orbit.second) !=
                decrypted_PRs.end()) {  // if the orbit of the node is already decrypted
                std::cout << orbit.first << " is " << decrypted_PRs.at(orbit.second) << std::endl;
            } else {
                Plaintext res;
                std::cout.precision(8);
                this->cc->Decrypt(keys.secretKey, PRs_deg0.at(orbit.second), &res);
                res->SetLength(s);
                auto vect_res = res->GetRealPackedValue();
                std::cout << orbit.first << " is " << vect_res[0] << " (dangling)" << std::endl;

                decrypted_PRs.insert({orbit.second, vect_res[0]});
            }
        } else
            std::cout << orbit.first << " is " << 0.15 / this->max_deg_in << " (constant)" << std::endl;
#endif
    }
}

std::map<int, double> NaiveMultiDC::export_results() {
    // We store the decrypted values in the map to decrypt only one per orbit
    std::map<Node, double> decrypted_PRs;
    // We store each result in the second map
    std::map<int, double> results;

    for (auto orbit : orbits) {
        if (partition[orbit.first] != datacenter_id)
            continue;

        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbit.second);
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        auto it_dang = std::find(dangling_nodes.begin(), dangling_nodes.end(), orbit.second);
#endif
        if (it != active_nodes.end()) {
            if (decrypted_PRs.find(orbit.second) !=
                decrypted_PRs.end()) {  // if the orbit of the node is already decrypted
                results.insert({orbit.first, decrypted_PRs.at(orbit.second)});
            } else {
                Plaintext res;
                std::cout.precision(8);
                this->cc->Decrypt(keys.secretKey, PRs.at(orbit.second), &res);
                res->SetLength(s);
                auto vect_res = res->GetRealPackedValue();
                results.insert({orbit.first, vect_res[0]});

                decrypted_PRs.insert({orbit.second, vect_res[0]});
            }
        }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        else if (it_dang != dangling_nodes.end()) {
            if (decrypted_PRs.find(orbit.second) !=
                decrypted_PRs.end()) {  // if the orbit of the node is already decrypted
                results.insert({orbit.first, decrypted_PRs.at(orbit.second)});
            } else {
                Plaintext res;
                std::cout.precision(8);
                this->cc->Decrypt(keys.secretKey, PRs_deg0.at(orbit.second), &res);
                res->SetLength(s);
                auto vect_res = res->GetRealPackedValue();
                results.insert({orbit.first, vect_res[0]});

                decrypted_PRs.insert({orbit.second, vect_res[0]});
            }
        } else
            results.insert({orbit.first, 0.15 / this->max_deg_in});
#endif
    }
    return results;
}

/**
 * @brief Main function to compute the PRs of the datacenters in the Naive multi DC scenario
 * First we parse graph of the datacenters to obtain multiple graphs
 * We let the datacenters initialize their local graph
 * We generate a common crypto context for all the datacenters
 * We run the PR algorithm
 */
void run_multi_dc_naive(int nb_iteration) {
#ifdef EXPE
    number_add_scalar = 0;
    number_add = 0;
    number_mult = 0;
    number_mult_scalar = 0;
    number_mult_plain = 0;
    number_encoding = 0;
    number_encryption = 0;
#endif
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC("../dot/students-partition1.dot", partition);

    std::cout << "start multi DC naive initialization" << std::endl;

    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    std::cout << "start gen_crypto_context" << std::endl;
    int multDepth = naive_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {1, 1};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = naive_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif

    // Initialize the datacenters
    std::vector<std::shared_ptr<NaiveMultiDC>> datacenters;
    int datacenter_num = 0;
#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    std::vector<int> max_deg_in;
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<NaiveMultiDC>(datacenter_num));
        datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci);
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }

#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time to init all datacenters",
        Experiment::Operation::INIT_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    int s = 1;  // size vectors
    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);

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
    expe_naive.setRingsize(cc->GetRingDimension());
    expe_naive.setCurrentRun(0);
    start = std::chrono::high_resolution_clock::now();
#endif
    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time for Preprocessing step all datacenters",
        Experiment::Operation::PREPROC_DC,
        Experiment::Unit::MILLISECONDS);
    expe_naive.addValue(
        number_add_scalar,
        "Number of scalar additions after the preprocessing",
        Experiment::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_add, "Number of additions after the preprocessing", Experiment::ADD, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult, "Number of multiplications after the preprocessing", Experiment::MUL, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult_scalar,
        "Number of scalar multiplications after the preprocessing",
        Experiment::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encoding, "Number of encoding after the preprocessing", Experiment::ENCODE, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encryption,
        "Number of encryptions after the preprocessing",
        Experiment::ENCRYPT,
        Experiment::Unit::UNIT);
#endif

    // COMPUTATION STEP
    std::cout << "start multi DC naive computation" << std::endl;
#ifdef EXPE
    start = std::chrono::high_resolution_clock::now();
#endif
    run(nb_iteration, datacenters, partition, multDepth);
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count(),
        "Time to compute " + std::to_string(nb_iteration) + " iterations",
        Experiment::Operation::EVAL_DC,
        Experiment::Unit::SECONDS);
#endif

    // COMPLETION STEP
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

#ifdef EXPE
    expe_naive.disable_crypto_context();
#endif
    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalSumKeys();
}

std::map<int, double> run_multi_dc_naive(int nb_iteration, string filename) {
#ifdef EXPE
    number_add_scalar = 0;
    number_add = 0;
    number_mult = 0;
    number_mult_scalar = 0;
    number_mult_plain = 0;
    number_encoding = 0;
    number_encryption = 0;
#endif
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC(filename, partition);

    std::cout << "start gen_crypto_context" << std::endl;
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    int multDepth = naive_MULT_DEPTH(nb_iteration);
#ifdef BOOTSTRAPPING
    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {1, 1};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = naive_LEVELS_CONSUMPTION + 1;  // Need a last level for decryption
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, UNIFORM_TERNARY);

    multDepth = levelsNeeded + levelsAvailableAfterBootstrap;
#endif

    std::cout << "start multi DC naive initialization" << std::endl;

    // Initialize the datacenters
    std::vector<std::shared_ptr<NaiveMultiDC>> datacenters;
    int datacenter_num = 0;
#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    std::vector<int> max_deg_in;
    for (Graph& graph : graphs) {
        int max_deg_in_dci = 0;
        datacenters.emplace_back(std::make_shared<NaiveMultiDC>(datacenter_num));
        datacenters.back()->init_local_graph(
            graph, partition.size(), compute_sub_partition(partition, graph.arcs, datacenter_num), max_deg_in_dci);
        datacenter_num++;
        max_deg_in.push_back(max_deg_in_dci);
    }

#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time to init all datacenters",
        Experiment::Operation::INIT_DC,
        Experiment::Unit::MILLISECONDS);
#endif

    int s = 1;  // size vectors
    if (!deserialize_cc(&cc, &keys, get_folder_by_params(multDepth, s))) {
        std::tie(cc, keys) = gen_crypto_context(multDepth, s);

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
    expe_naive.setRingsize(cc->GetRingDimension());
    expe_naive.setCurrentRun(0);
    start = std::chrono::high_resolution_clock::now();
#endif
    int max = *std::max_element(max_deg_in.begin(), max_deg_in.end());
    std::cout << "Max deg_in: " << max << std::endl;

    for (auto datacenter : datacenters) {
        datacenter->set_max_dg_in(max);
        datacenter->set_crypto_context(cc, keys, s);
        datacenter->pre_processing();
    }
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
        "Time for Preprocessing step all datacenters",
        Experiment::Operation::PREPROC_DC,
        Experiment::Unit::MILLISECONDS);
    expe_naive.addValue(
        number_add_scalar,
        "Number of scalar additions after the preprocessing",
        Experiment::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_add, "Number of additions after the preprocessing", Experiment::ADD, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult, "Number of multiplications after the preprocessing", Experiment::MUL, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult_scalar,
        "Number of scalar multiplications after the preprocessing",
        Experiment::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encoding, "Number of encoding after the preprocessing", Experiment::ENCODE, Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encryption,
        "Number of encryptions after the preprocessing",
        Experiment::ENCRYPT,
        Experiment::Unit::UNIT);
#endif

    // COMPUTATION STEP
    std::cout << "start multi DC naive computation" << std::endl;
#ifdef EXPE
    start = std::chrono::high_resolution_clock::now();
#endif
    run(nb_iteration, datacenters, partition, multDepth);
#ifdef EXPE
    end = std::chrono::high_resolution_clock::now();
    expe_naive.addValue(
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count(),
        "Time to compute " + std::to_string(nb_iteration) + " iterations",
        Experiment::Operation::EVAL_DC,
        Experiment::Unit::SECONDS);
    expe_naive.addValue(
        number_add_scalar,
        "Number of scalar addition after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::ADD_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_add,
        "Number of addition after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::ADD,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encoding,
        "Number of encoding after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::ENCODE,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_encryption,
        "Number of encryption after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::ENCRYPT,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult,
        "Number of multiplication after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::MUL,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult_scalar,
        "Number of multiplication scalar after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::MUL_SCALAR,
        Experiment::Unit::UNIT);
    expe_naive.addValue(
        number_mult_plain,
        "Number of multiplication plain after last iteration " + std::to_string(nb_iteration),
        Experiment::Operation::MUL_PLAIN,
        Experiment::Unit::UNIT);
#endif

    // COMPLETION STEP
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
    expe_naive.disable_crypto_context();
#endif
    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalMultKeys();
    cc->ClearEvalSumKeys();

    return results;
}
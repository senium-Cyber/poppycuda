#include "poppy_uj.cuh"

/**
 * @brief Pre-processing of the graph
 * Reads the graph from a csv file (init)
 * Load crypto context and keys
 * Initialize u_j
 * @param csv_name name of the csv file
 * @param option 'g' for generation, 'w' for generate and write, 'r' for read
 * @param folder folder where the keys are stored
 * @param separator separator of the csv file
 * @param begin_by_zero if the nodes begin by 0
 * @param directed if the graph is directed
 */
void POPPYuj::pre_processing(std::string csv_name, std::string separator, bool directed) {
    init(csv_name, separator, directed);

    int multDepth = 8;
    gen_crypto_context(multDepth);
    gen_rotation_keys();

    initialize_uj();
}

/**
 * @brief Read the graph from a csv file and compute the orbits with nauty
 * At the end create the maps for the u_j
 *
 * @param csv_name
 * @param separator
 * @param begin_by_zero
 * @param directed
 */
void POPPYuj::init(std::string csv_name, std::string separator, bool directed) {
    // Read the csv file
    Graph nodes_and_arcs = read_graph(csv_name, separator, directed);

    // time here
    auto start = std::chrono::high_resolution_clock::now();
    // Creation of the graph G
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, local_orbits, orbits_sz);
    static DEFAULTOPTIONS_SPARSEDIGRAPH(options);
    statsblk stats;
    SG_DECL(sg);  // Declare and initialize sparse graph structure

    // We set the options of the graph.
    options.digraph = TRUE;
    options.writemarkers = TRUE;  // Write levels on terminal
    options.writeautoms = TRUE;
    options.defaultptn = TRUE;  // No coloration

    SG_INIT(sg);

    int n = nodes_and_arcs.nodes.size();
    int nb_arcs = nodes_and_arcs.arcs.size();
    int m = SETWORDSNEEDED(n);

    std::map<int, int> trad_tab;
    int index = 0;
    for (auto node : nodes_and_arcs.nodes) {
        trad_tab[node] = index;
        index++;
    }
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, local_orbits, orbits_sz, n, "malloc");

    // Initialisation of the graph
    SG_ALLOC(sg, n, nb_arcs, "malloc");

    sg.nv = n;         // Number of vertices
    sg.nde = nb_arcs;  // Number of arcs

    // Creation of map_in and map_out of G
    int edge_index = 0;
    for (auto edge : nodes_and_arcs.arcs) {
        int i = trad_tab.at(edge.first);
        int j = trad_tab.at(edge.second);
        map_out.insert(std::pair<int, int>(i, j));
        map_in.insert(std::pair<int, int>(j, i));
        assert(edge_index < nb_arcs);
        sg.e[edge_index++] = j;  // Add each neighbor consecutively
    }

    // Out degree of each node of sg
    int i = 0;
    for (auto node : nodes_and_arcs.nodes) {
        auto deg = map_out.count(node);
        sg.v[trad_tab.at(node)] = i;    // table of starting indices for neighbors of each vertex in sg.e
        sg.d[trad_tab.at(node)] = deg;  // table of outdegrees
        i += deg;
    }
    assert(sg.dlen == n);
    assert(sg.vlen == n);

    std::cout << "Computing the orbits " << std::endl;
    sparsenauty(&sg, lab, ptn, local_orbits, &options, &stats, NULL);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken to compute the orbits: " << duration.count() << " seconds" << std::endl;

    this->number_of_nodes = n;

    this->uj_size = set_u_j_size(n, local_orbits, map_in, map_out);

    // Now we creat the map_u_j.
    this->nb_of_uj = creation_map_u_j(map_uj, n, local_orbits, map_in, map_out);

    this->nb_u_deg0 = creation_map_deg_out_0(map_deg0, n, local_orbits, map_in, map_out);

    this->orbits = local_orbits;

    // We delete g because we don't need him anymore
    SG_FREE(sg);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
}

/**
 * @brief Generate the crypto context
 *
 * @param multDepth
 */
void POPPYuj::gen_crypto_context(uint32_t multDepth) {
    uint32_t nbSlots = nextPowerOfTwo(uj_size);

    std::cout << "s : " << nbSlots << "(the next power of two of " << this->uj_size << ", the number of active nodes)"
              << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetBatchSize(nbSlots);

    parameters.SetSecurityLevel(HEStd_128_classic);

// Scaling parameters
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits = 78;
    usint firstMod = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

    // Keys generation
    keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
}

/**
 * @brief Generate the rotation keys
 *
 */
void POPPYuj::gen_rotation_keys() {
    std::vector<int> steps;
    steps.push_back(0);
    for (int i = 1; i < nextPowerOfTwo(this->uj_size); i++) {
        steps.push_back(i);
        steps.push_back(-i);
    }

    this->cc->EvalRotateKeyGen(keys.secretKey, steps);
}

/**
 * @brief Compute the PRs for each iteration
 * @param nb_iteration
 */
void POPPYuj::run(int nb_iteration) {
    creation_vect_PR_const();

    for (int i = 1; i < nb_iteration; i++) {
        std::cout << "iteration " << i << std::endl;
        compute_PR();
        reconstruction_u_j();
    }

    // Last iteration
    std::cout << "iteration " << nb_iteration << std::endl;
    if (this->nb_u_deg0 > 0) {
        // PR computation for dangling nodes
        std::cout << "compute_PR_dangling_nodes" << std::endl;
        compute_PR_dangling_node(nb_iteration);
        compute_PR();
    } else {
        compute_PR();
    }
}

/**
 * @brief Compute the PRs for standard nodes
 *
 */
void POPPYuj::compute_PR() {
    // this->print_data_center();
    PRs = u[0];
    for (int j = 1; j < nb_of_uj; j++) {
        cc->EvalAddInPlace(PRs, u[j]);
    }

    // End of PR algorithm
    cc->EvalMultInPlace(PRs, 0.85);
    cc->EvalAddInPlace(PRs, 0.15);
}

Ciphertext<DCRTPoly> POPPYuj::get_PRs() {
    return this->PRs;
}

std::map<int, int> POPPYuj::get_orbits() {
    std::map<int, int> orbits_map;
    for (int i = 0; i < this->number_of_nodes; i++) {
        orbits_map[i] = this->orbits[i];
    }
    return orbits_map;
}

/**
 * @brief Decrypt and print the PRs
 *
 */
void POPPYuj::print_PRs() {
    auto keys = get_keys();

    Plaintext res;
    std::cout.precision(8);
    cc->Decrypt(keys.secretKey, PRs, &res);
    res->SetLength(uj_size);
    auto vect_PRs = res->GetRealPackedValue();

    std::vector<double> vect_PRs_deg0;

    if (nb_u_deg0 > 0) {
        Plaintext res_deg0;
        std::cout.precision(8);
        cc->Decrypt(keys.secretKey, PRs_deg0, &res_deg0);
        res_deg0->SetLength(uj_size);
        vect_PRs_deg0 = res_deg0->GetRealPackedValue();
    }

    for (int i = 0; i < number_of_nodes; i++) {
        if (map_uj.count(orbits[i]) != 0) {
            int indice = map_uj.find(orbits[i])->second.indice;
            std::cout << i << " is " << vect_PRs[indice] << std::endl;
        } else if (map_deg0.count(orbits[i]) != 0) {
            int indice = map_deg0.find(orbits[i])->second.indice;
            std::cout << i << " is " << vect_PRs_deg0[indice] << std::endl;
        } else
            std::cout << i << " is 0.15" << std::endl;
    }
}

/**
 * @brief Generate the crypto context and serialize it
 *
 * @param multDepth
 * @param folder
 */
void POPPYuj::gen_crypto_context_and_serialize(uint32_t multDepth, std::string folder) {
    gen_crypto_context(multDepth);

    gen_rotation_keys();

    serialize_cc(cc, keys, multDepth);
}

/**
 * @brief Serialize the PRs
 *
 * @param PRs
 * @param PRs_deg0
 * @param folder
 */
void POPPYuj::serialize_ciphers_PRs(std::string folder) {
    serialize_cipher(PRs, folder, "PRs");
    serialize_cipher(PRs_deg0, folder, "PRs_deg0");
}

/**
 * @brief Deserialize the crypto context and keys
 *
 * @param folder
 */
void POPPYuj::read_crypto_context_and_keys(std::string folder) {
    deserialize_cc(&cc, &keys, folder);
}

/**
 * @brief Print the data center
 *
 */
void POPPYuj::print_data_center() {
    std::cout << std::setw(50) << std::setfill('*') << '\n' << std::endl;
    std::cout << "number_of_nodes : " << number_of_nodes << std::endl;
    std::cout << "map_in size : " << map_in.size() << std::endl;
    std::cout << "map_out size : " << map_out.size() << std::endl;
    std::cout << "uj_size : " << uj_size << std::endl;
    std::cout << "nb_of_uj : " << nb_of_uj << std::endl;
    std::cout << "u size : " << u.size() << std::endl;
    std::cout << "nb_u_deg0 : " << nb_u_deg0 << std::endl;
    std::cout << "map_uj size : " << map_uj.size() << std::endl;
    std::cout << std::setw(50) << std::setfill('*') << '\n' << std::endl;
}

/**
 * @brief Decrypt and print a ciphertext
 *
 * @param ct_name
 * @param ct
 */
void POPPYuj::decrypt_and_print(std::string ct_name, Ciphertext<DCRTPoly> ct) {
    Plaintext res;
    std::cout.precision(8);
    cc->Decrypt(keys.secretKey, ct, &res);
    res->SetLength(uj_size);
    auto vect = res->GetRealPackedValue();

    std::cout << ct_name << "={";
    for (double elt : vect) {
        std::cout << elt << ", ";
    }
    std::cout << "}" << std::endl;
}

std::map<int, double> POPPYuj::export_results() {
    std::map<int, double> results;

    Plaintext res;
    cc->Decrypt(keys.secretKey, PRs, &res);
    res->SetLength(uj_size);
    auto vect_PRs = res->GetRealPackedValue();

    std::vector<double> vect_PRs_deg0;

    if (nb_u_deg0 > 0) {
        Plaintext res_deg0;
        cc->Decrypt(keys.secretKey, PRs_deg0, &res_deg0);
        res_deg0->SetLength(uj_size);
        vect_PRs_deg0 = res_deg0->GetRealPackedValue();
    }

    for (int i = 0; i < number_of_nodes; i++) {
        if (map_uj.count(orbits[i]) != 0) {
            int indice = map_uj.find(orbits[i])->second.indice;
            results[i] = vect_PRs[indice];
        } else if (map_deg0.count(orbits[i]) != 0) {
            int indice = map_deg0.find(orbits[i])->second.indice;
            results[i] = vect_PRs_deg0[indice];
        } else
            results[i] = 0.15;
    }
    return results;
}

/**
 * @brief Initialize the u_j and encrypt them
 *
 */
void POPPYuj::initialize_uj() {
    std::vector<double> vect(uj_size);
    for (int i = 0; i < uj_size; i++)
        vect[i] = 0.0;

    std::vector<std::vector<double>> iv(nb_of_uj);
    for (size_t i = 0; i < iv.size(); i++)
        iv[i] = vect;

    int index_u_j = 0;  // Index to know each u_j to use

    auto it = map_uj.begin();
    auto it2 = map_uj.begin();
    while (it != map_uj.end()) {
        while (it2 != map_uj.end() && it2->first == it->first) {
            triplet value = it2->second;
            iv[index_u_j][value.indice] = 1.0 / (number_of_nodes * value.deg_out);

            index_u_j++;
            it2++;
        }

        index_u_j = 0;
        it = it2;
    }
    // instanciate u with a vector of nb_of_uj ciphertexts
    this->u = std::vector<Ciphertext<DCRTPoly>>(nb_of_uj);
    for (int j = 0; j < nb_of_uj; j++) {
        Plaintext plain = cc->MakeCKKSPackedPlaintext(iv[j]);
        this->u[j] = cc->Encrypt(keys.publicKey, plain);
    }
}

/**
 * @brief Compute the PRs for the constant nodes, this is called only once at the beginning
 *
 */
void POPPYuj::creation_vect_PR_const() {
    std::vector<double> vect_zero(uj_size);
    for (int i = 0; i < uj_size; i++)
        vect_zero[i] = 0.0;

    for (auto& it : map_uj) {
        triplet value = it.second;
        if (map_uj.count(value.vertex) == 0) {  // value.vertex is a constant node
            std::vector<double> vect_tmp = vect_zero;
            vect_tmp[value.indice] = 0.15 / value.deg_out;

            Plaintext cst_vect = cc->MakeCKKSPackedPlaintext(vect_tmp);
            auto cst_cipher = cc->Encrypt(keys.publicKey, cst_vect);

            doublet d = {it.first, cst_cipher};
            this->PRs_cst.insert(std::pair<int, doublet>(value.vertex, d));
        }
    }

    for (auto& it : map_deg0) {
        triplet value = it.second;
        if (map_uj.count(value.vertex) == 0) {  // value.vertex is a constant node
            std::vector<double> vect_tmp = vect_zero;
            vect_tmp[value.indice] = 0.15 / value.deg_out;

            Plaintext cst_vect = cc->MakeCKKSPackedPlaintext(vect_tmp);
            auto cst_cipher = cc->Encrypt(keys.publicKey, cst_vect);

            doublet d = {it.first, cst_cipher};
            this->PRs_cst.insert(std::pair<int, doublet>(value.vertex, d));
        }
    }
}

/**
 * @brief Reconstruct the u_j, called at the end of each iteration
 *
 */
void POPPYuj::reconstruction_u_j() {
    // First we reset the u_j to zero ciphertext
    std::vector<double> zero(uj_size);
    for (int i = 0; i < uj_size; i++) {
        zero[i] = 0.0;
    }

    Plaintext pl_zero = cc->MakeCKKSPackedPlaintext(zero);

    for (int j = 0; j < nb_of_uj; j++)
        this->u[j] = cc->Encrypt(keys.publicKey, pl_zero);

    // Secondly we fill them with PRs values
    int index_in_uj = 0;  // Index to know each u_j to use

    auto it = map_uj.begin();
    auto it2 = map_uj.begin();
    while (it != map_uj.end()) {
        while (it2 != map_uj.end() && it2->first == it->first) {
            triplet value = it2->second;
            if (map_uj.count(value.vertex) == 0) {
                Ciphertext<DCRTPoly> vect_cst;
                auto range = PRs_cst.equal_range(value.vertex);
                for (auto it_range = range.first; it_range != range.second; ++it_range) {
                    if ((it_range->second).node == it->first)
                        vect_cst = (it_range->second).cipher;
                }

                cc->EvalAddInPlace(u[index_in_uj], vect_cst);
            } else {
                int i = map_uj.find(value.vertex)->second.indice;
                int step_rot = value.indice - i;
                Ciphertext<DCRTPoly> PR_inter;  // A ciphertext to put the modification without modifying PR

                // Mask creation
                std::vector<double> mask(uj_size);
                for (int indice = 0; indice < uj_size; indice++) {
                    if (indice == i)
                        mask[i] = 1.0 / value.deg_out;
                    else
                        mask[indice] = 0.0;
                }
                Plaintext pl_mask = cc->MakeCKKSPackedPlaintext(mask);

                if (step_rot == 0) {
                    PR_inter = cc->EvalMult(PRs, pl_mask);
                    cc->EvalAddInPlace(u[index_in_uj], PR_inter);
                } else {
                    PR_inter = cc->EvalMult(PRs, pl_mask);
                    PR_inter = cc->EvalRotate(
                        PR_inter,
                        -step_rot);  // - step_rot > 0 : rotation to the left or - step_rot < 0 : rotation to the right
                    cc->EvalAddInPlace(u[index_in_uj], PR_inter);
                }
            }

            index_in_uj++;
            it2++;
        }

        index_in_uj = 0;
        it = it2;
    }
}

/**
 * @brief Compute the PRs for the dangling nodes, called once at the last iteration
 *
 * @param nb_iteration Used to know if we are in the first iteration
 */
void POPPYuj::compute_PR_dangling_node(int nb_iteration) {
    Ciphertext<DCRTPoly> last_PRs;  // Storage of last iteration PR

    if (nb_iteration == 1) {
        std::vector<double> vect_init(uj_size);
        for (int i = 0; i < uj_size; i++)
            vect_init[i] = 1.0 / number_of_nodes;
        Plaintext plain_init = cc->MakeCKKSPackedPlaintext(vect_init);
        last_PRs = cc->Encrypt(keys.publicKey, plain_init);
    } else
        last_PRs = PRs;

    std::vector<Ciphertext<DCRTPoly>> u_dangling =
        construction_u_deg0(last_PRs, uj_size, nb_u_deg0, map_deg0, map_uj, PRs_cst, keys, cc);

    this->PRs_deg0 = u_dangling[0];
    for (int j = 1; j < nb_u_deg0; j++) {
        cc->EvalAddInPlace(PRs_deg0, u_dangling[j]);
    }

    // End of PR algorithm
    cc->EvalMultInPlace(PRs_deg0, 0.85);
    cc->EvalAddInPlace(PRs_deg0, 0.15);
}

/**
 * @brief Main function of poppy_uj algorithm
 *
 */
void run_poppy_uj(int nb_iteration) {
    POPPYuj data_center;

    string filename = "../Graphs/students.txt";
    data_center.pre_processing(GRAPH_DIR + filename, " ", true);  // for student graph

    data_center.run(nb_iteration);

    std::cout << "\nAfter " << nb_iteration << " iterations, the PR : " << std::endl;
    data_center.print_PRs();

    // data_center.print_PRs();
    std::map<int, double> results = data_center.export_results();
    // Write the results in a file
    std::ofstream file;
    file.open("../results/" + filename + "_uj-" + std::to_string(nb_iteration) + ".csv");
    file << "Node,PR" << std::endl;
    for (auto result : results) {
        file << result.first << "," << result.second << std::endl;
    }
}
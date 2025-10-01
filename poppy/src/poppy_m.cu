#include "poppy_m.cuh"
#include <fstream>
#include <iomanip>
#include <iostream>

/**
 * @brief From a csv file, init the graph and launch the pre processing phase
 * Compute the adjacency matrix (init phase) and derive submatrices from it
 * Then prepare the sub diag vectors
 * We end the phase by initializing the PRs
 *
 * @param csv_name File name
 * @param option g: generate crypto context, w: generate and write crypto context, r: read crypto context
 * @param folder folder where the crypto context is stored
 * @param separator csv separator
 * @param begin_by_zero node in the csv file begin by 0 or 1
 * @param directed is the graph directed
 */
void POPPYm::pre_processing(std::string csv_name, std::string separator, bool directed) {
    if (!this->is_initialized)
        init(csv_name, separator, directed);

    int multDepth = 8;
    gen_crypto_context(multDepth);
    gen_rotation_keys();

    int n = mask_matrix.size();
    // find the least common multiple of n (size of matrix) and s (size of
    // vectors)
    this->m = this->s;
    while (m < n) {
        this->m += this->s;
    }
    std::vector<std::vector<std::vector<double>>> subMatrices = prepareMatrix(mask_matrix, this->s, this->m);
#ifdef DEBUG
    std::cout << "Number of submatrix = " << subMatrices.size() << std::endl;
#endif
    std::vector<Diags> subDiags((m / s) * (m / s));
    for (int l = 0; l < m / s; l++) {
        for (int c = 0; c < m / s; c++) {
            subDiags[l * m / s + c] = MakeCKKSPackedDiagonals(cc, subMatrices[l * m / s + c]);
        }
    }

    this->subDiags = subDiags;

    if (nb_dangling_nodes > 0) {
        std::vector<std::vector<std::vector<double>>> subMatrices_dangling_nodes =
            prepareMatrix(mask_matrix_dangling_nodes, this->s, this->m);

        std::vector<Diags> subDiags_dangling_nodes((m / s) * (m / s));
        for (int l = 0; l < m / s; l++) {
            for (int c = 0; c < m / s; c++) {
                subDiags_dangling_nodes[l * m / s + c] =
                    MakeCKKSPackedDiagonals(cc, subMatrices_dangling_nodes[l * m / s + c]);
            }
        }

        this->subDiags_dangling_nodes = subDiags_dangling_nodes;
    }
#ifdef DEBUG
    std::cout << "m / s = " << m / s << std::endl;
#endif
    initialize_PRs();
}

/**
 * @brief We use nauty to compute the orbits of the graph
 * Also in this function we create the maps for the incoming and outgoing nodes of G
 *
 * @param nodes_and_arcs Graph to be processed
 */
void POPPYm::nauty_routine(Graph nodes_and_arcs) {
    // Creation of the graph G
    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, local_orbits, orbits_sz);

    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    // We set the options of the graph.
    options.writeautoms = TRUE;
    options.writemarkers = TRUE;
    options.defaultptn = TRUE;  // No coloration
    options.getcanon = FALSE;
    options.digraph = TRUE;
    options.maxinvarlevel = 999;

    int n = nodes_and_arcs.nodes.size();
    int m = SETWORDSNEEDED(n);

    std::map<int, int> trad_tab;
    int index = 0;
    for (auto node : nodes_and_arcs.nodes) {
        trad_tab[node] = index;
        index++;
    }
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    DYNALLOC2(graph, g, g_sz, m, n, "malloc");
    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, local_orbits, orbits_sz, n, "malloc");

    // Initialisation of the graph
    EMPTYGRAPH(g, m, n);

    // Créez les maps pour les sommets entrants et sortants de G
    for (auto edge : nodes_and_arcs.arcs) {
        map_out.insert(std::pair<int, int>(trad_tab[edge.first], trad_tab[edge.second]));
        map_in.insert(std::pair<int, int>(trad_tab[edge.second], trad_tab[edge.first]));
        ADDONEARC(g, trad_tab[edge.first], trad_tab[edge.second], m);
    }

    std::cout << "Computing the orbits" << std::endl;
    densenauty(g, lab, ptn, local_orbits, &options, &stats, m, n, nullptr);
    std::cout << "Finished to compute the orbits" << std::endl;
    this->orbits = local_orbits;
    this->number_of_nodes = n;

    // We delete g because we don't need him anymore
    DYNFREE(g, g_sz);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
}

/**
 * @brief From the map_in and map_out, we identify the active and dangling nodes
 *
 */
void POPPYm::sort_active_and_dangling_nodes() {
    std::set<int> active_nodes_set;  // noeuds nécéssaires pour le calcul de PR
    std::set<int> dangling_nodes;    // noeuds inutiles
    std::multimap<int /*dangling node*/, int /*active node*/>
        dangling_actives_nodes;  // noeuds actifs connectés à des noeuds inutiles
#ifdef DEBUG
    for (auto edge : map_out) {
        std::cout << "Edge : " << edge.first << " -> " << edge.second << std::endl;
    }
#endif
    int nb_dangling_actives_nodes = 0;
    // Boucle qui permet de remplir actives_nodes et actives_dangling_nodes
    for (int i = 0; i < this->number_of_nodes; i++) {
        if (map_out.count(i) > 0)  // Not a dangling node
            active_nodes_set.insert(this->orbits[i]);
        else  // A dangling node
        {
            dangling_nodes.insert(this->orbits[i]);
            // Add also the neigthbors of the dangling node
            auto values = map_in.equal_range(i);
            if (values.first != map_in.end())
                nb_dangling_actives_nodes++;
            for (auto it = values.first; it != values.second; ++it) {
                dangling_actives_nodes.insert(std::pair<int, int>(this->orbits[i], this->orbits[it->second]));
                nb_dangling_actives_nodes++;
            }
        }
    }
    this->nb_dangling_nodes = dangling_nodes.size();

#ifdef DEBUG
    // List all actives nodes
    for (auto active_node : active_nodes_set) {
        std::cout << "Active node : " << active_node << std::endl;
    }
    // List all dangling nodes
    for (auto dangling_node : dangling_nodes) {
        std::cout << "dangling node : " << dangling_node << std::endl;
    }
    // List all actives nodes connected to dangling nodes
    for (auto dangling_active_node : dangling_actives_nodes) {
        std::cout << "Active node connected to dangling node : " << dangling_active_node.second << " -> "
                  << dangling_active_node.first << std::endl;
    }
#endif

    std::cout << "Number of orbits of active nodes : " << active_nodes_set.size() << std::endl;
    std::cout << "Number of dangling nodes : " << this->nb_dangling_nodes << std::endl;
    std::cout << "Number of active nodes connected to dangling nodes : " << dangling_actives_nodes.size() << std::endl;

    // We transform sets into vectors
    this->active_nodes = std::vector<int>(active_nodes_set.begin(), active_nodes_set.end());

    // Set nb_active_nodes equal to the max between active and dangling size
    this->number_active_nodes = std::max(this->active_nodes.size(), (std::size_t)nb_dangling_actives_nodes);

    this->active_dangling_nodes = std::vector<int>(this->number_active_nodes);
    std::vector<bool> completed_cells(this->number_active_nodes, false);

    // We order the active nodes (connected to dangling one) in a manner that the nodes take the same place as in
    // actives_nodes And we place between them the dangling nodes
    for (auto actives_dangling_node : dangling_actives_nodes) {
        auto it = std::find(active_nodes.begin(), active_nodes.end(), actives_dangling_node.second);
        if (it != active_nodes.end()) {
            int index_in_actives = it - active_nodes.begin();
            active_dangling_nodes[index_in_actives] = active_nodes[index_in_actives];
            completed_cells[index_in_actives] = true;
            index_nodes_in_active_nodes.push_back(index_in_actives);
        }
    }

    // We place the dangling nodes in the remaining cells
    int index = 0;
    for (auto dangling_node : dangling_nodes) {
        if (index == this->number_active_nodes)
            break;
        for (int i = index; i < this->number_active_nodes; i++) {
            if (!completed_cells[i]) {
                active_dangling_nodes[i] = dangling_node;
                completed_cells[i] = true;
                index = i;
                break;
            }
        }
    }
#ifdef DEBUG
    std::cout << "Active dangling nodes : ";
    for (auto active_dangling_node : active_dangling_nodes) {
        std::cout << active_dangling_node << ", ";
    }
    std::cout << std::endl;
#endif
}

std::map<int, double> POPPYm::export_results() {
    std::map<int, double> results;

    std::vector<Plaintext> res(m / s);
    std::vector<double> vect_PRs;

    // Decrypt page rank of standard nodes
    for (int i = 0; i < m / s; i++) {
        cc->Decrypt(keys.secretKey, PRs[i], &res[i]);
        res[i]->SetLength(s);
        auto subVect = res[i]->GetRealPackedValue();
        vect_PRs.insert(vect_PRs.end(), subVect.begin(), subVect.end());
    }

    // Decrypt page rank of dangling nodes
    std::vector<double> vect_PRs_deg0;
    if (this->nb_dangling_nodes > 0) {
        std::vector<Plaintext> res_deg0(m / s);

        for (int i = 0; i < m / s; i++) {
            cc->Decrypt(keys.secretKey, PRs_deg0[i], &res_deg0[i]);
            res_deg0[i]->SetLength(s);
            auto subVect = res_deg0[i]->GetRealPackedValue();
            vect_PRs_deg0.insert(vect_PRs_deg0.end(), subVect.begin(), subVect.end());
        }
    }

    // Print the PRs in the order of the nodes
    for (int i = 0; i < number_of_nodes; i++) {
        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbits[i]);
        if (it != active_nodes.end()) {
            int indice = it - active_nodes.begin();
            results[i] = vect_PRs[indice];
        } else {
            auto it_dangling_node = std::find(active_dangling_nodes.begin(), active_dangling_nodes.end(), orbits[i]);
            int indice = it_dangling_node - active_dangling_nodes.begin();
            results[i] = vect_PRs_deg0[indice];
        }
    }
    return results;
}

/**
 * @brief Initialize the coeus data structure to perform the computation
 * We read the graph from a csv file
 * We compute the orbits of the graph
 * We sort the active and dangling nodes
 * We create the adjacency matrix
 *
 * @param csv_name The file containing the graph
 * @param separator The file format separator
 * @param begin_by_zero
 * @param directed
 */
void POPPYm::init(std::string csv_name, std::string separator, bool directed) {
    // Read the csv file
    Graph nodes_and_arcs = read_graph(csv_name, separator, directed);
    nauty_routine(nodes_and_arcs);

    // PREPROCESSING
    sort_active_and_dangling_nodes();

    std::vector<std::vector<double>> mask_matrix = create_adjacency_matrix(active_nodes, orbits, map_in, map_out);
    this->mask_matrix = mask_matrix;

    if (nb_dangling_nodes > 0) {
        this->mask_matrix_dangling_nodes =
            create_adjacency_matrix(this->active_dangling_nodes, orbits, map_in, map_out);
    }

    this->is_initialized = true;
}

/**
 * @brief Initialize the PRs
 * We create a vector of size s with 1/n values (n is the number of nodes in the graph)
 *
 */
void POPPYm::initialize_PRs() {
    int n = this->number_of_nodes;
    auto public_key = this->keys.publicKey;

    std::vector<double> x(this->s, 1.0 / n);

    // We encode and enrypt the vector in CKKS to creat PRs
    Plaintext iv = this->cc->MakeCKKSPackedPlaintext(x);
    auto ct_iv = this->cc->Encrypt(public_key, iv);

    this->PRs = std::vector<Ciphertext<DCRTPoly>>(this->m / this->s, ct_iv);
}

/**
 * @brief Generate the crypto context
 * In this function, we set the vector size to the next power of two of the number of active nodes
 * or to the ring dimension divided by 2 if the ring dimension is smaller than the power of two of the number of active
 * nodes
 *
 * @param multDepth
 */
void POPPYm::gen_crypto_context(uint32_t multDepth) {
    printf("gen_crypto_context\n");

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

// Scaling parameters for crypto context
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits = 78;
    usint firstMod = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits = 59;
    usint firstMod = 60;
#endif

    // ringDims vector which give ringDim depending on the multiplicative depth (ringDims[3] is the ringDim needed for a
    // multDepth of 4)
    std::vector<int> ringDims = {8192,   16384,  16384,  16384,  32768,  32768,  32768,  32768,  32768,  65536,  65536,
                                 65536,  65536,  65536,  65536,  65536,  65536,  65536,  65536,  65536,  131072, 131072,
                                 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072,
                                 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072,
                                 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072,
                                 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072};

    uint32_t nbSlots = nextPowerOfTwo(number_active_nodes);

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecurityLevel(HEStd_128_classic);

    if (nbSlots < ringDims[multDepth - 1]) {
        parameters.SetBatchSize(nbSlots);
        this->s = nbSlots;
    } else
        this->s = ringDims[multDepth - 1] / 2;

    std::cout << "s : " << s << "(the next power of two of " << this->number_active_nodes
              << ", the number of active nodes)" << std::endl;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    this->cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    // cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    assert(cc->GetRingDimension() == ringDims[multDepth - 1]);
    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

    keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
}

/**
 * @brief Generate the crypto context and serialize it
 *
 * @param multDepth
 * @param folder
 */
void POPPYm::gen_crypto_context_and_serialize(uint32_t multDepth, std::string folder) {
    gen_crypto_context(multDepth);

    gen_rotation_keys();

    serialize_cc(cc, keys, multDepth);
}

/**
 * @brief Read the crypto context and keys from a folder
 *
 * @param folder
 */
void POPPYm::read_crypto_context_and_keys(std::string folder) {
    deserialize_cc(&cc, &keys, folder);
}

/**
 * @brief Generate the rotation keys
 *
 */
void POPPYm::gen_rotation_keys() {
    std::vector<int> steps;
    steps.push_back(0);
    for (int i = 1; i < s; i *= 2) {
        steps.push_back(i);
    }

    this->cc->EvalRotateKeyGen(keys.secretKey, steps);
}

/**
 * @brief Compute the PageRank
 *
 * @param nb_iteration
 */
void POPPYm::run(int nb_iteration) {
    for (int i = 1; i < nb_iteration; i++) {
        std::cout << "iteration " << i << std::endl;
        compute_PR();
    }

    // Last iteration
    std::cout << "iteration " << nb_iteration << std::endl;
    if (this->nb_dangling_nodes > 0) {
        std::cout << "compute_PR_dangling_nodes" << std::endl;
        compute_PR_dangling_node(nb_iteration);
    }
    compute_PR();
}

/**
 * @brief Compute the PageRank of the current iteration
 *
 * @param nb_iteration
 */
void POPPYm::compute_PR() {
    PRs = matmul_coeus_splited(cc, keys, subDiags, PRs);

    // End of PR algorithm
    for (int i = 0; i < m / s; i++) {
        cc->EvalMultInPlace(PRs[i], 0.85);
        cc->EvalAddInPlace(PRs[i], 0.15);
    }
}

/**
 * @brief Compute the PageRank for the dangling nodes
 * This function is called once at the last iteration
 *
 * @param nb_iteration
 */
void POPPYm::compute_PR_dangling_node(int nb_iteration) {
    PRs_deg0 = matmul_coeus_splited(cc, keys, subDiags_dangling_nodes, PRs);

    // End of PR algorithm
    for (int i = 0; i < m / s; i++) {
        cc->EvalMultInPlace(PRs_deg0[i], 0.85);
        cc->EvalAddInPlace(PRs_deg0[i], 0.15);
    }
}

/**
 * @brief Print the data center
 *
 */
void POPPYm::print_data_center() {
    std::cout << std::setw(50) << std::setfill('*') << '\n' << std::endl;
    std::cout << "number_of_nodes : " << number_of_nodes << std::endl;
    std::cout << "map_in size : " << map_in.size() << std::endl;
    std::cout << "map_out size : " << map_out.size() << std::endl;
    std::cout << "number_of_active_nodes : " << number_active_nodes << std::endl;
    std::cout << "mask_matrix size : " << mask_matrix.size() << std::endl;
    std::cout << std::setw(50) << std::setfill('*') << '\n' << std::endl;
}

/**
 * @brief Decrypt and print the PRs
 *
 */
void POPPYm::print_PRs() {
    std::vector<Plaintext> res(m / s);
    std::vector<double> vect_PRs;

    // Decrypt page rank of standard nodes
    for (int i = 0; i < m / s; i++) {
        cc->Decrypt(keys.secretKey, PRs[i], &res[i]);
        res[i]->SetLength(s);
        auto subVect = res[i]->GetRealPackedValue();
        vect_PRs.insert(vect_PRs.end(), subVect.begin(), subVect.end());
    }

    // Decrypt page rank of dangling nodes
    std::vector<double> vect_PRs_deg0;
    if (this->nb_dangling_nodes > 0) {
        std::vector<Plaintext> res_deg0(m / s);

        for (int i = 0; i < m / s; i++) {
            cc->Decrypt(keys.secretKey, PRs_deg0[i], &res_deg0[i]);
            res_deg0[i]->SetLength(s);
            auto subVect = res_deg0[i]->GetRealPackedValue();
            vect_PRs_deg0.insert(vect_PRs_deg0.end(), subVect.begin(), subVect.end());
        }
    }

    // Print the PRs in the order of the nodes
    for (int i = 0; i < number_of_nodes; i++) {
        auto it = std::find(active_nodes.begin(), active_nodes.end(), orbits[i]);
        if (it != active_nodes.end()) {
            int indice = it - active_nodes.begin();
            std::cout << i << " is " << vect_PRs[indice] << std::endl;
        } else {
            auto it_dangling_node = std::find(active_dangling_nodes.begin(), active_dangling_nodes.end(), orbits[i]);
            int indice = it_dangling_node - active_dangling_nodes.begin();
            std::cout << i << " is " << vect_PRs_deg0[indice] << std::endl;
        }
    }
}

/**
 * @brief Return the first ciphertext of the m/s ciphertexts
 *
 * @return Ciphertext<DCRTPoly>
 */
Ciphertext<DCRTPoly> POPPYm::get_PRs() {
    return PRs[0];
}

/**
 * @brief Main function Coeus algorithm
 *
 */
void run_poppy_m(int nb_iteration) {
    POPPYm data_center;

    string filename = "../Graphs/students.txt";
    data_center.pre_processing(GRAPH_DIR + filename, " ", true);  // for student graph

    data_center.run(nb_iteration);

    std::cout << "\nAfter " << nb_iteration << " iterations, the PR of" << std::endl;
    data_center.print_PRs();
    // Trim the extension
    // string dot_name = filename.substr(0, filename.find_last_of("."));
    // data_center.export_to_dot(DOT_DIR + dot_name + "_" + std::to_string(nb_iteration) + ".dot");
    std::map<int, double> results = data_center.export_results();
    // Write the results to a file
    std::ofstream file;
    file.open("../results/" + filename + "_coeus-" + std::to_string(nb_iteration) + ".csv");
    file << "Node,PR" << std::endl;
    for (auto result : results) {
        file << result.first << "," << result.second << std::endl;
    }
}

std::map<int, double> run_poppy_m(int nb_iteration, std::string filename, std::string delimitor) {
    POPPYm data_center;
    data_center.pre_processing(GRAPH_DIR + filename, delimitor, true);

    data_center.run(nb_iteration);

    std::map<int, double> results = data_center.export_results();
    data_center.get_cc()->ClearEvalMultKeys();
    data_center.get_cc()->ClearEvalAutomorphismKeys();
    data_center.get_cc()->ClearEvalSumKeys();
    return results;
}
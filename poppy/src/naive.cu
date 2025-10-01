#include "naive.cuh"

/**
 * @brief Pre-processing of the graph
 * Read the csv file (init)
 * Generate the crypto context (gen_crypto_context)
 * We end by initializing the PRs (initialize_PRs)
 *
 * @param csv_name file name
 * @param option
 * @param folder
 * @param separator
 * @param begin_by_zero
 * @param directed
 */
void Naive::pre_processing(std::string csv_name, std::string separator, bool directed) {
    init(csv_name, separator, directed);

    int multDepth = 8; 
    gen_crypto_context(multDepth);

    initialize_PRs();
}

/**
 * @brief Read the csv file and return the graph in the form of two maps of in edges and out edges
 *
 * @param csv_name file name
 * @param separator csv separator
 * @param begin_by_zero
 * @param directed
 */
void Naive::init(std::string csv_name, std::string separator, bool directed) {
    Graph nodes_and_arcs = read_graph(csv_name, separator, directed);

    // Creation du graph G
    DYNALLSTAT(graph, g, g_sz);
    static DEFAULTOPTIONS_GRAPH(options);

    // Initialisez les options du graph
    options.writeautoms = TRUE;
    options.writemarkers = TRUE;
    options.defaultptn = TRUE;  // Pas de coloration
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
    // Vérifiez la version de nauty
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    DYNALLOC2(graph, g, g_sz, m, n, "malloc");

    // Initialisez le graph
    EMPTYGRAPH(g, m, n);

    // Créez les maps pour les sommets entrants et sortants de G
    for (auto edge : nodes_and_arcs.arcs) {
        map_out.insert(std::pair<int, int>(trad_tab[edge.first], trad_tab[edge.second]));
        map_in.insert(std::pair<int, int>(trad_tab[edge.second], trad_tab[edge.first]));
        ADDONEARC(g, trad_tab[edge.first], trad_tab[edge.second], m);
    }
    this->number_of_nodes = n;
    DYNFREE(g, g_sz);
}

/**
 * @brief Generate the crypto context
 *
 * @param multDepth
 */
void Naive::gen_crypto_context(uint32_t multDepth) {
    // Crypto material

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);

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

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl
              << std::endl;

    keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
}

/**
 * @brief Compute the PageRank for a given number of iteration
 *
 * @param nb_iteration
 */
void Naive::run(int nb_iteration) {
    for (int i = 1; i <= nb_iteration; i++) {
        std::cout << "iteration " << i << std::endl;
        compute_PR();
    }
}

/**
 * @brief Compute the PageRank for the current iteration
 *
 */
void Naive::compute_PR() {
    int n = this->number_of_nodes;
    auto cc = this->cc;
    auto keys = this->keys;


    // We create ciphers of 0

    std::vector<double> x_0 = {0.0};
    Plaintext zero = cc->MakeCKKSPackedPlaintext(x_0);
    auto ct_zero = cc->Encrypt(keys.publicKey, zero);

    std::map<int, Ciphertext<DCRTPoly>> PR_end_iteration;

    for (int i = 0; i < n; i++) {
        auto nodes_in = map_in.equal_range(i);  // To get all the nodes getting in i
        auto sum = ct_zero;

        for (auto it = nodes_in.first; it != nodes_in.second; ++it) {  // For each nodes_in
            int j = it->second;                                        // node j -> i

            // If j as deg_out > 1
            if (map_out.count(j) > 1) {
                auto search = PRs.find(j);
                Ciphertext<DCRTPoly> val = cc->EvalMult(search->second, 1.0 / map_out.count(j));  // PR(j)/deg_out(j)
                sum = cc->EvalAdd(sum, val);
            } else {
                auto search = PRs.find(j);
                sum = cc->EvalAdd(sum, search->second);
            }
        }

        sum = cc->EvalMult(sum, 0.85);
        sum = cc->EvalAdd(sum, 0.15);

        PR_end_iteration.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, sum));
    }

    this->PRs = PR_end_iteration;
}

std::map<int, Ciphertext<DCRTPoly>> Naive::get_PRs() {
    return PRs;
}

/**
 * @brief Decrypt and print the PageRank
 *
 */
void Naive::print_PRs() {
    for (int i = 0; i < this->number_of_nodes; i++) {
        Plaintext res;
        std::cout.precision(8);
        this->cc->Decrypt(keys.secretKey, PRs[i], &res);
        res->SetLength(1);
        auto vect_res = res->GetRealPackedValue();
        std::cout << i << " is " << vect_res[0] << std::endl;
    }
}

/**
 * @brief Initialize the PageRank, all the nodes have the same PR equal to 1/n, n the number of nodes
 *
 */
void Naive::initialize_PRs() {
    int n = this->number_of_nodes;
    auto cc = this->cc;
    auto public_key = this->keys.publicKey;
    std::vector<double> x = {1.0 / n};

    // We create the initialisation
    Plaintext iv = cc->MakeCKKSPackedPlaintext(x);
    auto ct_iv = cc->Encrypt(public_key, iv);

    // We fill the map with the initialisation ciphertext
    for (int i = 0; i < n; i++) {
        this->PRs.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, ct_iv));
    }
}

void Naive::gen_rotation_keys() {}

/**
 * @brief Main function to run the naive algorithm
 *
 */
void run_naive(int nb_iteration) {
    Naive data_center;

    std::string filename = "../Graphs/students.txt";
    data_center.pre_processing(GRAPH_DIR + filename, " ", true);  // for student graph
    data_center.print_data_center();

    data_center.run(nb_iteration);

    std::cout << "\nAfter " << nb_iteration << " iterations, the PR of" << std::endl;
    data_center.print_PRs();
}
#include "poppy.cuh"

using std::string;
using namespace lbcrypto;

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

class Clear : public DataCenter<std::map<int, double>> {
public:
    void pre_processing(std::string csv_name, std::string separator = ",", bool directed = true);
    void init(std::string csv_name, std::string separator = ",", bool directed = true);
    void run(int nb_iteration);

    // getters
    std::map<int, double> get_PRs();

    // printers
    void print_PRs() override;

    // generation of the crypto material
    void gen_crypto_context(uint32_t multDepth);
    void gen_rotation_keys();

    std::map<int, double> export_results() { return get_PRs(); };

private:
    int max_deg_in = 0; 
    void initialize_PRs();
    void compute_PR();

    std::map<int, double> PRs;
};

void Clear::pre_processing(std::string csv_name, std::string separator, bool directed) {
    init(csv_name, separator, directed);

    gen_crypto_context(16);

    initialize_PRs();
}

void Clear::init(std::string csv_name, std::string separator, bool directed) {
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

    // Compue the deg_in max in the sub graph
    for (auto it = map_in.begin(); it != map_in.end(); ) {
        auto range = map_in.equal_range(it->first);
        int count = std::distance(range.first, range.second);

        if (count > this->max_deg_in)
            this->max_deg_in = count;

        it = range.second; // move to the next key
    }
}

void Clear::run(int nb_iteration) {
    std::cout << "Start clear computation" << std::endl;
    for (int i = 1; i <= nb_iteration; i++) {
        std::cout << "iteration " << i << std::endl;
        compute_PR();
    }
}

/**
 * @brief Page rank computation
 *
 */
void Clear::compute_PR() {
    int n = this->number_of_nodes;

    std::map<int, double> PR_end_iteration;  // Storage of the pagerank values

    // Initialization of the value 0
    double zero = 0.0;

    // Boucle sur chaque nœud pour calculer le PageRank
    for (int i = 0; i < n; i++) {
        auto nodes_in = map_in.equal_range(i);  // To get all nodes with a connexion to i
        double sum = zero;

        // Calculation of the sum of the pagerank of the incoming nodes
        for (auto it = nodes_in.first; it != nodes_in.second; ++it) {
            int j = it->second;  // nœud j -> i

            // If out-degree if j > 1
            if (map_out.count(j) > 1) {
                double val = PRs[j] / map_out.count(j);  // PR(j)/deg_out(j)
                sum += val;
            } else {
                sum += PRs[j];
            }
        }

        // Application de la formule de PageRank
        sum = 0.85 * sum + (0.15 / this->max_deg_in);

        // Stockage du PageRank calculé pour le nœud i
        PR_end_iteration.insert(std::pair<int, double>(i, sum));
    }

    // Stockage des PageRank calculés
    this->PRs = PR_end_iteration;
}

std::map<int, double> Clear::get_PRs() {
    return PRs;
}

/**
 * @brief At first iteration, all nodes have the same PageRank value equal to 1/n (n = number of nodes)
 *
 */
void Clear::initialize_PRs() {
    int n = this->number_of_nodes;

    // Calcul de la valeur initiale du PageRank
    double initial_value = 1.0 / n;

    // Stockage de la valeur initiale du PageRank
    std::map<int, double> PRs;

    // Remplissage de la map avec la valeur initiale
    for (int i = 0; i < n; i++) {
        PRs.insert(std::pair<int, double>(i, initial_value));
    }

    // Stockage des PageRanks initialisés
    this->PRs = PRs;
}

void Clear::print_PRs() {
    std::map<int, double> PRs = get_PRs();

    // Affichage des PageRanks
    for (auto PR : PRs) {
        std::cout << PR.first << " is " << PR.second << std::endl;
    }
}
void Clear::gen_crypto_context(uint32_t multDepth) {}
void Clear::gen_rotation_keys() {}

std::map<Node, double> clear_PR(int nb_iteration, std::string filename, std::string delimitor) {
    Clear data_center;

    data_center.pre_processing(GRAPH_DIR + filename, delimitor, true);
    // data_center.print_data_center();

    data_center.run(nb_iteration);
    // data_center.print_PRs();

    return data_center.export_results();
}
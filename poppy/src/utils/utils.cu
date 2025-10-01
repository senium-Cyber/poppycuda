#include "include/utils.cuh"

#include <fstream>
#include <iostream>
#include <vector>

using namespace lbcrypto;

uint32_t nextPowerOfTwo(uint32_t n) {
    --n;

    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    return n + 1;
}

void delete_vector_of_ciphers(std::vector<Ciphertext<DCRTPoly>> vector_of_ciphers) {
    for (auto& ct : vector_of_ciphers) {
        ct.reset();
    }
    vector_of_ciphers.clear();
    vector_of_ciphers.shrink_to_fit();
}

void delete_vector_of_plains(std::vector<Plaintext> vector_of_plains) {
    for (auto& ct : vector_of_plains) {
        ct.reset();
    }
    vector_of_plains.clear();
    vector_of_plains.shrink_to_fit();
}

void delete_map_of_plains(std::map<int, std::vector<Plaintext>> map_of_plains) {
    for (auto& ct : map_of_plains) {
        delete_vector_of_plains(ct.second);
    }
    map_of_plains.clear();
}

void print_map_in_out(std::multimap<int, int> map, int n) {
    for (int i = 0; i < n; i++) {
        if (map.count(i) != 0) {
            std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> values;
            values = map.equal_range(i);
            std::cout << i << " => ";

            for (auto it = values.first; it != values.second; ++it)
                std::cout << it->second << " ";
            std::cout << std::endl;
        }
    }
}

void print_map_u_j(std::multimap<int, triplet> map) {
    auto it = map.begin();
    auto it2 = map.begin();
    while (it != map.end()) {
        std::cout << it->first << " => ";
        // auto it2 = it;
        while (it2 != map.end() && it2->first == it->first) {
            it2->second.print_triplet();
            std::cout << '\n';
            it2++;
        }
        it = it2;
    }
}

void print_map_reconstruction_u_j(std::vector<std::multimap<int, quadruplet>> map) {
    for (int i = 0; i < map.size(); i++) {
        auto it = map[i].begin();
        auto it2 = map[i].begin();
        while (it != map[i].end()) {
            std::cout << it->first << " => ";
            // auto it2 = it;
            while (it2 != map[i].end() && it2->first == it->first) {
                it2->second.print_quadruplet();
                std::cout << '\n';
                it2++;
            }
            it = it2;
        }
    }
}

void print_PR(
    int n,
    std::multimap<int, triplet> map_u_j,
    std::multimap<int, triplet> map_deg0,
    std::vector<double> PR,
    std::vector<double> PR_deg0,
    int* orbits) {
    for (int i = 0; i < n; i++) {
        if (map_u_j.count(orbits[i]) != 0) {
            int indice = map_u_j.find(orbits[i])->second.indice;
            std::cout << i << " is " << PR[indice] << std::endl;
        } else if (map_deg0.count(orbits[i]) != 0) {
            int indice = map_deg0.find(orbits[i])->second.indice;
            std::cout << i << " is " << PR_deg0[indice] << std::endl;
        } else
            std::cout << i << " is 0.15" << std::endl;
    }
}

void print_cipher(Ciphertext<DCRTPoly> ct, CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, int s) {
    Plaintext res;
    std::vector<double> plain;
    std::cout.precision(8);
    cc->Decrypt(keys.secretKey, ct, &res);
    res->SetLength(s);
    plain = res->GetRealPackedValue();
    for (int i = 0; i < s; i++)
        std::cout << plain[i] << " ";
    std::cout << std::endl;
}

bool is_number(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(), ::isdigit);
}

Graph read_graph(string name, string sep, bool directed) {
    std::fstream newfile;

    Graph graph;  // nodes and arcs

    newfile.open(name, std::ios::in);
    if (newfile.is_open()) {
        string line;
        int line_nb = 0;
        while (getline(newfile, line)) {
            if (line_nb == 0 && line[0] == '#') {
                // Check if it starts with '#'

                // Trim the first character
                line = line.substr(1);
                // Check if first line = "undirected"
                // Erase first spaces
                line.erase(0, line.find_first_not_of(" "));
                // Get first word
                line = line.substr(0, line.find(" "));
                // Transform to lowercase
                std::transform(line.begin(), line.end(), line.begin(), ::tolower);
                if (line.compare("undirected") > 0) {
                    directed = false;
                } else if (line.compare("directed") == 0) {
                    directed = true;
                }
            }
            if (line[0] == '#')
                continue;
            if (line.find(sep) != string::npos) {
                size_t pos = line.find(sep);
                string node_in = line.substr(0, pos);
                string node_out = line.substr(pos + 1, line.find_first_not_of("0123456789", pos + 1) - pos - 1);

                if (is_number(node_in) && is_number(node_out)) {
                    int i = stoi(node_in);
                    int j = stoi(node_out);

                    graph.nodes.insert(i);
                    graph.nodes.insert(j);

                    if (!directed) {
                        if (line_nb % 2 != 0) {
                            if (graph.arcs.find({i, j}) == graph.arcs.end()) {
                                std::pair<int, int> arc = {i, j};
                                graph.arcs.insert(arc);
                            } else if (graph.arcs.find({j, i}) == graph.arcs.end()) {
                                std::pair<int, int> arc = {j, i};
                                graph.arcs.insert(arc);
                            }
                        } else {
                            if (graph.arcs.find({j, i}) == graph.arcs.end()) {
                                std::pair<int, int> arc = {j, i};
                                graph.arcs.insert(arc);
                            } else if (graph.arcs.find({i, j}) == graph.arcs.end()) {
                                std::pair<int, int> arc = {i, j};
                                graph.arcs.insert(arc);
                            }
                        }
                    } else {
                        std::pair<int, int> arc = {i, j};
                        graph.arcs.insert(arc);
                    }
                    line_nb++;
                }
            }
        }
        newfile.close();
        int min = *std::min_element(graph.nodes.begin(), graph.nodes.end());
        int max = *std::max_element(graph.nodes.begin(), graph.nodes.end());
        if (graph.nodes.size() - 1 != max && min != 0) {
            std::map<int, int> trad_tab;
            int index = 0;
            Graph g_bis;
            for (auto node : graph.nodes) {
                trad_tab[node] = index;
                g_bis.nodes.insert(index);
                index++;
            }
            for (auto arc : graph.arcs) {
                std::pair<int, int> arc_bis = {trad_tab[arc.first], trad_tab[arc.second]};
                g_bis.arcs.insert(arc_bis);
            }
            return g_bis;
        }
        return graph;
    } else {
        std::cerr << "Unable to open csv file" << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Parse a directed graph in a dot file and return the graph and the partition of the nodes in datacenters.
 *
 * @param filename
 * @param partition
 * @return Graph
 */
std::vector<Graph> parse_graph_multi_DC(string filename, std::map<int /*node*/, int /*datacenter*/>& partition) {
    std::fstream newfile;

    std::vector<Graph> graphs;

    newfile.open(filename, std::ios::in);
    if (newfile.is_open()) {
        bool is_subgraph = false;
        string line;
        while (getline(newfile, line)) {
            // Parse dot line
            if (line.find("subgraph") != std::string::npos) {
                is_subgraph = true;
                graphs.emplace_back();
            } else if (size_t pos = line.find("->"); pos != std::string::npos) {
                string node_in = line.substr(0, pos);
                string node_out = line.substr(pos + 2);

                // Nettoyage des espaces et du point-virgule
                node_in.erase(std::remove(node_in.begin(), node_in.end(), ' '), node_in.end());
                node_out.erase(std::remove(node_out.begin(), node_out.end(), ' '), node_out.end());
                node_out.erase(std::remove(node_out.begin(), node_out.end(), ';'), node_out.end());

                if (is_number(node_in) && is_number(node_out)) {
                    int i = stoi(node_in);
                    int j = stoi(node_out);

                    if (is_subgraph) {
                        graphs.back().nodes.insert(i);
                        graphs.back().nodes.insert(j);
                        graphs.back().arcs.insert({i, j});

                        if (partition.find(i) == partition.end())
                            partition.insert({i, graphs.size() - 1});
                        if (partition.find(j) == partition.end())
                            partition.insert({j, graphs.size() - 1});
                    } else {
                        // Case it is an inter-DC edge
                        int graph_index = partition.at(i);
                        int graph_index2 = partition.at(j);

                        graphs[graph_index].arcs.insert({i, j});
                        graphs[graph_index].nodes.insert(i);
                        graphs[graph_index].nodes.insert(j);

                        graphs[graph_index2].arcs.insert({i, j});
                        graphs[graph_index2].nodes.insert(i);
                        graphs[graph_index2].nodes.insert(j);
                    }
                }
            } else if (line.find("}") != std::string::npos) {
                // Case it is the end of a subgraph
                is_subgraph = false;
            } else {
                line.erase(std::remove(line.begin(), line.end(), ','), line.end());
                line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
                line.erase(std::remove(line.begin(), line.end(), ';'), line.end());
                if (is_number(line)) {
                    int i = stoi(line);
                    if (partition.find(i) == partition.end())
                        partition.insert({i, graphs.size() - 1});
                }
            }
        }
        newfile.close();
        return graphs;
    } else {
        std::cerr << "Unable to open csv file" << std::endl;
        exit(EXIT_FAILURE);
    }
}

string choose_graph() {
    // Open directory and read file name
    DIR* dr;
    struct dirent* en;
    std::vector<string> filenames;
    dr = opendir(GRAPH_DIR);  // open all directory
    if (dr) {
        while ((en = readdir(dr)) != nullptr) {
            // Store filenames in vect
            if (en->d_type == DT_REG)
                filenames.emplace_back(en->d_name);
        }
        closedir(dr);  // close all directory
    }

    // Print filenames
    std::cout << "Choose a graph among the following:" << std::endl;
    for (int i = 0; i < filenames.size(); i++) {
        std::cout << i << ": " << filenames[i] << std::endl;
    }

    // Choose a file
    int choice;
    std::cin >> choice;
    // Verify the choice
    while (choice < 0 || choice >= filenames.size()) {
        std::cout << "Please choose a number between 0 and " << filenames.size() - 1 << std::endl;
        std::cin >> choice;
    }

    return filenames[choice];
}

/**
 * @brief Compute the sub partition of the datacenter
 * We keep only the needed information for the datacenter
 *
 * @param partition
 * @param map_in
 * @return std::map<int, int>
 */
std::map<int, int>
compute_sub_partition(std::map<int, int> partition, std::set<std::pair<int, int>> arcs, int datacenter_id) {
    std::map<int, int> sub_partition;
    for (auto arc : arcs) {
        if (partition[arc.first] == datacenter_id || partition[arc.second] == datacenter_id) {
            sub_partition[arc.first] = partition[arc.first];
            sub_partition[arc.second] = partition[arc.second];
        }
    }
    return sub_partition;
}

std::pair<std::multimap<int, int>, std::multimap<int, int>> process_maps(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    std::map<int, int>& orbits,
    int datacenter_id) {
    std::map<int, int> trad_tab;  // Reindexation des noeuds
    std::set<int> distant_nodes;  // Ensemble des noeuds distants
    int index = 0;
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id) {
            trad_tab[node] = index;
            index++;
        } else {
            distant_nodes.insert(node);
        }
    }
    number_of_nodes = index;  // Local number of nodes

    std::multimap<int, int> map_in;
    std::multimap<int, int> map_out;

    for (auto edge : nodes_and_arcs->arcs) {
        map_out.insert(std::pair<int, int>(edge.first, edge.second));
        map_in.insert(std::pair<int, int>(edge.second, edge.first));
    }

    for (auto node : distant_nodes) {
        trad_tab[node] = index;
        index++;
    }

    int index_orbits = 0;
    for (auto node_trad : trad_tab) {
        orbits[node_trad.first] = index_orbits;
        index_orbits++;
    }

    return {map_in, map_out};
}

/**
 * @brief Compute the orbits for the local graph of the datacenter
 *  using sparse nauty (for sparse graphs)
 * Returns the map of the in and out edges of the datacenter
 *
 * @param nodes_and_arcs The graph of the datacenter
 * @param partition Partition of the nodes
 * @param number_of_nodes Number of local nodes
 * @param orbits Orbits
 * @return std::pair<std::multimap<int, int>, std::multimap<int, int>>
 */
std::pair<std::multimap<int, int>, std::multimap<int, int>> sparse_nauty_routine(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    int& nb_local_nodes_after_pruning,
    std::map<int, int>& orbits,
    int datacenter_id,
    int& duration,
    std::set<Node>& set_dangling_nodes) {
    // Creation of the graph G
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, local_orbits, orbits_sz);
    static DEFAULTOPTIONS_SPARSEDIGRAPH(options);
    statsblk stats;
    SG_DECL(sg);  // Declare and initialize sparse graph structure

    // We set the options of the graph.
#ifdef DEBUG
    options.writemarkers = TRUE;  // Write levels on terminal
#endif
    options.defaultptn = FALSE;  // Activate coloration

    std::map<int, int> trad_tab;  // Reindexation of the nodes
    std::set<int> distant_nodes;  // Set of distant nodes
    int index = 0;
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id) {
            trad_tab[node] = index;
            index++;
        } else
            distant_nodes.insert(node);
    }

    int n = index;
    number_of_nodes = n;
    int m = SETWORDSNEEDED(n);

    SG_INIT(sg);

#ifdef DEBUG
    std::cout << "Trad tab" << std::endl;
    for (auto it : trad_tab) {
        std::cout << it.first << " -> " << it.second << std::endl;
    }
#endif

    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");

    local_orbits = NULL;  // Initialiser à NULL
    DYNALLOC1(int, local_orbits, orbits_sz, n, "malloc");
    SG_ALLOC(sg, n, nodes_and_arcs->arcs.size(), "malloc");

    int nb_arcs = 0;
    std::set<int> border_nodes;
    for (auto arc : nodes_and_arcs->arcs) {
        if (partition->at(arc.first) == datacenter_id && partition->at(arc.second) == datacenter_id)
            nb_arcs += 1;  // Number of inter edges
        if (partition->at(arc.first) != datacenter_id && partition->at(arc.second) == datacenter_id)
            border_nodes.insert(arc.second);
        // Remove for opti (not in same orbit but same PR)
        if (partition->at(arc.first) == datacenter_id && partition->at(arc.second) != datacenter_id)
            border_nodes.insert(arc.first);
    }

    sg.nv = n;         // Number of vertices
    sg.nde = nb_arcs;  // Number of directed inter edges

    // Color the graph
    // Attribute for each local node connected to a distant node a unique color,
    // This is the same as having a unique orbits for each distant nodes
    // A colouring (partition) of the vertices is specified by a pair of arrays, usually called lab and ptn
    // The array lab contains a list of the vertices in some order
    // The array ptn indicates the division into colours : if ptn[i] = 0, then a cell(colour class) ends at position
    // i.
    std::fill_n(ptn, n, 1);
    // In lab write first all the intern node and then the others
    int lab_index = 0;

    // Add all local nodes that are node border nodes
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id && border_nodes.find(node) == border_nodes.end()) {
            assert(lab_index < n);
            lab[lab_index] = trad_tab.at(node);
            lab_index++;
        }
    }
    if (lab_index > 0) {
        ptn[lab_index - 1] = 0;
    }
    // Then add the border nodes, each with a different color
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id && border_nodes.find(node) != border_nodes.end()) {
            assert(lab_index < n);
            lab[lab_index] = trad_tab.at(node);
            ptn[lab_index] = 0;
            lab_index++;
        }
    }
    assert(lab_index == n);

    std::multimap<int, int> map_in;
    std::multimap<int, int> map_out;
    std::multimap<int, int> map_out_without_inter_edges;
    int edge_index = 0;
    // Create the map of the in and out edges, and the nauty graph
    for (auto edge : nodes_and_arcs->arcs) {
        map_out.insert(std::pair<int, int>(edge.first, edge.second));
        map_in.insert(std::pair<int, int>(edge.second, edge.first));

        if (partition->at(edge.first) == datacenter_id && partition->at(edge.second) == datacenter_id)
            map_out_without_inter_edges.insert(std::pair<int, int>(edge.first, edge.second));

        if (partition->at(edge.first) == datacenter_id && partition->at(edge.second) == datacenter_id) {
            assert(edge_index < nb_arcs);
            auto j = trad_tab.at(edge.second);
            sg.e[edge_index++] = j;  // Add each neighbor consecutively
        }
    }

    // Out degree of each node of sg
    int i = 0;
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id) {
            auto deg = map_out_without_inter_edges.count(node);
            sg.v[trad_tab.at(node)] = i;    // table of starting indices for neighbors of each vertex in sg.e
            sg.d[trad_tab.at(node)] = deg;  // table of outdegrees
            i += deg;
        }
    }
    assert(sg.dlen == n);
    assert(sg.vlen == n);

#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    sparsenauty(&sg, lab, ptn, local_orbits, &options, &stats, NULL);
#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
#endif

    // We add distant nodes to trad_tab
    for (auto node : distant_nodes) {
        trad_tab[node] = index;
        index++;
    }

    // Copy Allocated memory for orbits in a map
    // Create set of orbits representants
    std::set<int> orbits_representants;
    for (auto node_trad : trad_tab) {
        if (partition->at(node_trad.first) == datacenter_id) {
            orbits[node_trad.first] = local_orbits[node_trad.second];
            orbits_representants.insert(local_orbits[node_trad.second]);

            // Create set of dangling nodes
            if (map_out.count(node_trad.first) == 0)
                set_dangling_nodes.insert(orbits[node_trad.first]);
        } else
            orbits[node_trad.first] = node_trad.second;

#ifdef DEBUG
        std::cout << "Node " << node_trad.first << " is in orbit " << orbits[node_trad.first] << std::endl;
#endif
    }

    nb_local_nodes_after_pruning = orbits_representants.size();
    std::cout << "Finished to compute the orbits" << std::endl;

    // We delete sg
    SG_FREE(sg);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    DYNFREE(local_orbits, orbits_sz);
    return {map_in, map_out};
}

void printPatternGroups(const std::map<std::map<int, float>, std::vector<int>>& patternGroups) {
    std::cout << "-----patternGroups-----" << std::endl;
    // Print pattern groups
    for (const auto& pg : patternGroups) {
        std::cout << "Pattern: ";
        // Print the pattern (incoming edges)
        for (const auto& edge : pg.first) {
            std::cout << "(" << edge.first << ":" << edge.second << ") ";
        }
        std::cout << "\nNodes: ";
        // Print nodes having this pattern
        for (const auto& node : pg.second) {
            std::cout << node << " ";
        }
        std::cout << "\n-------------------\n";
    }
}

// // Fonction auxiliaire pour normaliser un pattern
std::map<int, float> normalizePattern(const std::map<int, float>& pattern, int nodeId) {
    std::map<int, float> normalized = pattern;
    // Si le pattern contient une référence à soi-même
    auto selfRef = normalized.find(nodeId);
    if (selfRef != normalized.end()) {
        // Remplacer l'ID du nœud par une valeur spéciale (par exemple -1)
        float weight = selfRef->second;
        normalized.erase(selfRef);
        normalized[-1] = weight;
    }
    return normalized;
}

/**
 * Cette fonction, findEquivalentNodes_with_dangling, prend en paramètre une multimap et un entier.
 *
 * @param map_in : une multimap où chaque clé représente un nœud cible et chaque valeur représente un nœud source.
 * @param n : un entier représentant le nombre total de nœuds.
 *
 * @return une map où chaque clé est un nœud et chaque valeur est le représentant équivalent de ce nœud
 * et la liste des dangling nodes sur laquelle faire les calculs à la derniere iteration seulement.
 *
 */

std::pair<std::map<int, int>, std::set<int>> findEquivalentNodes_with_dangling(
    const std::multimap<int, int>& map_in,
    int n,
    const std::set<int>& protected_nodes_out = std::set<int>(),
    const std::set<int>& protected_nodes_in = std::set<int>()) {
    // Calculer le degré sortant pour chaque nœud
    std::map<int, int> outDegrees;
    for (const auto& entry : map_in) {
        int sourceNode = entry.second;
        outDegrees[sourceNode]++;
    }

    // nodeToIncomingEdges[node] = map<sourceNode, weight>
    // où weight = 1.0f/outDegrees[sourceNode]
    std::vector<std::map<int, float>> nodeToIncomingEdges(n);

    // Initialisation des arêtes entrantes
    for (const auto& entry : map_in) {
        int targetNode = entry.first;
        int sourceNode = entry.second;
        nodeToIncomingEdges[targetNode][sourceNode] = 1.0f / outDegrees[sourceNode];
    }
    DisjointSet ds(n);

    bool changement_made = true;
    while (changement_made) {
        changement_made = false;

        // Étape 1 : Regrouper les noeuds par pattern d'arêtes entrantes
        std::map<std::map<int, float>, std::vector<int>> patternGroups;
        for (int node = 0; node < n; ++node) {
            // Si le noeud n'est pas un représentant ou est protégé, on ne le considère pas
            if (ds.find(node) == node && protected_nodes_out.find(node) == protected_nodes_out.end()) {
                auto normalizedPattern = normalizePattern(nodeToIncomingEdges[node], node);
                patternGroups[normalizedPattern].push_back(node);
            }
        }

        // Étape 2 : Unifier les groupes ayant plus d'un noeud
        for (auto& pg : patternGroups) {
            const auto& group = pg.second;
            if (group.size() > 1) {
                // Unir tous les noeuds de ce groupe, sauf si le premier noeud est protégé
                int first = group[0];
                if (protected_nodes_out.find(first) != protected_nodes_out.end()) {
                    continue;  // Skip si le premier noeud est protégé
                }
                for (size_t i = 1; i < group.size(); ++i) {
                    if (protected_nodes_out.find(group[i]) == protected_nodes_out.end()) {
                        ds.unite(first, group[i]);
                        changement_made = true;
                    }
                }
            }
        }
        // Si aucun changement n'a été fait, on peut arrêter
        if (!changement_made)
            break;

        // Étape 3 : Mettre à jour les patterns après union

        for (int node = 0; node < n; ++node) {
            if (ds.find(node) == node) {
                // On doit regrouper les sources maintenant unifiés
                std::map<int, float> newEdges;
                for (auto& src_w : nodeToIncomingEdges[node]) {
                    int src = src_w.first;
                    float w = src_w.second;
                    int srcRep = ds.find(src);
                    newEdges[srcRep] += w;
                }
                nodeToIncomingEdges[node] = std::move(newEdges);
            } else {
                // Les nœuds non-représentants n'ont plus besoin d'un pattern distinct.
                nodeToIncomingEdges[node].clear();
            }
        }
    }

    // Construire la map des équivalences finales
    std::map<int, int> equivalent;
    std::set<int> dangling;
    for (int i = 0; i < n; i++) {
        int rep = ds.find(i);
        if (outDegrees[i] != 0) {
            // Cas facile : ce n'est pas un dangling node
            equivalent[i] = rep;
        } else {
            // Cas dangling node
            // 1er cas : i est un dangling node représenté par un autre noeud
            if (i != rep) {
                equivalent[i] = rep;
            }
            // 2eme cas : i est un dangling node représenté par lui même
            else {
                for (int j = 0; j < n; j++) {
                    // On cherche tous les noeuds représentés par i
                    if (ds.find(j) == i) {
                        equivalent[i] = j;
                    } else {
                        // si i n'est pas protégé, on l'ajoute aux dangling nodes
                        if (protected_nodes_in.find(i) == protected_nodes_in.end()) {
                            equivalent[i] = i;
                            dangling.insert(i);
                        }
                        // si i est protégé, on ne l'ajoute pas aux dangling nodes et on lui assigne le représentant
                        else {
                            equivalent[i] = i;
                        }
                    }
                }
            }
        }
    }

    return {equivalent, dangling};
}

std::pair<std::multimap<int, int>, std::multimap<int, int>> equivalent_nodes_routine(
    Graph* nodes_and_arcs,
    const std::map<int, int>* partition,
    int& number_of_nodes,
    int& nb_local_nodes_after_pruning,
    std::map<int, int>& orbits,
    int datacenter_id,
    int& duration,
    std::set<int>& set_dangling_nodes) {
    // Creation of the graph G
    std::map<int, int> trad_tab;  // Reindexation of the nodes
    std::set<int> distant_nodes;  // Set of distant nodes
    int index = 0;
    for (auto node : nodes_and_arcs->nodes) {
        if (partition->at(node) == datacenter_id) {
            trad_tab[node] = index;
            index++;
        } else
            distant_nodes.insert(node);
    }
    number_of_nodes = index;  // Local number of nodes

    std::multimap<int, int> map_in;
    std::multimap<int, int> map_out;
    std::multimap<int, int> map_in_without_inter_edges;
    std::set<int> protected_nodes_out;
    std::set<int> protected_nodes_in;

    // Create the map of the in and out edges
    for (auto edge : nodes_and_arcs->arcs) {
        map_out.insert(std::pair<int, int>(edge.first, edge.second));
        map_in.insert(std::pair<int, int>(edge.second, edge.first));

        if (partition->at(edge.first) == datacenter_id && partition->at(edge.second) == datacenter_id)
            map_in_without_inter_edges.insert(std::pair<int, int>(trad_tab[edge.second], trad_tab[edge.first]));

        // If the edge is an incoming edge from a distant node to the datacenter, we protect the distant node
        if (partition->at(edge.first) != datacenter_id && partition->at(edge.second) == datacenter_id)
            protected_nodes_out.insert(trad_tab[edge.second]);

        // If the edge is an outgoing edge from the datacenter to a distant node, we protect the datacenter node
        if (partition->at(edge.first) == datacenter_id && partition->at(edge.second) != datacenter_id)
            protected_nodes_out.insert(trad_tab[edge.first]);

        // If the edge is an outgoing edge from the datacenter to a distant node, we protect the datacenter node
        if (partition->at(edge.first) == datacenter_id && partition->at(edge.second) != datacenter_id)
            protected_nodes_in.insert(trad_tab[edge.first]);
    }

#ifdef EXPE
    auto start = std::chrono::high_resolution_clock::now();
#endif
    // Use findEquivalentNodes to instantiate orbits
    auto set = findEquivalentNodes_with_dangling(
        map_in_without_inter_edges, number_of_nodes, protected_nodes_out, protected_nodes_in);
    auto local_orbits = set.first;
    set_dangling_nodes = set.second;
#ifdef EXPE
    auto end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
#endif

    // We add distant nodes to trad_tab
    for (auto node : distant_nodes) {
        trad_tab[node] = index;
        index++;
    }

    // Copy Allocated memory for orbits in a map
    // Create set of equivalent nodes representants
    std::set<int> equiv_nodes_representants;
    for (auto node_trad : trad_tab) {
        if (partition->at(node_trad.first) == datacenter_id) {
            orbits[node_trad.first] = local_orbits[node_trad.second];
            equiv_nodes_representants.insert(local_orbits[node_trad.second]);
        } else
            orbits[node_trad.first] = node_trad.second;
#ifdef DEBUG
        std::cout << "Node " << node_trad.first << " is in the SPEC-equiv of " << orbits[node_trad.first] << std::endl;
#endif
    }

    nb_local_nodes_after_pruning = equiv_nodes_representants.size();  //
    std::cout << "Finished to compute the equivalent nodes" << std::endl;

    return {map_in, map_out};
}

void export_to_dot(std::vector<Graph> graphs, std::string filename, std::map<int, int> partition) {
    std::ofstream file;
    // Open or create the file
    file.open(filename);
    if (file.is_open()) {
        file << "digraph G {" << std::endl;
        for (int i = 0; i < graphs.size(); i++) {
            file << "   "
                 << "subgraph cluster_" << i << " {" << std::endl;
            file << "       "
                 << "label = \"DC " << i << "\";" << std::endl;
            for (auto node : graphs[i].nodes) {
                if (partition[node] == i)
                    file << "       " << node << ";" << std::endl;
            }
            for (auto arc : graphs[i].arcs) {
                if (partition[arc.first] == partition[arc.second])
                    file << "        " << arc.first << " -> " << arc.second << ";" << std::endl;
            }
            file << "    "
                 << "}" << std::endl;
        }
        for (int i = 0; i < graphs.size(); i++) {
            for (auto arc : graphs[i].arcs) {
                if (partition[arc.first] != partition[arc.second] && partition[arc.first] == i)
                    file << arc.first << " -> " << arc.second << ";" << std::endl;
            }
        }
        file << "}" << std::endl;
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
        exit(EXIT_FAILURE);
    }
}

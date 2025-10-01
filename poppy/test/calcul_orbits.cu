#include <iostream>
#include <vector>
#include "../../nauty/nauty.h"
#include "../src/utils/include/utils.cuh"
#include "openfhe.h"

using namespace std::chrono;
using std::string;
using namespace lbcrypto;

void calcul_orbits() {
    // Creation of the graph G
    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);

    static DEFAULTOPTIONS_GRAPH(options);

    statsblk stats;
    int n, m;

    // We set the options of the graph.
    options.writeautoms = TRUE;
    options.writemarkers = TRUE;
    options.defaultptn = TRUE;  // No coloration
    options.getcanon = FALSE;
    options.digraph = TRUE;
    options.maxinvarlevel = 999;

    std::string filename = choose_graph();

    Graph graph_parsed = read_graph(GRAPH_DIR + filename, string(1, '\t'), true);

    n = graph_parsed.nodes.size();
    m = SETWORDSNEEDED(n);

    std::map<int, int> trad_tab;
    int index = 0;
    for (auto node : graph_parsed.nodes) {
        trad_tab[node] = index;
        index++;
    }

    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    DYNALLOC2(graph, g, g_sz, m, n, "malloc");
    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

    EMPTYGRAPH(g, m, n);
    for (auto edge : graph_parsed.arcs) {
        ADDONEARC(g, trad_tab[edge.first], trad_tab[edge.second], m);
    }

    std::cout << "Compute orbits for graph of size " << graph_parsed.nodes.size() << " nodes and "
              << graph_parsed.arcs.size() << " edges" << std::endl;

    auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, nullptr);

    std::cout << "\nTime for finding orbits when n = " << n << " : "
              << duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin << " ms.\n"
              << std::endl;
}
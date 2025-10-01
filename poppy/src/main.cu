#include <argp.h>
#include <iostream>
#include "calcul_PR.cuh"
#include "poppy_uj.cuh"
#include "utils/include/utils_stats.cuh"

using namespace std;

#ifdef EXPE
#include "../src/utils/include/utils_export.cuh"
std::string folder_expe = "expe-" + currentDateTimeExpe();
int number_encoding = 0;
int number_encryption = 0;
int number_mult = 0;
int number_add_scalar = 0;
int number_add = 0;
int number_rot = 0;
int number_mult_scalar = 0;
int number_mult_plain = 0;
#endif

struct arguments {
    int nb_iteration;
    std::string filename;
    std::string reference_file;
    std::string delimitor;
    int selection;
    bool test_uj;
    bool test_m;
    bool test_naive;
};

static error_t parse_opt(int key, char* arg, struct argp_state* state);
void launch_selection(int selection, arguments arguments);

static void printGroupedNodes(const std::map<int, int>& equivalent_map) {
    // Group nodes with the same image, excluding nodes that map to themselves
    std::map<int, std::vector<int>> grouped_nodes;
    for (const auto& pair : equivalent_map) {
        if (pair.first != pair.second) {
            grouped_nodes[pair.second].push_back(pair.first);
        }
    }

    // print the size of grouped_nodes
    std::cout << "Size of grouped_nodes: " << grouped_nodes.size() << std::endl;

    // Print groups of equivalent nodes
    int group_num = 1;
    for (const auto& group : grouped_nodes) {
        std::cout << "groupe[" << group_num << "] = (";
        std::cout << group.first << ", ";
        for (size_t i = 0; i < group.second.size(); ++i) {
            std::cout << group.second[i];
            if (i < group.second.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ")" << std::endl;
        group_num++;
    }

    // Calculate and print the total number of nodes in all groups
    int total_nodes = 0;
    for (const auto& group : grouped_nodes) {
        // Add 1 for the key node plus size of the equivalent nodes vector
        total_nodes += 1 + group.second.size();
    }
    std::cout << "Total number of nodes in all groups: " << total_nodes << std::endl;
}

void compute_orbits(std::string filename) {
    const string folder_serial = "../keys_graph_185";
    std::cout << "Compute orbits" << std::endl;
    POPPYuj data_center;
    if (filename == "") {
        filename = choose_graph();
    }
    data_center.init(GRAPH_DIR + filename, " ", true);

    auto orbits_map = data_center.get_orbits();
    printGroupedNodes(orbits_map);
}

int main(int argc, char* argv[]) {
    arguments arguments;
    arguments.nb_iteration = 1;
    arguments.filename = "";
    arguments.reference_file = "";
    arguments.delimitor = "";
    arguments.selection = 0;
    arguments.test_uj = true;
    arguments.test_m = true;
    arguments.test_naive = true;

    // Parse positional arguments with argp and add a --help option
    static struct argp_option options[] = {
        {"iteration", 'i', "ITERATION", 0, "Number of iterations to perform", 0},
        {"filename", 'f', "FILENAME", 0, "File to test", 1},
        {"reference_file", 'r', "REFFILE", 0, "Reference file", 2},
        {"delimitor", 'd', "DELIMITOR", 0, "Delimitor used in reference file", 3},
        {"selection", 's', "SELECTION", 0, "Selection of the version to run", 4},
        {"prog", 'p', "PROG", 0, "Program to run", 5},
        {"help", 'h', nullptr, 0, "Print usage information", 0},
        {nullptr}};

    // Call argp_parse to parse the arguments
    struct argp argp = {options, parse_opt, nullptr, nullptr, nullptr, nullptr, nullptr};
    argp_parse(&argp, argc, argv, 0, nullptr, &arguments);

    if (arguments.selection != 0) {
        launch_selection(arguments.selection, arguments);
        return 0;
    } else {
        while (true) {
            cout << "+--------------------------------------------------------------+" << endl;
            cout << "|                 PageRank calculation using                   |" << endl;
            cout << "|                  homomorphique encryption.                   |" << endl;
            cout << "+--------------------------------------------------------------+" << endl;
            cout << "| Versions                        | Source Files               |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 1. Naive PageRank               | naive.cpp                  |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 2. POPPYuj PageRank             | poppy_uj.cpp               |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 3. POPPYm PageRank              | poppy_m.cpp                |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 4. Clear PageRank               | clear_PR.cpp               |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 5. Unit Test On results         | test.cpp                   |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 6. Operation time execution     | op_compar.cpp              |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 7. Multi DC POPPYm              | multi_DC_poppy_m.cpp       |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 8. Multi DC POPPYuj             | multi_DC_poppy_uj.cpp      |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 9. Multi DC Naive               | multi_DC_naive.cpp         |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 10. Running expes               | expes.cpp                  |" << endl;
            cout << "+---------------------------------+----------------------------+" << endl;
            cout << "| 11. Show Multi DC graph stats                                |" << endl;
            cout << "+--------------------------------------------------------------+" << endl;
            cout << "| 12. Compute orbits                                           |" << endl;
            cout << "+--------------------------------------------------------------+" << endl;

            int selection = 0;
            bool valid = true;
            do {
                cout << endl << "> Run example (1 ~ 12) or exit (0) : ";
                if (!(cin >> selection)) {
                    valid = false;
                } else if (selection < 0 || selection > 12) {
                    valid = false;
                } else {
                    valid = true;
                }
                if (!valid) {
                    cout << "  [Beep~~] valid option: type 0 ~ 12" << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
            } while (!valid);
            launch_selection(selection, arguments);
        }
    }

    return 0;
}

/**
 * @brief Parser function for argp
 *
 * @param key
 * @param arg
 * @param state
 * @return error_t
 */
static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    auto args = (arguments*)(state->input);

    switch (key) {
        case 'i':
            args->nb_iteration = stoi(arg);
            break;
        case 'f':
            args->filename = arg;
            break;
        case 'r':
            args->reference_file = arg;
            break;
        case 'd':
            args->delimitor = arg;
            break;
        case 's':
            args->selection = stoi(arg);
            break;
        case 'h':
            argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
            break;
        case 'p':
            if (stoi(arg) == 1) {  // run POPPYuj
                args->test_m = false;
                args->test_naive = false;
            } else if (stoi(arg) == 2) {  // run POPPYm
                args->test_uj = false;
                args->test_naive = false;
            } else if (stoi(arg) == 3) {  // run Naive
                args->test_uj = false;
                args->test_m = false;
            } else if (stoi(arg) == 12) {  // run POPPYuj + POPPYm
                args->test_naive = false;
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

void launch_selection(int selection, arguments arguments) {
    switch (selection) {
        case 1:
            run_naive(arguments.nb_iteration);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
            break;

        case 2:
            run_poppy_uj(arguments.nb_iteration);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
            break;

        case 3:
            run_poppy_m(arguments.nb_iteration);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
            break;

        case 4:
            if (arguments.filename.empty()) {
                cout << "Please provide a filename with the --filename option" << endl;
            } else if (arguments.delimitor.empty()) {
                cout << "Please provide a delimitor with the --delimitor option" << endl;
            } else {
                auto results = clear_PR(arguments.nb_iteration, arguments.filename, arguments.delimitor);
                for (const auto& [node, pr] : results) {
                    std::cout << "Node " << node << ": " << pr << std::endl;
                }
            }
            break;

        case 5: {
            if (arguments.filename.empty()) {
                cout << "Please provide a filename with the --filename option" << endl;
            } else if (arguments.reference_file.empty()) {
                cout << "Please provide a reference file with the --reference_file option" << endl;
            } else if (arguments.delimitor.empty()) {
                cout << "Please provide a delimitor with the --delimitor option" << endl;
            } else {
                test(
                    arguments.nb_iteration,
                    arguments.filename,
                    arguments.reference_file,
                    arguments.delimitor,
                    arguments.test_uj,
                    arguments.test_m,
                    arguments.test_naive);
            }
        } break;

        case 6:
            operations_comparison();
            break;

        case 7: {
            if (!arguments.filename.empty()) {
                run_multi_dc_poppy_m(arguments.nb_iteration, arguments.filename);
            } else {
                run_multi_dc_poppy_m(arguments.nb_iteration);
            }
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        } break;
        case 8: {
            if (!arguments.filename.empty()) {
                run_multi_dc_poppy_uj(arguments.nb_iteration, arguments.filename);
            } else {
                run_multi_dc_poppy_uj(arguments.nb_iteration);
            }
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        } break;
        case 9: {
            if (!arguments.filename.empty()) {
                run_multi_dc_naive(arguments.nb_iteration, arguments.filename);
            } else {
                run_multi_dc_naive(arguments.nb_iteration);
            }
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        } break;

        case 10: {
            if (arguments.filename.empty()) {
                cout << "Please provide a filename with the --filename option" << endl;
            } else if (arguments.reference_file.empty()) {
                cout << "Please provide a reference file with the --reference_file option" << endl;
            } else if (arguments.delimitor.empty()) {
                cout << "Please provide a delimitor with the --delimitor option" << endl;
            } else {
                run_expes(
                    arguments.nb_iteration,
                    arguments.filename,
                    arguments.reference_file,
                    arguments.delimitor,
                    arguments.test_uj,
                    arguments.test_m,
                    arguments.test_naive);
            }
        } break;

        case 11:
            if (!arguments.filename.empty()) {
                multi_dc_graph_stats(arguments.filename);
            } else {
                cout << "Please provide a filename with the --filename option" << endl;
            }
            break;
        case 12:
            compute_orbits(arguments.filename);
            break;

        case 0:
            exit(0);
    }
}

#include "calcul_PR.cuh"
#ifdef EXPE
#include "../src/utils/include/utils_export.cuh"
#endif

void test(int max_iter, std::string to_test_filename, std::string reference_file, std::string delimitor, bool test_uj, bool test_m, bool test_naive) {
    std::vector<std::string> filenames = {to_test_filename};
    double threshold = 0.0001;

#ifdef EXPE
    std::string filename = to_test_filename;
    // SPlit to test file by "/" and get last element
    while (filename.find("/") != std::string::npos) {
        filename = filename.substr(filename.find("/") + 1);
    }

    expe_poppy_uj.saveToFile(folder_expe + "/expe_poppy_uj-" + filename + ".csv");
    expe_poppy_m.saveToFile(folder_expe + "/expe_poppy_m-" + filename + ".csv");
    expe_naive.saveToFile(folder_expe + "/expe_naive-" + filename + ".csv");
    // Duplicate stdout to file

    std::streambuf* cout_buf = std::cout.rdbuf();
    std::streambuf* cerr_buf = std::cerr.rdbuf();

    std::ofstream outf(folder_expe + "/log-" + filename + ".txt");
    std::ofstream errf(folder_expe + "/err-" + filename + ".txt");
    std::cout.rdbuf(outf.rdbuf());
    std::cerr.rdbuf(errf.rdbuf());
#endif

    // Running POPPY algorithm for 1 iteration until max_iter and compare result with clear PR
    for (int nb_iteration = 1; nb_iteration <= max_iter; nb_iteration++) {
#ifdef EXPE
        expe_poppy_m.setCurrentRun(nb_iteration);
        expe_poppy_uj.setCurrentRun(nb_iteration);
        expe_naive.setCurrentRun(nb_iteration);
#endif
        std::cerr << "Test iteration " << nb_iteration << std::endl;
        std::map<int, double> ref_results;
        ref_results = clear_PR(nb_iteration, reference_file, delimitor);
        for (auto filename : filenames) {
            std::cerr << "Test on " << filename << std::endl;
            std::map<int, double> results1;
            std::map<int, double> results2;
            std::map<int, double> results3;
            if (test_m) {
                results1 = run_multi_dc_poppy_m(nb_iteration, filename);
                CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
                for (auto ref_result : ref_results) {
                    if (results1[ref_result.first] > ref_result.second + threshold ||
                        results1[ref_result.first] < ref_result.second - threshold) {
                        std::cerr << "Error in poppym_version for iteration " << nb_iteration << " on " << filename
                                  << std::endl;
                        // Show diff
                        std::cerr << "Reference Value (" << ref_result.first << ") : " << ref_result.second
                                  << std::endl;
                        std::cerr << "Multi DC Coeus : " << results1[ref_result.first] << std::endl;
                    }
                }
            }
            if (test_uj) {
                results2 = run_multi_dc_poppy_uj(nb_iteration, filename);
                CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
                for (auto ref_result : ref_results) {
                    if (results2[ref_result.first] > ref_result.second + threshold ||
                        results2[ref_result.first] < ref_result.second - threshold) {
                        std::cerr << "Error in uj_version for iteration " << nb_iteration << " on " << filename
                                  << std::endl;
                        // Show diff
                        std::cerr << "Reference Value (" << ref_result.first << ") : " << ref_result.second
                                  << std::endl;
                        std::cerr << "Multi DC Uj : " << results2[ref_result.first] << std::endl;
                    }
                }
            }
            if (test_naive) {
                results3 = run_multi_dc_naive(nb_iteration, filename);
                CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
                for (auto ref_result : ref_results) {
                    if (results3[ref_result.first] > ref_result.second + threshold ||
                        results3[ref_result.first] < ref_result.second - threshold) {
                        std::cerr << "Error in naive_version for iteration " << nb_iteration << " on " << filename
                                  << std::endl;
                        // Show diff
                        std::cerr << "Reference Value (" << ref_result.first << ") : " << ref_result.second
                                  << std::endl;
                        std::cerr << "Multi DC naive : " << results3[ref_result.first] << std::endl;
                    }
                }
            }
        }
    }
#ifdef EXPE
    std::cout.rdbuf(cout_buf);
    std::cerr.rdbuf(cerr_buf);
#endif
}
#include "calcul_PR.cuh"
#ifdef EXPE
#include "../src/utils/include/utils_export.cuh"
#endif

void run_expes(
    int max_iter,
    std::string to_test_filename,
    std::string reference_file,
    std::string delimitor,
    bool test_uj,
    bool test_m,
    bool test_naive) {
    std::vector<std::string> filenames = {to_test_filename};

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
    std::ofstream errf_m(folder_expe + "/err-poppy_m-" + filename + ".txt");
    std::ofstream errf_uj(folder_expe + "/err-poppy_uj-" + filename + ".txt");
    std::ofstream errf_naive(folder_expe + "/err-naive-" + filename + ".txt");

    std::cout.rdbuf(outf.rdbuf());
#endif
#ifdef EXPE
    errf_m << "Test iteration " << max_iter << std::endl;
    errf_uj << "Test iteration " << max_iter << std::endl;
    errf_naive << "Test iteration " << max_iter << std::endl;
#endif
    std::map<int, double> ref_results;
    ref_results = clear_PR(max_iter, reference_file, delimitor);
    for (auto filename : filenames) {
#ifdef EXPE
        errf_m << "Test on " << filename << std::endl;
        errf_uj << "Test on " << filename << std::endl;
        errf_naive << "Test on " << filename << std::endl;
#endif
        std::map<int, double> results1;
        std::map<int, double> results2;
        std::map<int, double> results3;
        if (test_m) {
            results1 = run_multi_dc_poppy_m(max_iter, filename);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
#ifdef EXPE
            for (auto ref_result : ref_results) {
                errf_m << "Reference Value (" << ref_result.first << ") : " << ref_result.second << std::endl;
                errf_m << "Multi DC Coeus : " << results1[ref_result.first] << std::endl;
            }
#endif
        }
        if (test_uj) {
            results2 = run_multi_dc_poppy_uj(max_iter, filename);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
#ifdef EXPE
            for (auto ref_result : ref_results) {
                errf_uj << "Reference Value (" << ref_result.first << ") : " << ref_result.second << std::endl;
                errf_uj << "Multi DC Uj : " << results2[ref_result.first] << std::endl;
            }
#endif
        }
        if (test_naive) {
            results3 = run_multi_dc_naive(max_iter, filename);
            CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
#ifdef EXPE
            for (auto ref_result : ref_results) {
                errf_naive << "Reference Value (" << ref_result.first << ") : " << ref_result.second << std::endl;
                errf_naive << "Multi DC Naive : " << results3[ref_result.first] << std::endl;
            }
#endif
        }
    }
#ifdef EXPE
    std::cout.rdbuf(cout_buf);
    std::cerr.rdbuf(cerr_buf);

    // Close error files
    errf_m.close();
    errf_uj.close();
    errf_naive.close();
#endif
}
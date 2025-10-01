#include <poppy.cuh>
#include "openfhe.h"

std::map<int, double> clear_PR(int nb_iteration, std::string filename, std::string delimitor);
void run_naive(int nb_iteration);
void run_poppy_uj(int nb_iteration);
void run_poppy_m(int nb_iteration);
void run_poppy_m(int nb_iteration, std::string to_test_filename, std::string delimitor);
void run_multi_dc_poppy_uj(int nb_iteration);
void run_multi_dc_poppy_m(int nb_iteration);
void run_multi_dc_naive(int nb_iteration);
std::map<int, double> run_multi_dc_poppy_uj(int nb_iteration, std::string to_test_filename);
std::map<int, double> run_multi_dc_poppy_m(int nb_iteration, std::string to_test_filename);
std::map<int, double> run_multi_dc_naive(int nb_iteration, std::string to_test_filename);
void test(int max_iter, std::string to_test_filename, std::string reference_file, std::string delimitor, bool test_uj, bool test_m, bool test_naive);
void run_expes(int max_iter, std::string to_test_filename, std::string reference_file, std::string delimitor, bool test_uj, bool test_m, bool test_naive);
void calcul_orbits();
void operations_comparison();

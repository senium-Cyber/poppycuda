#pragma once

#include <vector>
#include "../../nauty/nausparse.h"
#include "openfhe.h"

using std::string;
using namespace lbcrypto;

/* Return the set of initialised ciphers.
   Each cipher contains 1/n in each cells where n is the number of nodes in the graph. */
std::map<int, Ciphertext<DCRTPoly>> initialisation_naif_PR(int n, int DC, KeyPair<DCRTPoly> keys, CryptoContext<DCRTPoly> cc);

/* For each nodes with an edges from DC_i to DC_j, add is PR divided by deg_out in inter_DC_PRs[DC_j]. */
void inter_DC_com(std::vector<std::map<int, Ciphertext<DCRTPoly>>>& inter_DC_PRs, std::map<int, Ciphertext<DCRTPoly>> PRs, int n, int DC, std::multimap<int, int> map_out, std::vector<int> partition, CryptoContext<DCRTPoly> cc);

/* Performe an iteration of PR algorithm using last PR vectors PRs. */
// std::map<int, Ciphertext<DCRTPoly>> naif_PR_calcul(std::map<int, Ciphertext<DCRTPoly>> PRs, std::multimap<int, int> map_in,
// std::multimap<int, int> map_out, int n, KeyPair<DCRTPoly> keys, CryptoContext<DCRTPoly> cc);


std::map<int, Ciphertext<DCRTPoly>> multi_DC_naif_PR(std::map<int, Ciphertext<DCRTPoly>> PRs, std::map<int, Ciphertext<DCRTPoly>> inter_DC_PRs, std::multimap<int, int> map_in, std::multimap<int, int> map_out, int n, std::vector<int> partition, int DC, KeyPair<DCRTPoly> keys, CryptoContext<DCRTPoly> cc);

/* Return the set of initialised ciphers in the multi DC setting.
   Each cipher contains 1/n in each cells where n is the number of nodes in the graph.*/
std::map<int, Ciphertext<DCRTPoly>> initialisation_multi_DC_naif_PR(int n, int DC, std::vector<int> partition, KeyPair<DCRTPoly> keys, CryptoContext<DCRTPoly> cc);

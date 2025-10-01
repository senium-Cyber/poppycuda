#pragma once

#include <vector>
#include "../../nauty/nauty.h"
#include "../include/poppy_uj.cuh"
#include "utils.cuh"

using std::string;
using namespace lbcrypto;

/* Function to calculate the size needed by the u_js.
   This size is equal to the number of orbits - number of nodes without in edges - number of nodes without out edges.
   We set n_deg0 = number of nodes without out edges before returning size of u_js.
*/
int set_u_j_size(int n, int* orbits, std::multimap<int, int> map_in, std::multimap<int, int> map_out);
int set_u_j_size(
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<Node /*node*/, int /*datacenter*/> partition,
    int datacenter_id,
    int& nb_constant_nodes,
    int nb_dangling_nodes);
/*  Function to creat the map_u_j.
    Parameters : - the empty map_u_j we will fill,
                 - n : the number of vertices of the graph,
                 - orbits the orbits of the graph given by nauty,
                 - map_in and map_out : maps for the in and out vertics avec the graph.
    Return : The number of u_j we will need.
*/
int creation_map_u_j(
    std::multimap<int, triplet>& map_u_j,
    int n,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out);
int creation_map_u_j(
    std::multimap<int, triplet>& map_u_j,
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> partition,
    int datacenter_id,
    std::map<int /*orbit*/, int /*repr*/> representative_orbits,
    std::set<Node> set_dangling_nodes);

/** For reconstruction uj give map_uj for the two first maps.
 *  For construction map deg0 give map_uj and map_deg0 for the second map.
 */
void creation_map_reconstruction_u(
    std::multimap<int, triplet> map_uj,
    std::multimap<int, triplet> map_deg0,
    std::vector<std::multimap<int, quadruplet>>& map_reconstruction_uj,
    int s,
    int m,
    std::map<int, int> representative_orbits,
    std::map<Node /*node*/, int /*data center*/> partition,
    int datacenter_id);

/*  Function to creat the map of u_j for the calcul of PR of nodes with deg_out = 0.
    Parameters : - the empty map_deg_out_0 we will fill,
                 - n : the number of vertices of the graph,
                 - orbits : the orbits of the graph given by nauty,
                 - map_in and map_out : maps for the in and out vertics avec the graph.
    Return : The number of u_j we will need.
*/
int creation_map_deg_out_0(
    std::multimap<int, triplet>& map_deg_out_0,
    int n,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out);
int creation_map_deg_out_0(
    std::multimap<int, triplet>& map_deg_out_0,
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> partition,
    int datacenter_id,
    std::map<int /*orbit*/, int /*repr*/> representative_orbits,
    std::set<Node> set_dangling_nodes);

/*  Function to initialise the u_j, this is the iteration 0 of the PR calculation.
    Parameters : - n : the number of vertices of the graph,
                 - nb_u_j the number of u_j nwe need,
                 - map_u_j the map for the creation of the u_j,
    Return : A vector of size the number of u_j. Each cells represent an encrypted u_j.
*/
std::vector<Ciphertext<DCRTPoly>> initialisation_u_j(
    int n,
    int n_r,
    int nb_u_j,
    std::multimap<int, triplet> map_u_j,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc);

std::vector<Ciphertext<DCRTPoly>> construction_u_deg0(
    Ciphertext<DCRTPoly> PR,
    int n_r,
    int nb_u_j,
    std::multimap<int, triplet> map_deg_out_0,
    std::multimap<int, triplet> map_u_j,
    std::multimap<int, doublet> PR_cst,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc);

std::vector<Ciphertext<DCRTPoly>> construction_uj(
    std::vector<Ciphertext<DCRTPoly>> PRs,
    int m,
    int s,
    std::vector<std::multimap<int, quadruplet>> map_reconstruction_uj,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc,
    std::map<int /*orbit*/, int /*repr*/> representative_orbits,
    std::map<int, int> partition,
    int datacenter_id);

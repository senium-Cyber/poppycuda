#include "include/utils_poppy_uj.cuh"

#include <fstream>
#include <iostream>
#include <vector>
#ifdef EXPE
#include "include/utils_export.cuh"
#endif
using namespace lbcrypto;

// For a unique DC
int set_u_j_size(int n, int* orbits, std::multimap<int, int> map_in, std::multimap<int, int> map_out) {
    int cpt_in = 0;
    int cpt_out = 0;

    int orbs = 0;
    for (int i = 0; i < n; i++) {
        if (orbits[i] == i) {  // orbit
            orbs++;
            if (map_in.count(i) == 0)  // constancy
                cpt_in++;
            if (map_out.count(i) == 0)  // inutility
                cpt_out++;
        }
    }

    int n_deg0 = cpt_out;
    int size_uj = orbs - cpt_in - cpt_out;

    if (n_deg0 < size_uj)
        return size_uj;
    else
        return n_deg0;
}

// For multi DC
int set_u_j_size(
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<Node /*node*/, int /*datacenter*/> partition,
    int datacenter_id,
    int& nb_constant_nodes,
    int nb_dangling_nodes) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
    std::set<int> orbs;
    for (auto orbit : orbits) {
        if (partition[orbit.first] == datacenter_id) {
            if (orbs.count(orbit.second) == 0) {     // new orbit
                if (map_in.count(orbit.first) == 0)  // Constant node
                    nb_constant_nodes++;
            }
            orbs.insert(orbit.second);
        }
    }
    if (nb_constant_nodes > 0)
        return orbs.size() - 1 - nb_dangling_nodes; // All constant nodes are in the same orbit
    else
        return orbs.size() - nb_dangling_nodes;

#else
    int size_uj = 0;
    for (auto orbit : orbits) {
        if (partition[orbit.first] == datacenter_id) {
            if (map_in.count(orbit.first) != 0)  // Not a constant node
                size_uj++;
            else
                nb_constant_nodes++;
        }
    }
    return size_uj;
#endif
}

// For a unique DC
int creation_map_u_j(
    std::multimap<int, triplet>& map_u_j,
    int n,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out) {
    int nb_u_j = 0;
    unsigned int indice = 0;

    for (int i = 0; i < n; i++) {
        if (map_u_j.count(orbits[i]) == 0) {
            if (map_in.count(i) != 0 && map_out.count(i) != 0) {  // If deg_in and deg_out of vertex i is not null
                std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> values;
                values = map_in.equal_range(i);  // To get all the values of the key i

                int calcul_nb_u_j = 0;

                for (auto it = values.first; it != values.second; ++it) {
                    triplet t = {indice, orbits[it->second], (int)map_out.count(it->second)};
                    map_u_j.insert(std::pair<int, triplet>(i, t));
                    calcul_nb_u_j++;
                }
                indice++;

                if (calcul_nb_u_j > nb_u_j)
                    nb_u_j = calcul_nb_u_j;
            }
        }
    }
    return nb_u_j;
}

// For multi DC
int creation_map_u_j(
    std::multimap<int, triplet>& map_u_j,
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> partition,
    int datacenter_id,
    std::map<int /*orbit*/, int /*repr*/> representative_orbits,
    std::set<Node> set_dangling_nodes) {
    // The max number of in edges in the graph (among the effective nodes)
    int nb_u_j = 0;
    unsigned int indice = 0;

    for (auto orbit : orbits) {
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
        if (map_u_j.count(representative_orbits.at(orbit.second)) == 0) {  // Do once by orbit
            if (set_dangling_nodes.find(orbit.second) == set_dangling_nodes.end()) {
                // If the orbit is not a dangling node (effective)
#endif
                if (map_in.count(orbit.first) != 0) {
                    // If the orbit is not a constant node
                    if (partition.at(orbit.first) == datacenter_id) {   // If the node belongs to the datacenter
                        auto values = map_in.equal_range(orbit.first);  // To get all nodes going to orbit

                        int calcul_nb_u_j = 0;

                        for (auto it = values.first; it != values.second; ++it) {
                            triplet t;
                            if (partition.at(it->second) == datacenter_id) {
                                t = {indice, orbits.at(it->second), (int)map_out.count(it->second)};

                                if (map_in.count(it->second) >
                                    0) {  // Constant nodes are added in u_const and not in u_j
                                    calcul_nb_u_j++;
                                }

                            } else
                                t = {
                                    indice,
                                    orbits.at(it->second),
                                    1};  // We don't know the real value of deg_out of distant nodes
                            map_u_j.insert({representative_orbits.at(orbit.second), t});
                        }
                        indice++;

                        if (calcul_nb_u_j > nb_u_j)
                            nb_u_j = calcul_nb_u_j;
                    }
                }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
            }
        }
#endif
    }
    return nb_u_j;
}

// For multi DC
void creation_map_reconstruction_u(
    std::multimap<int, triplet> map_uj,
    std::multimap<int, triplet> map_deg0,
    std::vector<std::multimap<int, quadruplet>>& map_reconstruction_uj,
    int s,
    int m,
    std::map<int, int> representative_orbits,
    std::map<Node /*node*/, int /*data center*/> partition,
    int datacenter_id) {
    // Initialize map_reconstruction_uj
    map_reconstruction_uj = std::vector<std::multimap<int, quadruplet>>(m / s);
    auto it = map_deg0.begin();
    auto it2 = map_deg0.begin();
    while (it != map_deg0.end()) {
        while (it2 != map_deg0.end() && it2->first == it->first) {
            triplet value = it2->second;
            Node node = it2->first;

            if (partition.at(representative_orbits.at(value.vertex)) == datacenter_id) {
                // Node is from our datacenter
                if (map_uj.count(representative_orbits.at(value.vertex)) != 0) {
                    // Node is not constant
                    int i = map_uj.find(representative_orbits.at(value.vertex))->second.indice;
                    int rot_value = (i + s - (value.indice % s)) % s;
                    quadruplet infos = {node, value.indice, value.vertex, value.deg_out};
                    map_reconstruction_uj[value.indice / s].insert({rot_value, infos});
                }
            }
            it2++;
        }
        it = it2;
    }
}

// For a unique DC
int creation_map_deg_out_0(
    std::multimap<int, triplet>& map_deg_out_0,
    int n,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out) {
    int nb_u_deg0 = 0;
    unsigned int indice = 0;

    for (int i = 0; i < n; i++) {
        if (map_deg_out_0.count(orbits[i]) == 0) {
            if (map_out.count(i) == 0) {
                std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> values;
                values = map_in.equal_range(i);  // To get all the values of the key i

                int calcul_nb_u_deg0 = 0;

                for (auto it = values.first; it != values.second; ++it) {
                    triplet t = {indice, orbits[it->second], (int)map_out.count(it->second)};
                    map_deg_out_0.insert(std::pair<int, triplet>(i, t));
                    calcul_nb_u_deg0++;
                }
                indice++;

                if (calcul_nb_u_deg0 > nb_u_deg0)
                    nb_u_deg0 = calcul_nb_u_deg0;
            }
        }
    }
    return nb_u_deg0;
}

// For multi DC
int creation_map_deg_out_0(
    std::multimap<int, triplet>& map_deg_out_0,
    std::map<int, int> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> partition,
    int datacenter_id,
    std::map<int /*orbit*/, int /*repr*/> representative_orbits,
    std::set<Node> set_dangling_nodes) {
    // The max number of in edges in the graph (among the effective nodes)
    int nb_u_j = 0;
    unsigned int indice = 0;

    for (auto orbit : orbits) {
        if (set_dangling_nodes.find(orbit.second) != set_dangling_nodes.end()) {  // Orbit is a dangling node
            if (map_in.count(orbit.first) != 0) {
                // If the node is not a constant node
                if (partition.at(orbit.first) == datacenter_id) {   // If the node belongs to the datacenter
                    auto values = map_in.equal_range(orbit.first);  // To get all nodes going to orbit

                    int calcul_nb_u_j = 0;

                    for (auto it = values.first; it != values.second; ++it) {
                        triplet t;
                        if (partition.at(it->second) == datacenter_id) {
                            t = {indice, orbits.at(it->second), (int)map_out.count(it->second)};

                            if (map_in.count(it->second) >
                                0) {  // Constant nodes are not added in u_deg0
                                calcul_nb_u_j++;
                            }

                        } else
                            t = {
                                indice,
                                orbits.at(it->second),
                                1};  // We don't know the real value of deg_out of distant nodes
                        map_deg_out_0.insert({representative_orbits.at(orbit.second), t});
                    }
                    indice++;

                    if (calcul_nb_u_j > nb_u_j)
                        nb_u_j = calcul_nb_u_j;
                }
            }
        }
    }
    return nb_u_j;
}

std::vector<Ciphertext<DCRTPoly>> initialisation_u_j(
    int n,
    int n_r,
    int nb_u_j,
    std::multimap<int, triplet> map_u_j,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    std::vector<double> vect(n_r);
    for (int i = 0; i < n_r; i++)
        vect[i] = 0.0;

    std::vector<std::vector<double>> iv(nb_u_j);
    for (int i = 0; i < nb_u_j; i++)
        iv[i] = vect;

    int index_u_j = 0;  // Index to know each u_j to use

    auto it = map_u_j.begin();
    auto it2 = map_u_j.begin();
    while (it != map_u_j.end()) {
        while (it2 != map_u_j.end() && it2->first == it->first) {
            triplet value = it2->second;
            iv[index_u_j][value.indice] = 1.0 / (n * value.deg_out);

            index_u_j++;
            it2++;
        }

        index_u_j = 0;
        it = it2;
    }

    std::vector<Ciphertext<DCRTPoly>> ciph_iv(nb_u_j);
    for (int i = 0; i < nb_u_j; i++) {
        Plaintext plain = cc->MakeCKKSPackedPlaintext(iv[i]);
        ciph_iv[i] = cc->Encrypt(keys.publicKey, plain);
#ifdef EXPE
        number_encryption++;
#endif
    }

    return ciph_iv;
}

// For a unique DC
std::vector<Plaintext> creation_vect_PR_const(
    int n_r,
    int n_deg0,
    std::multimap<int, triplet> map_u_j,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    int index = 0;
    std::vector<double> vect(n_r);
    for (int i = 0; i < n_r; i++)
        vect[i] = 0.0;

    std::vector<std::vector<double>> cst_PRs(n_deg0);
    for (int i = 0; i < n_deg0; i++)
        cst_PRs[i] = vect;

    for (auto it = map_u_j.begin(); it != map_u_j.end(); it++) {
        triplet value = it->second;
        if (map_u_j.count(value.vertex) == 0) {
            cst_PRs[index][value.indice] = 0.15 / value.deg_out;
            index++;
        }
    }

    std::vector<Plaintext> pl_cst_PR(n_deg0);
    for (int i = 0; i < n_deg0; i++) {
        pl_cst_PR[i] = cc->MakeCKKSPackedPlaintext(cst_PRs[i]);
    }

    return pl_cst_PR;
}

void reconstruction_u_j(
    std::vector<Ciphertext<DCRTPoly>>& u_j,
    Ciphertext<DCRTPoly> PR,
    int nb_u_j,
    std::multimap<int, triplet> map_u_j,
    std::vector<Plaintext> PR_cst,
    int n_r,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    // First we reset the u_j to zero ciphertext
    std::vector<double> zero(n_r, 0.0);
    Plaintext pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encryption++;
#endif

    for (int j = 0; j < nb_u_j; j++)
        u_j[j] = ct_zero->Clone();

    // Secondly we fill them with PRs values
    int index_u_j = 0;       // Index to know each u_j to use
    int index_PR_const = 0;  // Index to know each vector of PR_const to use

    auto it = map_u_j.begin();
    auto it2 = map_u_j.begin();
    while (it != map_u_j.end()) {
        while (it2 != map_u_j.end() && it2->first == it->first) {
            triplet value = it2->second;
            if (map_u_j.count(value.vertex) == 0) {
                cc->EvalAddInPlace(u_j[index_u_j], PR_cst[index_PR_const]);
#ifdef EXPE
                number_add++;
#endif
                index_PR_const++;
            } else {
                int i = map_u_j.find(value.vertex)->second.indice;
                int step_rot = value.indice - i;  // roration associé de PRs à uj
                Ciphertext<DCRTPoly> PR_inter;    // A ciphertext to put the modification without modifie PR

                // Mask creation
                std::vector<double> mask(n_r);
                for (int indice = 0; indice < n_r; indice++) {
                    if (indice == i)
                        mask[i] = 1.0 / value.deg_out;
                    else
                        mask[indice] = 0.0;
                }
                Plaintext pl_mask = cc->MakeCKKSPackedPlaintext(mask);

                if (step_rot == 0) {
                    PR_inter = cc->EvalMult(PR, pl_mask);
                    cc->EvalAddInPlace(u_j[index_u_j], PR_inter);
#ifdef EXPE
                    number_add++;
                    number_mult_plain++;
#endif
                } else {
                    PR_inter = cc->EvalMult(PR, pl_mask);
                    PR_inter = cc->EvalRotate(
                        PR_inter,
                        -step_rot);  // - step_rot > 0 : rotation to the left or - step_rot < 0 : rotation to
                                     // the right
                    cc->EvalAddInPlace(u_j[index_u_j], PR_inter);
#ifdef EXPE
                    number_add++;
                    number_mult_plain++;
                    number_rot++;
#endif
                }
            }

            index_u_j++;  // balader dans les values
            it2++;
        }

        index_u_j = 0;
        it = it2;
    }
}

// For a unique DC
std::vector<Ciphertext<DCRTPoly>> construction_u_deg0(
    Ciphertext<DCRTPoly> PR,
    int n_r,
    int nb_u_j,
    std::multimap<int, triplet> map_deg_out_0,
    std::multimap<int, triplet> map_uj,
    std::multimap<int, doublet> PR_cst,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    std::vector<double> zeros(n_r);
    Plaintext pl_zeros = cc->MakeCKKSPackedPlaintext(zeros);
    auto ct_zeros = cc->Encrypt(keys.publicKey, pl_zeros);
#ifdef EXPE
    number_encryption++;
#endif
    std::vector<Ciphertext<DCRTPoly>> u(nb_u_j);
    for (int i = 0; i < nb_u_j; i++)
        u[i] = ct_zeros;

    int j = 0;  // Index to know each u_j to use

    auto it = map_deg_out_0.begin();
    auto it2 = map_deg_out_0.begin();
    while (it != map_deg_out_0.end()) {
        while (it2 != map_deg_out_0.end() && it2->first == it->first) {
            triplet value = it2->second;

            if (map_uj.count(value.vertex) == 0) {  // case the node attached to the uj is a constant node
                Ciphertext<DCRTPoly> vect_cst;

                auto range = PR_cst.equal_range(value.vertex);
                for (auto it_range = range.first; it_range != range.second; ++it_range) {
                    if ((it_range->second).node == it->first) {
                        vect_cst = (it_range->second).cipher;
                        break;
                    }
                }

                cc->EvalAddInPlace(u[j], vect_cst);
#ifdef EXPE
                number_add++;
#endif
            } else {
                int i = map_uj.find(value.vertex)->second.indice;
                int step_rot = value.indice - i;
                Ciphertext<DCRTPoly> PR_inter;  // A ciphertext to put the modification without modifying PR

                // Mask creation
                std::vector<double> mask(n_r);
                for (int indice = 0; indice < n_r; indice++) {
                    if (indice == i)
                        mask[i] = 1.0 / value.deg_out;
                    else
                        mask[indice] = 0.0;
                }
                Plaintext pl_mask = cc->MakeCKKSPackedPlaintext(mask);

                if (step_rot == 0) {
                    PR_inter = cc->EvalMult(PR, pl_mask);
                    cc->EvalAddInPlace(u[j], PR_inter);
#ifdef EXPE
                    number_add++;
                    number_mult_plain++;
#endif
                } else {
                    PR_inter = cc->EvalMult(PR, pl_mask);
                    PR_inter = cc->EvalRotate(
                        PR_inter,
                        -step_rot);  // - step_rot > 0 : rotation to the left or - step_rot < 0 : rotation to
                                     // the right
                    cc->EvalAddInPlace(u[j], PR_inter);
#ifdef EXPE
                    number_add++;
                    number_mult_plain++;
                    number_rot++;
#endif
                }
            }

            j++;
            it2++;
        }

        j = 0;
        it = it2;
    }
    return u;
}

// For multi DC

std::vector<Ciphertext<DCRTPoly>> construction_uj(
    std::vector<Ciphertext<DCRTPoly>> PRs,
    int m,
    int s,
    std::vector<std::multimap<int, quadruplet>> map_reconstruction_uj,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc,
    std::map<int, int> representative_orbits,
    std::map<int, int> partition,
    int datacenter_id) {
    std::vector<double> zero(s, 0.0);
    Plaintext pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encryption++;
#endif

    std::vector<Ciphertext<DCRTPoly>> u_deg0(m / s);

    // To save memory, we fill each sub u_deg0 one by one.
    for (int i = 0; i < m / s; i++) {
        u_deg0[i] = ct_zero->Clone();  // Set each sub u_deg0 to a vector of zeros

        auto map_reconstruction = map_reconstruction_uj[i];
        std::shared_ptr<std::vector<DCRTPoly>> precomp = cc->EvalFastRotationPrecompute(PRs[i]);
        Ciphertext<DCRTPoly> rotated_PR;

        auto it = map_reconstruction.begin();
        auto it2 = map_reconstruction.begin();

        while (it != map_reconstruction.end()) {
            // for each rotation : save the rotated PRs
            int rot = it->first;

#ifdef COEUS
            // Check if rot is a power of 2
            if ((rot & (rot - 1)) == 0) {
                rotated_PR = cc->EvalFastRotation(PRs[i], rot, 2 * cc->GetRingDimension(), precomp);
#ifdef EXPE
                number_rot++;
#endif
            } else {
                int rotated = 0;
                // Find the largest power of 2 that is smaller than rot
                int k = 1;
                while (!(rot & k)) {
                    k = k << 1;
                }
                rotated = k;
                rotated_PR = cc->EvalFastRotation(PRs[i], k, 2 * cc->GetRingDimension(), precomp);
#ifdef EXPE
                number_rot++;
#endif

                // Compute the remaining rotation
                while (rotated != rot) {
                    // Find the largest power of 2 that is smaller than remaining rot
                    int k = 1;
                    while (!((rot - rotated) & k)) {
                        k = k << 1;
                    }
                    rotated += k;
                    // Rotate the PR
                    rotated_PR = cc->EvalRotate(rotated_PR, k);
#ifdef EXPE
                    number_rot++;
#endif
                }
            }
#else
            rotated_PR = cc->EvalFastRotation(PRs[i], rot, 2 * cc->GetRingDimension(), precomp);
#ifdef EXPE
            number_rot++;
#endif
#endif

            // Now we mask the rotated cipher
            // We creat a unique mask per sub_uj
            std::vector<double> mask_of_zeros(s, 0.0);
            std::vector<std::vector<double>> masks(m / s, mask_of_zeros);
            while (it2 != map_reconstruction.end() && it2->first == it->first) {
                auto quadruplet = it2->second;

                assert(quadruplet.deg_out != 0);
                masks[quadruplet.indice / s][quadruplet.indice % s] += 1.0 / quadruplet.deg_out;
                it2++;
            }
            for (int sub_uj = 0; sub_uj < m / s; sub_uj++) {
                Plaintext pl_mask = cc->MakeCKKSPackedPlaintext(masks[i]);
                auto PR_inter = cc->EvalMult(rotated_PR, pl_mask);
                cc->EvalAddInPlace(u_deg0[sub_uj], PR_inter);

#ifdef EXPE
                number_add++;
                number_mult_plain++;
#endif
            }
            it = it2;
        }
    }
    return u_deg0;
}
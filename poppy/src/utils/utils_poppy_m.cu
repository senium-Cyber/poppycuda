#include "include/utils_poppy_m.cuh"
#include <fstream>
#include <iostream>
#include <queue>
#include "include/utils.cuh"  //  print ciphers for debug

using namespace lbcrypto;
#ifdef EXPE
#include "include/utils_export.cuh"
#endif
/**
 * @brief Create the adjacency matrix such as (i, j) is equal to the number of "in
 * edges" of node i from the orbit j divided by the deg out of node i
 *
 * @param actives_nodes_vect Loop over the first size_vect elements of this
 * array
 * @param size_vect
 * @param orbits
 * @param map_in
 * @param map_out
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> create_adjacency_matrix(
    std::vector<int> actives_nodes_vect,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out) {
    std::vector<std::vector<double>> mask_matrix(actives_nodes_vect.size());

    for (int i = 0; i < actives_nodes_vect.size(); i++) {
        int indice_i = actives_nodes_vect[i];  // indice_i is the orbit of the node i
        mask_matrix[i] = std::vector<double>(actives_nodes_vect.size(), 0.0);

        auto nodes_in_i = map_in.equal_range(indice_i);

        for (auto it_range = nodes_in_i.first; it_range != nodes_in_i.second; ++it_range) {
            auto it = std::find(actives_nodes_vect.begin(), actives_nodes_vect.end(), orbits[it_range->second]);
            int j = it - actives_nodes_vect.begin();
            if (it == actives_nodes_vect.end()) {
                // We are in the case that we get a node connected to an active node
                continue;
            }
            mask_matrix[i][j] += 1.0 / (double)map_out.count(it_range->second);
        }
    }

    return mask_matrix;
}

/**
 * @brief Create the adjacency matrix such as (i, j) is equal to the number of "in
 * edges" of node i from the orbit j divided by the deg out of node i
 *
 * @param actives_nodes_vect
 * @param size_vect
 * @param orbits
 * @param map_in
 * @param map_out
 * @param representative_orbits
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> create_adjacency_matrix(
    std::vector<Orbit> actives_nodes_vect,
    std::map<int /*node*/, int /*orbit*/> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> representative_orbits,
    std::map<Node, int> partition,
    int datacenter_id) {
    std::vector<std::vector<double>> mask_matrix;

    for (Orbit orbit_active_node : actives_nodes_vect) {
        std::vector<double> mask_line(actives_nodes_vect.size(), 0.0);

        // If the active node is one of our node
        if (partition[representative_orbits[orbit_active_node]] == datacenter_id) {
            auto nodes_in_i = map_in.equal_range(representative_orbits[orbit_active_node]);

            for (auto edge_in = nodes_in_i.first; edge_in != nodes_in_i.second; ++edge_in) {
                auto it = std::find(actives_nodes_vect.begin(), actives_nodes_vect.end(), orbits[edge_in->second]);

                if (it == actives_nodes_vect.end()) {
                    continue;
                }
#if defined(COMPUTE_ORBITS) ^ defined(COMPUTE_EQUIVALENT_NODES)
                if (partition[edge_in->second] == datacenter_id && map_in.count(edge_in->second) == 0) {
                    // The node is our node and a constant node
                    continue;
                }
#endif
                int j = it - actives_nodes_vect.begin();
                if (partition[edge_in->second] == datacenter_id)
                    mask_line[j] += 1.0 / map_out.count(edge_in->second);
                else
                    mask_line[j] += 1.0;
            }
        }
        mask_matrix.push_back(mask_line);
    }

    return mask_matrix;
}

void print_matrix(std::vector<std::vector<double>> matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::vector<Diags> extract_and_encode_diags(
    CryptoContext<DCRTPoly> cc,
    std::vector<std::vector<double>> mask_matrix,
    int batchSize,
    int& m) {
    auto sub_matrices = prepareMatrix(mask_matrix, batchSize, m);

    std::vector<Diags> sub_diags((m / batchSize) * (m / batchSize));
    for (int l = 0; l < m / batchSize; l++) {
        for (int c = 0; c < m / batchSize; c++) {
            sub_diags[l * m / batchSize + c] = MakeCKKSPackedDiagonals(cc, sub_matrices[l * m / batchSize + c]);
        }
    }

    return sub_diags;
}

std::vector<std::vector<double>> padMatrixWithZeros(std::vector<std::vector<double>> matrix, int batchSize) {
    int n = matrix.size();
    // add zeros to each lines of the matrix
    for (int i = 0; i < n; i++) {
        for (int j = n; j < batchSize; j++) {
            matrix[i].push_back(0.0);
        }
    }

    // add lines of zeros to the matrix
    std::vector<double> line_zeros(batchSize, 0.0);
    for (int i = n; i < batchSize; i++) {
        matrix.push_back(line_zeros);
    }

    return matrix;
}

void padMatrixWithZerosInPlace(std::vector<std::vector<double>>& matrix, int batchSize) {
    int n = matrix.size();
    // add zeros to each lines of the matrix
    for (int i = 0; i < n; i++) {
        for (int j = n; j < batchSize; j++) {
            matrix[i].push_back(0.0);
        }
    }

    // add lines of zeros to the matrix
    std::vector<double> line_zeros(batchSize, 0.0);
    for (int i = n; i < batchSize; i++) {
        matrix.push_back(line_zeros);
    }
}

std::vector<double> padVectorWithZeros(std::vector<double> vector, int batchSize) {
    int n = vector.size();
    // add zeros at the end of the vector
    for (int i = n; i < batchSize; i++) {
        vector.push_back(0.0);
    }

    return vector;
}

Diags MakeCKKSPackedDiagonals(CryptoContext<DCRTPoly> cc, std::vector<std::vector<double>>& matrix) {
    Diags diags;

    int n = matrix.size();

    // extract the diagonals and encode them as Plaintexts
    std::vector<double> zero(n, 0.0);
    for (int i = 0; i < n; i++) {
        std::vector<double> diag;
        for (int j = 0; j < n; j++) {
            diag.push_back(matrix[j][(i + j) % n]);
        }

        if (diag == zero) {
            diags.push_back({nullptr, true});
        } else {
            Plaintext ptxt = cc->MakeCKKSPackedPlaintext(diag);
            diags.push_back({ptxt, false});
        }
    }

    return diags;
}

Ciphertext<DCRTPoly> matmul_a_la_Coeus(
    CryptoContext<DCRTPoly> cc,
    std::vector<Plaintext> diags,
    Ciphertext<DCRTPoly> ct1,
    std::vector<int> diags_null) {
    // Get the diagonals of the matrix as Plaintexts
    std::stack<Ciphertext<DCRTPoly>> st;
    Ciphertext<DCRTPoly> ct_result;
    Ciphertext<DCRTPoly> ct_temp;
    int cpt_diag_not_null = 0;

    for (size_t i = 0; i < diags_null.size(); i++) {
        Ciphertext<DCRTPoly> ct_rot;
        if (st.empty()) {
            ct_rot = cc->EvalRotate(ct1, i);
        } else {
            // find the largest power of 2 that is smaller than i
            int k = 1;
            while (!(i & k)) {
                k = k << 1;
            }
            // rotate the ciphertext
            Ciphertext<DCRTPoly> cached_ct = st.top();
            ct_rot = cc->EvalRotate(cached_ct, k);

            // pop the stack if i is a power of 2
            if (i & (k << 1)) {
                st.pop();
            }
        }

        // push the rotated ciphertext to the stack if i is even
        if (i && !(i % 2)) {
            st.push(ct_rot);
        }

        // multiply the rotated ciphertext with the diagonal
        if (diags_null[i] == 0) {
            if (cpt_diag_not_null == 0) {
                ct_result = cc->EvalMult(diags[cpt_diag_not_null], ct_rot);
            } else {
                ct_temp = cc->EvalMult(diags[cpt_diag_not_null], ct_rot);
                ct_result = cc->EvalAdd(ct_result, ct_temp);
            }

            cpt_diag_not_null += 1;
        }
    }

    return ct_result;
}

std::vector<Ciphertext<DCRTPoly>> matmul_coeus_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<Diags> subDiags,
    std::vector<Ciphertext<DCRTPoly>> subVectors) {
    int nbSubMatrices = subVectors.size();
    int s = subDiags[0].size();
    std::vector<double> zero(s, 0.0);
    auto pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
    std::vector<Ciphertext<DCRTPoly>> result(nbSubMatrices);
    for (int i = 0; i < nbSubMatrices; i++) {
        result[i] = ct_zero->Clone();
    }

    for (int c = 0; c < nbSubMatrices; c++) {
        std::stack<Ciphertext<DCRTPoly>> st;
        Ciphertext<DCRTPoly> RotatedVect = subVectors[c];
        auto precomp = cc->EvalFastRotationPrecompute(RotatedVect);

        for (int i = 0; i < s; i++) {
            // rotation
            Ciphertext<DCRTPoly> ct_rot;
            if (st.empty()) {
                ct_rot = cc->EvalFastRotation(RotatedVect, i, cc->GetRingDimension() * 2, precomp);
#ifdef EXPE
                number_rot++;
#endif
            } else {
                // find the largest power of 2 that is smaller than i
                int k = 1;
                while (!(i & k)) {
                    k = k << 1;
                }
                // rotate the ciphertext
                Ciphertext<DCRTPoly> cached_ct = st.top();
                ct_rot = cc->EvalRotate(cached_ct, k);
#ifdef EXPE
                number_rot++;
#endif

                // pop the stack if i is a power of 2
                if (i & (k << 1)) {
                    st.pop();
                }
            }

            // push the rotated ciphertext to the stack if i is even
            if (i && !(i % 2)) {
                st.push(ct_rot);
            }

            // For each diagonal not nul we perform mult + add with current result
            for (int l = 0; l < nbSubMatrices; l++) {
                if (!subDiags[l * nbSubMatrices + c][i].is_zero) {
                    auto mult = cc->EvalMult(subDiags[l * nbSubMatrices + c].at(i).packedDiag, ct_rot);
                    result[l] = cc->EvalAdd(mult, result[l]);
#ifdef EXPE
                    number_mult_plain++;
                    number_add++;
#endif
                }
            }
        }
    }

    return result;
}

Diag encode_diag(
    CryptoContext<DCRTPoly> cc,
    std::vector<std::vector<std::vector<double>>> mask_matrices,
    int subMat,
    int i) {
    Diag diag;
    auto matrix = mask_matrices[subMat];
    std::vector<double> diag_num;
    for (int j = 0; j < matrix.size(); j++) {
        diag_num.push_back(matrix[j][(i + j) % matrix.size()]);
    }

    if (std::accumulate(diag_num.begin(), diag_num.end(), 0.0) == 0) {
        diag.is_zero = true;
        diag.packedDiag = nullptr;
    } else {
        diag.packedDiag = cc->MakeCKKSPackedPlaintext(diag_num);
#ifdef EXPE
        number_encoding++;
#endif
        diag.is_zero = false;
    }

    return diag;
}

std::vector<Ciphertext<DCRTPoly>> matmul_coeus_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<std::vector<std::vector<double>>> mask_matrices,
    std::vector<Ciphertext<DCRTPoly>> subVectors) {
    int nbSubMatrices = subVectors.size();
    int s = cc->GetEncodingParams()->GetBatchSize();
    std::vector<double> zero(s, 0.0);
    auto pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encoding++;
    number_encryption++;
#endif
    std::vector<Ciphertext<DCRTPoly>> result(nbSubMatrices);
    for (int c = 0; c < nbSubMatrices; c++) {
        result[c] = ct_zero->Clone();
    }

    for (int c = 0; c < nbSubMatrices; c++) {
        std::stack<Ciphertext<DCRTPoly>> st;
        Ciphertext<DCRTPoly> RotatedVect = subVectors[c];
        auto precomp = cc->EvalFastRotationPrecompute(RotatedVect);

         for (int i = 0; i < s; i++) {
            // rotation
            Ciphertext<DCRTPoly> ct_rot;
            if (st.empty()) {
                ct_rot = cc->EvalFastRotation(RotatedVect, i, cc->GetRingDimension() * 2, precomp);
#ifdef EXPE
                number_rot++;
#endif
            } else {
                // find the largest power of 2 that is smaller than i
                int k = 1;
                while (!(i & k)) {
                    k = k << 1;
                }
                // rotate the ciphertext
                Ciphertext<DCRTPoly> cached_ct = st.top();
                ct_rot = cc->EvalRotate(cached_ct, k);
#ifdef EXPE
                number_rot++;
#endif

                // pop the stack if i is a power of 2
                if (i & (k << 1)) {
                    st.pop();
                }
            }

            // push the rotated ciphertext to the stack if i is even
            if (i && !(i % 2)) {
                st.push(ct_rot);
            }

            // For each diagonal not nul we perform mult + add with current result
            for (int l = 0; l < nbSubMatrices; l++) {
                Diag diag = encode_diag(cc, mask_matrices, l * nbSubMatrices + c, i);
                if (!diag.is_zero) {
                    auto mult = cc->EvalMult(diag.packedDiag, ct_rot);
                    result[l] = cc->EvalAdd(mult, result[l]);
#ifdef EXPE
                    number_mult_plain++;
                    number_add++;
#endif
                }
            }
        }
    }


    return result;
}

std::vector<Ciphertext<DCRTPoly>> matmul_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<std::vector<std::vector<double>>> mask_matrices,
    std::vector<Ciphertext<DCRTPoly>> subVectors) {
    int nbSubMatrices = subVectors.size();
    int s = cc->GetEncodingParams()->GetBatchSize();
    std::vector<double> zero(s, 0.0);
    auto pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encoding++;
    number_encryption++;
#endif
    std::vector<Ciphertext<DCRTPoly>> result(nbSubMatrices);
    for (int i = 0; i < nbSubMatrices; i++) {
        result[i] = ct_zero->Clone();
    }

    for (int c = 0; c < nbSubMatrices; c++) {
        Ciphertext<DCRTPoly> RotatedVect = subVectors[c];
        auto precomp = cc->EvalFastRotationPrecompute(RotatedVect);

        for (int i = 0; i < s; i++) {
            // rotation
            Ciphertext<DCRTPoly> ct_rot;
            ct_rot = cc->EvalFastRotation(RotatedVect, i, cc->GetRingDimension() * 2, precomp);
#ifdef EXPE
            number_rot++;
#endif
            // For each diagonal not nul we perform mult + add with current result
            for (int l = 0; l < nbSubMatrices; l++) {
                Diag diag = encode_diag(cc, mask_matrices, l * nbSubMatrices + c, i);
                if (!diag.is_zero) {
                    auto mult = cc->EvalMult(diag.packedDiag, ct_rot);
                    result[l] = cc->EvalAdd(mult, result[l]);
#ifdef EXPE
                    number_mult_plain++;
                    number_add++;
#endif
                }
            }
        }
    }

    return result;
}

std::vector<Ciphertext<DCRTPoly>> matmul_Halevi_Shoup_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<Diags> subDiags,
    std::vector<Ciphertext<DCRTPoly>> subVectors) {
    int nbSubMatrices = subVectors.size();
    int s = subDiags[0].size();
    std::vector<double> zero(s, 0.0);
    auto pl_zero = cc->MakeCKKSPackedPlaintext(zero);
    auto ct_zero = cc->Encrypt(keys.publicKey, pl_zero);
#ifdef EXPE
    number_encryption++;
#endif
    std::vector<Ciphertext<DCRTPoly>> result(nbSubMatrices, ct_zero);

    for (int c = 0; c < nbSubMatrices; c++) {
        std::vector<Ciphertext<DCRTPoly>> orderedSubVectors(nbSubMatrices);
        for (int l = 0; l < nbSubMatrices; l++) {
            // Rotation of vectors with respect to each other
            orderedSubVectors[l] = subVectors[(c + l) % nbSubMatrices];
        }

#ifdef EXPE
        int local_mult = 0;
#endif
        for (int i = 0; i < s; i++) {
            // For each diagonal not nul we perform mult + add with current result
            for (int l = 0; l < nbSubMatrices; l++) {
                if (!subDiags[l * nbSubMatrices + c].at(i).is_zero) {
                    auto mult = cc->EvalMult(subDiags[l * nbSubMatrices + c].at(i).packedDiag, orderedSubVectors[l]);
#ifdef EXPE
                    number_mult_plain++;
                    number_add++;
                    local_mult++;
#endif
                    result[l] = cc->EvalAdd(mult, result[l]);
                }
            }

            // Then we perform rotation for the next diagonal
            if (i < s - 1)
                orderedSubVectors = SplitRotation(cc, orderedSubVectors, s);
        }
#ifdef EXPE
        expe_poppy_m.addValue(
            local_mult, "Mat mult Halevi Shoup local mult number", Experiment::Operation::MUL, Experiment::UNIT);
#endif
    }

    return result;
}

/**
 * @brief Create submatrix such as we split the matrix in submatrix of size s*s
 *
 * @param matrix
 * @param s
 * @param m
 * @return std::vector<std::vector<std::vector<double>>>
 */
std::vector<std::vector<std::vector<double>>> prepareMatrix(std::vector<std::vector<double>> matrix, int s, int m) {
    padMatrixWithZerosInPlace(matrix, s);
    std::vector<std::vector<std::vector<double>>> subMatrices(
        (m / s) * (m / s), std::vector<std::vector<double>>(s, std::vector<double>(s)));

    for (int i = 0; i < m / s; i++) {
        for (int j = 0; j < m / s; j++) {
            // Copier les éléments de la matrice principale dans les sous-matrices
            for (int k = 0; k < s; ++k) {
                for (int l = 0; l < s; ++l) {
                    subMatrices[i * (m / s) + j][k][l] = matrix[i * s + k][j * s + l];
                }
            }
        }
    }
    return subMatrices;
}

/**
 * @brief From vector, create subvectors of size s and encrypt these
 *
 * @param cc
 * @param keys
 * @param vector
 * @param s Size of the vectors
 * @param m Size of the matrix
 * @return std::vector<Ciphertext<DCRTPoly>>
 */
std::vector<Ciphertext<DCRTPoly>>
prepareVector(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, std::vector<double> vector, int s, int m) {
    vector = padVectorWithZeros(vector, m);

    std::vector<Ciphertext<DCRTPoly>> subVectors(m / s);

    for (int l = 0; l < m / s; l++) {
        // Copier les éléments du vecteur dans les sous-vecteurs
        std::vector<double> subVector(s);
        for (int i = 0; i < s; ++i) {
            subVector[i] = vector[l * s + i];
        }
        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(subVector);
        Ciphertext<DCRTPoly> ct_vect = cc->Encrypt(keys.publicKey, ptxt);
#ifdef EXPE
        number_encryption++;
#endif

        subVectors[l] = ct_vect;
    }

    return subVectors;
}

/**
 * @brief From the matrix, extract the lower left triangle
 *
 * @param matrix a square matrix
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> lowerTriangle(std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));  // Lower triangle matrix

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            L[i][j] = matrix[i][j];
        }
    }

    return L;
}

/**
 * @brief From the matrix, extract the upper right triangle
 *
 * @param matrix a square matrix
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> upperTriangle(std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));  // Upper triangle matrix

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            U[i][j] = matrix[i][j];
        }
    }

    return U;
}

// Fonction pour créer les nouvelles sous-matrices Ti
std::vector<std::vector<std::vector<double>>> permuteSubMatrices(
    std::vector<std::vector<std::vector<double>>>& subMatrices,
    int s) {
    int nbSubMatrices = subMatrices.size();  // Nombre de sous-matrices

    // Sub matrix S_{l,c} at position l * nbSubMatricesPerLine + c of subMatrices
    // and of size s
    std::vector<std::vector<std::vector<double>>> newSubMatrices(
        nbSubMatrices,
        std::vector<std::vector<double>>(subMatrices[0].size(), std::vector<double>(subMatrices[0][0].size(), 0.0)));

    int nbSubMatricesPerLine = sqrt(nbSubMatrices);
    for (int l = 0; l < nbSubMatricesPerLine; l++) {
        for (int c = 0; c < nbSubMatricesPerLine; c++) {
            int i = l * nbSubMatricesPerLine + c;
            int nextIndex = l * nbSubMatricesPerLine + (c + 1) % nbSubMatricesPerLine;

            std::vector<std::vector<double>> upperTri = upperTriangle(subMatrices[i]);  // U(T_i) = U(S_i)
            std::vector<std::vector<double>> lowerTri =
                lowerTriangle(subMatrices[nextIndex]);  // L(T_i) = L(S_{nextIndex}})

            int i2 = l * nbSubMatricesPerLine + (c - l + nbSubMatricesPerLine) % nbSubMatricesPerLine;
            for (int row = 0; row < s; ++row) {
                for (int col = 0; col < s; ++col) {
                    if (row <= col)  // Si on est dans le triangle supérieur
                    {
                        newSubMatrices[i2][row][col] += upperTri[row][col];
                    } else {
                        newSubMatrices[i2][row][col] += lowerTri[row][col];
                    }
                }
            }
        }
    }
    return newSubMatrices;
}

std::vector<Ciphertext<DCRTPoly>>
SplitRotation(CryptoContext<DCRTPoly> cc, std::vector<Ciphertext<DCRTPoly>> vects, int s) {
    int vectSize = vects.size();
    std::vector<Ciphertext<DCRTPoly>> rotVectors(vectSize);
    for (int i = 0; i < vectSize; i++) {
        rotVectors[i] = cc->EvalRotate(vects[i], 1);
#ifdef EXPE
        number_rot++;
#endif
    }

    std::vector<double> mask1(s, 1.0);
    std::vector<double> maskNo1(s, 0.0);
    mask1[s - 1] = 0.0;
    maskNo1[s - 1] = 1.0;
    Plaintext pl_mask1 = cc->MakeCKKSPackedPlaintext(mask1);
    Plaintext pl_maskNo1 = cc->MakeCKKSPackedPlaintext(maskNo1);

    std::vector<Ciphertext<DCRTPoly>> splitVectors(vectSize);
    if (vectSize > 1) {
        for (int i = 0; i < vectSize; i++) {
            splitVectors[i] = cc->EvalAdd(
                cc->EvalMult(rotVectors[i], pl_mask1),
                cc->EvalMult(rotVectors[(i + 1 + vectSize) % vectSize], pl_maskNo1));
#ifdef EXPE
            number_mult_plain++;
            number_mult_plain++;
            number_add++;
#endif
        }
    } else {
        splitVectors[0] = rotVectors[0];
    }

    return splitVectors;
}

// For stats on Coeus/HaleviShoup with/without Fast Rotation
std::vector<Plaintext>
PackedDiagonals(CryptoContext<DCRTPoly> cc, std::vector<std::vector<double>>& matrix, int batchSize) {
    int n = matrix.size();

    // extract the diagonals and encode them as Plaintexts
    std::vector<Plaintext> result;
    for (int i = 0; i < n; i++) {
        std::vector<double> diag;
        for (int j = 0; j < n; j++) {
            diag.push_back(matrix[j][(i + j) % n]);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(diag, 1, 0, nullptr, batchSize);
        result.push_back(ptxt);
    }

    return result;
}

long matmult_HS(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect) {
    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> ct_result = cc->EvalMult(diags[0], ct_vect);
    Ciphertext<DCRTPoly> ct_rot;

    for (size_t i = 1; i < diags.size(); i++) {
        ct_rot = cc->EvalRotate(ct_vect, i);
        ct_result = cc->EvalAdd(ct_result, cc->EvalMult(diags[i], ct_rot));
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    return elapsed.count();
}

long matmult_HS_FastRot(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect) {
    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> ct_result = cc->EvalMult(diags[0], ct_vect);
    Ciphertext<DCRTPoly> ct_rot;
    auto precomp = cc->EvalFastRotationPrecompute(ct_vect);

    for (size_t i = 1; i < diags.size(); i++) {
        ct_rot = cc->EvalFastRotation(ct_vect, i, cc->GetRingDimension() * 2, precomp);
        ct_result = cc->EvalAdd(ct_result, cc->EvalMult(diags[i], ct_rot));
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    return elapsed.count();
}

long matmult_Coeus(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect) {
    auto start = std::chrono::high_resolution_clock::now();

    std::stack<Ciphertext<DCRTPoly>> st;
    Ciphertext<DCRTPoly> ct_result = cc->EvalMult(diags[0], ct_vect);
    Ciphertext<DCRTPoly> ct_temp;

    for (size_t i = 1; i < diags.size(); i++) {
        Ciphertext<DCRTPoly> ct_rot;
        if (st.empty()) {
            // i is a power of two
            assert((i == 0) || ((i & (i - 1)) == 0));
            ct_rot = cc->EvalRotate(ct_vect, i);
        } else {
            // find the largest power of 2 that is smaller than i
            int k = 1;
            while (!(i & k)) {
                k = k << 1;
            }
            // rotate the ciphertext
            Ciphertext<DCRTPoly> cached_ct = st.top();
            ct_rot = cc->EvalRotate(cached_ct, k);

            // pop the stack if i is a power of 2
            if ((i != 0) && ((i & (i - 1)) == 0)) {
                st.pop();
            }
        }

        // push the rotated ciphertext to the stack if i is even
        if (i && !(i % 2)) {
            st.push(ct_rot);
        }

        ct_temp = cc->EvalMult(diags[i], ct_rot);
        ct_result = cc->EvalAdd(ct_result, ct_temp);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    return elapsed.count();
}

long matmult_Coeus_FastRot(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect) {
    auto start = std::chrono::high_resolution_clock::now();

    std::stack<Ciphertext<DCRTPoly>> st;
    Ciphertext<DCRTPoly> ct_result = cc->EvalMult(diags[0], ct_vect);
    Ciphertext<DCRTPoly> ct_temp;
    auto precomp = cc->EvalFastRotationPrecompute(ct_vect);

    for (size_t i = 1; i < diags.size(); i++) {
        Ciphertext<DCRTPoly> ct_rot;
        if (st.empty()) {
            // i is a power of two
            assert((i == 0) || ((i & (i - 1)) == 0));
            ct_rot = cc->EvalFastRotation(ct_vect, i, cc->GetRingDimension() * 2, precomp);
        } else {
            // find the largest power of 2 that is smaller than i
            int k = 1;
            while (!(i & k)) {
                k = k << 1;
            }
            // rotate the ciphertext
            Ciphertext<DCRTPoly> cached_ct = st.top();
            ct_rot = cc->EvalRotate(cached_ct, k);

            // pop the stack if i is a power of 2
            if ((i != 0) && ((i & (i - 1)) == 0)) {
                st.pop();
            }
        }

        // push the rotated ciphertext to the stack if i is even
        if (i && !(i % 2)) {
            st.push(ct_rot);
        }

        ct_temp = cc->EvalMult(diags[i], ct_rot);
        ct_result = cc->EvalAdd(ct_result, ct_temp);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    return elapsed.count();
}
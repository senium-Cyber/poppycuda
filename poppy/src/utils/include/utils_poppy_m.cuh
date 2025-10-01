#pragma once
#include <cassert>
#include <map>
#include <vector>
#include "../../nauty/nausparse.h"
#include "utils.cuh"

using std::string;
using namespace lbcrypto;

struct Diag {
    Plaintext packedDiag;
    bool is_zero;
};

typedef std::vector<Diag> Diags;

std::vector<std::vector<double>> create_adjacency_matrix(
    std::vector<int> actives_nodes_vect,
    int* orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out);
std::vector<std::vector<double>> create_adjacency_matrix(
    std::vector<int> actives_nodes_vect,
    std::map<int /*node*/, int /*orbit*/> orbits,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    std::map<int, int> representative_orbits,
    std::map<Node, int> partition,
    int datacenter_id);
void print_matrix(std::vector<std::vector<double>> matrix, int size);

std::vector<Diags> extract_and_encode_diags(
    CryptoContext<DCRTPoly> cc,
    std::vector<std::vector<double>> mask_matrix,
    int batchSize,
    int& m);
/** Function to pad matrix with zeros to creat a square matrix of size batchSize
 */
std::vector<std::vector<double>> padMatrixWithZeros(std::vector<std::vector<double>> matrix, int batchSize);
void padMatrixWithZerosInPlace(std::vector<std::vector<double>>& matrix, int batchSize);

/** Function to pad vector with zeros to creat a vector of size batchSize
 */
std::vector<double> padVectorWithZeros(std::vector<double> vector, int batchSize);

/**
 * Creat diagonals of Matrix from left to right and encode them with CKKS Crypto Context. This function don't packed
 * null diagonals.
 */
Diags MakeCKKSPackedDiagonals(CryptoContext<DCRTPoly> cc, std::vector<std::vector<double>>& matrix);

/**
 * Perform vector (ct1) and matrix (given by diags) multiplication with Coeus algorithm.
 * The vector diags_null have size the number of diagonals in the matrix with 1 if the diag is null, 0 otherwise.
 */
Ciphertext<DCRTPoly> matmul_a_la_Coeus(
    CryptoContext<DCRTPoly> cc,
    std::vector<Plaintext> diags,
    Ciphertext<DCRTPoly> ct1,
    std::vector<int> diags_null);
std::vector<Ciphertext<DCRTPoly>> matmul_Halevi_Shoup_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<Diags> subDiags,
    std::vector<Ciphertext<DCRTPoly>> subVectors);
std::vector<Ciphertext<DCRTPoly>> matmul_coeus_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<Diags> subDiags,
    std::vector<Ciphertext<DCRTPoly>> subVectors);

std::vector<std::vector<std::vector<double>>> prepareMatrix(std::vector<std::vector<double>> matrix, int s, int m);

std::vector<Ciphertext<DCRTPoly>>
prepareVector(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, std::vector<double> vector, int s, int m);

std::vector<std::vector<double>> lowerTriangle(std::vector<std::vector<double>>& matrix);

std::vector<std::vector<double>> upperTriangle(std::vector<std::vector<double>>& matrix);

std::vector<std::vector<std::vector<double>>> permuteSubMatrices(
    std::vector<std::vector<std::vector<double>>>& subMatrices,
    int s);

std::vector<Ciphertext<DCRTPoly>>
SplitRotation(CryptoContext<DCRTPoly> cc, std::vector<Ciphertext<DCRTPoly>> vects, int s);

// For stats on Coeus/HaleviShoup with/without Fast Rotation
std::vector<Plaintext>
PackedDiagonals(CryptoContext<DCRTPoly> cc, std::vector<std::vector<double>>& matrix, int batchSize);

std::vector<Ciphertext<DCRTPoly>> matmul_coeus_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<std::vector<std::vector<double>>> mask_matrices,
    std::vector<Ciphertext<DCRTPoly>> subVectors);

std::vector<Ciphertext<DCRTPoly>> matmul_splited(
    CryptoContext<DCRTPoly> cc,
    KeyPair<DCRTPoly> keys,
    std::vector<std::vector<std::vector<double>>> mask_matrices,
    std::vector<Ciphertext<DCRTPoly>> subVectors);

// The four next functions return the time to compute the matrix multiplication in milliseconds
long matmult_HS(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect);
long matmult_HS_FastRot(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect);
long matmult_Coeus(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect);
long matmult_Coeus_FastRot(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct_vect);

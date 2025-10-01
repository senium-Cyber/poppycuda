#include <chrono>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "openfhe.h"
using namespace lbcrypto;

using namespace std;

uint32_t nextPowerOfTwo(uint32_t n) {
    --n;

    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    return n + 1;
}

struct Diags {
    std::vector<Plaintext> PackedDiags;
    std::vector<int> diags_null;
};

vector<double> globalDegree;

int findIndex(vector<pair<int, double>> a, int x) {
    for (int i = 0; i < a.size(); i++)
        if (a[i].first == x)
            return i;
    return -1;
}

bool compareDegree(int i, int j) {
    return ::globalDegree[i] < ::globalDegree[j];
}

template <typename T>
ostream& operator<<(ostream& out, vector<T> const& v) {
    for (int i = 0; i < v.size(); i++)
        out << v[i] << ' ';
    return out;
}

class AdjacencyMatrix {
   private:
    vector<vector<double>> _matrix;

   public:
    // Constructor and Destructor
    AdjacencyMatrix(vector<vector<double>> m) { _matrix = m; }

    AdjacencyMatrix() {}
    ~AdjacencyMatrix() {}

    // class methods

    // Function to get the matrix
    vector<vector<double>> getMatrix() { return _matrix; }

    // Function to get the size of the matrix
    int size() { return _matrix.size(); }

    // Function to compute degree of all the nodes
    vector<double> degree_from_matrix() {
        vector<double> degrees;

        for (int i = 0; i < _matrix.size(); i++) {
            double count = 0;
            for (int j = 0; j < _matrix[0].size(); j++) {
                if (_matrix[i][j] != 0)
                    count++;
            }

            degrees.push_back(count);
        }

        return degrees;
    }

    // Implementation of Cuthill-Mckee algorithm
    vector<int> CuthillMckee() {
        vector<double> degrees = degree_from_matrix();

        ::globalDegree = degrees;

        queue<int> Q;
        vector<int> R;
        vector<pair<int, double>> notVisited;

        for (int i = 0; i < degrees.size(); i++)
            notVisited.push_back(make_pair(i, degrees[i]));

        // Vector notVisited helps in running BFS
        // even when there are dijoind graphs
        while (notVisited.size()) {
            int minNodeIndex = 0;

            for (int i = 0; i < notVisited.size(); i++)
                if (notVisited[i].second < notVisited[minNodeIndex].second)
                    minNodeIndex = i;

            Q.push(notVisited[minNodeIndex].first);

            notVisited.erase(notVisited.begin() + findIndex(notVisited, notVisited[Q.front()].first));

            // Simple BFS
            while (!Q.empty()) {
                vector<int> toSort;

                for (int i = 0; i < _matrix[0].size(); i++) {
                    if (i != Q.front() && _matrix[Q.front()][i] != 0 && findIndex(notVisited, i) != -1) {
                        toSort.push_back(i);
                        notVisited.erase(notVisited.begin() + findIndex(notVisited, i));
                    }
                }

                sort(toSort.begin(), toSort.end(), compareDegree);

                for (int i = 0; i < toSort.size(); i++)
                    Q.push(toSort[i]);

                R.push_back(Q.front());
                Q.pop();
            }
        }

        return R;
    }

    // Implementation of reverse Cuthill-Mckee algorithm
    vector<int> ReverseCuthillMckee() {
        vector<int> cuthill = CuthillMckee();

        int n = cuthill.size();

        if (n % 2 == 0)
            n -= 1;

        n = n / 2;

        for (int i = 0; i <= n; i++) {
            int j = cuthill[cuthill.size() - 1 - i];
            cuthill[cuthill.size() - 1 - i] = cuthill[i];
            cuthill[i] = j;
        }

        return cuthill;
    }

    // Function to permute a matrix based on permutation vector
    AdjacencyMatrix permute(const vector<int>& permutation) {
        vector<vector<double>> permuted(_matrix.size(), vector<double>(_matrix.size()));
        for (int i = 0; i < _matrix.size(); ++i) {
            for (int j = 0; j < _matrix[i].size(); ++j) {
                permuted[i][j] = _matrix[permutation[i]][permutation[j]];
            }
        }
        return AdjacencyMatrix(permuted);
    }

    AdjacencyMatrix transpose() {
        // Vérifier si la matrice est vide
        if (_matrix.empty() || _matrix[0].empty())
            return {};

        int rows = _matrix.size();
        int cols = _matrix[0].size();

        // Initialiser la transposée avec les dimensions inversées
        vector<vector<double>> transposed(cols, std::vector<double>(rows));

        // Remplir la transposée
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                transposed[j][i] = _matrix[i][j];
            }
        }

        return AdjacencyMatrix(transposed);
    }

    // Fonction pour additionner deux matrices
    AdjacencyMatrix add_to(AdjacencyMatrix matrix) {
        auto matrix2 = matrix.getMatrix();

        // Vérifier si les matrices ont la même taille
        if (_matrix.size() != matrix2.size() || _matrix.empty() || _matrix[0].size() != matrix2[0].size()) {
            std::cerr << "Erreur: Les matrices doivent avoir la même taille pour l'addition." << std::endl;
            return {};
        }

        int rows = _matrix.size();
        int cols = _matrix[0].size();

        // Initialiser la matrice résultante avec les dimensions de l'une des matrices
        std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

        // Additionner les éléments des deux matrices
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = _matrix[i][j] + matrix2[i][j];
            }
        }

        return AdjacencyMatrix(result);
    }

    // Fonction pour soustraire deux matrices
    AdjacencyMatrix subtract_to(AdjacencyMatrix matrix) {
        auto matrix2 = matrix.getMatrix();

        // Vérifier si les matrices ont la même taille
        if (_matrix.size() != matrix2.size() || _matrix.empty() || _matrix[0].size() != matrix2[0].size()) {
            std::cerr << "Erreur: Les matrices doivent avoir la même taille pour la soustraction." << std::endl;
            return {};
        }

        int rows = _matrix.size();
        int cols = _matrix[0].size();

        // Initialiser la matrice résultante avec les dimensions de l'une des matrices
        std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

        // Soustraire les éléments des deux matrices
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i][j] = _matrix[i][j] - matrix2[i][j];
            }
        }

        return AdjacencyMatrix(result);
    }

    // Fonction pour imprimer la matrice
    void printMatrix() {
        for (int i = 0; i < _matrix.size(); ++i) {
            for (int j = 0; j < _matrix[i].size(); ++j) {
                cout << _matrix[i][j] << " ";
            }
            cout << endl;
        }
        printf("***********\n");
    }
};

AdjacencyMatrix UnsymetricReverseCuthillMcKee(AdjacencyMatrix matrix) {
    auto m_t = matrix.transpose();                                              // At
    auto sum = matrix.add_to(m_t);                                              // A + At
    auto permutation = sum.ReverseCuthillMckee();                               // P
    auto permuted_matrix_sum = sum.permute(permutation);                        // Pt(A + At)P
    auto permuted_matrix_t = m_t.permute(permutation);                          // PtAtP
    auto permuted_matrix = permuted_matrix_sum.subtract_to(permuted_matrix_t);  // Pt(A + At)P - PtAtP = Pt A P
    return permuted_matrix;
}

std::vector<std::vector<double>> padMatrixWithZeros(std::vector<std::vector<double>> matrix, int batchSize) {
    int n = matrix.size();
    // add zeros to each lines of the matrix
    for (int i = 0; i < n; i++) {
        for (int j = n; j < batchSize; j++) {
            matrix[i].push_back(0);
        }
    }

    // add lines of zeros to the matrix
    std::vector<double> line_zeros(batchSize, 0);
    for (int i = n; i < batchSize; i++) {
        matrix.push_back(line_zeros);
    }

    return matrix;
}

/* Extract diagonals form the matrix and identify diagonals null*/
Diags MakeCKKSPackedDiagonals(CryptoContext<DCRTPoly> cc, std::vector<std::vector<double>>& matrix, int batchSize) {
    Diags diags;

    int n = matrix.size();
    if (n < batchSize) {
        matrix = padMatrixWithZeros(matrix, batchSize);
        n = matrix.size();  // batchSize
    }

    // extract the diagonals and encode them as Plaintexts
    std::vector<double> zero(n, 0.0);
    // std::vector<int> diags_null(n,0);
    int nb_diag_zero = 0;
    for (int i = 0; i < n; i++) {
        std::vector<double> diag;
        for (int j = 0; j < n; j++) {
            diag.push_back(matrix[j][(i + j) % n]);
        }

        if (diag == zero) {
            diags.diags_null.push_back(1);
            nb_diag_zero += 1;
        } else {
            diags.diags_null.push_back(0);
            Plaintext ptxt = cc->MakeCKKSPackedPlaintext(diag);
            diags.PackedDiags.push_back(ptxt);
        }
    }

    std::cout << "nb diags : " << diags.diags_null.size() << ", nb diag null : " << nb_diag_zero << std::endl;

    return diags;
}

Ciphertext<DCRTPoly> matmul_a_la_Coeus_diag_null_opt(
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

        // multiply the rotated ciphertext with the diagonal if it is not null
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
    // print cpt_diag_not_null
    std::cout << "cpt_diag_not_null : " << cpt_diag_not_null << std::endl;

    return ct_result;
}

Ciphertext<DCRTPoly>
matmul_a_la_Coeus_wo_any_opt(CryptoContext<DCRTPoly> cc, std::vector<Plaintext> diags, Ciphertext<DCRTPoly> ct1) {
    // Get the diagonals of the matrix as Plaintexts
    std::stack<Ciphertext<DCRTPoly>> st;
    Ciphertext<DCRTPoly> ct_result;
    Ciphertext<DCRTPoly> ct_temp;

    for (size_t i = 0; i < diags.size(); i++) {
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
        if (i == 0) {
            ct_result = cc->EvalMult(diags[i], ct_rot);
        } else {
            ct_temp = cc->EvalMult(diags[i], ct_rot);
            ct_result = cc->EvalAdd(ct_result, ct_temp);
        }
    }

    return ct_result;
}

// TODO : a modifier pour ne pas faire les rotations pas necessaire (i.e quand diag null)
Ciphertext<DCRTPoly> mat_mul_a_la_halevi_shoup(
    CryptoContext<DCRTPoly> cc,
    std::vector<Plaintext>& diag,
    Ciphertext<DCRTPoly>& v,
    int batchSize,
    std::vector<int> diags_null) {
    Ciphertext<DCRTPoly> ct_result = cc->EvalMult(diag[0], v);

    Ciphertext<DCRTPoly> ct_rot = cc->EvalRotate(v, 1);
    for (size_t i = 1; i < diag.size(); i++) {
        if (diags_null[i] == 0) {
            ct_result = cc->EvalAdd(ct_result, cc->EvalMult(diag[i], ct_rot));
            ct_rot = cc->EvalRotate(ct_rot, 1);
        }
    }

    return ct_result;
}

std::vector<std::vector<double>> readGraph(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::unordered_map<int, int> nodeIndexMap;

    if (file.is_open()) {
        // Skip the lines until we reach the actual data
        while (std::getline(file, line) && line.find("#") != std::string::npos)
            ;

        // Read the data and build the adjacency matrix
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int fromNode, toNode;

            if (!(iss >> fromNode >> toNode)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }

            // If the fromNode index doesn't exist in the map, add it
            if (nodeIndexMap.find(fromNode) == nodeIndexMap.end()) {
                int newIndex = nodeIndexMap.size();
                nodeIndexMap[fromNode] = newIndex;
            }

            // If the toNode index doesn't exist in the map, add it
            if (nodeIndexMap.find(toNode) == nodeIndexMap.end()) {
                int newIndex = nodeIndexMap.size();
                nodeIndexMap[toNode] = newIndex;
            }
        }

        // Initialize the adjacency matrix with the correct size
        int numNodes = nodeIndexMap.size();

        // Initialize the adjacency matrix with the correct size
        std::vector<std::vector<double>> adjacencyMatrix(numNodes, std::vector<double>(numNodes, 0.0));

        printf("numNodes : %d\n", numNodes);

        // Rewind to the beginning of the file
        file.clear();
        file.seekg(0);

        // Skip the lines until we reach the actual data
        while (std::getline(file, line) && line.find("#") != std::string::npos)
            ;

        // Read the data and fill in the adjacency matrix
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int fromNode, toNode;

            if (!(iss >> fromNode >> toNode)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            // printf("edge : %d -> %d\n", fromNode, nodeIndexMap[fromNode]);
            // Update the adjacency matrix
            adjacencyMatrix[nodeIndexMap[fromNode]][nodeIndexMap[toNode]] = 1.0;
        }

        file.close();

        return adjacencyMatrix;
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }
}

int main() {
    // declare start end and elapsed withour initializing them
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed;

    /* Instantiation of the adjacency matrix */
    std::cout << "------ Reading the graph ------- " << std::endl;

    // std::string graphName = "Wiki-Vote";
    std::string graphName = "p2p-Gnutella06";
    // std::string graphName = "Cit-HepTh";

    std::string filename = "../Graphs/" + graphName + ".txt";
    std::vector<std::vector<double>> matrix = readGraph(filename);

    // std::string filename = "./Graphs/Email-EuAll.txt";
    // std::vector<std::vector<double>> matrix = readGraph(filename);

    AdjacencyMatrix m(matrix);

    /* Instantiation of the crypto material */
    uint32_t multDepth = 10;
    uint32_t scaleModSize = 50;
    uint32_t batchSize = nextPowerOfTwo(m.getMatrix().size());

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetSecurityLevel(HEStd_NotSet);  // disable security
    int N = 1 << 15;
    parameters.SetRingDim(N);  // a power of 2 integer
    parameters.SetBatchSize(batchSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    // int N = cc->GetRingDimension();

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    std::cout << "CKKS scheme is using ring dimension " << N << std::endl << std::endl;

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    std::vector<int> steps;
    steps.push_back(0);
    for (int i = 1; i < N / 2; i *= 2) {
        steps.push_back(i);
    }

    cc->EvalRotateKeyGen(keys.secretKey, steps);

    // Encrypt the encoded vectors
    std::vector<double> v(batchSize, 1.0);
    Plaintext ptxt_v = cc->MakeCKKSPackedPlaintext(v);
    auto ct_v = cc->Encrypt(keys.publicKey, ptxt_v);

    std::string log_N = std::to_string(int(log2(N)));

    /* Test MatMul with the Reverse Cuthill-McKee optimization */

    // Ouverture du fichier pour écrire les résultats au format CSV
    std::string outputName_rcm = "../results-rcm/" + log_N + graphName + ".txt";
    std::ofstream outputFile_rcm(outputName_rcm);
    if (!outputFile_rcm.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier pour écrire les résultats." << std::endl;
        return 1;
    }

    // Écrire l'en-tête du fichier CSV
    outputFile_rcm << "Test, Reverse Cuthill McKee (s), Multiplication matrice-vecteur (s)\n";

    // Nombre de tests à effectuer
    int numTests_rcm = 10;
    std::cout << "------ Matmul with RCM ------- " << std::endl;
    for (int test = 0; test < numTests_rcm; ++test) {
        std::cout << "Test " << test + 1 << "..." << std::endl;
        // Effectuer le test Reverse Cuthill McKee
        start = std::chrono::high_resolution_clock::now();
        auto permuted_matrix = UnsymetricReverseCuthillMcKee(m);
        auto A = permuted_matrix.getMatrix();
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;  // ce temps est le temps de RCM + permuté la matrice

        auto diags = MakeCKKSPackedDiagonals(cc, A, batchSize);
        // Écrire le résultat dans le fichier CSV
        outputFile_rcm << test + 1 << ", " << elapsed.count() << ", ";

        // Effectuer le test de multiplication matrice-vecteur homomorphique
        start = std::chrono::high_resolution_clock::now();
        auto ct_result = matmul_a_la_Coeus_diag_null_opt(cc, diags.PackedDiags, ct_v, diags.diags_null);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        // Écrire le résultat dans le fichier CSV
        outputFile_rcm << elapsed.count() << "\n";

        // Forcer l'écriture immédiate dans le fichier
        outputFile_rcm.flush();
    }

    // Fermeture du fichier
    outputFile_rcm.close();

    /* Test MatMul without the Reverse Cuthill-McKee optimization */

    // Ouverture du fichier pour écrire les résultats au format CSV
    std::string outputName_worcm = "../results-rcm/" + log_N + graphName + "-worcm.txt";
    std::ofstream outputFile_worcm(outputName_worcm);
    if (!outputFile_worcm.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier pour écrire les résultats (without rcm)." << std::endl;
        return 1;
    }

    // Écrire l'en-tête du fichier CSV
    outputFile_worcm << "Test, Multiplication matrice-vecteur (s)\n";

    // Nombre de tests à effectuer
    int numTests = 10;
    std::cout << "------ Matmul without RCM ------- " << std::endl;
    for (int test = 0; test < numTests; ++test) {
        std::cout << "Test " << test + 1 << "..." << std::endl;
        auto A = m.getMatrix();
        auto diags = MakeCKKSPackedDiagonals(cc, A, batchSize);

        // Effectuer le test de multiplication matrice-vecteur homomorphique
        start = std::chrono::high_resolution_clock::now();
        auto ct_result = matmul_a_la_Coeus_diag_null_opt(cc, diags.PackedDiags, ct_v, diags.diags_null);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        // Écrire le résultat dans le fichier CSV
        outputFile_worcm << test + 1 << ", " << elapsed.count() << "\n";

        // Forcer l'écriture immédiate dans le fichier
        outputFile_worcm.flush();
    }

    // Fermeture du fichier
    outputFile_worcm.close();

    /* Test MatMul without any optimization (i.e neither RCM nor ignoring the diags null) */

    // Ouverture du fichier pour écrire les résultats au format CSV
    std::string outputName_wo_any_opt = "../results-rcm/" + log_N + graphName + "-wo_any_opt.txt";
    std::ofstream outputFile_wo_any(outputName_wo_any_opt);
    if (!outputFile_wo_any.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier pour écrire les résultats (without rcm)." << std::endl;
        return 1;
    }

    // Écrire l'en-tête du fichier CSV
    outputFile_wo_any << "Test, Multiplication matrice-vecteur (s)\n";

    // Nombre de tests à effectuer
    int numTests_wo_opt = 10;
    std::cout << "------ Matmul without any optimization ------- " << std::endl;
    for (int test = 0; test < numTests_wo_opt; ++test) {
        std::cout << "Test " << test + 1 << "..." << std::endl;
        auto A = m.getMatrix();
        auto diags = MakeCKKSPackedDiagonals(cc, A, batchSize);

        // Effectuer le test de multiplication matrice-vecteur homomorphique
        start = std::chrono::high_resolution_clock::now();
        auto ct_result = matmul_a_la_Coeus_wo_any_opt(cc, diags.PackedDiags, ct_v);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        // Écrire le résultat dans le fichier CSV
        outputFile_wo_any << test + 1 << ", " << elapsed.count() << "\n";

        // Forcer l'écriture immédiate dans le fichier
        outputFile_wo_any.flush();
    }

    // Fermeture du fichier
    outputFile_wo_any.close();

    return 0;
}

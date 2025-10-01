#include "include/utils_crypto.cuh"

#ifdef EXPE
#include "include/utils_export.cuh"
#endif

CryptoContext<DCRTPoly> gen_crypto_context_test(uint32_t multDepth) {
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    // First crypto context generation for learning ringDim
    CCParams<CryptoContextCKKSRNS> paramTest;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    // Currently, only FIXEDMANUAL and FIXEDAUTO modes are supported for 128-bit CKKS bootstrapping.
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    // All modes are supported for 64-bit CKKS bootstrapping.
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits               = 59;
    usint firstMod               = 60;
#endif

#ifdef BOOTSTRAPPING
    uint32_t dnum = 6; // OpenFHE recommandation for large precision
    if (multDepth >= 20)
        dnum = 6; // OpenFHE recommandation for large precision
    else
       dnum = 4;
#endif

    paramTest.SetSecretKeyDist(secretKeyDist);
    paramTest.SetMultiplicativeDepth(multDepth);
    paramTest.SetSecurityLevel(HEStd_128_classic);

    paramTest.SetScalingModSize(dcrtBits);
    paramTest.SetScalingTechnique(rescaleTech);
    paramTest.SetFirstModSize(firstMod);
#ifdef BOOTSTRAPPING
    paramTest.SetNumLargeDigits(dnum);
#endif

    auto ccTest = GenCryptoContext(paramTest);
    return ccTest;
}

/**
 * @brief Generate the crypto context for the PR computation
 *
 * @param multDepth
 * @param number_active_nodes
 * @return std::pair<CryptoContext<DCRTPoly>, KeyPair<DCRTPoly>>
 */
std::pair<CryptoContext<DCRTPoly>, KeyPair<DCRTPoly>> gen_crypto_context(uint32_t multDepth, int batchSize) {
    std::cout << "Generating crypto context" << std::endl;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    // Currently, only FIXEDMANUAL and FIXEDAUTO modes are supported for 128-bit CKKS bootstrapping.
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    // All modes are supported for 64-bit CKKS bootstrapping.
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits               = 59;
    usint firstMod               = 60;
#endif

#ifdef BOOTSTRAPPING
    uint32_t dnum = 6; // OpenFHE recommandation for large precision
    if (multDepth >= 25)
        dnum = 6; // OpenFHE recommandation for large precision
    else
       dnum = 4;
#endif

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecurityLevel(HEStd_128_classic);

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    parameters.SetBatchSize(batchSize);
#ifdef BOOTSTRAPPING
    parameters.SetNumLargeDigits(dnum);
#endif

#ifdef EXPE
    expe_poppy_m.setCryptoParams(parameters);
    expe_poppy_uj.setCryptoParams(parameters);
    expe_naive.setCryptoParams(parameters);
#endif

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

    KeyPair<DCRTPoly> keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    return {cc, keys};
}

/**
 * @brief Generate the rotation keys
 *
 */
void gen_rotation_keys(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, uint32_t nbSlots) {
    std::vector<int> steps;
    steps.push_back(0);

    for (size_t i = 1; i < nbSlots; i++) {
        steps.push_back(i);
    }

    cc->EvalRotateKeyGen(keys.secretKey, steps);
}

/**
 * @brief Generate the rotation keys for Coeus
 * In coeus we only use powers of 2 rotations
 *
 */
void gen_rotation_keys_coeus(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, uint32_t nbSlots) {
    std::vector<int> steps;
    steps.push_back(0);
    for (int i = 1; i < nbSlots; i *= 2) {
        steps.push_back(i);
    }

    cc->EvalRotateKeyGen(keys.secretKey, steps);
}

/**
 * @brief We assume that the crypto context is already generated, multKey and rotKey
 *
 * @param cc
 */
void serialize_cc(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, int multDepth) {
    auto batchSize = cc->GetCryptoParameters()->GetEncodingParams()->GetBatchSize();
    std::string dir = get_folder_by_params(multDepth, batchSize);
    struct stat sb;
    if (stat(dir.data(), &sb) != 0) {
        if (stat(SERIAL_FOLDER, &sb) != 0) {
            if (mkdir(SERIAL_FOLDER, 0777) != 0) {
                std::cerr << "Error creating directory " << SERIAL_FOLDER << std::endl;
                std::exit(1);
            }
        }
        if (mkdir(dir.data(), 0777) != 0) {
            std::cerr << "Error creating directory " << dir << std::endl;
            std::exit(1);
        }
    } else {
        std::cerr << "Directory " << dir << " already exists" << std::endl;
    }
    if (!Serial::SerializeToFile(dir + "/" SERIAL_CC, cc, SerType::BINARY)) {
        std::cerr << "Error writing serialization of the crypto context to "
                     "cryptocontext.txt"
                  << std::endl;
        std::exit(1);
    }

    std::cout << "Cryptocontext serialized" << std::endl;

    if (!Serial::SerializeToFile(dir + "/" SERIAL_PUB, keys.publicKey, SerType::BINARY)) {
        std::cerr << "Exception writing public key to " SERIAL_PUB << std::endl;
        std::exit(1);
    }
    std::cout << "Public key serialized" << std::endl;

    if (!Serial::SerializeToFile(dir + "/" SERIAL_SECRET, keys.secretKey, SerType::BINARY)) {
        std::cerr << "Exception writing private key to " SERIAL_SECRET << std::endl;
        std::exit(1);
    }
    std::cout << "Private key serialized" << std::endl;

    std::ofstream multKeyFile(dir + "/" SERIAL_MULT, std::ios::out | std::ios::binary);
    if (multKeyFile.is_open()) {
        if (!cc->SerializeEvalMultKey(multKeyFile, SerType::BINARY)) {
            std::cerr << "Error writing eval mult keys" << std::endl;
            std::exit(1);
        }
        std::cout << "EvalMult/ relinearization keys have been serialized" << std::endl;
        multKeyFile.close();
    } else {
        std::cerr << "Error serializing EvalMult keys" << std::endl;
        std::exit(1);
    }

    std::ofstream rotationKeyFile(dir + "/" SERIAL_ROT, std::ios::out | std::ios::binary);
    if (rotationKeyFile.is_open()) {
        if (!cc->SerializeEvalAutomorphismKey(rotationKeyFile, SerType::BINARY)) {
            std::cerr << "Error writing rotation keys" << std::endl;
            std::exit(1);
        }
        std::cout << "Rotation keys have been serialized" << std::endl;
    } else {
        std::cerr << "Error serializing Rotation keys" << std::endl;
        std::exit(1);
    }
}

bool deserialize_cc(CryptoContext<DCRTPoly>* cc, KeyPair<DCRTPoly>* keys, std::string folder_name) {
    struct stat sb;
    if (stat(folder_name.data(), &sb) != 0) {
        std::cerr << "directory " << folder_name << " does not exist" << std::endl;
        return false;
    }

    if (!Serial::DeserializeFromFile(folder_name + SERIAL_CC, *cc, SerType::BINARY)) {
        std::cerr << "I cannot read serialized data from: " << folder_name << "/cryptocontext.txt" << std::endl;
        std::exit(1);
    }
    std::cout << "Client CC deserialized" << std::endl;

    if (!Serial::DeserializeFromFile(folder_name + SERIAL_PUB, keys->publicKey, SerType::BINARY)) {
        std::cerr << "I cannot read serialized data from: " << folder_name << "/key_public.txt" << std::endl;
        std::exit(1);
    }
    std::cout << "Client Public Key deserialized" << std::endl;

    if (!Serial::DeserializeFromFile(folder_name + SERIAL_SECRET, keys->secretKey, SerType::BINARY)) {
        std::cerr << "I cannot read serialized data from: " << folder_name << "/key_secret.txt" << std::endl;
        std::exit(1);
    }
    std::cout << "Client Secret Key deserialized" << std::endl;

    std::ifstream multKeyIStream(folder_name + SERIAL_MULT, std::ios::in | std::ios::binary);
    if (!multKeyIStream.is_open()) {
        std::cerr << "Cannot read serialization from " << folder_name + SERIAL_MULT << std::endl;
        std::exit(1);
    }
    if (!(*cc)->DeserializeEvalMultKey(multKeyIStream, SerType::BINARY)) {
        std::cerr << "Could not deserialize eval mult key file" << std::endl;
        std::exit(1);
    }
    std::cout << "Eval mult keys deserialized" << std::endl;

    std::ifstream rotKeyIStream(folder_name + SERIAL_ROT, std::ios::in | std::ios::binary);
    if (!rotKeyIStream.is_open()) {
        std::cerr << "Cannot read serialization from " << folder_name + SERIAL_ROT << std::endl;
        std::exit(1);
    }
    if (!(*cc)->DeserializeEvalAutomorphismKey(rotKeyIStream, SerType::BINARY)) {
        std::cerr << "Could not deserialize eval rot key file" << std::endl;
        std::exit(1);
    }
    std::cout << "Eval rotation keys deserialized" << std::endl;
    return true;
}

void serialize_cipher(Ciphertext<DCRTPoly> ct, std::string folder_name, std::string filename) {
    struct stat sb;
    if (stat(folder_name.data(), &sb) != 0) {
        if (stat(SERIAL_FOLDER, &sb) != 0) {
            if (mkdir(SERIAL_FOLDER, 0777) != 0) {
                std::cerr << "Error creating directory " << SERIAL_FOLDER << std::endl;
                std::exit(1);
            }
        }
        if (mkdir(folder_name.data(), 0777) != 0) {
            std::cerr << "Error creating directory " << folder_name << std::endl;
            std::exit(1);
        }
    } else {
        std::cerr << "Directory " << folder_name << " already exists" << std::endl;
    }

    if (!Serial::SerializeToFile( folder_name + "/" + folder_name + filename , ct, SerType::BINARY)) {
        std::cerr << "Error writing cipher to " << folder_name + "/" + folder_name + filename << std::endl;
        std::exit(1);
    }
    std::cout << "Cipher serialized" << std::endl;  
}

std::string get_folder_by_params(int multDepth, int batchSize) {
    return SERIAL_FOLDER "m-" + std::to_string(multDepth) + "_b-" + std::to_string(batchSize) + "/";
}
#include <iostream>
#include "../src/utils/include/utils.cuh"
#include "../src/utils/include/utils_export.cuh"
#include "../src/utils/include/utils_graph.cuh"
#include "../src/utils/include/utils_poppy_m.cuh"
#include "../src/utils/include/utils_stats.cuh"
#include "openfhe.h"
using namespace std::chrono;
using namespace lbcrypto;
#define NB_TRY 100

uint32_t nextPowerOf2(uint32_t n) {
    --n;

    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    return n + 1;
}

CryptoContext<DCRTPoly> generate_common_crypto_context(
    int nbSlots,
    CCParams<CryptoContextCKKSRNS>& parameters,
    uint32_t multDepth = 16,
    int ringDim = 0) {
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecurityLevel(HEStd_128_classic);
    if (ringDim != 0) {
        parameters.SetRingDim(ringDim);
    }

// Scaling parameters for first and second crypto context
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits = 78;
    usint firstMod = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // usint dcrtBits = 59;
    // usint firstMod = 60;
    usint dcrtBits = 30;
    usint firstMod = 35;
    // usint dcrtBits = 50;
    // usint firstMod = dcrtBits + 5;
    // usint dcrtBits = 50;
    // usint firstMod = 60;
#endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    if (nbSlots != 0)
        parameters.SetBatchSize(nbSlots);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    return cc;
}

CryptoContext<DCRTPoly> generate_common_crypto_context(CCParams<CryptoContextCKKSRNS>& parameters) {
    uint32_t multDepth = 8;

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecurityLevel(HEStd_128_classic);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // usint dcrtBits = 59;
    // usint firstMod = 60;
    usint dcrtBits = 30;
    usint firstMod = 35;
    // usint dcrtBits = 50;
    // usint firstMod = dcrtBits + 5;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << ", with batching." << std::endl
              << std::endl;
    return cc;
}

#define NB_TESTS 1000
void compare_rotation_key_gen() {
    std::cout << "Comparison of the different operations proposed by openFHE CKKS scheme for " << NB_TESTS
              << " iterations.\n"
              << std::endl;

    int number_active_nodes = 165;
    uint32_t nbSlots = nextPowerOf2(number_active_nodes);

    // Creation of vectors, Plaintexts and Ciphertexts:

    std::vector<double> v_1(nbSlots);
    std::vector<double> v_2(nbSlots);
    std::vector<double> v_zero(nbSlots);
    for (int i = 0; i < (int)nbSlots; i++) {
        v_1[i] = 3.1415;
        v_2[i] = 1.234;
        v_zero[i] = 0;
    }

    std::vector<double> v_N_1(nbSlots);
    std::vector<double> v_N_2(nbSlots);
    for (int i = 0; i < (int)nbSlots; i++) {
        v_N_1[i] = i;
        v_N_2[i] = 42;
    }

    {  // Setting up CKKS scheme with batching:

        std::vector<int> steps;
        steps.push_back(0);
        for (size_t i = 1; i < nbSlots; i++) {
            steps.push_back(i);
        }

        std::cout << "Start generating crypto context" << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(nbSlots, parameters);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        cc->EvalRotateKeyGen(keys.secretKey, steps);
        auto timeGenRotKey = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to generate rotation keys with batching (" << nbSlots << " slots) : " << timeGenRotKey
                  << " ms." << std::endl;

        // We create the initialisation
        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        Plaintext plain_v_2 = cc->MakeCKKSPackedPlaintext(v_2);
        Plaintext plain_v_zero = cc->MakeCKKSPackedPlaintext(v_zero);

        Plaintext plain_v_N_1 = cc->MakeCKKSPackedPlaintext(v_N_1);
        Plaintext plain_v_N_2 = cc->MakeCKKSPackedPlaintext(v_N_2);

        auto ct_v_2 = cc->Encrypt(keys.publicKey, plain_v_2);
        auto ct_v_N_1 = cc->Encrypt(keys.publicKey, plain_v_N_1);
        // auto ct_v_N_2 = cc->Encrypt(keys.publicKey, plain_v_N_2);
        auto ct_v_zero = cc->Encrypt(keys.publicKey, plain_v_zero);

        // Computation:

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            ct_v_zero = cc->EvalAdd(plain_v_1, ct_v_zero);

        auto timeAddCst = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " additions between 1 cipher and 1 plain : " << timeAddCst
                  << " ms." << std::endl;
        std::cout << "Avarage timeAddCst : " << timeAddCst / NB_TESTS << std::endl;

        ct_v_zero = cc->Encrypt(keys.publicKey, plain_v_zero);

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            ct_v_zero = cc->EvalAdd(ct_v_zero, ct_v_2);

        auto timeAdd = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " additions between 2 ciphers : " << timeAdd << " ms."
                  << std::endl;
        std::cout << "Avarage timeAdd : " << timeAdd / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            auto sol_rotat = cc->EvalRotate(ct_v_2, i % 10 + 1);

        auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " rotations : " << timeRot << " ms." << std::endl;
        std::cout << "Avarage timeRot : " << timeRot / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = cc->EvalMult(ct_v_2, 0.85);
        }

        auto timeMultCst = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in ]0,1[) : " << timeMultCst << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst 0.85 : " << timeMultCst / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = cc->EvalMult(ct_v_N_1, 42);
        }

        auto timeMultCstN = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in N) : " << timeMultCstN << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst 42 : " << timeMultCstN / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = cc->EvalMult(ct_v_N_1, plain_v_N_2);
        }

        auto timeMultPlain = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in Plaintext) : " << timeMultPlain << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst ckks : " << timeMultPlain / NB_TESTS << std::endl;
    }
    // Setting up CKKS scheme without batching:
    std::cout << "\n\n" << std::endl;

    {
        std::cout << "Start generating crypto context" << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        CryptoContext<DCRTPoly> ccWO = generate_common_crypto_context(parameters);
        auto ringDim = ccWO->GetRingDimension();

        auto keys = ccWO->KeyGen();
        ccWO->EvalMultKeyGen(keys.secretKey);

        std::vector<int> stepsWO;
        stepsWO.push_back(0);
        for (size_t i = 1; i < ringDim / 2; i++) {
            stepsWO.push_back(i);
        }

        auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        ccWO->EvalRotateKeyGen(keys.secretKey, stepsWO);
        auto timeGenRotKey = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to generate rotation keys without batching (" << ringDim << " slots) : " << timeGenRotKey
                  << " ms." << std::endl;

        // Creation of vectors, Plaintexts and Ciphertexts:

        // We create the initialisation
        auto plain_v_1 = ccWO->MakeCKKSPackedPlaintext(v_1);
        auto plain_v_2 = ccWO->MakeCKKSPackedPlaintext(v_2);
        auto plain_v_zero = ccWO->MakeCKKSPackedPlaintext(v_zero);

        auto plain_v_N_1 = ccWO->MakeCKKSPackedPlaintext(v_N_1);
        auto plain_v_N_2 = ccWO->MakeCKKSPackedPlaintext(v_N_2);

        auto ct_v_2 = ccWO->Encrypt(keys.publicKey, plain_v_2);
        auto ct_v_N_1 = ccWO->Encrypt(keys.publicKey, plain_v_N_1);
        auto ct_v_zero = ccWO->Encrypt(keys.publicKey, plain_v_zero);

        // Computation:

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            ct_v_zero = ccWO->EvalAdd(plain_v_1, ct_v_zero);

        auto timeAddCst = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " additions between 1 cipher and 1 plain : " << timeAddCst
                  << " ms." << std::endl;
        std::cout << "Avarage timeAddCst : " << timeAddCst / NB_TESTS << std::endl;

        ct_v_zero = ccWO->Encrypt(keys.publicKey, plain_v_zero);

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            ct_v_zero = ccWO->EvalAdd(ct_v_zero, ct_v_2);

        auto timeAdd = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " additions between 2 ciphers : " << timeAdd << " ms."
                  << std::endl;
        std::cout << "Avarage timeAdd : " << timeAdd / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++)
            auto sol_rotat = ccWO->EvalRotate(ct_v_2, i % 10 + 1);

        auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS << " rotations : " << timeRot << " ms." << std::endl;
        std::cout << "Avarage timeRot : " << timeRot / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = ccWO->EvalMult(ct_v_2, 0.85);
        }

        auto timeMultCst = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in ]0,1[) : " << timeMultCst << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst : " << timeMultCst / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = ccWO->EvalMult(ct_v_N_1, 42);
        }

        auto timeMultCstN = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in N) : " << timeMultCstN << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst : " << timeMultCstN / NB_TESTS << std::endl;

        begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        for (int i = 0; i < NB_TESTS; i++) {
            auto sol_mult = ccWO->EvalMult(ct_v_N_1, plain_v_N_2);
        }

        auto timeMultPlain = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
        std::cout << "  * Time to compute " << NB_TESTS
                  << " multiplications between 1 cipher and 1 plain (plain in Plaintext) : " << timeMultPlain << " ms."
                  << std::endl;
        std::cout << "Avarage timeMultCst : " << timeMultPlain / NB_TESTS << std::endl;
    }
}

void get_cipher_size() {
    std::vector<uint32_t> slots;
    for (int i = 0; i < 7; i++) {
        slots.push_back(256 << i);
    }
    for (auto nbSlots : slots) {
        std::vector<double> v_1(nbSlots);
        // Fill with random values
        for (int i = 0; i < (int)nbSlots; i++) {
            v_1[i] = rand() % 100;
        }
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(nbSlots, parameters);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);

        // Serialization
        std::stringstream ss;
        Serial::Serialize(ct_v_1, ss, SerType::BINARY);
        ss.seekp(0, std::ios::end);

        std::cout << "Size of the cipher, batchSize = " << nbSlots << ": " << ss.tellp() << " octets" << std::endl;
    }
}

void compare_first_and_following() {
    // Compare the time taken by the first addition and the following ones
    CCParams<CryptoContextCKKSRNS> parameters;
    auto cc = generate_common_crypto_context(256, parameters);
    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    std::vector<double> v_1(256);
    std::vector<double> v_2(256);
    std::vector<double> v_zero(256);
    for (int i = 0; i < 256; i++) {
        // Set random values
        v_1[i] = rand() % 100;
        v_2[i] = rand() % 100;
    }

    // Make plaintext
    auto plt_v1 = cc->MakeCKKSPackedPlaintext(v_1);
    auto plt_v2 = cc->MakeCKKSPackedPlaintext(v_2);
    auto ct_v1 = cc->Encrypt(keys.publicKey, plt_v1);
    auto ct_v2 = cc->Encrypt(keys.publicKey, plt_v2);

    Experiment expe;
    expe.setCryptoParams(parameters);
    expe.setRingsize(cc->GetRingDimension());
    {
        std::vector<long> times_first;
        std::vector<long> times_second;
        for (int i = 0; i < NB_TRY; i++) {
            // For each computation take times
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalAdd(ct_v1, ct_v2);
                auto timeAdd = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_first.push_back(timeAdd);
            }
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalAdd(ct_v1, ct_v2);
                auto timeAdd = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_second.push_back(timeAdd);
            }
        }
        expe.addValue(
            {std::accumulate(times_first.begin(), times_first.end(), 0) / NB_TRY},
            "First addition mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ADD);
        expe.addValue(
            {std::accumulate(times_second.begin(), times_second.end(), 0) / NB_TRY},
            "Following additions mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ADD);
    }
    std::cout << "End of addition\n";
    {
        std::vector<long> times_first;
        std::vector<long> times_second;
        for (int i = 0; i < NB_TRY; i++) {
            // For each computation take times
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalMult(ct_v1, 0.85);
                auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_first.push_back(timeMult);
            }
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalMult(ct_v1, 0.85);
                auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_second.push_back(timeMult);
            }
        }
        expe.addValue(
            {std::accumulate(times_first.begin(), times_first.end(), 0) / NB_TRY},
            "First multiplication mean on " + std::to_string(NB_TRY),
            Experiment::Operation::MUL);
        expe.addValue(
            {std::accumulate(times_second.begin(), times_second.end(), 0) / NB_TRY},
            "Following multiplications mean on " + std::to_string(NB_TRY),
            Experiment::Operation::MUL);
    }
    std::cout << "End of multiplication\n";
    cc->EvalRotateKeyGen(keys.secretKey, {1, 10});
    {
        std::vector<long> times_first;
        std::vector<long> times_second;
        for (int i = 0; i < NB_TRY; i++) {
            // For each computation take times
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalRotate(ct_v1, 1);
                auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_first.push_back(timeRot);
            }
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalRotate(ct_v1, 1);
                auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_second.push_back(timeRot);
            }
        }
        expe.addValue(
            {std::accumulate(times_first.begin(), times_first.end(), 0) / NB_TRY},
            "First rotation by 1 mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ROT);
        expe.addValue(
            {std::accumulate(times_second.begin(), times_second.end(), 0) / NB_TRY},
            "Following rotations by 1 mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ROT);
    }
    std::cout << "End of rotation\n";
    {
        std::vector<long> times_first;
        std::vector<long> times_second;
        for (int i = 0; i < NB_TRY; i++) {
            // For each computation take times
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalRotate(ct_v1, 10);
                auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_first.push_back(timeRot);
            }
            {
                auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
                auto ct_res = cc->EvalRotate(ct_v1, 10);
                auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
                times_second.push_back(timeRot);
            }
        }
        expe.addValue(
            {std::accumulate(times_first.begin(), times_first.end(), 0) / NB_TRY},
            "First rotation by 10 mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ROT);
        expe.addValue(
            {std::accumulate(times_second.begin(), times_second.end(), 0) / NB_TRY},
            "Following rotations by 10 mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ROT);
    }
    std::cout << "End of rotation\n";
    std::cout << expe << std::endl;
}

void test_all_operation(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<lbcrypto::DCRTPoly>& keys,
    Experiment& expe,
    Plaintext& plain_v_1) {
    {  // Time to encrypt
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);
            auto timeEncrypt = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeEncrypt, "Encrypt", Experiment::Operation::ENCRYPT, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of encryption\n";
    auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);

    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalMult(ct_v_1, 0.85);
            auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeMult, "Mult scalar", Experiment::Operation::MUL_SCALAR, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of multiplication by scalar\n";

    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalMult(ct_v_1, ct_v_1);
            auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeMult, "Mult", Experiment::Operation::MUL, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of multiplication\n";

    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalMult(plain_v_1, ct_v_1);
            auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeMult, "Mult plain", Experiment::Operation::MUL_PLAIN, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of plain multiplication\n";

    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalAdd(ct_v_1, ct_v_1);
            auto timeAdd = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeAdd, "Add", Experiment::Operation::ADD, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of addition\n";

    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalAdd(ct_v_1, 0.85);
            auto timeMult = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeMult, "Add scalar", Experiment::Operation::ADD_SCALAR, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of addition by scalar\n";

    cc->EvalRotateKeyGen(keys.secretKey, {1, 10});
    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalRotate(ct_v_1, 1);
            auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeRot, "Rotate of 1", Experiment::Operation::ROT, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of rotation by 1\n";
    {
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalRotate(ct_v_1, 10);
            auto timeRot = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeRot, "Rotate of 10", Experiment::Operation::ROT, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of rotation by 10\n";
    {  // Time to decrypt
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            Plaintext res;
            cc->Decrypt(ct_v_1, keys.secretKey, &res);
            auto timeDecrypt = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeDecrypt, "Decrypt mean", Experiment::Operation::DECRYPT, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of decryption\n";
    {
        auto precomp = cc->EvalFastRotationPrecompute(ct_v_1);
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalFastRotation(ct_v_1, 1, cc->GetRingDimension() * 2, precomp);
            auto timeSum = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin;
            expe.addValue(timeSum, "Fast rotation", Experiment::Operation::FAST_ROT, Experiment::Unit::MILLISECONDS);
        }
    }
    std::cout << "End of fast rotation by 1\n";
}

void compare_multiple_batch() {
    std::vector<uint32_t> slots;
    for (int i = 0; i < 7; i++) {
        slots.push_back(256 << i);
    }
    Experiment expe("multiple-batch.csv");
    for (auto nbSlots : slots) {
        std::vector<double> v_1(nbSlots);
        // Fill with random values
        for (int i = 0; i < (int)nbSlots; i++) {
            v_1[i] = rand() % 100;
        }
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(nbSlots, parameters);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        expe.setCryptoParams(parameters);
        expe.setRingsize(cc->GetRingDimension());

        test_all_operation(cc, keys, expe, plain_v_1);
    }
}

void compare_multiple_ring() {
    Experiment expe("multiple-ring-same-depth.csv");
    for (int ringDim = 2 << 13; ringDim <= 131072; ringDim *= 2) {
        std::vector<double> v_1(256);
        // Fill with random values
        for (int i = 0; i < (int)256; i++) {
            v_1[i] = rand() % 100;
        }
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(0, parameters, 4, ringDim);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        expe.setCryptoParams(parameters);
        expe.setRingsize(cc->GetRingDimension());

        test_all_operation(cc, keys, expe, plain_v_1);
    }
}

void compare_multiple_depth() {
    Experiment expe("multiple-depth-same-ring-131k-10it-bits30-first35.csv");
    for (int multDepth = 25; multDepth < 65; multDepth += 5) {
        std::vector<double> v_1(256);
        // Fill with random values
        for (int i = 0; i < (int)256; i++) {
            v_1[i] = rand() % 100;
        }
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(256, parameters, multDepth, 1 << 17);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        expe.setCryptoParams(parameters);
        expe.setRingsize(cc->GetRingDimension());

        test_all_operation(cc, keys, expe, plain_v_1);
        cc->ClearEvalAutomorphismKeys();
        cc->ClearEvalMultKeys();
        cc->ClearEvalSumKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
    }
}

void compare_ring_dim() {
    Experiment expe("ring-size-bits30-first35.csv");
    for (int multDepth = 8; multDepth < 64; multDepth++) {
        std::cout << multDepth << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;
        auto duration = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        CryptoContext<DCRTPoly> cc;
        cc = generate_common_crypto_context(256, parameters, multDepth);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        auto timeGen = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - duration;
        expe.setCryptoParams(parameters);
        expe.setRingsize(cc->GetRingDimension());
        expe.addValue(
            timeGen, "Time to generate crypto context", Experiment::Operation::GEN, Experiment::Unit::MILLISECONDS);

        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
    }
}

void compare_rotation(int nbSlots, int ringDim) {
    Experiment expe("rotation" + std::to_string(nbSlots) + "-" + std::to_string(ringDim) + ".csv");
    // Generate crypto context
    CCParams<CryptoContextCKKSRNS> parameters;
    auto cc = generate_common_crypto_context(nbSlots, parameters, 16, ringDim);
    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    expe.setCryptoParams(parameters);
    expe.setRingsize(cc->GetRingDimension());
    std::vector<int> slots;
    for (int i = 1; i <= cc->GetRingDimension() / 2; i *= 2) {
        slots.push_back(i);
    }
    cc->EvalRotateKeyGen(keys.secretKey, slots);

    std::vector<double> v_1(nbSlots);
    // Fill with random values
    for (int i = 0; i < (int)(nbSlots); i++) {
        v_1[i] = rand() % 100;
    }
    Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
    auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);
    for (auto slot : slots) {
        std::vector<long> times;
        for (int i = 0; i < NB_TRY; i++) {
            auto begin = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
            auto ct_res = cc->EvalRotate(ct_v_1, slot);
            auto timeRot = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count() - begin;
            times.push_back(timeRot);
        }
        expe.addValue(
            std::accumulate(times.begin(), times.end(), 0) / NB_TRY,
            "Rotate of " + std::to_string(slot) + " mean on " + std::to_string(NB_TRY),
            Experiment::Operation::ROT,
            Experiment::Unit::MICROSECONDS);
    }
}

void try_to_explode() {
    Experiment expe("time-explode.csv");

    for (int multDepth = 8; multDepth <= 32; multDepth *= 2) {
        CCParams<CryptoContextCKKSRNS> parameters;
        auto cc = generate_common_crypto_context(512, parameters, multDepth);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        expe.setCryptoParams(parameters);
        expe.setRingsize(cc->GetRingDimension());
        int number = 0;
        std::vector<double> v_1(512);
        // Fill with random values
        for (int i = 0; i < (int)512; i++) {
            v_1[i] = rand() % 100;
        }
        Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
        // Encrypt
        auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);
        std::vector<double> v_2(512, 1.0);
        auto ct_v_2 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(v_2));
        std::cout << "Original value " << v_1 << std::endl;
        try {
            while (true) {
                auto a = cc->EvalMult(ct_v_1, ct_v_2);
                auto b = cc->EvalMult(ct_v_1, ct_v_2);
                ct_v_1 = cc->EvalAdd(a, b);
                number++;
            }
        } catch (const std::exception& e) {
            std::cout << "Number of multiplications : " << number << std::endl;
            Plaintext res;
            cc->Decrypt(keys.secretKey, ct_v_1, &res);
            std::cout << "value " << res->GetCKKSPackedValue() << std::endl;
            std::cerr << e.what() << '\n';
        }
    }
}

void try_to_explode_rot() {
    Experiment expe("time-explode.csv");
    int multDepth = 2;
    CCParams<CryptoContextCKKSRNS> parameters;
#define BATCH_SIZE 8
    auto cc = generate_common_crypto_context(BATCH_SIZE, parameters, multDepth);
    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    cc->EvalRotateKeyGen(keys.secretKey, {1});
    expe.setCryptoParams(parameters);
    expe.setRingsize(cc->GetRingDimension());
    int number = 0;
    std::vector<double> v_1(BATCH_SIZE);
    // Fill with random values
    for (int i = 0; i < (int)BATCH_SIZE; i++) {
        v_1[i] = i;
    }
    Plaintext plain_v_1 = cc->MakeCKKSPackedPlaintext(v_1);
    // Encrypt
    auto ct_v_1 = cc->Encrypt(keys.publicKey, plain_v_1);
    std::cout << "Original value " << v_1 << std::endl;
    try {
        while (true) {
            ct_v_1 = cc->EvalRotate(ct_v_1, 1);
            number++;
            if (number % 100 == 0) {
                Plaintext res;
                cc->Decrypt(keys.secretKey, ct_v_1, &res);
                std::cout << "value: " << res->GetRealPackedValue() << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cout << "Number of Rotation : " << number << std::endl;
        std::cerr << e.what() << '\n';
    }
}

void generate_and_compare_partition() {
    auto g = read_graph("../expes_papier/synthetic-10000-nodes_2.txt", " ", true);  // 185 nodes
    // auto g = read_graph("../Graphs/email-Eu-core.txt", " ", true);                 // 1005 nodes
    // auto g = read_graph("../Graphs/p2p-Gnutella04.txt", "	", true);              // 10876 nodes
    // auto g = read_graph("../Graphs/Slashdot0811.txt", "	", true);                  // 77360 nodes
    // auto g = read_graph("../Graphs/Amazon0302.txt", "	", true);                  // 262111 nodes
    // auto g = read_graph("../Graphs/elliptic_bitcoin_transaction.csv", ",", true);  // 203,769 nodes

    // auto g = read_graph("../Graphs/Wiki-Vote.txt", " ", true);

    std::map<int, int> partition1 = split_uniformaly_ascendant(g, 5);
    std::map<int, int> partition2 = split_uniformaly_random(g, 5);
    std::map<int, int> partition3 = split_following_stats_uniformaly_random(g, {0.43, 0.31, 0.11, 0.1, 0.05});
    std::map<int, int> partition4 = split_following_stats_ascendant(g, {0.43, 0.31, 0.11, 0.1, 0.05});
    /*
        std::cout << "Partition 1" << std::endl;
        auto stats1 = get_edges_stats(g, partition1, false);
        show_stats(stats1);
        auto graphs1 = split_following_partition(g, partition1);
        export_to_dot(graphs1, "../dot/email-Eu-core-partition1.dot", partition1);

        std::cout << "Partition 2" << std::endl;
        auto stats2 = get_edges_stats(g, partition2, false);
        show_stats(stats2);
        auto graphs2 = split_following_partition(g, partition2);
        export_to_dot(graphs2, "../dot/email-Eu-core-partition2.dot", partition2);

        std::cout << "Partition 3" << std::endl;
        auto stats3 = get_edges_stats(g, partition3, false);
        show_stats(stats3);
        auto graphs3 = split_following_partition(g, partition3);
        export_to_dot(graphs3, "../dot/email-Eu-core-partition3.dot", partition3);
    */
    std::cout << "Partition 4" << std::endl;
    auto stats4 = get_edges_stats(g, partition4, false);
    show_stats(stats4);
    auto graphs4 = split_following_partition(g, partition4);
    export_to_dot(graphs4, "../expes_papier/synthetic-10000-nodes_2-partition4.dot", partition4);
}

void compare_metis_part() {
    Graph g = read_graph("../Graphs/students.csv", ",", false);
    std::vector<int> parts = {
        0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
        0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0,
        0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1,
        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0};
    std::map<int, int> partition;
    for (int i = 0; i < (int)parts.size(); i++) {
        partition[i] = parts[i];
    }
    auto stats = get_edges_stats(g, partition, false);
    show_stats(stats);
    auto graphs = split_following_partition(g, partition);
    export_to_dot(graphs, "../dot/partitionMetis.dot", partition);
}
/*
Partition 1
[
[ 15, 15, 14, 4, 14, ],
[ 15, 6, 20, 14, 12, ],
[ 13, 11, 19, 19, 19, ],
[ 17, 7, 15, 12, 18, ],
[ 11, 16, 17, 16, 12, ],
]
Partition 2
[
[ 14, 17, 15, 12, 12, ],
[ 10, 28, 13, 10, 25, ],
[ 16, 9, 13, 17, 16, ],
[ 12, 7, 16, 13, 9, ],
[ 13, 18, 14, 7, 15, ],
]
Partition 3
[
[ 0, 14, 11, 4, 3, ],
[ 18, 73, 48, 16, 12, ],
[ 10, 39, 22, 7, 9, ],
[ 3, 15, 6, 3, 4, ],
[ 1, 20, 10, 2, 1, ],
]
*/

void test_partition() {
    Graph g = read_graph("../Graphs/students.csv", ",", false);
    Graph gbis;
    for (auto node : g.nodes) {
        gbis.nodes.insert(node - 1);
    }
    for (auto arc : g.arcs) {
        gbis.arcs.insert({arc.first - 1, arc.second - 1});
    }
    std::map<int, int> partition;
    std::vector<Graph> graphs = parse_graph_multi_DC("../Graphs/multi_DC_students_intercom-5.dot", partition);

    for (auto subG : graphs) {
        for (auto node : subG.nodes) {
            assert(gbis.nodes.find(node) != gbis.nodes.end());
        }
        for (auto arc : subG.arcs) {
            assert(gbis.arcs.find(arc) != gbis.arcs.end());
        }
    }
}

void test_bootstrap() {
    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_128_classic);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits = 59;
    usint firstMod = 60;
    int numSlots = 8;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    parameters.SetBatchSize(numSlots);

    // < ceil(log2(slots))
    std::vector<uint32_t> levelBudget = {3, 3};
    // Auto params
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = 15;
    auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    usint depth = levelsAvailableAfterBootstrap + levelsNeeded;
    std::cout << "Level needed for bootstrap " << levelsNeeded << std::endl;
    parameters.SetMultiplicativeDepth(depth);
    std::cout << "Mult Depth " << depth << std::endl;

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    usint ringDim = cryptoContext->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

    cryptoContext->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots);

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    // Generate bootstrapping keys.
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    std::vector<double> v_1;
    for (int i = 0; i < numSlots; i++) {
        v_1.push_back(i);
    }
    Plaintext plain_v_1 = cryptoContext->MakeCKKSPackedPlaintext(v_1);

    auto ct_1 = cryptoContext->Encrypt(keyPair.publicKey, plain_v_1);

    auto ct_res = cryptoContext->EvalMult(ct_1, ct_1);
    auto ct_res2 = cryptoContext->EvalMult(ct_res, ct_res);
    auto ct_res3 = cryptoContext->EvalMult(ct_res2, ct_res2);  // 21
    auto ct_res4 = cryptoContext->EvalMult(ct_res3, ct_res3);
    auto ct_res5 = cryptoContext->EvalMult(ct_res4, ct_res4);
    auto ct_res6 = cryptoContext->EvalMult(ct_res5, ct_res5);  // 18
    auto ct_res7 = cryptoContext->EvalMult(ct_res6, ct_res6);
    auto ct_res8 = cryptoContext->EvalMult(ct_res7, ct_res7);  // 16
    auto ct_res9 = cryptoContext->EvalMult(ct_res8, ct_res8);
    auto ct_res10 = cryptoContext->EvalMult(ct_res9, ct_res9);
    auto ct_res11 = cryptoContext->EvalMult(ct_res10, ct_res10);
    auto ct_res12 = cryptoContext->EvalMult(ct_res11, ct_res11);
    auto ct_res13 = cryptoContext->EvalMult(ct_res12, ct_res12);
    auto ct_res14 = cryptoContext->EvalMult(ct_res13, ct_res13);
    auto ct_res15 = cryptoContext->EvalMult(ct_res14, ct_res14);
    auto ct_res16 = cryptoContext->EvalMult(ct_res15, ct_res15);
    auto ct_res17 = cryptoContext->EvalMult(ct_res16, ct_res16);
    auto ct_res18 = cryptoContext->EvalMult(ct_res17, ct_res17);
    auto ct_res19 = cryptoContext->EvalMult(ct_res18, ct_res18);  // 5
    auto ct_res20 = cryptoContext->EvalMult(ct_res19, ct_res19);
    auto ct_res21 = cryptoContext->EvalMult(ct_res20, ct_res20);
    auto ct_res22 = cryptoContext->EvalMult(ct_res21, ct_res21);
    auto ct_res23 = cryptoContext->EvalMult(ct_res22, ct_res22);

    std::cout << "Before bootstrap level remaining " << parameters.GetMultiplicativeDepth() - ct_res22->GetLevel()
              << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto bootstrapped = cryptoContext->EvalBootstrap(ct_res22);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time to bootstrap " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
    std::cout << "After bootstrap level remaining " << parameters.GetMultiplicativeDepth() - bootstrapped->GetLevel()
              << std::endl;

    auto bootsrapped1 = cryptoContext->EvalMult(bootstrapped, bootstrapped);
    auto bootsrapped2 = cryptoContext->EvalMult(bootsrapped1, bootsrapped1);
    auto bootsrapped3 = cryptoContext->EvalMult(bootsrapped2, bootsrapped2);
    auto bootsrapped4 = cryptoContext->EvalMult(bootsrapped3, bootsrapped3);

    std::cout << "Before bootstrap2 level remaining " << parameters.GetMultiplicativeDepth() - bootsrapped4->GetLevel()
              << std::endl;
    auto bootstrapped2times = cryptoContext->EvalBootstrap(bootsrapped4);
    std::cout << "After bootstrap2 level remaining "
              << parameters.GetMultiplicativeDepth() - bootstrapped2times->GetLevel() << std::endl;
}

void compare_bootstrap() {
    Experiment expe("compare-bootstrap.csv");
    for (int i = 5; i < 45; i++) {
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetSecurityLevel(HEStd_128_classic);

        ScalingTechnique rescaleTech = FLEXIBLEAUTO;
        usint dcrtBits = 59;
        usint firstMod = 60;
        int numSlots = 8;

        parameters.SetScalingModSize(dcrtBits);
        parameters.SetScalingTechnique(rescaleTech);
        parameters.SetFirstModSize(firstMod);
        parameters.SetBatchSize(numSlots);

        // < ceil(log2(slots))
        std::vector<uint32_t> levelBudget = {3, 3};
        // Auto params
        std::vector<uint32_t> bsgsDim = {0, 0};

        uint32_t levelsAvailableAfterBootstrap = i;
        auto levelsNeeded = FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        usint depth = levelsAvailableAfterBootstrap + levelsNeeded;
        std::cout << "Level needed for bootstrap " << levelsNeeded << std::endl;
        parameters.SetMultiplicativeDepth(depth);
        std::cout << "Mult Depth " << depth << std::endl;

        CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
        expe.setCryptoParams(parameters);
        expe.setRingsize(cryptoContext->GetRingDimension());

        cryptoContext->Enable(PKE);
        cryptoContext->Enable(KEYSWITCH);
        cryptoContext->Enable(LEVELEDSHE);
        cryptoContext->Enable(ADVANCEDSHE);
        cryptoContext->Enable(FHE);

        usint ringDim = cryptoContext->GetRingDimension();
        std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

        cryptoContext->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots);

        auto keyPair = cryptoContext->KeyGen();
        cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        // Generate bootstrapping keys.
        cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

        std::vector<double> v_1;
        for (int i = 0; i < numSlots; i++) {
            v_1.push_back(i);
        }
        Plaintext plain_v_1 = cryptoContext->MakeCKKSPackedPlaintext(v_1, 1, parameters.GetMultiplicativeDepth() - 1);
        auto ct_1 = cryptoContext->Encrypt(keyPair.publicKey, plain_v_1);
        auto start = std::chrono::high_resolution_clock::now();
        auto ciphertextAfter = cryptoContext->EvalBootstrap(ct_1);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        expe.addValue(
            duration,
            "Time to bootstrap for " + std::to_string(i) + " levels",
            Experiment::Operation::BOOTSTRAP,
            Experiment::Unit::SECONDS);
        cryptoContext->ClearEvalAutomorphismKeys();
        cryptoContext->ClearEvalMultKeys();
        cryptoContext->ClearEvalSumKeys();
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
    }
}

void test_fast_rotation() {
    CCParams<CryptoContextCKKSRNS> parameters;
    auto cc = generate_common_crypto_context(256, parameters, 8);
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);

    cc->EvalRotateKeyGen(keyPair.secretKey, {10, -10});

    std::vector<double> vect1;
    for (int i = 0; i < 256; i++) {
        vect1.push_back(i);
    }
    auto plt_1 = cc->MakeCKKSPackedPlaintext(vect1);
    auto ct_1 = cc->Encrypt(keyPair.publicKey, plt_1);

    auto precomp = cc->EvalFastRotationPrecompute(ct_1);
    std::vector<Ciphertext<DCRTPoly>> res1;
    auto begin1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    for (int i = 0; i < NB_TRY; i++) {
        res1.push_back(cc->EvalFastRotation(ct_1, 10, cc->GetRingDimension() * 2, precomp));
    }
    // Sum all res
    Ciphertext<DCRTPoly> sum = res1[0];
    for (int i = 1; i < NB_TRY; i++) {
        sum = cc->EvalAdd(sum, res1[i]);
    }
    auto timeRot1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin1;
    std::cout << "Time to rotate by 10 " << timeRot1 / NB_TRY << " ms" << std::endl;

    auto begin2 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    std::vector<Ciphertext<DCRTPoly>> res2;
    for (int i = 0; i < NB_TRY; i++) {
        res2.push_back(cc->EvalFastRotation(ct_1, -10, cc->GetRingDimension() * 2, precomp));
    }
    // Sum all res
    Ciphertext<DCRTPoly> sum2 = res2[0];
    for (int i = 1; i < NB_TRY; i++) {
        sum2 = cc->EvalAdd(sum2, res2[i]);
    }
    auto timeRot2 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin2;
    std::cout << "Time to rotate by -10 " << timeRot2 / NB_TRY << " ms" << std::endl;

    auto begin3 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    std::vector<Ciphertext<DCRTPoly>> res3;
    for (int i = 0; i < NB_TRY; i++) {
        res3.push_back(cc->EvalFastRotationExt(ct_1, 10, precomp, true));
    }
    // Sum all res
    Ciphertext<DCRTPoly> sum3 = res3[0];
    for (int i = 1; i < NB_TRY; i++) {
        sum3 = cc->EvalAdd(sum3, res3[i]);
    }
    sum3 = cc->KeySwitchDown(sum3);
    auto timeRot3 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin3;
    std::cout << "Time to rotate by 10 with keyswitch " << timeRot3 / NB_TRY << " ms" << std::endl;

    auto begin4 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    std::vector<Ciphertext<DCRTPoly>> res4;
    for (int i = 0; i < NB_TRY; i++) {
        res4.push_back(cc->EvalFastRotationExt(ct_1, -10, precomp, true));
    }
    // Sum all res
    Ciphertext<DCRTPoly> sum4 = res4[0];
    for (int i = 1; i < NB_TRY; i++) {
        sum4 = cc->EvalAdd(sum4, res4[i]);
    }
    sum4 = cc->KeySwitchDown(sum4);
    auto timeRot4 = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - begin4;
    std::cout << "Time to rotate by -10 with keyswitch " << timeRot4 / NB_TRY << " ms" << std::endl;
    Plaintext plt_res1;
    cc->Decrypt(keyPair.secretKey, sum, &plt_res1);
    Plaintext plt_res2;
    cc->Decrypt(keyPair.secretKey, sum2, &plt_res2);
    Plaintext plt_res3;
    cc->Decrypt(keyPair.secretKey, sum3, &plt_res3);
    Plaintext plt_res4;
    cc->Decrypt(keyPair.secretKey, sum4, &plt_res4);
    plt_res1->SetLength(256);
    plt_res2->SetLength(256);
    plt_res3->SetLength(256);
    plt_res4->SetLength(256);
    std::cout << "Result of rotation by 10 " << plt_res1->GetCKKSPackedValue() << std::endl;
    std::cout << "Result of rotation by 10 with keyswitch " << plt_res3->GetCKKSPackedValue() << std::endl;
    std::cout << "Result of rotation by -10 " << plt_res2->GetCKKSPackedValue() << std::endl;
    std::cout << "Result of rotation by -10 with keyswitch " << plt_res4->GetCKKSPackedValue() << std::endl;
}

void test_nauty() {
    // auto graph = read_graph("../Graphs/synthetic-graph-2.txt", ",", true);
    auto graph = read_graph("../Graphs/students.csv", ",", false);  // 185 nodes
    // auto graph = read_graph("../Graphs/email-Eu-core.txt", " ", true);               // 1005 nodes
    //  auto graph = read_graph("../Graphs/Slashdot0811.txt", "	", true);               // 77360 nodes
    // auto graph = read_graph("../Graphs/Amazon0302.txt", "	", true);  // 262111 nodes
    // auto graph = read_graph("../Graphs/web-NotreDame.txt", " ", true);               // 325729 nodes
    // auto graph = read_graph("../Graphs/soc-pokec-relationships.txt", "	", true);   // 1632803 nodes

    std::map<int, int> oneg_part;
    for (auto node : graph.nodes) {
        oneg_part[node] = 0;
    }
    int number_of_nodes;
    int nb_local_nodes_after_pruning;
    std::map<int, int> orbits_one;
    int duration;
    std::set<Node> set_dangling_nodes;
    sparse_nauty_routine(
        &graph, &oneg_part, number_of_nodes, nb_local_nodes_after_pruning, orbits_one, 0, duration, set_dangling_nodes);
    std::cout << "Time to compute orbits for graph " << duration << " s" << std::endl;
    std::set<Orbit> unique_orbits_one;
    for (auto orbit : orbits_one) {
        unique_orbits_one.insert(orbit.second);
    }

    std::cout << "Number of nodes for graph " << graph.nodes.size() << std::endl;
    std::cout << "Number of orbits for graph " << unique_orbits_one.size() << std::endl;
    std::cout << "Pruned nodes " << graph.nodes.size() - unique_orbits_one.size() << std::endl;
}

void test_nauty_with_partition() {
    // auto graph = read_graph("../Graphs/students.csv", ",", false);
    // auto graph = read_graph("../Graphs/Slashdot0811.txt", "	", true);
    auto graph = read_graph("../Graphs/email-Eu-core.txt", " ", false);
    // auto graph = read_graph("../Graphs/Amazon0302.txt", "	", true);
    // auto graph = read_graph("../Graphs/web-NotreDame.txt", " ", true);
    // auto graph = read_graph("../Graphs/elliptic_bitcoin_transaction.csv", ",", true);
    std::map<int, int> partition;
    // auto graphs = parse_graph_multi_DC("../Graphs/partition4.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/Slashdot0811-partition4.dot", partition);
    auto graphs = parse_graph_multi_DC("../dot/email-Eu-core-partition4.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/Amazon0302-partition4.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/web-NotreDame-partition3.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/elliptic_bitcoin_transaction-partition4.dot", partition);

    std::map<int, int> full_g_part;
    for (auto node : graph.nodes) {
        // For the full graph, each node is in the same orbit
        full_g_part[node] = 0;
    }
    int number_of_nodes;
    int nb_local_nodes_after_pruning;
    std::map<int, int> orbits_one;
    int time_full_graph;
    std::set<Node> set_dangling_nodes;
    sparse_nauty_routine(
        &graph,
        &full_g_part,
        number_of_nodes,
        nb_local_nodes_after_pruning,
        orbits_one,
        0,
        time_full_graph,
        set_dangling_nodes);

    std::set<Orbit> unique_orbits_one;
    for (auto orbit : orbits_one) {
        std::cout << orbit.first << " - > " << orbit.second << std::endl;
        unique_orbits_one.insert(orbit.second);
    }

    std::vector<std::map<Node, Orbit>> orbits_part;
    std::vector<int> times;
    int i = 0;
    for (auto part : graphs) {
        std::map<Node, Orbit> orbits;
        int number_of_nodes;
        int nb_local_nodes_after_pruning;
        int duration;
        std::set<Node> subset_dangling_nodes;
        sparse_nauty_routine(&part, &partition, number_of_nodes, nb_local_nodes_after_pruning, orbits, i, duration, subset_dangling_nodes);
        times.push_back(duration);
        orbits_part.push_back(orbits);

        i++;
    }
    std::cout << "nb de noeuds " << graph.nodes.size() << std::endl;

    std::cout << "Time to compute orbits for full graph " << time_full_graph << " s" << std::endl;
    std::cout << "Number of orbits for full graph " << unique_orbits_one.size() << std::endl;
    std::cout << "=> Pruned nodes " << graph.nodes.size() - unique_orbits_one.size() << std::endl;
    int part_orbits_number = 0;
    std::vector<int> nb_orbits_per_orbits;
    for (int i = 0; i < orbits_part.size(); i++) {
        std::set<Orbit> unique_orbits;
        for (auto orbit : orbits_part[i]) {
            if (partition[orbit.first] == i) {
                unique_orbits.insert(orbit.second);
            }
        }
        part_orbits_number += unique_orbits.size();
        nb_orbits_per_orbits.push_back(unique_orbits.size());
    }
    std::cout << "Number of orbits for all the partitions " << part_orbits_number << std::endl;
    std::cout << "=> Pruned nodes " << graph.nodes.size() - part_orbits_number << std::endl;
    std::cout << "times: " << times << std::endl;

    std::vector<int> partitions_size;
    std::vector<int> local_nodes;
    std::vector<int> border_nodes;
    int datacenter_id = 0;
    for (auto part : graphs) {
        partitions_size.push_back(part.nodes.size());

        int index_local = 0;
        for (auto node : part.nodes) {
            if (partition[node] == datacenter_id)
                index_local++;
        }
        local_nodes.push_back(index_local);

        std::set<int> index_border;
        for (auto arc : part.arcs) {
            if (partition[arc.first] != datacenter_id && partition[arc.second] == datacenter_id)
                index_border.insert(arc.second);
            if (partition[arc.first] == datacenter_id && partition[arc.second] != datacenter_id)
                index_border.insert(arc.first);
        }
        border_nodes.push_back(index_border.size());

        datacenter_id++;
    }

    std::cout << "Partitions_size: " << partitions_size << std::endl;
    std::cout << "Only local nodes: " << local_nodes << std::endl;
    std::cout << "Border local nodes: " << border_nodes << std::endl;

    std::cout << "Number nodes after pruning by orbits: " << nb_orbits_per_orbits << std::endl;
}

void test_constant_dangling_nodes() {
    // auto graph = read_graph("../Graphs/students.csv", ",", false);
    // auto graph = read_graph("../Graphs/Slashdot0811.txt", "	", true);
    auto graph = read_graph("../Graphs/email-Eu-core.txt", " ", false);
    // auto graph = read_graph("../Graphs/Amazon0302.txt", "	", true);
    std::map<int, int> partition;
    // auto graphs = parse_graph_multi_DC("../Graphs/partition3.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/Slashdot0811-partition3.dot", partition);
    auto graphs = parse_graph_multi_DC("../dot/email-Eu-core-partition3.dot", partition);
    // auto graphs = parse_graph_multi_DC("../dot/Amazon0302-partition3.dot", partition);

    std::multimap<Node, Node> map_out;
    std::multimap<Node, Node> map_in;
    for (auto arc : graph.arcs) {
        map_in.insert({arc.second, arc.first});
        map_out.insert(arc);
    }
    int nb_cst_nodes = 0;
    int nb_dgl_nodes = 0;
    for (auto node : graph.nodes) {
        if (map_out.count(node) == 0)
            nb_dgl_nodes += 1;
        if (map_in.count(node) == 0)
            nb_cst_nodes += 1;
    }
    std::cout << "Nb nodes : " << graph.nodes.size() << std::endl;
    std::cout << "Nb constant nodes : " << nb_cst_nodes << std::endl;
    std::cout << "Nb dangling nodes : " << nb_dgl_nodes << std::endl;

    auto cst_nodes = number_constant_nodes(graphs, partition);
    auto dgl_nodes = number_dangling_nodes(graphs, partition);

    std::vector<int> partitions_size;
    for (auto part : graphs) {
        partitions_size.push_back(part.nodes.size());
    }

    std::cout << "Partitions_size : " << partitions_size << std::endl;
    std::cout << "Nb constant nodes : " << cst_nodes << std::endl;
    std::cout << "Nb dangling nodes : " << dgl_nodes << std::endl;
}

void test_fast_rotation_matmult() {
    Experiment expe("compare-matmult-fastRot.csv");

    CCParams<CryptoContextCKKSRNS> parameters;
    auto cc = generate_common_crypto_context(0, parameters, 1);
    expe.setRingsize(cc->GetRingDimension());
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::vector<int> steps;
    steps.push_back(0);
    for (int i = 1; i < 2048; i++) {
        steps.push_back(i);
    }
    cc->EvalRotateKeyGen(keyPair.secretKey, steps);

    std::cout << "cryptoContext generated" << std::endl;

    for (int batchSize = 16; batchSize <= 2048; batchSize = batchSize << 1) {
        parameters.SetBatchSize(batchSize);
        expe.setCryptoParams(parameters);

        std::cout << parameters.GetBatchSize() << std::endl;
        std::cout << "batchSize: " << batchSize << std::endl;

        std::vector<double> vect;
        for (int i = 0; i < batchSize; i++) {
            vect.push_back(i);
        }
        auto plt_vect = cc->MakeCKKSPackedPlaintext(vect, 1, 0, nullptr, batchSize);
        auto ct_vect = cc->Encrypt(keyPair.publicKey, plt_vect);

        std::vector<std::vector<double>> mat(batchSize);
        for (int i = 0; i < batchSize; i++) {
            for (int j = 0; j < batchSize; j++) {
                mat[i].push_back(i * batchSize + j);
            }
        }
        std::vector<Plaintext> diags = PackedDiagonals(cc, mat, batchSize);

        std::vector<std::vector<long>> result(4);
        for (int nb_try = 0; nb_try < NB_TRY; nb_try++) {
            result[0].push_back(matmult_HS(cc, diags, ct_vect));
            result[1].push_back(matmult_HS_FastRot(cc, diags, ct_vect));
            result[2].push_back(matmult_Coeus(cc, diags, ct_vect));
            result[3].push_back(matmult_Coeus_FastRot(cc, diags, ct_vect));
        }

        expe.addValue(
            std::accumulate(result[0].begin(), result[0].end(), 0) / NB_TRY,
            "Matrix multiplication with HaleviShoup mean on " + std::to_string(NB_TRY),
            Experiment::Operation::EVAL_OP,
            Experiment::Unit::MILLISECONDS);
        expe.addValue(
            std::accumulate(result[1].begin(), result[1].end(), 0) / NB_TRY,
            "Matrix multiplication with HaleviShoup and FastRot mean on " + std::to_string(NB_TRY),
            Experiment::Operation::EVAL_OP,
            Experiment::Unit::MILLISECONDS);
        expe.addValue(
            std::accumulate(result[2].begin(), result[2].end(), 0) / NB_TRY,
            "Matrix multiplication with Coeus mean on " + std::to_string(NB_TRY),
            Experiment::Operation::EVAL_OP,
            Experiment::Unit::MILLISECONDS);
        expe.addValue(
            std::accumulate(result[3].begin(), result[3].end(), 0) / NB_TRY,
            "Matrix multiplication with Coeus and FastRot mean on " + std::to_string(NB_TRY),
            Experiment::Operation::EVAL_OP,
            Experiment::Unit::MILLISECONDS);
    }
}
#include "rcm.hpp"

void create_adj_unsym(
    Graph g,
    std::multimap<Node, Node> map_out,
    std::vector<int>& adj_row,
    std::vector<int>& adj,
    int node_num,
    int adj_max,
    int& adj_num,
    std::set<Node> activeNodes) {
    for (auto node : g.nodes) {
        if (activeNodes.count(node) == 0) {
            continue;
        }
        auto range = map_out.equal_range(node);
        for (auto it = range.first; it != range.second; ++it) {
            if (activeNodes.count(it->second) != 0) {
                adj_set(node_num, adj_max, &adj_num, adj_row.data(), adj.data(), node + 1, it->second + 1);
            }
        }
    }
}

void create_adj_sym(
    Graph g,
    std::multimap<Node, Node> map_out,
    std::multimap<Node, Node> map_in,
    std::vector<int>& adj_row,
    std::vector<int>& adj,
    int node_num,
    int adj_max,
    int& adj_num,
    std::set<Node> activeNodes) {
    for (auto node : g.nodes) {
        if (activeNodes.count(node) == 0) {
            continue;
        }
        auto range = map_out.equal_range(node);
        for (auto it = range.first; it != range.second; ++it) {
            if (activeNodes.count(it->second) != 0) {
                adj_set(node_num, adj_max, &adj_num, adj_row.data(), adj.data(), node + 1, it->second + 1);
            }
        }
        auto range_in = map_in.equal_range(node);
        for (auto it = range_in.first; it != range_in.second; ++it) {
            if (activeNodes.count(it->second) != 0) {
                adj_set(node_num, adj_max, &adj_num, adj_row.data(), adj.data(), node + 1, it->second + 1);
            }
        }
    }
}

int number_null_diags(int node_num, std::vector<int> adj_row, std::vector<int> adj) {
    // Count number of null diags
    int diag_null = 0;
    for (int i = 0; i < node_num; i++) {
        bool found = false;
        for (int j = 0; j < node_num; j++) {
            for (int k = adj_row[(i + j) % node_num]; k < adj_row[(i + j) % node_num + 1]; k++) {
                if (adj[k] == i) {
                    found = true;
                    break;
                }
            }
            if (found) {
                diag_null++;
                break;
            }
        }
    }
    return diag_null;
}

void test_mat(Graph g) {
    std::vector<std::vector<int>> mat;
    for (int i = 0; i < g.nodes.size(); i++) {
        std::vector<int> row(g.nodes.size(), 0);
        mat.push_back(row);
    }
    for (auto arc : g.arcs) {
        mat[arc.first][arc.second] = 1;
    }
    int diag_null = 0;
    for (int i = 0; i < mat.size(); i++) {
        std::vector<double> diag;
        for (int j = 0; j < mat.size(); j++) {
            diag.push_back(mat[j][(i + j) % mat.size()]);
        }
        if (std::accumulate(diag.begin(), diag.end(), 0) == 0) {
            diag_null++;
        }
    }
    std::cout << "Number of null diags " << diag_null << std::endl;
    std::cout << "Percent of null diags " << (double)diag_null / (double)(2 * g.nodes.size()) * 100 << " %"
              << std::endl;
}

/* Cuthill-Mckee, an Algorithm for Reducing the Bandwidth of the matrix */
void test_rcm() {
    // Graph g = read_graph("../Graphs/email-Eu-core.txt", " ", true);
    Graph g = read_graph("../Graphs/Email-EuAll.txt", "	", true);
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, int ADJ_NUM, the number of adjacency entries.
    //
    //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    //    Input, int ADJ[ADJ_NUM], the adjacency structure.
    //    For each row, it contains the column indices of the nonzero entries.
    std::multimap<Node, Node> map_in;
    std::multimap<Node, Node> map_out;
    std::set<Node> activeNodes;
    for (auto arc : g.arcs) {
        map_out.insert(arc);
        map_in.insert({arc.second, arc.first});
    }
    for (auto node : g.nodes) {
        if (map_out.count(node) != 0)
            activeNodes.insert(node);
    }
    test_mat(g);

    int node_num;
    node_num = g.nodes.size();
    std::cout << "Number of nodes in matrix : " << node_num << std::endl;

    // Symmetric Matrix (A + A^T)
    int adj_num_sym;
    int adj_max = node_num * (node_num - 1);
    std::vector<int> adj_row_sym(node_num + 1);
    std::vector<int> adj_sym(adj_max);
    adj_set(node_num, adj_max, &adj_num_sym, adj_row_sym.data(), adj_sym.data(), -1, -1);
    create_adj_sym(g, map_out, map_in, adj_row_sym, adj_sym, node_num, adj_max, adj_num_sym, activeNodes);

    // adj_print(node_num, adj_num_sym, adj_row_sym.data(), adj_sym.data(), "  Symmetric adjacency matrix:");
    // adj_show(node_num, adj_num_sym, adj_row_sym.data(), adj_sym.data());

    // Call to rcm
    std::vector<int> perm_sym(node_num);
    genrcm(node_num, adj_num_sym, adj_row_sym.data(), adj_sym.data(), perm_sym.data());
    // Create a new graph with perm
    // std::cout << "Permutation : " << perm_sym << std::endl;
    Graph g_perm;
    g_perm.nodes = g.nodes;
    for (auto arc : g.arcs) {
        g_perm.arcs.insert({perm_sym[arc.first] - 1, perm_sym[arc.second] - 1});
    }

    std::multimap<Node, Node> map_in_perm;
    std::multimap<Node, Node> map_out_perm;
    for (auto arc : g_perm.arcs) {
        map_out_perm.insert(arc);
        map_in_perm.insert({arc.second, arc.first});
    }

    test_mat(g_perm);

    // Unsymmetric Matrix (original)
    std::vector<int> adj_row_unsym(node_num + 1);
    std::vector<int> adj_unsym(adj_max);
    int adj_num_unsym;
    adj_set(node_num, adj_max, &adj_num_unsym, adj_row_unsym.data(), adj_unsym.data(), -1, -1);
    create_adj_unsym(g, map_out, adj_row_unsym, adj_unsym, node_num, adj_max, adj_num_unsym, activeNodes);

    // adj_print(node_num, adj_num_unsym, adj_row_unsym.data(), adj_unsym.data(), "  Unsymmetric adjacency
    // matrix:"); adj_show(node_num, adj_num_unsym, adj_row_unsym.data(), adj_unsym.data());

    std::vector<int> perm_unsym(node_num);
    genrcm(node_num, adj_num_unsym, adj_row_unsym.data(), adj_unsym.data(), perm_unsym.data());

    // Create a new graph with perm
    Graph g_perm_unsym;
    g_perm_unsym.nodes = g.nodes;
    for (auto arc : g.arcs) {
        g_perm_unsym.arcs.insert({perm_unsym[arc.first] - 1, perm_unsym[arc.second] - 1});
    }

    std::multimap<Node, Node> map_in_perm_unsym;
    std::multimap<Node, Node> map_out_perm_unsym;
    for (auto arc : g_perm_unsym.arcs) {
        map_out_perm_unsym.insert(arc);
        map_in_perm_unsym.insert({arc.second, arc.first});
    }
    test_mat(g_perm_unsym);
}

void operations_comparison() {
    // test_rcm();  // Cuthill-Mckee
    // test_fast_rotation_matmult();
    generate_and_compare_partition();
    // test_nauty();
    // test_nauty_with_partition();
    // test_constant_dangling_nodes();
    // test_fast_rotation();
    // compare_bootstrap();
    // try_to_explode_rot();
    // compare_ring_dim();
    // compare_multiple_ring();
    // compare_multiple_batch();
    // compare_multiple_depth();
}
#pragma once

#include <openfhe.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <ostream>

using namespace lbcrypto;

class Experiment {
public:
    Experiment() { this->save_to_file = false; };
    Experiment(std::string filename);
    ~Experiment() = default;

    void saveToFile(std::string filename);

    void save(std::string filename);
    void setCryptoParams(CCParams<CryptoContextCKKSRNS> params);
    void setRingsize(int ringsize);
    void setCurrentRun(int run) { this->current_run = run; };
    enum Operation {
        ADD_SCALAR, //unit op
        ADD,
        MUL,
        MUL_SCALAR,
        MUL_PLAIN,
        ROT,
        FAST_ROT,
        ENCODE,
        ENCRYPT,
        DECRYPT,
        GEN,
        INIT_DC,
        SYNC,
        PREPROC_DC,
        EVAL_DC,
        EVAL_OP,
        PRUNING,
        BOOTSTRAP,
        ADD_SCALAR_SYNC, //synch
        ADD_SYNC,
        ENCODE_SYNC,
        ENCRYPT_SYNC,
        MUL_SYNC,
        MUL_SCALAR_SYNC,
        MUL_PLAIN_SYNC,
        ROT_SYNC,
        CIPHER_SYNC,
        ADD_SCALAR_COMPUTE, //Compute
        ADD_COMPUTE,
        MUL_COMPUTE,
        MUL_SCALAR_COMPUTE,
        MUL_PLAIN_COMPUTE,
        ROT_COMPUTE,
        ENCODE_COMPUTE,
        ENCRYPT_COMPUTE,
        ADD_SCALAR_PREP, // Prepar
        ADD_PREP,
        MUL_PREP,
        MUL_SCALAR_PREP,
        MUL_PLAIN_PREP,
        ROT_PREP,
        ENCODE_PREP,
        ENCRYPT_PREP,
        PREP
    };
    enum Unit {
        MICROSECONDS,
        MILLISECONDS,
        SECONDS,
        UNIT
    };
    void addValue(std::vector<long> values, std::string title, Operation operation);
    void addValue(long value, std::string title, Operation operation, Unit unit);
    void addValue(std::vector<long> values, std::string title, std::vector<Operation> operations, Unit unit);
    friend std::ostream& operator<<(std::ostream& os, const Experiment& dt);

    void disable_crypto_context() { this->crypto_context_available = false; };

private:
    bool save_to_file = false;
    bool crypto_context_available = false;
    std::ofstream* file;
    struct CryptoParams {
        int multDepth;
        int batchSize;
    } current_params;

    struct Statement {
        CryptoParams crypto_params;
        long value;
        std::string title;
        Operation operation;
        Unit unit;
        int ring_size;
        int run;
    };
    int current_run = 0;
    int current_ring_size;

    std::vector<Statement> times;
};

extern Experiment expe_poppy_m;
extern Experiment expe_poppy_uj;
extern Experiment expe_naive;
extern int number_encoding;
extern int number_encryption;
extern int number_mult;
extern int number_add_scalar;
extern int number_add;
extern int number_mult_scalar;
extern int number_mult_plain;
extern int number_rot;
extern std::string folder_expe;

const std::string currentDateTimeExpe();
#include <sys/stat.h>
#include <iostream>
#include "openfhe.h"
// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
using namespace lbcrypto;

#ifdef BOOTSTRAPPING
#define SERIAL_FOLDER "serial-bootstrap/"
#else
#define SERIAL_FOLDER "serial/"
#endif
#define SERIAL_CC "cryptocontext"
#define SERIAL_PUB "pub_key"
#define SERIAL_SECRET "secr_key"
#define SERIAL_MULT "key_mult"
#define SERIAL_ROT "rot_keys"

CryptoContext<DCRTPoly> gen_crypto_context_test(uint32_t multDepth);
std::pair<CryptoContext<DCRTPoly>, KeyPair<DCRTPoly>> gen_crypto_context(uint32_t multDepth, int batchSize);
void gen_rotation_keys(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, uint32_t nbSlots);
void gen_rotation_keys_coeus(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, uint32_t nbSlots);
void serialize_cc(CryptoContext<DCRTPoly> cc, KeyPair<DCRTPoly> keys, int multDepth);
bool deserialize_cc(CryptoContext<DCRTPoly>* cc, KeyPair<DCRTPoly>* keys, std::string folder_name);
void serialize_cipher(Ciphertext<DCRTPoly> ct, std::string folder_name, std::string filename);
std::string get_folder_by_params(int multDepth, int batchSize);

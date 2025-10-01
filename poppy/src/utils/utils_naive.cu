#include "include/utils_naive.cuh"
#include "../include/naive.cuh"
#include "../include/poppy.cuh"
//include files from fideslib
// #include "gtest/gtest.h"
#include "/home/comp/csxtchen/cbct/poppy/poppy-icde/openfhe/openfhe-development/src/pke/include/openfhe.h"

#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Ciphertext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Context.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/KeySwitchingKey.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Limb.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Plaintext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/openfhe-interface/RawCiphertext.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib//include/ConstantsGPU.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/Math.cuh"
// #include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/ParametrizedTest.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/cpuNTT.hpp"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/cpuNTT_nega.hpp"

//#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/test/hook.h"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/ApproxModEval.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/Bootstrap.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/BootstrapPrecomputation.cuh"
#include "/home/comp/csxtchen/cbct/FIDE/FIDESlib/FIDESlib/include/CKKS/CoeffsToSlots.cuh"

//include normal files
#include <fstream>
#include <iostream>
#include <vector>
using namespace FIDESlib;
using namespace lbcrypto;
std::vector<FIDESlib::PrimeRecord> p64{
    {.p = 2305843009218281473}, {.p = 2251799661248513}, {.p = 2251799661641729}, {.p = 2251799665180673},
    {.p = 2251799682088961},    {.p = 2251799678943233}, {.p = 2251799717609473}, {.p = 2251799710138369},
    {.p = 2251799708827649},    {.p = 2251799707385857}, {.p = 2251799713677313}, {.p = 2251799712366593},
    {.p = 2251799716691969},    {.p = 2251799714856961}, {.p = 2251799726522369}, {.p = 2251799726129153},
    {.p = 2251799747493889},    {.p = 2251799741857793}, {.p = 2251799740416001}, {.p = 2251799746707457},
    {.p = 2251799756013569},    {.p = 2251799775805441}, {.p = 2251799763091457}, {.p = 2251799767154689},
    {.p = 2251799765975041},    {.p = 2251799770562561}, {.p = 2251799769776129}, {.p = 2251799772266497},
    {.p = 2251799775281153},    {.p = 2251799774887937}, {.p = 2251799797432321}, {.p = 2251799787995137},
    {.p = 2251799787601921},    {.p = 2251799791403009}, {.p = 2251799789568001}, {.p = 2251799795466241},
    {.p = 2251799807131649},    {.p = 2251799806345217}, {.p = 2251799805165569}, {.p = 2251799813554177},
    {.p = 2251799809884161},    {.p = 2251799810670593}, {.p = 2251799818928129}, {.p = 2251799816568833},
    {.p = 2251799815520257}};

std::vector<FIDESlib::PrimeRecord> sp64{
    {.p = 2305843009218936833}, {.p = 2305843009220116481}, {.p = 2305843009221820417}, {.p = 2305843009224179713},
    {.p = 2305843009225228289}, {.p = 2305843009227980801}, {.p = 2305843009229160449}, {.p = 2305843009229946881},
    {.p = 2305843009231650817}, {.p = 2305843009235189761}, {.p = 2305843009240301569}, {.p = 2305843009242923009},
    {.p = 2305843009244889089}, {.p = 2305843009245413377}, {.p = 2305843009247641601}};

FIDESlib::CKKS::Parameters params{.logN = 13, .L = 5, .dnum = 2, .primes = p64, .Sprimes = sp64};
std::map<int, Ciphertext<DCRTPoly>> initialisation_naif_PR(Naive& data_center) {
    int n = data_center.get_number_of_nodes();
    auto cc = data_center.get_cc();
    auto public_key = data_center.get_public_key();
    std::vector<double> x = {1.0 / n};

    // We create the initialisation
    Plaintext iv = cc->MakeCKKSPackedPlaintext(x);
    auto ct_iv = cc->Encrypt(public_key, iv);

    // We fill the map with the initialisation ciphertext
    std::map<int, Ciphertext<DCRTPoly>> map;
    for (int i = 0; i < n; i++) {
        map.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, ct_iv));
    }
    return map;
}

std::map<int, Ciphertext<DCRTPoly>> initialisation_multi_DC_naif_PR(
    int n,
    int DC,
    std::vector<int> partition,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    std::vector<double> x = {1.0 / n};

    // We create the initialisation
    Plaintext iv = cc->MakeCKKSPackedPlaintext(x);
    auto ct_iv = cc->Encrypt(keys.publicKey, iv);

    // We fill the map with the initialisation ciphertext
    std::map<int, Ciphertext<DCRTPoly>> map;
    for (int i = 0; i < n; i++) {
        if (partition[i] == DC)
            map.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, ct_iv));
    }
    return map;
}

void inter_DC_com(
    std::vector<std::map<int, Ciphertext<DCRTPoly>>>& inter_DC_PRs,
    std::map<int, Ciphertext<DCRTPoly>> PRs,
    int n,
    int DC,
    std::multimap<int, int> map_out,
    std::vector<int> partition,
    CryptoContext<DCRTPoly> cc) {
        
    FIDESlib::CKKS::RawParams raw_param = FIDESlib::CKKS::GetRawParams(cc);
    FIDESlib::CKKS::Context GPUcc(params.adaptTo(raw_param), {0});
    std::size_t size_partition = partition.size();
    for (std::size_t i = 0; i < size_partition; i++) {
        if (partition[i] == DC) {
            if (map_out.count(i) > 0) {
                std::pair<std::multimap<int, int>::iterator,
                          std::multimap<int, int>::iterator>
                    values;
                values = map_out.equal_range(i);  // To get all the values of the key i

                for (auto it = values.first;
                     it != values.second; ++it) {  // For each values :
                    int j = it->second;            // node j -> i

                    if (partition[i] != partition[j]) {  // We add PR(i)/deg_out(i) to the
                                                         // knowledge of the DC of j
                        if (map_out.count(i) > 1) {
                            auto search = PRs.find(i);
                            //GPU CKKS
                            FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                            FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                            GPUct1.multScalar(1.0 / map_out.count(i));
                            FIDESlib::CKKS::RawCipherText raw_res;
                            GPUct1.store(GPUcc, raw_res);
                            Ciphertext<DCRTPoly> val;
                            GetOpenFHECipherText(val, raw_res);
                            /*Ciphertext<DCRTPoly>*/ val =
                                cc->EvalMult(search->second, 1.0 / map_out.count(i));
                            inter_DC_PRs[partition[j]].insert(
                                std::pair<int, Ciphertext<DCRTPoly>>(i, val));
                        }

                        else {
                            auto search = PRs.find(i);
                            inter_DC_PRs[partition[j]].insert(
                                std::pair<int, Ciphertext<DCRTPoly>>(i, search->second));
                        }
                    }
                }
            }
        }
    }
}

// PRs = naif_PR_calcul(PRs, data_center.get_map_in(),
// data_center.get_map_out(), data_center.get_number_of_nodes(),
// data_center.get_keys(), data_center.get_cc());

std::map<int, Ciphertext<DCRTPoly>> naif_PR_calcul(
    Naive data_center,
    std::map<int, Ciphertext<DCRTPoly>> PRs) {
    int n = data_center.get_number_of_nodes();
    auto cc = data_center.get_cc();
    auto keys = data_center.get_keys();
    FIDESlib::CKKS::RawParams raw_param = FIDESlib::CKKS::GetRawParams(cc);
    
    FIDESlib::CKKS::Context GPUcc(params.adaptTo(raw_param), {0});
    std::vector<double> x_0 = {0.0};
    std::vector<double> x_15 = {0.15};

    // We create ciphers of 0 and 0.15

    Plaintext zero = cc->MakeCKKSPackedPlaintext(x_0);
    Plaintext zero_15 = cc->MakeCKKSPackedPlaintext(x_15);
    auto ct_zero = cc->Encrypt(keys.publicKey, zero);
    auto ct_zero_15 = cc->Encrypt(keys.publicKey, zero_15);

    std::map<int, Ciphertext<DCRTPoly>> PR_end_iteration;

    for (int i = 0; i < n; i++) {
        if (data_center.get_map_in().count(i) > 0) {
            auto values = data_center.get_map_in().equal_range(
                i);  // To get all the values of the key i

            data_center.print_data_center();

            auto sum = ct_zero;

            for (auto it = values.first; it != values.second;
                 ++it) {  // For each values :
                std::cout << " (debug1)" << std::endl;
                int j = it->second;  // node j -> i

                // If j as deg_out > 1
                if (data_center.get_map_out().count(j) > 1) {
                    auto search = PRs.find(j);
                     //GPU CKKS
                    FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                    FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                    GPUct1.multScalar(1.0 / data_center.get_map_out().count(j));
                    Ciphertext<DCRTPoly> val = cc->EvalMult(
                        search->second,
                        1.0 / data_center.get_map_out().count(j));  // PR(j)/deg_out(j)
                    FIDESlib::CKKS::RawCipherText raw2 = FIDESlib::CKKS::GetRawCipherText(cc, sum);
                    FIDESlib::CKKS::Ciphertext GPUct2(GPUcc, raw2);
                    GPUct2.add(GPUct1);                 
                    FIDESlib::CKKS::RawCipherText raw_res;
                    GPUct2.store(GPUcc, raw_res);
                    GetOpenFHECipherText(sum, raw_res);
                    sum = cc->EvalAdd(sum, val);
                } else {
                    auto search = PRs.find(j);
                    FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                    FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                    FIDESlib::CKKS::RawCipherText raw2 = FIDESlib::CKKS::GetRawCipherText(cc, sum);
                    FIDESlib::CKKS::Ciphertext GPUct2(GPUcc, raw2);
                    GPUct2.add(GPUct1);
                    FIDESlib::CKKS::RawCipherText raw_res;
                    GPUct2.store(GPUcc, raw_res);
                    GetOpenFHECipherText(sum, raw_res);
                    sum = cc->EvalAdd(sum, search->second);
                }
                std::cout << " (debug2)" << std::endl;
            }
            FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, sum );
            FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
            GPUct1.multScalar(0.85);
            GPUct1.addScalar(0.15);
            FIDESlib::CKKS::RawCipherText raw_res;
            GPUct1.store(GPUcc, raw_res);
            GetOpenFHECipherText(sum, raw_res);
            std::cout << " (debug3)" << std::endl;
            sum = cc->EvalMult(sum, 0.85);
            std::cout << " (debug4)" << std::endl;
            sum = cc->EvalAdd(sum, 0.15);
            std::cout << " (debug5)" << std::endl;

            PR_end_iteration.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, sum));
            std::cout << " (debug6)" << std::endl;
        } else {  // If i is without in_edges PR_i = 0.15
            PR_end_iteration.insert(
                std::pair<int, Ciphertext<DCRTPoly>>(i, ct_zero_15));
        }
        std::cout << " (debug7)" << std::endl;
    }
    std::cout << " (debug8)" << std::endl;
    return PR_end_iteration;
}

std::map<int, Ciphertext<DCRTPoly>> multi_DC_naif_PR(
    std::map<int, Ciphertext<DCRTPoly>> PRs,
    std::map<int, Ciphertext<DCRTPoly>> inter_DC_PRs,
    std::multimap<int, int> map_in,
    std::multimap<int, int> map_out,
    int n,
    std::vector<int> partition,
    int DC,
    KeyPair<DCRTPoly> keys,
    CryptoContext<DCRTPoly> cc) {
    std::vector<double> x_0(n);
    std::vector<double> x_15(n);
    for (int i = 0; i < n; i++) {
        x_0[i] = 0.0;
        x_15[i] = 0.15;
    }
    FIDESlib::CKKS::RawParams raw_param = FIDESlib::CKKS::GetRawParams(cc);
    FIDESlib::CKKS::Context GPUcc(params.adaptTo(raw_param), {0});
    // We create ciphers of 0 and 0.15
    Plaintext zero = cc->MakeCKKSPackedPlaintext(x_0);
    Plaintext zero_15 = cc->MakeCKKSPackedPlaintext(x_15);
    auto ct_zero = cc->Encrypt(keys.publicKey, zero);
    auto ct_zero_15 = cc->Encrypt(keys.publicKey, zero_15);

    std::map<int, Ciphertext<DCRTPoly>> PR_end_iteration;

    // À remettre lorsque l'on fera la version multi-DC
    // int index_PR_inter_DC = 0;
    // int index_node_in_DC = 0;

    std::size_t size_partition = partition.size();
    for (std::size_t i = 0; i < size_partition; i++) {
        if (partition[i] == DC) {  // Check for nodes in DC

            if (map_in.count(i) > 0) {
                std::pair<std::multimap<int, int>::iterator,
                          std::multimap<int, int>::iterator>
                    values;
                values = map_in.equal_range(i);  // To get all the values of the key i

                auto sum = ct_zero;

                for (auto it = values.first;
                     it != values.second; ++it) {  // For each values :
                    int j = it->second;            // node j -> i

                    if (partition[j] == DC) {        // If node j is in DC
                        if (map_out.count(j) > 1) {  // If j as deg_out > 1
                            auto search = PRs.find(j);
                            FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                            FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                            GPUct1.multScalar(1.0 / map_out.count(i));
                            FIDESlib::CKKS::RawCipherText raw_res;
                            GPUct1.store(GPUcc, raw_res);
                            Ciphertext<DCRTPoly> val;
                            GetOpenFHECipherText(val, raw_res);
                            FIDESlib::CKKS::RawCipherText raw2 = FIDESlib::CKKS::GetRawCipherText(cc, sum);
                            FIDESlib::CKKS::Ciphertext GPUct2(GPUcc, raw2);
                            GPUct2.add(GPUct1);                 
                           // FIDESlib::CKKS::RawCipherText raw_res;
                            GPUct2.store(GPUcc, raw_res);
                            GetOpenFHECipherText(sum, raw_res);
                            /*Ciphertext<DCRTPoly>*/ val =
                                cc->EvalMult(search->second, 1.0 / map_out.count(j));                           
                            cc->EvalAddInPlace(sum, val);
                        }

                        else {
                            auto search = PRs.find(j);
                            FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                            FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                            FIDESlib::CKKS::RawCipherText raw2 = FIDESlib::CKKS::GetRawCipherText(cc, sum);
                            FIDESlib::CKKS::Ciphertext GPUct2(GPUcc, raw2);
                            GPUct2.add(GPUct1);
                            FIDESlib::CKKS::RawCipherText raw_res;
                            GPUct2.store(GPUcc, raw_res);
                            GetOpenFHECipherText(sum, raw_res);
                            sum = cc->EvalAdd(sum, search->second);
                        }
                    } else {  // If node j is not in DC we add inter-DC communication
                        auto search = inter_DC_PRs.find(j);
                        FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, search->second);
                        FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                        FIDESlib::CKKS::RawCipherText raw2 = FIDESlib::CKKS::GetRawCipherText(cc, sum);
                        FIDESlib::CKKS::Ciphertext GPUct2(GPUcc, raw2);
                        GPUct2.add(GPUct1);
                        FIDESlib::CKKS::RawCipherText raw_res;
                        GPUct2.store(GPUcc, raw_res);
                        GetOpenFHECipherText(sum, raw_res);
                        sum = cc->EvalAdd(sum, search->second);
                        // À remettre lorsque l'on fera la version multi-DC
                        // index_PR_inter_DC++;
                    }
                }
                FIDESlib::CKKS::RawCipherText raw1 = FIDESlib::CKKS::GetRawCipherText(cc, sum );
                FIDESlib::CKKS::Ciphertext GPUct1(GPUcc, raw1);
                GPUct1.multScalar(0.85);
                GPUct1.addScalar(0.15);
                FIDESlib::CKKS::RawCipherText raw_res;
                GPUct1.store(GPUcc, raw_res);
                GetOpenFHECipherText(sum, raw_res);
                sum = cc->EvalMult(sum, 0.85);
                sum = cc->EvalAdd(sum, 0.15);

                PR_end_iteration.insert(std::pair<int, Ciphertext<DCRTPoly>>(i, sum));
            } else  // If i is without in_edges PR_i = 0.15
                PR_end_iteration.insert(
                    std::pair<int, Ciphertext<DCRTPoly>>(i, ct_zero_15));

            // À remettre lorsque l'on fera la version multi-DC
            // index_node_in_DC++;
        }
    }
    return PR_end_iteration;
}
/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 19:59:21
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-11 20:30:51
 * @FilePath: \openfhe-development\src\pke\include\homoencrypt-compute.h
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#ifndef SRC_LIB_HOMOENCRYPT_COMPUTE_H_
#define SRC_LIB_HOMOENCRYPT_COMPUTE_H_

#include "openfhe.h"
#include "ciphertext-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"

using namespace lbcrypto;
using namespace std;
using namespace std::chrono;

// 自定义类型简化书写
using Plain = Plaintext;
using Cipher = Ciphertext<DCRTPoly>;

class HomoEncryptCompute {
    CryptoContext<DCRTPoly> context;

public:
    // 好像不需要构造函数
    HomoEncryptCompute(){};

    /**
     * Generate the cryptocontext for the evaluation of the bitonic sorting network
     *
     * @param num_slots The number of slots in the ciphertext.
     * @param levels_required The required circuit depth
     * @param toy_parameters Choose whether to use toy parameters (true) or 128-bit security parameters (false)
     * @return the total depth of the circuit, including the bootstrapping operation
     */
    int generate_context_network(int num_slots, int levels_required, bool toy_parameters);

    /**
     * Generate the rotation keys required by the network-based sorting
     *
     * @param num_slots The number of slots in the ciphertext.
     */
    void generate_rotation_keys_network(int num_slots);

    /**
     * Generate the cryptocontext for the evaluation of the permutation sorting network
     *
     * @param num_slots The number of slots in the ciphertext.
     * @param levels_required The required circuit depth
     * @param toy_parameters Choose whether to use toy parameters (true) or 128-bit security parameters (false)
     */
    void generate_context_permutation(int num_slots, int levels_required, bool toy_parameters);

    /**
     * Generate a rotation key
     *
     * @param index The index of the rotation
     */
    void generate_rotation_key(int index);

    /**
      * Basic FHE operations
      */

    // Encode a vector of doubles into a plaintext
    Plain encode(const vector<double>& vec, int level, int num_slots);

    // Encodes a value in a plaintext (it will be repeated)
    Plain encode(double value, int level, int num_slots);

    // Encrypt a vector of doubles
    Cipher encrypt(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0);

    // Encrypt a vector in expanded encoding
    Cipher encrypt_expanded(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0, int repetitions = 1);

    // Encrypt a vector in repeated encoding
    Cipher encrypt_repeated(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0, int repetitions = 1);

    //Decodes a ciphertexxt
    vector<double> decode(const Plain& p);

    // Decrypt a ciphertext
    Plain decrypt(const Cipher& c);


    // Add two ciphertexts/plaintexts
    Cipher add(const Cipher& c1, const Cipher& c2);
    Cipher add(const Cipher& c, const Plain& p);
    Cipher add(const Cipher& c, double d);

    // Add multiple ciphertexts using a tree structure
    Cipher add_tree(vector<Cipher> v);

    // Subtract two ciphertexts/plaintexts
    Cipher sub(const Cipher& c1, const Cipher& c2);
    Cipher sub(const Cipher& c, const Plain& p);

    // Multiply two ciphertexts/plaintext
    Cipher mult(const Cipher& c, const Plain& p);
    Cipher mult(const Cipher& c, double d);
    Cipher mult(const Cipher& c1, const Cipher& c2);

    // Rotate a ciphertext by a specified index
    Cipher rot(const Cipher& c, int index);

    // Perform bootstrapping operation on a ciphertext
    Cipher bootstrap(const Cipher& c);

    /**
      * Permutation-based operations
      */

    // Approximation of the k-scaled Sigmoid
    Cipher sigmoid(const Cipher& in, int n, int degree, int scaling);

    // Rotate-and-sum elements at distance n
    Cipher rotsum(const Cipher& in, int n);

    // Approximation of sinc function
    Cipher sinc(const Cipher& in, int degree, int n);


    /**
      * Network-based operations
      */

    // Approximation of the ReLU function
    Cipher relu(const Cipher& in, int degree, int n);


    /**
      * Utilities
      */

    // Print the values of slots in a ciphertext
    void print(const Cipher& c, int slots = 0, string prefix = "");

private:
    KeyPair<DCRTPoly> key_pair; // Key pair for the FHE system

    void print_moduli_chain(const DCRTPoly& poly);

};

#endif  // SRC_LIB_HOMOENCRYPT_COMPUTE_H_
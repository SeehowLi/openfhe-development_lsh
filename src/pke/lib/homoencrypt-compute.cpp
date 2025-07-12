/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 20:23:55
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 15:26:05
 * @FilePath: \openfhe-development\src\pke\lib\homoencrypt-compute.cpp
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

#include "homoencrypt-compute.h"

// 用于设置context等加密参数--放在测试程序中，不用单独封装函数
int HomoEncryptCompute::generate_context_network(int num_slots, int levels_required, bool toy_parameters) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    vector<uint32_t> level_budget;

    level_budget = {2, 2};
    int dcrtBits = 54;
    int firstMod = 55;
    // 是否满足标准
    if (toy_parameters) {
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(1 << 12);
    } else {
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(1 << 16);
    }

    cout << "Levels required: " << levels_required << endl;
    // 密钥存储的大小
    parameters.SetNumLargeDigits(5);
    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    // 额外的一层？
    int levelsUsedBeforeBootstrap = levels_required + 1;

    int circuit_depth = levelsUsedBeforeBootstrap + FHECKKSRNS::GetBootstrapDepth(level_budget, SPARSE_TERNARY);

    parameters.SetMultiplicativeDepth(circuit_depth);

    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    key_pair = context->KeyGen();
    // 有意思的函数
    print_moduli_chain(key_pair.publicKey->GetPublicElements()[0]);

    cout << endl;

    context->EvalMultKeyGen(key_pair.secretKey);

    context->EvalBootstrapSetup(level_budget, {0, 0}, num_slots);
    context->EvalBootstrapKeyGen(key_pair.secretKey, num_slots);


    return circuit_depth;
}

// 用于设置context等加密参数--放在测试程序中，不用单独封装函数
void HomoEncryptCompute::generate_context_permutation(int num_slots, int levels_required, bool toy) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);

    int dcrtBits = 45;
    int firstMod = 48;

    if (toy) {
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);

        if (num_slots <= 1 << 14) parameters.SetRingDim(1 << 15);
        if (num_slots <= 1 << 13) parameters.SetRingDim(1 << 14);
        if (num_slots <= 1 << 12) parameters.SetRingDim(1 << 13);
        if (num_slots <= 1 << 11) parameters.SetRingDim(1 << 12);

        cout << "n: " << num_slots << endl;
    } else {
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(1 << 16);
    }

    cout << "N: " << parameters.GetRingDim() << ", ";

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(2);

    parameters.SetMultiplicativeDepth(levels_required);


    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    //context->Enable(FHE);

    key_pair = context->KeyGen();

    print_moduli_chain(key_pair.publicKey->GetPublicElements()[0]);

    cout << ", λ: 128 bits" << endl;

    context->EvalMultKeyGen(key_pair.secretKey);

}

// 可以当作测试函数的子函数
void HomoEncryptCompute::generate_rotation_keys_network(int num_slots) {
    vector<int> rotations;

    for (int i = 0; i < log2(num_slots); i++) {
        rotations.push_back(pow(2, i));
        rotations.push_back(-pow(2, i));
    }

    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

// 同样可以当作测试函数的子函数
void HomoEncryptCompute::generate_rotation_key(int index) {
    vector<int> rotations;

    rotations.push_back(index);

    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

// 放在主函数
Plain HomoEncryptCompute::encode(const vector<double> &vec, int level, int num_slots) {
    Plain p = context->MakeCKKSPackedPlaintext(vec, 1, level, nullptr, num_slots);
    p->SetLength(num_slots);

    return p;
}
// 放在主函数，上面那个encode的简易版本
Plain HomoEncryptCompute::encode(double value, int level, int num_slots) {
    if (!context) {
        throw std::runtime_error("Context not initialized in encode");
    }
    vector<double> repeated_value;
    for (int i = 0; i < num_slots; i++) repeated_value.push_back(value);

    return encode(repeated_value, level, num_slots);
}

// 加密，主函数
Cipher HomoEncryptCompute::encrypt(const vector<double> &vec, int level, int num_slots) {
    Plain p = encode(vec, level, num_slots);

    return context->Encrypt(p, key_pair.publicKey);
}

// 扩展加密，可以放在sorting函数中
Cipher HomoEncryptCompute::encrypt_expanded(const vector<double> &vec, int level, int num_slots, int repetitions) {
    vector<double> repeated;

    for (std::size_t i = 0; i < vec.size(); i++) {
        for (int j = 0; j < repetitions; j++) {
            repeated.push_back(vec[i]);
        }
    }

    Plain p = encode(repeated, level, num_slots);

    return context->Encrypt(p, key_pair.publicKey);
}

// 重复加密，可以放在sorting函数中
Cipher HomoEncryptCompute::encrypt_repeated(const vector<double> &vec, int level, int num_slots, int repetitions) {
    vector<double> repeated;

    for (int i = 0; i < repetitions; i++) {
        for (std::size_t j = 0; j < vec.size(); j++) {
            repeated.push_back(vec[j]);
        }
    }

    Plain p = encode(repeated, level, num_slots);

    return context->Encrypt(p, key_pair.publicKey);
}

// 解码，主函数
vector<double> HomoEncryptCompute::decode(const Plain& p) {
    return p->GetRealPackedValue();
}

// 解密，主函数
Plain HomoEncryptCompute::decrypt(const Cipher &c) {
    Plain p;
    context->Decrypt(key_pair.secretKey, c, &p);

    return p;
}

// 密文加法加法，主函数
Cipher HomoEncryptCompute::add(const Cipher &a, const Cipher &b) {
    return context->EvalAdd(a, b);
}

// 明密文加法，主函数
Cipher HomoEncryptCompute::add(const Cipher &a, const Plain &b) {
    return context->EvalAdd(a, b);
}

Cipher HomoEncryptCompute::add(const Cipher &a, double d) {
    return context->EvalAdd(a, encode(d, a->GetLevel(), a->GetSlots()));
}

Cipher HomoEncryptCompute::add_tree(vector<Cipher> v) {
    return context->EvalAddMany(v);
}

Cipher HomoEncryptCompute::sub(const Cipher &a, const Cipher &b) {
    return context->EvalSub(a, b);
}

Cipher HomoEncryptCompute::sub(const Cipher &c, const Plain &p) {
    return context->EvalSub(c, p);
}

Cipher HomoEncryptCompute::mult(const Cipher &c, const Plain& p) {
    return context->EvalMult(c, p);
}

Cipher HomoEncryptCompute::mult(const Cipher &c1, const Cipher &c2) {
    return context->EvalMult(c1, c2);
}

Cipher HomoEncryptCompute::mult(const Cipher &c, double v) {
    return context->EvalMult(c, encode(v, c->GetLevel(), c->GetSlots()));
}

Cipher HomoEncryptCompute::rot(const Cipher& c, int index) {
    return context->EvalRotate(c, index);
}

Cipher HomoEncryptCompute::bootstrap(const Cipher &c) {
    return context->EvalBootstrap(c);
}

// 旋转累加
Cipher HomoEncryptCompute::rotsum(const Cipher &in, int n) {
    Cipher result = add(in, rot(in, n));

    for (int i = 1; i < log2(n); i++) {
        result = add(result, rot(result, n * pow(2, i)));
    }

    return result;
}

// 重点的函数近似
Cipher HomoEncryptCompute::sigmoid(const Cipher &in, int n, int degree, int scaling) {
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }
    return context->EvalChebyshevFunction([scaling, n](double x) -> double {
        return 1/(n + n * pow(2.71828182846, -scaling*x));

    }, in, -1, 1, degree);
}


Cipher HomoEncryptCompute::sinc(const Cipher &in, int poly_degree, int n) {
    return context->EvalChebyshevFunction([n](double x) -> double { return sin(3.14159265358979323846 * x * n) / (3.14159265358979323846 * x * n); },
                                          in,
                                          -1,
                                          1, poly_degree);
}

Cipher HomoEncryptCompute::relu(const Cipher &in, int poly_degree, int n) {
    return context->EvalChebyshevFunction([](double x) -> double { if (x > 0) return x; return 0; },
                                          in,
                                          -1,
                                          1, poly_degree);
}

// 打印函数
void HomoEncryptCompute::print(const Cipher &c, int slots, string prefix) {
    if (slots == 0) {
        slots = c->GetSlots();
    }

    cout << prefix;

    Plain result;
    context->Decrypt(key_pair.secretKey, c, &result);
    result->SetSlots(slots);
    vector<double> v = result->GetRealPackedValue();

    cout << "[ ";

    for (int i = 0; i < slots; i += 1) {
        string segno = "";
        if (v[i] > 0) {
            segno = "";
        } else {
            segno = "-";
            v[i] = -v[i];
        }


        if (i == slots - 1) {
            cout << segno << v[i] << " ]";
        } else {
            if (abs(v[i]) <= 0.00001)
                cout << "0.0000" << " ";
            else
                cout << segno << v[i] << " ";
        }
    }

    cout << endl;
}

void HomoEncryptCompute::print_moduli_chain(const DCRTPoly& poly){
    int num_primes = poly.GetNumOfElements();
    double total_bit_len = 0.0;
    for (int i = 0; i < num_primes; i++) {
        auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
        total_bit_len += log(qi.ConvertToDouble()) / log(2);
    }
    std::cout << "log(QP): " << ((int)total_bit_len);
}
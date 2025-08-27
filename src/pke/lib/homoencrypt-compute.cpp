/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 20:23:55
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-08-02 10:58:06
 * @FilePath: \openfhe-131\src\pke\lib\homoencrypt-compute.cpp
 * @Description: Unified homomorphic encryption computation library implementation
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

#include "homoencrypt-compute.h"

// ========== Context Generation Functions ==========

void HomoEncryptCompute::generate_context_knn(int num_slots, int levels_required, uint32_t ring_dim, bool toy){
    CCParams<CryptoContextCKKSRNS> parameters;
    
    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);
    // parameters.SetSecretKeyDist(lbcrypto::SPARSE_TERNARY);

    int dcrtBits = 35;
    int firstMod = 40;

    if (toy) {
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(ring_dim);
        cout << "num_slots: " << num_slots << endl;
    } else {
        // 安全强度设置,1<<16环规模的最大的模数是1747,1<<15是881
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(ring_dim);
    }

    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(1);
    parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    parameters.SetMultiplicativeDepth(levels_required);
    
    // 生成加密上下文
    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    // context->Enable(FHE);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);

    // 预计算所有需要的值
    vector<int> rot_index_all;
    for (int i = 0; i < log2(128); i++) {
        rot_index_all.push_back(pow(2, i) * 128);
        rot_index_all.push_back(pow(2, i));//左移
        rot_index_all.push_back(-1*pow(2, i));//右移动

    }
    // for (int i = 1; i < log2(128)+1; i++) {
    //     rot_index_all.push_back((128 * 127) / pow(2,i));
    // }
    // rot_index_all.push_back(-128);
    // 并行生成所有密钥
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rot_index_all.size()); i++) {
        context->EvalRotateKeyGen(key_pair.secretKey, {rot_index_all[i]});
    }

    // const std::vector<DCRTPoly>& ckks_pk = key_pair.publicKey->GetPublicElements();
    // std::cout << "Moduli chain of pk: " << std::endl;
    // print_moduli_chain_detail(ckks_pk[0]);
}

void HomoEncryptCompute::generate_context_toy(int num_slots, int levels_required, uint32_t ring_dim){
    CCParams<CryptoContextCKKSRNS> parameters;
    
    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);
    // parameters.SetSecretKeyDist(lbcrypto::SPARSE_TERNARY);

    int dcrtBits = 59;
    int firstMod = 60;

    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(ring_dim);
    cout << "num_slots: " << num_slots << endl;
    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    parameters.SetMultiplicativeDepth(levels_required);
    
    // 生成加密上下文
    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, num_slots);

}

void HomoEncryptCompute::generate_context_ctsfirstbts_toy(int num_slots, int levels_extra, uint32_t ring_dim){
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // parameters.SetSecretKeyDist(secretKeyDist);

    int dcrtBits = 52;
    int firstMod = 60;

    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(ring_dim);
    cout << "num_slots: " << num_slots << endl;
    cout << "N: " << parameters.GetRingDim() << endl;

    // num_slots = ring_dim / 2;

    parameters.SetBatchSize(num_slots);

    // ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    // parameters.SetNumLargeDigits(3);
    // parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    std::vector<uint32_t> levelBudget = {3, 3};

    usint depth = levels_extra + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);
    std::cout << "Total depth is: " << depth << std::endl;
    
    // 生成加密上下文
    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);
    
    // context->EvalBootstrapSetup(levelBudget);
    EvalBTSSetup(levelBudget, {0,0}, num_slots , 0, false);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, num_slots);

}

void HomoEncryptCompute::generate_context_stcfirstbts_toy(int num_slots, int levels_extra, uint32_t ring_dim){
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // parameters.SetSecretKeyDist(secretKeyDist);

    int dcrtBits = 52;
    int firstMod = 60;

    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(ring_dim);
    cout << "num_slots: " << num_slots << endl;
    cout << "N: " << parameters.GetRingDim() << endl;

    // num_slots = ring_dim / 2;

    parameters.SetBatchSize(num_slots);

    // ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    // parameters.SetNumLargeDigits(3);
    // parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    std::vector<uint32_t> levelBudget = {3, 3};

    usint depth = levels_extra + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);
    std::cout << "total depth is: " << depth << std::endl;

    
    // 生成加密上下文
    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);
    
    // StC first BTS precompute
    EvalBTSSetup(levelBudget, {0,0}, num_slots, 0, true);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, num_slots);

}


void HomoEncryptCompute::generate_context_128bit(int num_slots, int levels_required, uint32_t ring_dim,
                                                 int dcrtBits, int firstMod) {
    CCParams<CryptoContextCKKSRNS> parameters;
    
    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);
    // parameters.SetSecretKeyDist(lbcrypto::SPARSE_TERNARY);

    // 安全强度设置
    parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
    parameters.SetRingDim(ring_dim);

    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    parameters.SetMultiplicativeDepth(levels_required);
    
    // 生成加密上下文
    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);

    const std::vector<DCRTPoly>& ckks_pk = key_pair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain_detail(ckks_pk[0]);
}

void HomoEncryptCompute::generate_context_bts(std::vector<uint32_t> levelBudget, int num_slots, int levels_required, uint32_t ring_dim,
                                                int lStC, int lCtS,
                                                 int dcrtBits, int firstMod) {
    CCParams<CryptoContextCKKSRNS> parameters;
    
    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);
    // parameters.SetSecretKeyDist(lbcrypto::SPARSE_TERNARY);

    // 安全强度设置
    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(ring_dim);

    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    // ScalingTechnique rescaleTech = FIXEDAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    // parameters.SetDigitSize(30);
    // parameters.SetKeySwitchTechnique(BV);

    parameters.SetMultiplicativeDepth(levels_required);
    
    // 生成加密上下文
    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    EvalStC_CtSSetup(levelBudget, num_slots, lStC, lCtS);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, num_slots);

    const std::vector<DCRTPoly>& ckks_pk = key_pair.publicKey->GetPublicElements();
    std::cout << "Moduli chain of pk: " << std::endl;
    print_moduli_chain_detail(ckks_pk[0]);
}


// ========== Key Generation Functions ==========

void HomoEncryptCompute::generate_rotation_keys_network(int num_slots) {
    vector<int> rotations;

    for (int i = 0; i < log2(num_slots); i++) {
        rotations.push_back(pow(2, i));
        rotations.push_back(-pow(2, i));
    }

    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

void HomoEncryptCompute::generate_rotation_key(int index) {
    vector<int> rotations;
    rotations.push_back(index);
    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

void HomoEncryptCompute::generate_rotation_key(vector<int> rotations) {
    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

// ========== Basic FHE Operations ==========
void HomoEncryptCompute::genbootkey(PrivateKey<DCRTPoly> secretkey, int num_slots) {
    if (!context) {
        throw std::runtime_error("Context not initialized in genbootkey");
    }
    context->EvalBootstrapKeyGen(secretkey, num_slots);
}

PrivateKey<DCRTPoly> HomoEncryptCompute::getsecretkey() {
    if (!context) {
        throw std::runtime_error("Context not initialized in getsecretkey");
    }
    return key_pair.secretKey;
}

Plain HomoEncryptCompute::encode(const vector<double> &vec, int level, int num_slots) {
    Plain p = context->MakeCKKSPackedPlaintext(vec, 1, level, nullptr, num_slots);
    p->SetLength(num_slots);
    return p;
}

Plain HomoEncryptCompute::encode(const vector<double> &vec, size_t scale_deg, int level, int num_slots) {
    Plain p = context->MakeCKKSPackedPlaintext(vec, scale_deg, level, nullptr, num_slots);
    p->SetLength(num_slots);
    return p;
}

Plain HomoEncryptCompute::encode(double value, int level, int num_slots) {
    if (!context) {
        throw std::runtime_error("Context not initialized in encode");
    }
    vector<double> repeated_value;
    for (int i = 0; i < num_slots; i++) repeated_value.push_back(value);
    return encode(repeated_value, 1, level, num_slots);
}

Plain HomoEncryptCompute::encode(const vector<double> &vec, int level) {
    Plain p = context->MakeCKKSPackedPlaintext(vec, 1, level);
    return p;
}

Cipher HomoEncryptCompute::encrypt(const Plain &p) {
    return context->Encrypt(key_pair.publicKey, p);
}

Cipher HomoEncryptCompute::encrypt(const vector<double> &vec, int level, int num_slots) {
    Plain p = encode(vec, 1, level, num_slots);
    return context->Encrypt(p, key_pair.publicKey);
}

Cipher HomoEncryptCompute::encrypt_expanded(const vector<double> &vec, int level, int num_slots, int repetitions) {
    vector<double> repeated;
    for (std::size_t i = 0; i < vec.size(); i++) {
        for (int j = 0; j < repetitions; j++) {
            repeated.push_back(vec[i]);
        }
    }
    Plain p = encode(repeated, 1, level, num_slots);
    return context->Encrypt(p, key_pair.publicKey);
}

Cipher HomoEncryptCompute::encrypt_repeated(const vector<double> &vec, int level, int num_slots, int repetitions) {
    vector<double> repeated;
    for (int i = 0; i < repetitions; i++) {
        for (std::size_t j = 0; j < vec.size(); j++) {
            repeated.push_back(vec[j]);
        }
    }
    Plain p = encode(repeated, 1, level, num_slots);
    return context->Encrypt(p, key_pair.publicKey);
}

vector<double> HomoEncryptCompute::decode(const Plain& p) {
    return p->GetRealPackedValue();
}

Plain HomoEncryptCompute::decrypt(const Cipher &c) {
    Plain p;
    context->Decrypt(key_pair.secretKey, c, &p);
    return p;
}

// Arithmetic operations
Cipher HomoEncryptCompute::add(const Cipher &ciphertext1, const Cipher &ciphertext2) {
    return context->EvalAdd(ciphertext1, ciphertext2);
}

Cipher HomoEncryptCompute::add(const Cipher& ciphertext, Plain& poly) {
    return context->EvalAdd(ciphertext, poly);
}

Cipher HomoEncryptCompute::add(const Cipher &ciphertext, double v) {
    Plain pt_temp =  encode(v, ciphertext->GetLevel(), ciphertext->GetSlots());
    return context->EvalAdd(ciphertext, pt_temp);
}

void HomoEncryptCompute::add_inplace(Cipher &ciphertext1, const Cipher &ciphertext2) {
    context->EvalAddInPlace(ciphertext1, ciphertext2);
}

Cipher HomoEncryptCompute::add_tree(vector<Cipher> v) {
    return context->EvalAddMany(v);
}

Cipher HomoEncryptCompute::sub(const Cipher &ciphertext1, const Cipher &ciphertext2) {
    return context->EvalSub(ciphertext1, ciphertext2);
}

Cipher HomoEncryptCompute::sub(const Cipher &ciphertext, Plain &poly) {
    return context->EvalSub(ciphertext, poly);
}

Cipher HomoEncryptCompute::mult(const Cipher &ciphertext, const Plain& poly) {
    return context->EvalMult(ciphertext, poly);
}

Cipher HomoEncryptCompute::mult(const Cipher &ciphertext1, const Cipher &ciphertext2) {
    return context->EvalMult(ciphertext1, ciphertext2);
}

Cipher HomoEncryptCompute::mult(const Cipher &ciphertext, double v) {
    return context->EvalMult(ciphertext, encode(v, ciphertext->GetLevel(), ciphertext->GetSlots()));
}

Cipher HomoEncryptCompute::square(const Cipher &ciphertext) {
    return context->EvalSquare(ciphertext);
}

Cipher HomoEncryptCompute::rot(const Cipher& ciphertext, int index) {
    return context->EvalRotate(ciphertext, index);
}

Cipher HomoEncryptCompute::rotsum(const Cipher &in, int n) {
    Cipher result = add(in, rot(in, n));
    for (int i = 1; i < log2(n); i++) {
        result = add(result, rot(result, n * pow(2, i)));
    }
    return result;
}

Cipher HomoEncryptCompute::modreduceinternal(Cipher& in, int degree) {
    auto algo = context->GetScheme();
    auto result = algo->ModReduceInternal(in, degree);   
    return result;
}

void HomoEncryptCompute::modreduceinternalinplace(Cipher& in, int degree) {
    auto algo = context->GetScheme();
    algo->ModReduceInternalInPlace(in, degree);
}

// ========== Polynomial Operations ==========

Poly HomoEncryptCompute::PolyFromDCRTPoly(const DCRTPoly& poly) {
    DCRTPoly polyCopy(poly);
    polyCopy.SetFormat(Format::COEFFICIENT);
    Poly result=polyCopy.CRTInterpolate();
    return result;
}

NativePoly HomoEncryptCompute::ShiftRight(const NativePoly& poly, uint32_t shift) {
    size_t N = poly.GetRingDimension();
    NativePoly shiftedPoly(poly);
    shiftedPoly.SetFormat(Format::COEFFICIENT);

    for (size_t j = 0; j < N; j++) {
        if (j >= shift ) {
            shiftedPoly[j] = poly[j - shift];
        } else {
            shiftedPoly[j] = poly.GetModulus() - poly[j - shift + N];       
        }
    }
    return shiftedPoly;
}

DCRTPoly HomoEncryptCompute::ShiftRight(const DCRTPoly& poly, uint32_t shift) {
    DCRTPoly shiftedPoly(poly);
    shiftedPoly.SetFormat(Format::COEFFICIENT);

    for (size_t i = 0; i < shiftedPoly.GetNumOfElements(); i++) {
        NativePoly element = shiftedPoly.GetElementAtIndex(i);
        element = ShiftRight(element, shift);
        shiftedPoly.SetElementAtIndex(i, element);
    }

    shiftedPoly.SetFormat(poly.GetFormat());
    return shiftedPoly;
}

// ========== Ciphertext Utilities ==========

Cipher HomoEncryptCompute::ShiftRight(const Cipher& ciphertext, uint32_t shift) {
    Cipher shiftedCiphertext = ciphertext->Clone();
    std::vector<DCRTPoly> elements = shiftedCiphertext->GetElements();
    for (size_t i = 0; i < elements.size(); i++) {
        elements[i] = ShiftRight(elements[i], shift);
    }
    shiftedCiphertext->SetElements(elements);
    return shiftedCiphertext;
}

Poly HomoEncryptCompute::myDecrypt(const Cipher& ciphertext, const PrivateKey<DCRTPoly> privateKey) {
    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s(privateKey->GetPrivateElement());
    s.DropLastElements(s.GetNumOfElements() - cv[0].GetNumOfElements());
    DCRTPoly b(cv[0]);
    b.SetFormat(Format::EVALUATION);
    DCRTPoly c11(cv[1]);
    c11.SetFormat(Format::EVALUATION);
    b += s * c11;
    return PolyFromDCRTPoly(b);
}

Poly HomoEncryptCompute::myDecrypt(const Cipher& ciphertext) {
    return myDecrypt(ciphertext, key_pair.secretKey);
}

BigVector HomoEncryptCompute::LWEfromCiph(const Cipher& ciphertext) {
    BigVector c1vtemp = PolyFromDCRTPoly(ciphertext->GetElements()[1]).GetValues();
    BigVector c1v(c1vtemp);
    BigInteger q = c1v.GetModulus();

    for (size_t i = 1; i < c1vtemp.GetLength(); i++) {
        c1v[i] = q - c1vtemp[c1vtemp.GetLength() - i]; // OpenFHE doesn't know negative integers !
    }

    c1v[0] = (c1vtemp[0] + PolyFromDCRTPoly(ciphertext->GetElements()[0])[0]) % q;

    return c1v;
}

BigInteger HomoEncryptCompute::DecryptLWE(const Cipher& ciphertext, const PrivateKey<DCRTPoly>& privateKey) {
    BigVector c1v = LWEfromCiph(ciphertext);
    BigVector skv = PolyFromDCRTPoly(privateKey->GetPrivateElement()).GetValues();

    BigInteger val = 0;
    for (size_t i = 0; i < c1v.GetLength(); i++) {
        val += c1v[i] * skv[i];
    }

    return val % c1v.GetModulus();
}

BigInteger HomoEncryptCompute::DecryptLWE(const Cipher& ciphertext) {
    return DecryptLWE(ciphertext, key_pair.secretKey);
}

BigVector HomoEncryptCompute::DecryptLWE(const Cipher& ciphertext, const PrivateKey<DCRTPoly>& privateKey, int n) {
    BigVector skv = PolyFromDCRTPoly(privateKey->GetPrivateElement()).GetValues();
    
    size_t N=ciphertext->GetCryptoContext()->GetRingDimension();
    BigVector v(n);

    BigInteger q=ciphertext->GetElements()[0].GetModulus();
    v.SetModulus(q);

    for (int j = 0; j < n; j++) {
        BigVector cv=LWEfromCiph(ShiftRight(ciphertext,N-j*N/n));
        BigInteger val=0;
        for (size_t i = 0; i < N; i++) {
            val +=(q-cv[i]) * skv[i];
        }
        v[j] = val % q;
    }

    return v;
}

BigVector HomoEncryptCompute::DecryptLWE(const Cipher& ciphertext, int n) {
    return DecryptLWE(ciphertext, key_pair.secretKey, n);
}

Cipher HomoEncryptCompute::CiphertextConjugate(const Cipher& ciphertext) {
    uint32_t indexConj = 2 * context->GetRingDimension() - 1;
    auto evalConjKeyMap = context->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());
    return context->EvalAutomorphism(ciphertext, indexConj, evalConjKeyMap);
}

// ========== SlotsToCoeffs Operations ==========

// Shadow classes for accessing private members
class FHECKKSRNSShadow : public FHERNS {
public:
    const uint32_t K_SPARSE  = 28;  
    const uint32_t K_UNIFORM = 512;  
    const uint32_t K_UNIFORMEXT = 768;
    static const uint32_t R_UNIFORM = 6; 
    static const uint32_t R_SPARSE = 3;  
    uint32_t m_correctionFactor = 0; 
    std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> m_bootPrecomMap;

    // Chebyshev series coefficients for the SPARSE case
    static const inline std::vector<double> g_coefficientsSparse{
        -0.18646470117093214,   0.036680543700430925,    -0.20323558926782626,     0.029327390306199311,
        -0.24346234149506416,   0.011710240188138248,    -0.27023281815251715,     -0.017621188001030602,
        -0.21383614034992021,   -0.048567932060728937,   -0.013982336571484519,    -0.051097367628344978,
        0.24300487324019346,    0.0016547743046161035,   0.23316923792642233,      0.060707936480887646,
        -0.18317928363421143,   0.0076878773048247966,   -0.24293447776635235,     -0.071417413140564698,
        0.37747441314067182,    0.065154496937795681,    -0.24810721693607704,     -0.033588418808958603,
        0.10510660697380972,    0.012045222815124426,    -0.032574751830745423,    -0.0032761730196023873,
        0.0078689491066424744,  0.00070965574480802061,  -0.0015405394287521192,   -0.00012640521062948649,
        0.00025108496615830787, 0.000018944629154033562, -0.000034753284216308228, -2.4309868106111825e-6,
        4.1486274737866247e-6,  2.7079833113674568e-7,   -4.3245388569898879e-7,   -2.6482744214856919e-8,
        3.9770028771436554e-8,  2.2951153557906580e-9,   -3.2556026220554990e-9,   -1.7691071323926939e-10,
        2.5459052150406730e-10};

    // Chebyshev series coefficients for the OPTIMIZED/uniform case
    static const inline std::vector<double> g_coefficientsUniform{
        0.15421426400235561,    -0.0037671538417132409,  0.16032011744533031,      -0.0034539657223742453,
        0.17711481926851286,    -0.0027619720033372291,  0.19949802549604084,      -0.0015928034845171929,
        0.21756948616367638,    0.00010729951647566607,  0.21600427371240055,      0.0022171399198851363,
        0.17647500259573556,    0.0042856217194480991,   0.086174491919472254,     0.0054640252312780444,
        -0.046667988130649173,  0.0047346914623733714,   -0.17712686172280406,     0.0016205080004247200,
        -0.22703114241338604,   -0.0028145845916205865,  -0.13123089730288540,     -0.0056345646688793190,
        0.078818395388692147,   -0.0037868875028868542,  0.23226434602675575,      0.0021116338645426574,
        0.13985510526186795,    0.0059365649669377071,   -0.13918475289368595,     0.0018580676740836374,
        -0.23254376365752788,   -0.0054103844866927788,  0.056840618403875359,     -0.0035227192748552472,
        0.25667909012207590,    0.0055029673963982112,   -0.073334392714092062,    0.0027810273357488265,
        -0.24912792167850559,   -0.0069524866497120566,  0.21288810409948347,      0.0017810057298691725,
        0.088760951809475269,   0.0055957188940032095,   -0.31937177676259115,     -0.0087539416335935556,
        0.34748800245527145,    0.0075378299617709235,   -0.25116537379803394,     -0.0047285674679876204,
        0.13970502851683486,    0.0023672533925155220,   -0.063649401080083698,    -0.00098993213448982727,
        0.024597838934816905,   0.00035553235917057483,  -0.0082485030307578155,   -0.00011176184313622549,
        0.0024390574829093264,  0.000031180384864488629, -0.00064373524734389861,  -7.8036008952377965e-6,
        0.00015310015145922058, 1.7670804180220134e-6,   -0.000033066844379476900, -3.6460909134279425e-7,
        6.5276969021754105e-6,  6.8957843666189918e-8,   -1.1842811187642386e-6,   -1.2015133285307312e-8,
        1.9839339947648331e-7,  1.9372045971100854e-9,   -3.0815418032523593e-8,   -2.9013806338735810e-10,
        4.4540904298173700e-9,  4.0505136697916078e-11,  -6.0104912807134771e-10,  -5.2873323696828491e-12,
        7.5943206779351725e-11, 6.4679566322060472e-13,  -9.0081200925539902e-12,  -7.4396949275292252e-14,
        1.0057423059167244e-12, 8.1701187638005194e-15,  -1.0611736208855373e-13,  -8.9597492970451533e-16,
        1.1421575296031385e-14};

    // Chebyshev series coefficients for the COMPOSITESCALING case where d > 2
    const std::vector<double> g_coefficientsUniformExt{
        // New Coefficients (K_UNIFORM = 768)
        0.12602195635248634,    -0.0030834928649740388,  0.1293538007310393,      -0.0029150296085609707,
        0.13880323885842225,    -0.0025534902415420128,  0.15259900956315636,     -0.0019572806381606537,
        0.16740348080390202,    -0.0010852123927167594,  0.17795704156012629,     7.3594791671716396e-05,
        0.17708229644467954,    0.0014573280941530976,   0.15661113656175465,     0.0028850600459592078,
        0.10984969661272398,    0.0040295575406054489,   0.035829873357113948,    0.004449523200499763,
        -0.055520186697616318,  0.0037264589074560098,   -0.14007871037019429,    0.001719720247528076,
        -0.18281801001428047,   -0.0011373848818829857,  -0.15209319897288492,    -0.0037123962122311092,
        -0.043785371196750272,  -0.0045107273507656552,  0.09756154430583093,     -0.002604845726688627,
        0.18481556762187912,    0.0012462519210521535,   0.1403768476069214,      0.0043541760219966428,
        -0.024293645826662724,  0.0037846793397644275,   -0.17560536795332429,    -0.0005605968506360667,
        -0.1519811728143392,    -0.0045192348096649545,  0.048231020943727741,    -0.0032001529516056853,
        0.19692074387699257,    0.0024419388214462485,   0.078182928643403107,    0.0047838249172446005,
        -0.16476594792427054,   -0.00036614509861925492, -0.14537982038722122,    -0.0050995116137312257,
        0.13564231010825495,    -0.00050653194386865278, 0.16465075644913021,     0.0052831338103145531,
        -0.1493249604350485,    -0.00016209880585104635, -0.13934114757550983,    -0.0054247353644288178,
        0.20649654831497111,    0.0026431561325639561,   0.032277990808412343,    0.0039463054621702767,
        -0.23636345040634044,   -0.0059041496654351176,  0.17831596275657194,     0.0017594032442182191,
        0.05094162125752931,    0.0040150842221901416,   -0.24841268578463685,    -0.0073080801617375155,
        0.3122522704364516,     0.0073316847629231194,   -0.26606798599442621,    -0.0054892692910619113,
        0.17878607636323862,    0.0033586935001791839,   -0.10066311654486482,    -0.001754132071278842,
        0.049074577561330504,   0.00080234886593034873,  -0.021150143470356698,   -0.0003269871328764949,
        0.0081757002802533667,  0.00012021127618051574,  -0.0028652357611661534,  -4.0244300629116574e-05,
        0.00091801734966694636, 1.2361006806444711e-05,  -0.0002707191913116332,  -3.504631720275642e-06,
        7.3888955616723944e-05, 9.2189772261859728e-07,  -1.8752943907614565e-05, -2.2597387576370175e-07,
        4.4436168671606267e-06, 5.1807959456553769e-08,  -9.8651004908533913e-07, -1.1146078152883018e-08,
        2.0582706963882007e-07, 2.2568126993711184e-09,  -4.0469622058265335e-08, -4.31163542777443e-10,
        7.517057515198321e-09,  7.7904840375183328e-11,  -1.3219720621636946e-09, -1.3342979848924908e-11,
        2.2055962238660182e-10, 2.1724065123826773e-12,  -3.4974624736954921e-11, -3.3609296485004418e-13,
        5.2789108285402917e-12, 4.9471164793087018e-14,  -7.5998777765849013e-13, -4.2492853307002972e-15,
        1.0768090434260388e-13, -2.1478500584069139e-15, -1.3891315735425435e-14};

};

class FHECKKSRNSDerived : public FHECKKSRNS {
public:
    void EvalSlotsToCoeffsSetup(const CryptoContextImpl<DCRTPoly>& cc, std::vector<uint32_t> levelBudget,
                                        std::vector<uint32_t> dim1, uint32_t numSlots,uint32_t lDec) {
        
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

        uint32_t M     = cc.GetCyclotomicOrder();
        uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;

        bool precompute = true;

        // This is a workaround because m_correctionFactor and m_bootPrecomMap are private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
    
        uint32_t &m_correctionFactor=shadow->m_correctionFactor;                                    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
    
        m_correctionFactor = 9;

        m_bootPrecomMap[slots]                      = std::make_shared<CKKSBootstrapPrecom>();
        std::shared_ptr<CKKSBootstrapPrecom> precom = m_bootPrecomMap[slots];

        precom->m_slots = slots;
        precom->m_dim1  = dim1[0];

        uint32_t logSlots = std::log2(slots);
        // even for the case of a single slot we need one level for rescaling
        if (logSlots == 0) {
            logSlots = 1;
        }

        // Perform some checks on the level budget and compute parameters
        std::vector<uint32_t> newBudget = levelBudget;

        if (newBudget[1] > logSlots) {
            std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
         newBudget[1] = logSlots;
        }
        if (newBudget[1] < 1) {
           std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
           newBudget[1] = 1;
        }

        // precom->m_paramsEnc = GetCollapsedFFTParams(slots, newBudget[0], dim1[0]);
        precom->m_paramsDec = GetCollapsedFFTParams(slots, newBudget[1], dim1[1]);

        if (precompute) {
            uint32_t m    = 4 * slots;
            //bool isSparse = (M != m) ? true : false;

            // computes indices for all primitive roots of unity
            std::vector<uint32_t> rotGroup(slots);
            uint32_t fivePows = 1;
            for (uint32_t i = 0; i < slots; ++i) {
                rotGroup[i] = fivePows;
                fivePows *= 5;
                fivePows %= m;
            }

            // computes all powers of a primitive root of unity exp(2 * M_PI/m)
            std::vector<std::complex<double>> ksiPows(m + 1);
            for (uint32_t j = 0; j < m; ++j) {
                double angle = 2.0 * M_PI * j / m;
               ksiPows[j].real(cos(angle));
                ksiPows[j].imag(sin(angle));
            }
             ksiPows[m] = ksiPows[0];

             // Extract the modulus prior to bootstrapping
             NativeInteger q = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
            double qDouble  = q.ConvertToDouble();

            uint128_t factor = ((uint128_t)1 << ((uint32_t)std::round(std::log2(qDouble))));
            double pre       = qDouble / factor;
            double scaleDec  = 1 / pre;

            precom->m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false, scaleDec, lDec);

        }
    }

    void EvalStC_CtSSetup(const CryptoContextImpl<DCRTPoly>& cc, std::vector<uint32_t> levelBudget,
                                        std::vector<uint32_t> dim1, uint32_t numSlots,uint32_t lStC, uint32_t lCtS) {
        
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

        uint32_t M     = cc.GetCyclotomicOrder();
        uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;

        // This is a workaround because m_correctionFactor and m_bootPrecomMap are private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
    
        uint32_t &m_correctionFactor=shadow->m_correctionFactor;                                    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
        uint32_t compositeDegree = cryptoParams->GetCompositeDegree();                                    
        // 防止 m_correctionFactor为0
        if(m_correctionFactor == 0){
            m_correctionFactor = 9;
        }

        m_bootPrecomMap[slots]                      = std::make_shared<CKKSBootstrapPrecom>();
        std::shared_ptr<CKKSBootstrapPrecom> precom = m_bootPrecomMap[slots];

        precom->m_slots = slots;
        precom->m_dim1  = dim1[0];

        uint32_t logSlots = std::log2(slots);
        // even for the case of a single slot we need one level for rescaling
        if (logSlots == 0) {
            logSlots = 1;
        }

        // Perform some checks on the level budget and compute parameters
        std::vector<uint32_t> newBudget = levelBudget;

        if (newBudget[0] > logSlots) {
            std::cerr << "\nWarning, the level budget for encoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[0] = logSlots;
        }
        if (newBudget[0] < 1) {
            std::cerr << "\nWarning, the level budget for encoding can not be zero. Setting it to 1" << std::endl;
            newBudget[0] = 1;
        }

        if (newBudget[1] > logSlots) {
            std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[1] = logSlots;
        }
        if (newBudget[1] < 1) {
            std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
            newBudget[1] = 1;
        }

        precom->m_paramsEnc = GetCollapsedFFTParams(slots, newBudget[0], dim1[0]);
        precom->m_paramsDec = GetCollapsedFFTParams(slots, newBudget[1], dim1[1]);

        uint32_t m    = 4 * slots;
        bool isSparse = (M != m) ? true : false;

        // computes indices for all primitive roots of unity
        std::vector<uint32_t> rotGroup(slots);
        uint32_t fivePows = 1;
        for (uint32_t i = 0; i < slots; ++i) {
            rotGroup[i] = fivePows;
            fivePows *= 5;
            fivePows %= m;
        }

        // computes all powers of a primitive root of unity exp(2 * M_PI/m)
        std::vector<std::complex<double>> ksiPows(m + 1);
        for (uint32_t j = 0; j < m; ++j) {
            double angle = 2.0 * M_PI * j / m;
            ksiPows[j].real(cos(angle));
            ksiPows[j].imag(sin(angle));
        }
            ksiPows[m] = ksiPows[0];

        // Extract the modulus prior to bootstrapping
        NativeInteger q = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
        double qDouble  = q.ConvertToDouble();

        uint128_t factor = ((uint128_t)1 << ((uint32_t)std::round(std::log2(qDouble))));
        double pre       = (compositeDegree > 1) ? 1.0 : qDouble / factor;
        double k         = (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) ? shadow->K_SPARSE : 1.0;
        // StC和CtS的缩放因子
        double scaleDec = (compositeDegree > 1) ? qDouble / cryptoParams->GetScalingFactorReal(0) : 1 / pre;
        std ::cout << "scaleStC: " << scaleDec << std::endl;
        double scaleEnc  = pre / k;
        std ::cout << "scaleCtS: " << scaleEnc << std::endl;

        bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                            (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);

        if (isLTBootstrap) {
            // allocate all vectors
            std::vector<std::vector<std::complex<double>>> U0(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U1(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U0hatT(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U1hatT(slots, std::vector<std::complex<double>>(slots));

            for (size_t i = 0; i < slots; i++) {
                for (size_t j = 0; j < slots; j++) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) % m];
                    U0hatT[j][i] = std::conj(U0[i][j]);
                    U1[i][j]     = std::complex<double>(0, 1) * U0[i][j];
                    U1hatT[j][i] = std::conj(U1[i][j]);
                }
            }

            if (!isSparse) {
                precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, scaleDec, lStC);
                precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, scaleEnc, lCtS);
            }
            else {
                precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, U1, 1, scaleDec, lStC);
                precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleEnc, lCtS);
            }
        }
        else {
            precom->m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false, scaleDec, lStC);
            precom->m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleEnc, lCtS);
        }
    }
    // 用于设置StC first和CtS first BTS的自举参数
    void EvalStC_CtSSetup_BTS(const CryptoContextImpl<DCRTPoly>& cc, std::vector<uint32_t> levelBudget,
                                        std::vector<uint32_t> dim1, uint32_t numSlots,bool isStCFirst) {
        
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

        uint32_t M     = cc.GetCyclotomicOrder();
        uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;
        uint32_t lStC = 1;
        uint32_t lCtS = 1;

        // This is a workaround because m_correctionFactor and m_bootPrecomMap are private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
    
        uint32_t &m_correctionFactor=shadow->m_correctionFactor;                                    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
        uint32_t compositeDegree = cryptoParams->GetCompositeDegree();                                    
        // 防止 m_correctionFactor为0
        if(m_correctionFactor == 0){
            m_correctionFactor = 9;
        }

        m_bootPrecomMap[slots]                      = std::make_shared<CKKSBootstrapPrecom>();
        std::shared_ptr<CKKSBootstrapPrecom> precom = m_bootPrecomMap[slots];

        precom->m_slots = slots;
        precom->m_dim1  = dim1[0];

        uint32_t logSlots = std::log2(slots);
        // even for the case of a single slot we need one level for rescaling
        if (logSlots == 0) {
            logSlots = 1;
        }

        // Perform some checks on the level budget and compute parameters
        std::vector<uint32_t> newBudget = levelBudget;

        if (newBudget[0] > logSlots) {
            std::cerr << "\nWarning, the level budget for encoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[0] = logSlots;
        }
        if (newBudget[0] < 1) {
            std::cerr << "\nWarning, the level budget for encoding can not be zero. Setting it to 1" << std::endl;
            newBudget[0] = 1;
        }

        if (newBudget[1] > logSlots) {
            std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[1] = logSlots;
        }
        if (newBudget[1] < 1) {
            std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
            newBudget[1] = 1;
        }

        precom->m_paramsEnc = GetCollapsedFFTParams(slots, newBudget[0], dim1[0]);
        precom->m_paramsDec = GetCollapsedFFTParams(slots, newBudget[1], dim1[1]);

        uint32_t m    = 4 * slots;
        bool isSparse = (M != m) ? true : false;

        // computes indices for all primitive roots of unity
        std::vector<uint32_t> rotGroup(slots);
        uint32_t fivePows = 1;
        for (uint32_t i = 0; i < slots; ++i) {
            rotGroup[i] = fivePows;
            fivePows *= 5;
            fivePows %= m;
        }

        // computes all powers of a primitive root of unity exp(2 * M_PI/m)
        std::vector<std::complex<double>> ksiPows(m + 1);
        for (uint32_t j = 0; j < m; ++j) {
            double angle = 2.0 * M_PI * j / m;
            ksiPows[j].real(cos(angle));
            ksiPows[j].imag(sin(angle));
        }
            ksiPows[m] = ksiPows[0];

        // Extract the modulus prior to bootstrapping
        NativeInteger q = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
        double qDouble  = q.ConvertToDouble();

        uint128_t factor = ((uint128_t)1 << ((uint32_t)std::round(std::log2(qDouble))));
        double pre       = (compositeDegree > 1) ? 1.0 : qDouble / factor;
        double k         = (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) ? shadow->K_SPARSE : 1.0;
        // StC和CtS的缩放因子
        double scaleDec = (compositeDegree > 1) ? qDouble / cryptoParams->GetScalingFactorReal(0) : 1 / pre;
        std ::cout << "scaleStC: " << scaleDec << std::endl;
        double scaleEnc  = pre / k;
        std ::cout << "scaleCtS: " << scaleEnc << std::endl;

        // 设置StC和CtS的层数
        uint32_t approxModDepth = GetModDepthInternal(cryptoParams->GetSecretKeyDist());
        uint32_t depthBT        = approxModDepth + precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] +
                           precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
        // compute # of levels to remain when encoding the coefficients
        uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();
        // for FLEXIBLEAUTOEXT we do not need extra modulus in auxiliary plaintexts
        if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
            L0 -= 1;
        // 判断自举类型--StC/CtS first
        if(isStCFirst) {
            // StC first -- 还需要仔细琢磨
            lStC = 2; // StC之后预留一层--因为有第一次乘法不缩放的操作，所以后续需要对这个地方进行拓展
            // lCtS = L0 - compositeDegree * depthBT + approxModDepth;
            lCtS = L0 - precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] - 1;
        }
        else {
            // CtS first
            lStC = L0 - compositeDegree * depthBT;
            lCtS = L0 - compositeDegree * (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] + 1);
        }
        std::cout << "SET---lStC: " << lStC << ", lCtS: " << lCtS << std::endl;

        bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                            (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);

        if (isLTBootstrap) {
            // allocate all vectors
            std::vector<std::vector<std::complex<double>>> U0(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U1(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U0hatT(slots, std::vector<std::complex<double>>(slots));
            std::vector<std::vector<std::complex<double>>> U1hatT(slots, std::vector<std::complex<double>>(slots));

            for (size_t i = 0; i < slots; i++) {
                for (size_t j = 0; j < slots; j++) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) % m];
                    U0hatT[j][i] = std::conj(U0[i][j]);
                    U1[i][j]     = std::complex<double>(0, 1) * U0[i][j];
                    U1hatT[j][i] = std::conj(U1[i][j]);
                }
            }

            if (!isSparse) {
                precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, scaleDec, lStC);
                precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, scaleEnc, lCtS);
            }
            else {
                precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, U1, 1, scaleDec, lStC);
                precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleEnc, lCtS);
            }
        }
        else {
            precom->m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false, scaleDec, lStC);
            precom->m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleEnc, lCtS);
        }
    }

    // ========= Bootstrapping Setup ==========
    void EvalBTSSetup(const CryptoContextImpl<DCRTPoly>& cc, std::vector<uint32_t> levelBudget,
                        std::vector<uint32_t> dim1, uint32_t numSlots,
                        uint32_t correctionFactor, bool isStCFirst) {
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

        if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
            OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
    #if NATIVEINT == 128
        if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
            OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
    #endif

        uint32_t M     = cc.GetCyclotomicOrder();
        uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;
        
        // This is a workaround because m_correctionFactor and m_bootPrecomMap are private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);                    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
        uint32_t &m_correctionFactor=shadow->m_correctionFactor;
        // Set correction factor by default, if it is not already set.
        if (correctionFactor == 0) {
            if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO ||
                cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT ||
                cryptoParams->GetScalingTechnique() == COMPOSITESCALINGAUTO ||
                cryptoParams->GetScalingTechnique() == COMPOSITESCALINGMANUAL) {
                // The default correction factors chosen yielded the best precision in our experiments.
                // We chose the best fit line from our experiments by running ckks-bootstrapping-precision.cpp.
                // The spreadsheet with our experiments is here:
                // https://docs.google.com/spreadsheets/d/1WqmwBUMNGlX6Uvs9qLXt5yeddtCyWPP55BbJPu5iPAM/edit?usp=sharing
                auto tmp = std::round(-0.265 * (2 * std::log2(M / 2) + std::log2(slots)) + 19.1);
                if (tmp < 7)
                    m_correctionFactor = 7;
                else if (tmp > 13)
                    m_correctionFactor = 13;
                else
                    m_correctionFactor = static_cast<uint32_t>(tmp);
            }
            else {
                m_correctionFactor = 9;
            }
        }
        else {
            m_correctionFactor = correctionFactor;
        }
        m_bootPrecomMap[slots]                      = std::make_shared<CKKSBootstrapPrecom>();
        std::shared_ptr<CKKSBootstrapPrecom> precom = m_bootPrecomMap[slots];
        // 计算CtS first的StC和CtS预计算参数
        EvalStC_CtSSetup_BTS(cc, levelBudget, dim1, numSlots, isStCFirst);
        
    }

    void EvalBTSSetupOrigin(const CryptoContextImpl<DCRTPoly>& cc, std::vector<uint32_t> levelBudget,
                        std::vector<uint32_t> dim1, uint32_t numSlots,
                        uint32_t correctionFactor, bool level_set) {
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

        if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
            OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
    #if NATIVEINT == 128
        if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
            OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
    #endif

        uint32_t M     = cc.GetCyclotomicOrder();
        uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;
        
        // This is a workaround because m_correctionFactor and m_bootPrecomMap are private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);                    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
        uint32_t &m_correctionFactor=shadow->m_correctionFactor;
        // Set correction factor by default, if it is not already set.
        if (correctionFactor == 0) {
            if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO ||
                cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT ||
                cryptoParams->GetScalingTechnique() == COMPOSITESCALINGAUTO ||
                cryptoParams->GetScalingTechnique() == COMPOSITESCALINGMANUAL) {
                // The default correction factors chosen yielded the best precision in our experiments.
                // We chose the best fit line from our experiments by running ckks-bootstrapping-precision.cpp.
                // The spreadsheet with our experiments is here:
                // https://docs.google.com/spreadsheets/d/1WqmwBUMNGlX6Uvs9qLXt5yeddtCyWPP55BbJPu5iPAM/edit?usp=sharing
                auto tmp = std::round(-0.265 * (2 * std::log2(M / 2) + std::log2(slots)) + 19.1);
                if (tmp < 7)
                    m_correctionFactor = 7;
                else if (tmp > 13)
                    m_correctionFactor = 13;
                else
                    m_correctionFactor = static_cast<uint32_t>(tmp);
            }
            else {
                m_correctionFactor = 9;
            }
        }
        else {
            m_correctionFactor = correctionFactor;
        }
        m_bootPrecomMap[slots]                      = std::make_shared<CKKSBootstrapPrecom>();
        std::shared_ptr<CKKSBootstrapPrecom> precom = m_bootPrecomMap[slots];

        precom->m_slots = slots;
        precom->m_dim1  = dim1[0];

        uint32_t logSlots = std::log2(slots);
        // even for the case of a single slot we need one level for rescaling
        if (logSlots == 0) {
            logSlots = 1;
        }

        // Perform some checks on the level budget and compute parameters
        std::vector<uint32_t> newBudget = levelBudget;

        if (newBudget[0] > logSlots) {
            std::cerr << "\nWarning, the level budget for encoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[0] = logSlots;
        }
        if (newBudget[0] < 1) {
            std::cerr << "\nWarning, the level budget for encoding can not be zero. Setting it to 1" << std::endl;
            newBudget[0] = 1;
        }

        if (newBudget[1] > logSlots) {
            std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
            newBudget[1] = logSlots;
        }
        if (newBudget[1] < 1) {
            std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
            newBudget[1] = 1;
        }

        precom->m_paramsEnc = GetCollapsedFFTParams(slots, newBudget[0], dim1[0]);
        precom->m_paramsDec = GetCollapsedFFTParams(slots, newBudget[1], dim1[1]);

        if (level_set) {
            uint32_t m    = 4 * slots;
            bool isSparse = (M != m) ? true : false;

            // computes indices for all primitive roots of unity
            std::vector<uint32_t> rotGroup(slots);
            uint32_t fivePows = 1;
            for (uint32_t i = 0; i < slots; ++i) {
                rotGroup[i] = fivePows;
                fivePows *= 5;
                fivePows %= m;
            }

            // computes all powers of a primitive root of unity exp(2 * M_PI/m)
            std::vector<std::complex<double>> ksiPows(m + 1);
            for (uint32_t j = 0; j < m; ++j) {
                double angle = 2.0 * M_PI * j / m;
                ksiPows[j].real(cos(angle));
                ksiPows[j].imag(sin(angle));
            }
            ksiPows[m] = ksiPows[0];

            uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

            // Extract the modulus prior to bootstrapping
            NativeInteger q = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
            double qDouble  = q.ConvertToDouble();

            uint128_t factor = ((uint128_t)1 << (static_cast<uint32_t>(std::round(std::log2(qDouble)))));
            double pre       = (compositeDegree > 1) ? 1.0 : qDouble / factor;
            double k         = (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) ? shadow->K_SPARSE : 1.0;
            double scaleEnc  = pre / k;
            // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
            double scaleDec = (compositeDegree > 1) ? qDouble / cryptoParams->GetScalingFactorReal(0) : 1 / pre;

            uint32_t approxModDepth = GetModDepthInternal(cryptoParams->GetSecretKeyDist());
            uint32_t depthBT        = approxModDepth + precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] +
                            precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];

            // compute # of levels to remain when encoding the coefficients
            uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();
            // for FLEXIBLEAUTOEXT we do not need extra modulus in auxiliary plaintexts
            if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
                L0 -= 1;
            std::cout << "L0: " << L0 << std::endl;
            uint32_t lCtS = L0 - compositeDegree * (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] + 1);
            uint32_t lStC = L0 - compositeDegree * depthBT;
            std::cout << "lStC: " << lStC << ", lCtS: " << lCtS << std::endl;
            bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                                (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);

            if (isLTBootstrap) {
                // allocate all vectors
                std::vector<std::vector<std::complex<double>>> U0(slots, std::vector<std::complex<double>>(slots));
                std::vector<std::vector<std::complex<double>>> U1(slots, std::vector<std::complex<double>>(slots));
                std::vector<std::vector<std::complex<double>>> U0hatT(slots, std::vector<std::complex<double>>(slots));
                std::vector<std::vector<std::complex<double>>> U1hatT(slots, std::vector<std::complex<double>>(slots));

                for (size_t i = 0; i < slots; i++) {
                    for (size_t j = 0; j < slots; j++) {
                        U0[i][j]     = ksiPows[(j * rotGroup[i]) % m];
                        U0hatT[j][i] = std::conj(U0[i][j]);
                        U1[i][j]     = std::complex<double>(0, 1) * U0[i][j];
                        U1hatT[j][i] = std::conj(U1[i][j]);
                    }
                }

                if (!isSparse) {
                    precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, scaleEnc, lCtS);
                    precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, scaleDec, lStC);
                }
                else {
                    precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleEnc, lCtS);
                    precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, U1, 1, scaleDec, lStC);
                }
            }
            else {
                precom->m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleEnc, lCtS);
                precom->m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false, scaleDec, lStC);
            }
        }

    }

    Cipher EvalSlotsToCoeffs(ConstCiphertext<DCRTPoly> ctxt)  {
        uint32_t slots = ctxt->GetSlots();

        // This is a workaround because m_bootPrecomMap is private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
        
        auto pair = m_bootPrecomMap.find(slots);
        if (pair == m_bootPrecomMap.end()) {
            std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                             std::string(" slots were not generated") +
                             std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
            OPENFHE_THROW(errorMsg);
        }
        const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

        return FHECKKSRNS::EvalSlotsToCoeffs(precom->m_U0PreFFT, ctxt);
    }

    Cipher EvalCoeffsToSlots(ConstCiphertext<DCRTPoly> ctxt)  {
        uint32_t slots = ctxt->GetSlots();

        // This is a workaround because m_bootPrecomMap is private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
    
        auto pair = m_bootPrecomMap.find(slots);
        if (pair == m_bootPrecomMap.end()) {
            std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                             std::string(" slots were not generated") +
                             std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
            OPENFHE_THROW(errorMsg);
        }
        const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

        return FHECKKSRNS::EvalCoeffsToSlots(precom->m_U0hatTPreFFT, ctxt);
    }

    Cipher EvalLinearTransformCtS(ConstCiphertext<DCRTPoly> ctxt)  {
        uint32_t slots = ctxt->GetSlots();

        // This is a workaround because m_bootPrecomMap is private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;

        auto pair = m_bootPrecomMap.find(slots);
        if (pair == m_bootPrecomMap.end()) {
            std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                                std::string(" slots were not generated") +
                                std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
            OPENFHE_THROW(errorMsg);
        }
        const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

        return FHECKKSRNS::EvalLinearTransform(precom->m_U0hatTPre, ctxt);
    }

    Cipher EvalLinearTransformStC(ConstCiphertext<DCRTPoly> ctxt)  {
        uint32_t slots = ctxt->GetSlots();

        // This is a workaround because m_bootPrecomMap is private in FHECKKSRNS
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);    
        std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;

        auto pair = m_bootPrecomMap.find(slots);
        if (pair == m_bootPrecomMap.end()) {
            std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                                std::string(" slots were not generated") +
                                std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
            OPENFHE_THROW(errorMsg);
        }
        const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

        return FHECKKSRNS::EvalLinearTransform(precom->m_U0Pre, ctxt);
    }

    Cipher Conjugate(ConstCiphertext<DCRTPoly> ciphertext,
                                   std::map<uint32_t, EvalKey<DCRTPoly>>& evalKeyMap){
        const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
        uint32_t N                      = cv[0].GetRingDimension();

        std::vector<uint32_t> vec(N);
        PrecomputeAutoMap(N, 2 * N - 1, &vec);

        auto algo = ciphertext->GetCryptoContext()->GetScheme();

        Ciphertext<DCRTPoly> result = ciphertext->Clone();

        algo->KeySwitchInPlace(result, evalKeyMap.at(2 * N - 1));

        std::vector<DCRTPoly>& rcv = result->GetElements();

        rcv[0] = rcv[0].AutomorphismTransform(2 * N - 1, vec);
        rcv[1] = rcv[1].AutomorphismTransform(2 * N - 1, vec);

        return result;
    }

private:
    uint32_t GetModDepthInternal(SecretKeyDist secretKeyDist) {
        auto shadow=reinterpret_cast<FHECKKSRNSShadow*>(this);
        if (secretKeyDist == UNIFORM_TERNARY) {
            return GetMultiplicativeDepthByCoeffVector(shadow->g_coefficientsUniform, false) + shadow->R_UNIFORM;
        }
        else {
            return GetMultiplicativeDepthByCoeffVector(shadow->g_coefficientsSparse, false) + shadow->R_SPARSE;
        }
    }

};

std::shared_ptr<FHECKKSRNS> HomoEncryptCompute::getFHEAlgorithm(const CryptoContext<DCRTPoly>& cryptoContext) {
    std::shared_ptr<SchemeBase<DCRTPoly>> scheme=cryptoContext->GetScheme();
    
    // This is because m_FHE is protected in SchemeBase<DCRTPoly>
    struct Shadow : public SchemeBase<DCRTPoly> {
        using SchemeBase<DCRTPoly>::m_FHE;
    };

    auto baseAlgo = static_cast<const Shadow&>(*scheme).m_FHE;

    if(!baseAlgo) {
        throw std::runtime_error("Failed to get FHEBase<DCRTPoly>");
    }

    std::shared_ptr<FHECKKSRNS> algo = std::dynamic_pointer_cast<FHECKKSRNS>(baseAlgo);
    if (!algo) {
        throw std::runtime_error("Failed to cast FHEBase<DCRTPoly> to FHECKKSRNS");
    }
    return algo;
}

void HomoEncryptCompute::EvalSlotsToCoeffsSetup(uint32_t levelBudget, uint32_t numSlots, uint32_t lDec) {
    std::vector<uint32_t> dim1 = {0, 0};
    std::vector<uint32_t> levelBudget2 = {0,levelBudget};
    
    std::shared_ptr<FHECKKSRNS> algo= getFHEAlgorithm(context);
    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    algo2.EvalSlotsToCoeffsSetup(*context, levelBudget2, dim1, numSlots,lDec);
}

void HomoEncryptCompute::EvalStC_CtSSetup(std::vector<uint32_t> levelBudget, uint32_t numSlots, uint32_t lStC, uint32_t lCtS) {
    std::vector<uint32_t> dim1 = {0, 0};
    
    std::shared_ptr<FHECKKSRNS> algo= getFHEAlgorithm(context);
    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    algo2.EvalStC_CtSSetup(*context, levelBudget, dim1, numSlots,lStC, lCtS);
}

void HomoEncryptCompute::EvalBTSSetup(vector<uint32_t> levelBudget, vector<uint32_t> dim1, 
                        uint32_t numSlots, uint32_t correctionFactor, bool isStCFirst) {
    std::shared_ptr<FHECKKSRNS> algo= getFHEAlgorithm(context);
    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    algo2.EvalBTSSetup(*context, levelBudget, dim1, numSlots, correctionFactor, isStCFirst);
    // algo2.EvalBTSSetupOrigin(*context, levelBudget, dim1, numSlots, correctionFactor, level_set);
}

Cipher HomoEncryptCompute::SlotsToCoeffs(const Cipher& ciph) {
    std::shared_ptr<FHECKKSRNS> algo = getFHEAlgorithm(context);

    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    Cipher ciphout = algo2.EvalSlotsToCoeffs(ciph);
    std::cout << "StC Finish!" << std::endl; 
    // 稀疏
    if(ciph->GetSlots() != context->GetCyclotomicOrder()/4){
        std::cout << "StC Rotate!" << std::endl;
        context->EvalAddInPlace(ciphout, context->EvalRotate(ciphout, ciphout->GetSlots()));
    } 
    
    return ciphout;
}

Cipher HomoEncryptCompute::CoeffsToSlots(const Cipher& ciph) {
    std::shared_ptr<FHECKKSRNS> algo = getFHEAlgorithm(context);

    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    Cipher ciphout = algo2.EvalCoeffsToSlots(ciph);
    std::cout << "CtS Finish!" << std::endl; 
    
    return ciphout;
}

Cipher HomoEncryptCompute::LinearTransformCtS(const Cipher& ciph){
    std::shared_ptr<FHECKKSRNS> algo = getFHEAlgorithm(context);

    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    std::cout << "CtS LinerTransformation Start!" << std::endl; 
    Cipher ciphout = algo2.EvalLinearTransformCtS(ciph);
    std::cout << "CtS LinerTransformation Finish!" << std::endl; 
    
    return ciphout;
}

Cipher HomoEncryptCompute::LinearTransformStC(const Cipher& ciph){
    std::shared_ptr<FHECKKSRNS> algo = getFHEAlgorithm(context);

    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    std::cout << "StC LinerTransformation Start!" << std::endl; 
    Cipher ciphout = algo2.EvalLinearTransformStC(ciph);
    std::cout << "StC LinerTransformation Finish!" << std::endl; 

    // 稀疏
    if(ciph->GetSlots() != context->GetCyclotomicOrder()/4){
        context->EvalAddInPlace(ciphout, context->EvalRotate(ciphout, ciphout->GetSlots()));
    } 
    
    return ciphout;
}

// ========== Bootstrapping Operations ==========
double HomoEncryptCompute::GetBigModulus(const std::shared_ptr<lbcrypto::CryptoParametersCKKSRNS> cryptoParams) {
    double qDouble           = 1.0;
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
    for (uint32_t j = 0; j < compositeDegree; ++j) {
        qDouble *= cryptoParams->GetElementParams()->GetParams()[j]->GetModulus().ConvertToDouble();
    }

    return qDouble;
}

void HomoEncryptCompute::ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ciphertext, uint32_t numIter) {
    int32_t r = numIter;
    for (int32_t j = 1; j < r + 1; j++) {
        context->EvalSquareInPlace(ciphertext);
        ciphertext    = context->EvalAdd(ciphertext, ciphertext);
        double scalar = -1.0 / std::pow((2.0 * M_PI), std::pow(2.0, j - r));
        context->EvalAddInPlace(ciphertext, scalar);
        context->ModReduceInPlace(ciphertext);
    }
}

void HomoEncryptCompute::AdjustCiphertext(Ciphertext<DCRTPoly>& ciphertext, double correction){
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto algo                = context->GetScheme();
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT ||
        cryptoParams->GetScalingTechnique() == COMPOSITESCALINGAUTO ||
        cryptoParams->GetScalingTechnique() == COMPOSITESCALINGMANUAL) {
        uint32_t lvl       = cryptoParams->GetScalingTechnique() != FLEXIBLEAUTOEXT ? 0 : 1;
        double targetSF    = cryptoParams->GetScalingFactorReal(lvl);
        double sourceSF    = ciphertext->GetScalingFactor();
        uint32_t numTowers = ciphertext->GetElements()[0].GetNumOfElements();
        double modToDrop = cryptoParams->GetElementParams()->GetParams()[numTowers - 1]->GetModulus().ConvertToDouble();
        for (uint32_t j = 2; j <= compositeDegree; ++j) {
            modToDrop *= cryptoParams->GetElementParams()->GetParams()[numTowers - j]->GetModulus().ConvertToDouble();
        }

        // in the case of FLEXIBLEAUTO, we need to bring the ciphertext to the right scale using a
        // a scaling multiplication. Note the at currently FLEXIBLEAUTO is only supported for NATIVEINT = 64.
        // So the other branch is for future purposes (in case we decide to add add the FLEXIBLEAUTO support
        // for NATIVEINT = 128.
#if NATIVEINT != 128
        // Scaling down the message by a correction factor to emulate using a larger q0.
        // This step is needed so we could use a scaling factor of up to 2^59 with q9 ~= 2^60.
        double adjustmentFactor = (targetSF / sourceSF) * (modToDrop / sourceSF) * std::pow(2, -correction);
#else
        double adjustmentFactor = (targetSF / sourceSF) * (modToDrop / sourceSF);
#endif
        context->EvalMultInPlace(ciphertext, adjustmentFactor);

        algo->ModReduceInternalInPlace(ciphertext, compositeDegree);
        ciphertext->SetScalingFactor(targetSF);
    }
    else {
#if NATIVEINT != 128
        // Scaling down the message by a correction factor to emulate using a larger q0.
        // This step is needed so we could use a scaling factor of up to 2^59 with q9 ~= 2^60.
        context->EvalMultInPlace(ciphertext, std::pow(2, -correction));
        algo->ModReduceInternalInPlace(ciphertext, compositeDegree);
#endif
    }
}

Cipher HomoEncryptCompute::conjugate(const Cipher& ciphertext) {
    std::shared_ptr<FHECKKSRNS> algo = getFHEAlgorithm(context);

    FHECKKSRNSDerived &algo2=static_cast<FHECKKSRNSDerived&>(*algo);

    auto evalKeyMap = context->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());

    return algo2.Conjugate(ciphertext, evalKeyMap);
}

void HomoEncryptCompute::ExtendCiphertext(std::vector<DCRTPoly>& ctxtDCRT, const CryptoContextImpl<DCRTPoly>& cc,
                                  const std::shared_ptr<DCRTPoly::Params> elementParamsRaisedPtr) {
    // TODO: YSP We should be able to use one of the DCRTPoly methods for this; If not, we can define a new method there and use it here

    // CompositeDegree = 2: [a]_q0q1     =     [a*q1^-1]_q0 *     q1 + [a*q0^-1]_q1 *q0
    // CompositeDegree = 3: [a]_q0q1q2   =   [a*q1q2^-1]_q0 *   q1q2 + [a*q0q2^-1]_q1 *q0q2 + [a*q0q1^-1]_q2 *q0q1
    // CompositeDegree = 4: [a]_q0q1q2q3 = [a*q1q2q3^-1]_q0 * q1q2q3 + [a*q0q2q3^-1]_q1 * q0q2q3 + [a*q0q1q3^-1]_q2 * q0q1q3 + [a*q0q1q2^-1]_q3 * q0q1q2

    const auto cryptoParams  = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    std::vector<NativeInteger> qj(compositeDegree);
    for (uint32_t j = 0; j < compositeDegree; ++j) {
        qj[j] = elementParamsRaisedPtr->GetParams()[j]->GetModulus().ConvertToInt();
    }

    std::vector<NativeInteger> qhat_modqj(compositeDegree);
    qhat_modqj[0] = qj[1].Mod(qj[0]);
    qhat_modqj[1] = qj[0].Mod(qj[1]);

    std::vector<NativeInteger> qhat_inv_modqj(compositeDegree);

    for (uint32_t d = 2; d < compositeDegree; d++) {
        for (uint32_t j = 0; j < d; ++j) {
            qhat_modqj[j] = qj[d].ModMul(qhat_modqj[j], qj[j]);
        }
        qhat_modqj[d] = qj[1].ModMul(qj[0], qj[d]);
        for (uint32_t j = 2; j < d; ++j) {
            qhat_modqj[d] = qj[j].ModMul(qhat_modqj[d], qj[d]);
        }
    }

    for (uint32_t j = 0; j < compositeDegree; ++j) {
        qhat_inv_modqj[j] = qhat_modqj[j].ModInverse(qj[j]);
    }

    NativeInteger qjProduct =
        std::accumulate(qj.begin() + 1, qj.end(), NativeInteger{1}, std::multiplies<NativeInteger>());
    uint32_t init_element_index = compositeDegree;
    for (size_t i = 0; i < ctxtDCRT.size(); i++) {
        std::vector<DCRTPoly> temp(compositeDegree + 1, DCRTPoly(elementParamsRaisedPtr, COEFFICIENT));
        std::vector<DCRTPoly> ctxtDCRT_modq(compositeDegree, DCRTPoly(elementParamsRaisedPtr, COEFFICIENT));

        ctxtDCRT[i].SetFormat(COEFFICIENT);
        for (size_t j = 0; j < ctxtDCRT[i].GetNumOfElements(); j++) {
            for (size_t k = 0; k < compositeDegree; k++)
                ctxtDCRT_modq[k].SetElementAtIndex(j, ctxtDCRT[i].GetElementAtIndex(j) * qhat_inv_modqj[k]);
        }
        //=========================================================================================================
        temp[0] = ctxtDCRT_modq[0].GetElementAtIndex(0);
        for (auto& el : temp[0].GetAllElements()) {
            el *= qjProduct;
        }
        //=========================================================================================================
        for (size_t d = 1; d < compositeDegree; d++) {
            temp[init_element_index] = ctxtDCRT_modq[d].GetElementAtIndex(d);

            for (size_t k = 0; k < compositeDegree; k++) {
                if (k != d) {
                    temp[d].SetElementAtIndex(k, temp[0].GetElementAtIndex(k) * qj[k]);
                }
            }
            //=========================================================================================================
            NativeInteger qjProductD{1};
            for (size_t k = 0; k < compositeDegree; k++) {
                if (k != d)
                    qjProductD *= qj[k];
            }

            for (size_t j = compositeDegree; j < elementParamsRaisedPtr->GetParams().size(); j++) {
                auto value = temp[init_element_index].GetElementAtIndex(j) * qjProductD;
                temp[d].SetElementAtIndex(j, value);
            }
            //=========================================================================================================
            {
                auto value = temp[init_element_index].GetElementAtIndex(d) * qjProductD;
                temp[d].SetElementAtIndex(d, value);
            }
            //=========================================================================================================
            temp[0] += temp[d];
        }

        temp[0].SetFormat(EVALUATION);
        ctxtDCRT[i] = temp[0];
    }
}



Cipher HomoEncryptCompute::bootstrap(Cipher &ciphertext) {
    return context->EvalBootstrap(ciphertext);
}
// 自举前的一些其他设置
void HomoEncryptCompute::bootothersetup() {
    // 这里可以添加一些其他的设置
    // 例如，设置一些参数或者预处理数据等
    // 目前没有具体的实现
}
// Mod Raise
Cipher HomoEncryptCompute::modraise(Cipher& ciphertext, const uint32_t& correction, const uint32_t& compositeDegree){
    // In FLEXIBLEAUTO, raising the ciphertext to a larger number
    // of towers is a bit more complex, because we need to adjust
    // it's scaling factor to the one that corresponds to the level
    // it's being raised to.
    // Increasing the modulus
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());
    
    uint32_t M = context->GetCyclotomicOrder();
    uint32_t L0              = cryptoParams->GetElementParams()->GetParams().size();

    auto elementParamsRaised = *(cryptoParams->GetElementParams());
    // For FLEXIBLEAUTOEXT we raised ciphertext does not include extra modulus
    // as it is multiplied by auxiliary plaintext
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT) {
        elementParamsRaised.PopLastParam();
    }
    auto paramsQ   = elementParamsRaised.GetParams();
    uint32_t sizeQ = paramsQ.size();

    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }
    auto elementParamsRaisedPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);

    auto raised = ciphertext->Clone();
    auto algo = context->GetScheme();
    algo->ModReduceInternalInPlace(raised, compositeDegree * (raised->GetNoiseScaleDeg() - 1));
    
    AdjustCiphertext(raised, correction);
    auto ctxtDCRT = raised->GetElements();

    if (compositeDegree > 1) {
        // RNS basis extension from level 0 RNS limbs to the raised RNS basis
        ExtendCiphertext(ctxtDCRT, *context, elementParamsRaisedPtr);
    }
    else {
        // We only use the level 0 ciphertext here. All other towers are automatically ignored to make
        // CKKS bootstrapping faster.
        for (size_t i = 0; i < ctxtDCRT.size(); i++) {
            DCRTPoly temp(elementParamsRaisedPtr, COEFFICIENT);
            ctxtDCRT[i].SetFormat(COEFFICIENT);
            temp = ctxtDCRT[i].GetElementAtIndex(0);
            temp.SetFormat(EVALUATION);
            ctxtDCRT[i] = temp;
        }
    }

    raised->SetLevel(L0 - ctxtDCRT[0].GetNumOfElements());
    raised->SetElements(std::move(ctxtDCRT));

    return raised;
}
// 部分和，用于稀疏打包
Cipher HomoEncryptCompute::partialsum(Cipher& ciphertext, const size_t& N, const size_t& n){
    for (uint32_t j = 1; j < N / (2 * n); j <<= 1) {
            auto temp = context->EvalRotate(ciphertext, j * n);
            context->EvalAddInPlace(ciphertext, temp);
        }
    return ciphertext;
}
// CtS first BTS
Cipher HomoEncryptCompute::bootstrapCtSfirst(Cipher&ciphertext, uint32_t numIterations, uint32_t precision){
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());
    auto algo = context->GetScheme();
    // 增加错误判断
    if (!cryptoParams) {
        OPENFHE_THROW("CryptoParametersCKKSRNS is not set for the ciphertext.");
    }
    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
#if NATIVEINT == 128
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
#endif
    if (numIterations != 1 && numIterations != 2) {
        OPENFHE_THROW("CKKS Iterative Bootstrapping is only supported for 1 or 2 iterations.");
    }

    uint32_t M               = context->GetCyclotomicOrder();
    uint32_t L0              = cryptoParams->GetElementParams()->GetParams().size();
    auto initSizeQ           = ciphertext->GetElements()[0].GetNumOfElements();
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    if (numIterations > 1) {
        // Step 1: Get the input.
        uint32_t powerOfTwoModulus = 1 << precision;

        // Step 2: Scale up by powerOfTwoModulus, and extend the modulus to powerOfTwoModulus * q.
        // Note that we extend the modulus implicitly without any code calls because the value always stays 0.
        Ciphertext<DCRTPoly> ctScaledUp = ciphertext->Clone();
        // We multiply by powerOfTwoModulus, and leave the last CRT value to be 0 (mod powerOfTwoModulus).
        context->GetScheme()->MultByIntegerInPlace(ctScaledUp, powerOfTwoModulus);
        ctScaledUp->SetLevel(L0 - ctScaledUp->GetElements()[0].GetNumOfElements());

        // Step 3: Bootstrap the initial ciphertext.
        auto ctInitialBootstrap = context->EvalBootstrap(ciphertext, numIterations - 1, precision);
        context->GetScheme()->ModReduceInternalInPlace(ctInitialBootstrap, compositeDegree);

        // Step 4: Scale up by powerOfTwoModulus.
        context->GetScheme()->MultByIntegerInPlace(ctInitialBootstrap, powerOfTwoModulus);

        // Step 5: Mod-down to powerOfTwoModulus * q
        // We mod down, and leave the last CRT value to be 0 because it's divisible by powerOfTwoModulus.
        auto ctBootstrappedScaledDown = ctInitialBootstrap->Clone();
        auto bootstrappingSizeQ       = ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements();

        // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
        if (bootstrappingSizeQ <= initSizeQ) {
            return ciphertext->Clone();
        }

        // TODO: YSP Can be removed for FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            for (auto& cv : ctBootstrappedScaledDown->GetElements()) {
                cv.DropLastElements(bootstrappingSizeQ - initSizeQ);
            }
            ctBootstrappedScaledDown->SetLevel(L0 - ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements());
        }

        // Step 6 and 7: Calculate the bootstrapping error by subtracting the original ciphertext from the bootstrapped ciphertext. Mod down to q is done implicitly.
        auto ctBootstrappingError = context->EvalSub(ctBootstrappedScaledDown, ctScaledUp);

        // Step 8: Bootstrap the error.
        auto ctBootstrappedError = context->EvalBootstrap(ctBootstrappingError, 1, 0);
        context->GetScheme()->ModReduceInternalInPlace(ctBootstrappedError, compositeDegree);

        // Step 9: Subtract the bootstrapped error from the initial bootstrap to get even lower error.
        auto finalCiphertext = context->EvalSub(ctInitialBootstrap, ctBootstrappedError);

        // Step 10: Scale back down by powerOfTwoModulus to get the original message.
        context->EvalMultInPlace(finalCiphertext, static_cast<double>(1) / powerOfTwoModulus);
        return finalCiphertext;
    }

    uint32_t slots = ciphertext->GetSlots();
    // 通过shadow赋值--注意这里的指针
    std::shared_ptr<FHECKKSRNS> algo_bridge = getFHEAlgorithm(context);
    auto shadow = reinterpret_cast<FHECKKSRNSShadow*>(algo_bridge.get());
    uint32_t &m_correctionFactor=shadow->m_correctionFactor;  
    std::cout << "m_correctionFactor: " << m_correctionFactor << std::endl;
    std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
    
    auto pair = m_bootPrecomMap.find(slots);
    if (pair == m_bootPrecomMap.end()) {
        std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                             std::string(" slots were not generated") +
                             std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
        OPENFHE_THROW(errorMsg);
    }
    const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;
    std::cout << "Shadow Find! " << std::endl;

    size_t N = context->GetRingDimension();
    // 这个q是q0，即初始模数，需要大于等于中间模数（缩放因子）；我们说的gap也是q0比缩放因子大的位数
    double qDouble = GetBigModulus(cryptoParams);
    // p是明文模数，在CKKS之中其实就是缩放因子也就是中间的模数
    const auto p = cryptoParams->GetPlaintextModulus();
    double powP  = pow(2, p);
    std::cout << "p: " << p << " powP: " << powP << std::endl;
    // deg其实就是gap
    int32_t deg = std::round(std::log2(qDouble / powP));
    std::cout << "degree: " << deg << std::endl;
#if NATIVEINT != 128
    if (deg > static_cast<int32_t>(m_correctionFactor) &&
        (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
         cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL)) {
        OPENFHE_THROW("Degree [" + std::to_string(deg) + "] must be less than or equal to the correction factor [" +
                      std::to_string(m_correctionFactor) + "].");
    }
#endif
    uint32_t correction = m_correctionFactor - deg;
    std::cout << "correction: " << correction << std::endl;
    double post         = std::pow(2, static_cast<double>(deg));
    
    // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
    double pre      = (compositeDegree > 1) ? cryptoParams->GetScalingFactorReal(0) / qDouble : 1. / post;
    uint64_t scalar = std::llround(post);
    
    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------
    
    // 先输出之前的模数
    auto paramsQ = context->GetElementParams()->GetParams();
    std::cout << "\nModuli in Q:" << std::endl;
    for (uint32_t i = 0; i < paramsQ.size(); i++) {
      std::cout << "q" << i << ": " << log2(paramsQ[i]->GetModulus().ConvertToDouble()) << std::endl;
    }
    std::cout << "Before BTS level:" << ciphertext->GetLevel() << " Encoding scale:" << log2(ciphertext->GetScalingFactor()) << std::endl;
    std::cout << "q: " << ciphertext->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << ciphertext->GetNoiseScaleDeg() << std::endl;

    // if(ciphertext->GetNoiseScaleDeg() == 2){
    //     algo->ModReduceInternalInPlace(ciphertext, 1);
    //     std::cout << "q: " << ciphertext->GetElements()[0].GetModulus().GetMSB() << std::endl;
    //     std::cout << "scale_deg: " << ciphertext->GetNoiseScaleDeg() << std::endl;
    // }
    
    // Mod Raise需要有一层的额外模数，也就是此时不能在最底层模数q0上--中间需要消耗一层模数，乘以correction
    // 乘以correction是为了缩放因子可以很接近甚至等于q0,虽然我没弄明白为什么
    // Mod Raise里面的correction是m_correctionFactor - deg，要大于等于0，这样modraise里面的除以1<<correction缩小密文
    // 这可能是就是为什么可以等于q0的原因（也还是不确定），不过整个BTS最后的corFactor是乘以1<<correction，相互抵消
    std::cout << "Number of levels remaining before modraise: " << 30 - ciphertext->GetLevel() << std::endl;
    std::cout << "Start Mod Raise..." << std::endl;
    auto raised = modraise(ciphertext, correction, compositeDegree);
    std::cout << "Mod Raise Finish!" << std::endl;

    paramsQ = context->GetElementParams()->GetParams();
    std::cout << "After Mod Raise level:" << raised->GetLevel() << " Encoding scale:" << log2(raised->GetScalingFactor()) << std::endl;
    std::cout << "q: " << raised->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << raised->GetNoiseScaleDeg() << std::endl;

    //------------------------------------------------------------------------------
    // SETTING PARAMETERS FOR APPROXIMATE MODULAR REDUCTION
    //------------------------------------------------------------------------------
    std::cout << "Start Approximate Modular Parameters Set..." << std::endl;
    // Coefficients of the Chebyshev series interpolating 1/(2 Pi) Sin(2 Pi K x)
    std::vector<double> coefficients;
    double k = 0;

    if (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) {
        coefficients = shadow->g_coefficientsSparse;
        // k = K_SPARSE;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else {
        // For larger composite degrees, larger K needs to be used to achieve a reasonable probability of failure
        if ((compositeDegree == 1) || ((compositeDegree == 2) && (N < (1 << 17)))) {
            coefficients = shadow->g_coefficientsUniform;
            k            = shadow->K_UNIFORM;
        }
        else {
            coefficients = shadow->g_coefficientsUniformExt;
            k            = shadow->K_UNIFORMEXT;
        }
    }
    std::cout << "Start Approximate Modular Parameters Set Finish!" << std::endl;
    // 缩放因子在这里了！
    std::cout << "pre: " << pre << " k: " << k << " N: " << N << std::endl;
    double constantEvalMult = pre * (1.0 / (k * N));

    context->EvalMultInPlace(raised, constantEvalMult);

    // no linear transformations are needed for Chebyshev series as the range has been normalized to [-1,1]
    double coeffLowerBound = -1;
    double coeffUpperBound = 1;

    Cipher ctxtDec;
    Cipher ct_test;

    // 判断是否使用最基础的线性变换--只消耗一层，在小有效槽数下很有用
    bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                         (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);
    
    if (slots == M / 4){
        //------------------------------------------------------------------------------
        // FULLY PACKED CASE
        //------------------------------------------------------------------------------
        std::cout << "Full Packed Case!" << std::endl;
        //------------------------------------------------------------------------------
        // Running CoeffToSlot
        //------------------------------------------------------------------------------
        std::cout << "Start CoeffToSlot..." << std::endl;
        // need to call internal modular reduction so it also works for FLEXIBLEAUTO
        algo->ModReduceInternalInPlace(raised, compositeDegree);
        
        // only one linear transform is needed as the other one can be derived
        auto ctxtEnc = (isLTBootstrap) ? LinearTransformCtS(raised) : CoeffsToSlots(raised);
        std::cout << "CoeffToSlot Finish!" << std::endl;
        std::cout << "Start Process Conjugate..." << std::endl;
        
        auto conj       = conjugate(ctxtEnc);
        auto ctxtEncI   = context->EvalSub(ctxtEnc, conj);
        context->EvalAddInPlace(ctxtEnc, conj);
        algo->MultByMonomialInPlace(ctxtEncI, 3 * M / 4); // 相当于除j

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                context->ModReduceInPlace(ctxtEnc);
                context->ModReduceInPlace(ctxtEncI);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);//回到单层
                algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
            }
        }
        std::cout << "Process Conjugate Finish!" << std::endl;

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------
        std::cout << "Start Approximate Mod Reduction..." << std::endl;
        // Evaluate Chebyshev series for the sine wave
        ctxtEnc  = context->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);
        ctxtEncI = context->EvalChebyshevSeries(ctxtEncI, coefficients, coeffLowerBound, coeffUpperBound);

        // Double-angle iterations
        if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
            (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
            if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);//回到单层
                algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
            }
            uint32_t numIter;
            if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                numIter = shadow->R_UNIFORM;
            else
                numIter = shadow->R_SPARSE;
            ApplyDoubleAngleIterations(ctxtEnc, numIter);
            ApplyDoubleAngleIterations(ctxtEncI, numIter);
        }

        algo->MultByMonomialInPlace(ctxtEncI, M / 4);
        context->EvalAddInPlace(ctxtEnc, ctxtEncI);

        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            // 与pre对消
            std::cout << "scalar: " << scalar << std::endl;
            algo->MultByIntegerInPlace(ctxtEnc, scalar);
        }
        std::cout << "Approximate Mod Reduction Finish!" << std::endl;

        //------------------------------------------------------------------------------
        // Running SlotToCoeff
        //------------------------------------------------------------------------------
        std::cout << "Start SlotToCoeff..." << std::endl;
        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        std::cout << "Number of levels remaining: " << 30 - ctxtEnc->GetLevel() << std::endl;
        std::cout << "scale_deg: " << ctxtEnc->GetNoiseScaleDeg() << std::endl;
        if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
        }
        std::cout << "Number of levels remaining: " << 30 - ctxtEnc->GetLevel() << std::endl;
        std::cout << "scale_deg: " << ctxtEnc->GetNoiseScaleDeg() << std::endl;

        // Only one linear transform is needed
        ctxtDec = (isLTBootstrap) ? LinearTransformStC(ctxtEnc) : SlotsToCoeffs(ctxtEnc);
        std::cout << "SlotToCoeff Finish!" << std::endl;
    }
    else{
        //------------------------------------------------------------------------------
        // SPARSELY PACKED CASE
        //------------------------------------------------------------------------------
        std::cout << "Sparse Packed Case!" << std::endl;
        //------------------------------------------------------------------------------
        // Running PartialSum
        //------------------------------------------------------------------------------
        std::cout << "Start PartialSum..." << std::endl;
        partialsum(raised, N, slots);
        std::cout << "PartialSum Finish!" << std::endl;
        //------------------------------------------------------------------------------
        // Running CoeffsToSlots
        //------------------------------------------------------------------------------
        std::cout << "Start CoeffsToSlots..." << std::endl;
        algo->ModReduceInternalInPlace(raised, compositeDegree);

        auto ctxtEnc = (isLTBootstrap) ? LinearTransformCtS(raised) : CoeffsToSlots(raised);
        std::cout << "CoeffsToSlots Finish!" << std::endl;
        std::cout << "Start Process Conjugate..." << std::endl;
        auto conj = conjugate(ctxtEnc);
        context->EvalAddInPlace(ctxtEnc, conj);

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                context->ModReduceInPlace(ctxtEnc);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            }
        }
        std::cout << "Process Conjugate Finish!" << std::endl;

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------
        std::cout << "Start Approximate Mod Reduction..." << std::endl;
        // Evaluate Chebyshev series for the sine wave
        ctxtEnc = context->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);

        // Double-angle iterations
        if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
            (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
            if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            }
            uint32_t numIter;
            if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                numIter = shadow->R_UNIFORM;
            else
                numIter = shadow->R_SPARSE;
            ApplyDoubleAngleIterations(ctxtEnc, numIter);
        }

        // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            algo->MultByIntegerInPlace(ctxtEnc, scalar);
        }
        std::cout << "Approximate Mod Reduction Finish!" << std::endl;

        //------------------------------------------------------------------------------
        // Running SlotsToCoeffs
        //------------------------------------------------------------------------------
        std::cout << "Start SlotsToCoeffs..." << std::endl;
        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
        }

        // linear transform for decoding
        ctxtDec = (isLTBootstrap) ? LinearTransformStC(ctxtEnc) : SlotsToCoeffs(ctxtEnc);
        std::cout << "SlotsToCoeffs Finish!" << std::endl;
    }

#if NATIVEINT != 128
    // 64-bit only: scale back the message to its original scale.
    uint64_t corFactor = static_cast<uint64_t>(1) << std::llround(correction);
    std::cout << "corFactor: " << corFactor << std::endl;
    algo->MultByIntegerInPlace(ctxtDec, corFactor);
#endif
    std::cout << "Number of levels remaining: " << 30 - ctxtDec->GetLevel() << std::endl;
    std::cout << "scale_deg: " << ctxtDec->GetNoiseScaleDeg() << std::endl;

    auto bootstrappingNumTowers = ctxtDec->GetElements()[0].GetNumOfElements();

    // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
    if (bootstrappingNumTowers <= initSizeQ) {
        return ciphertext->Clone();
    }

    // return ct_test;
    return ctxtDec;
}
// StC first BTS
Cipher HomoEncryptCompute::bootstrapStCfirst(Cipher&ciphertext, uint32_t numIterations, uint32_t precision){
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());
    auto algo = context->GetScheme();
    // 增加错误判断
    if (!cryptoParams) {
        OPENFHE_THROW("CryptoParametersCKKSRNS is not set for the ciphertext.");
    }
    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
#if NATIVEINT == 128
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
#endif
    if (numIterations != 1 && numIterations != 2) {
        OPENFHE_THROW("CKKS Iterative Bootstrapping is only supported for 1 or 2 iterations.");
    }

    uint32_t M               = context->GetCyclotomicOrder();
    uint32_t L0              = cryptoParams->GetElementParams()->GetParams().size();
    auto initSizeQ           = ciphertext->GetElements()[0].GetNumOfElements();
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    if (numIterations > 1) {
        // Step 1: Get the input.
        uint32_t powerOfTwoModulus = 1 << precision;

        // Step 2: Scale up by powerOfTwoModulus, and extend the modulus to powerOfTwoModulus * q.
        // Note that we extend the modulus implicitly without any code calls because the value always stays 0.
        Ciphertext<DCRTPoly> ctScaledUp = ciphertext->Clone();
        // We multiply by powerOfTwoModulus, and leave the last CRT value to be 0 (mod powerOfTwoModulus).
        context->GetScheme()->MultByIntegerInPlace(ctScaledUp, powerOfTwoModulus);
        ctScaledUp->SetLevel(L0 - ctScaledUp->GetElements()[0].GetNumOfElements());

        // Step 3: Bootstrap the initial ciphertext.
        auto ctInitialBootstrap = context->EvalBootstrap(ciphertext, numIterations - 1, precision);
        context->GetScheme()->ModReduceInternalInPlace(ctInitialBootstrap, compositeDegree);

        // Step 4: Scale up by powerOfTwoModulus.
        context->GetScheme()->MultByIntegerInPlace(ctInitialBootstrap, powerOfTwoModulus);

        // Step 5: Mod-down to powerOfTwoModulus * q
        // We mod down, and leave the last CRT value to be 0 because it's divisible by powerOfTwoModulus.
        auto ctBootstrappedScaledDown = ctInitialBootstrap->Clone();
        auto bootstrappingSizeQ       = ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements();

        // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
        if (bootstrappingSizeQ <= initSizeQ) {
            return ciphertext->Clone();
        }

        // TODO: YSP Can be removed for FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            for (auto& cv : ctBootstrappedScaledDown->GetElements()) {
                cv.DropLastElements(bootstrappingSizeQ - initSizeQ);
            }
            ctBootstrappedScaledDown->SetLevel(L0 - ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements());
        }

        // Step 6 and 7: Calculate the bootstrapping error by subtracting the original ciphertext from the bootstrapped ciphertext. Mod down to q is done implicitly.
        auto ctBootstrappingError = context->EvalSub(ctBootstrappedScaledDown, ctScaledUp);

        // Step 8: Bootstrap the error.
        auto ctBootstrappedError = context->EvalBootstrap(ctBootstrappingError, 1, 0);
        context->GetScheme()->ModReduceInternalInPlace(ctBootstrappedError, compositeDegree);

        // Step 9: Subtract the bootstrapped error from the initial bootstrap to get even lower error.
        auto finalCiphertext = context->EvalSub(ctInitialBootstrap, ctBootstrappedError);

        // Step 10: Scale back down by powerOfTwoModulus to get the original message.
        context->EvalMultInPlace(finalCiphertext, static_cast<double>(1) / powerOfTwoModulus);
        return finalCiphertext;
    }

    uint32_t slots = ciphertext->GetSlots();
    // 通过shadow赋值--注意这里的指针
    std::shared_ptr<FHECKKSRNS> algo_bridge = getFHEAlgorithm(context);
    auto shadow = reinterpret_cast<FHECKKSRNSShadow*>(algo_bridge.get());
    uint32_t &m_correctionFactor=shadow->m_correctionFactor;  
    std::cout << "m_correctionFactor: " << m_correctionFactor << std::endl;
    std::map<uint32_t, std::shared_ptr<CKKSBootstrapPrecom>> &m_bootPrecomMap=shadow->m_bootPrecomMap;
    
    auto pair = m_bootPrecomMap.find(slots);
    if (pair == m_bootPrecomMap.end()) {
        std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
                             std::string(" slots were not generated") +
                             std::string(" Need to call EvalBootstrapSetup and then EvalBootstrapKeyGen to proceed"));
        OPENFHE_THROW(errorMsg);
    }
    const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;
    std::cout << "Shadow Find! " << std::endl;

    size_t N = context->GetRingDimension();

    double qDouble = GetBigModulus(cryptoParams);

    const auto p = cryptoParams->GetPlaintextModulus();
    double powP  = pow(2, p);

    int32_t deg = std::round(std::log2(qDouble / powP));
#if NATIVEINT != 128
    if (deg > static_cast<int32_t>(m_correctionFactor) &&
        (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
         cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL)) {
        OPENFHE_THROW("Degree [" + std::to_string(deg) + "] must be less than or equal to the correction factor [" +
                      std::to_string(m_correctionFactor) + "].");
    }
#endif
    uint32_t correction = m_correctionFactor - deg;
    double post         = std::pow(2, static_cast<double>(deg));

    // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
    double pre      = (compositeDegree > 1) ? cryptoParams->GetScalingFactorReal(0) / qDouble : 1. / post;
    uint64_t scalar = std::llround(post);

    Cipher ctxtDec;
    auto ctxtAfterStC = ciphertext->Clone();
    Cipher ct_test;

    // 判断是否使用最基础的线性变换--只消耗一层，在小有效槽数下很有用
    bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                        (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);

    //------------------------------------------------------------------------------
    // Running SlotsToCoeffs
    //------------------------------------------------------------------------------
    std::cout << "Start SlotToCoeff..." << std::endl;
    // In the case of FLEXIBLEAUTO, we need one extra tower
    // TODO: See if we can remove the extra level in FLEXIBLEAUTO
    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
        // algo->ModReduceInternalInPlace(raised, compositeDegree);
        if(ctxtAfterStC->GetNoiseScaleDeg() == 2){
            ctxtAfterStC = algo->ModReduceInternal(ctxtAfterStC, compositeDegree);
        }
        else{
            std::cout << "No need to mod reduce before StC!" << std::endl;
        }
    }
    // Only one linear transform is needed
    ctxtAfterStC = (isLTBootstrap) ? LinearTransformStC(ctxtAfterStC) : SlotsToCoeffs(ctxtAfterStC);
    std::cout << "SlotToCoeff Finish!" << std::endl;
    // ct_test = ctxtAfterStC->Clone();
    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------
    std::cout << "Start Mod Raise..." << std::endl;
    auto raised = modraise(ctxtAfterStC, correction, compositeDegree);
    std::cout << "Mod Raise Finish!" << std::endl;
    //------------------------------------------------------------------------------
    // SETTING PARAMETERS FOR APPROXIMATE MODULAR REDUCTION
    //------------------------------------------------------------------------------
    std::cout << "Start Approximate Modular Parameters Set..." << std::endl;
    // Coefficients of the Chebyshev series interpolating 1/(2 Pi) Sin(2 Pi K x)
    std::vector<double> coefficients;
    double k = 0;

    if (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) {
        coefficients = shadow->g_coefficientsSparse;
        // k = K_SPARSE;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else {
        // For larger composite degrees, larger K needs to be used to achieve a reasonable probability of failure
        if ((compositeDegree == 1) || ((compositeDegree == 2) && (N < (1 << 17)))) {
            coefficients = shadow->g_coefficientsUniform;
            k            = shadow->K_UNIFORM;
        }
        else {
            coefficients = shadow->g_coefficientsUniformExt;
            k            = shadow->K_UNIFORMEXT;
        }
    }
    std::cout << "Start Approximate Modular Parameters Set Finish!" << std::endl;
    // 缩放因子在这里了！
    double constantEvalMult = pre * (1.0 / (k * N));

    context->EvalMultInPlace(raised, constantEvalMult);

    // no linear transformations are needed for Chebyshev series as the range has been normalized to [-1,1]
    double coeffLowerBound = -1;
    double coeffUpperBound = 1;

    if (slots == M / 4){
        //------------------------------------------------------------------------------
        // FULLY PACKED CASE
        //------------------------------------------------------------------------------
        std::cout << "Full Packed Case!" << std::endl;
        //------------------------------------------------------------------------------
        // Running CoeffToSlot
        //------------------------------------------------------------------------------
        std::cout << "Start CoeffToSlot..." << std::endl;
        algo->ModReduceInternalInPlace(raised, compositeDegree);

        auto ctxtEnc = (isLTBootstrap) ? LinearTransformCtS(raised) : CoeffsToSlots(raised);
        std::cout << "CoeffsToSlots Finish!" << std::endl;
        std::cout << "Start Process Conjugate..." << std::endl;
        // 应该还是两个密文，相同的处理
        auto conj = conjugate(ctxtEnc);
        auto ctxtEncI   = context->EvalSub(ctxtEnc, conj);
        context->EvalAddInPlace(ctxtEnc, conj);
        algo->MultByMonomialInPlace(ctxtEncI, 3 * M / 4);

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (raised->GetNoiseScaleDeg() > 1) {
                context->ModReduceInPlace(ctxtEnc);
                context->ModReduceInPlace(ctxtEncI);
            }
        }
        else {
            if (raised->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
                algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
            }
        }
        std::cout << "Extract Conjugate Finish!" << std::endl;        

        // ct_test = ctxtEncI->Clone();
        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------
        std::cout << "Start Approximate Mod Reduction..." << std::endl;
        // Evaluate Chebyshev series for the sine wave
        ctxtEnc  = context->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);
        ctxtEncI = context->EvalChebyshevSeries(ctxtEncI, coefficients, coeffLowerBound, coeffUpperBound);
        
        // Double-angle iterations
        if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
            (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
            if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
                algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
            }
            uint32_t numIter;
            if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                numIter = shadow->R_UNIFORM;
            else
                numIter = shadow->R_SPARSE;
            ApplyDoubleAngleIterations(ctxtEnc, numIter);
            ApplyDoubleAngleIterations(ctxtEncI, numIter);
        }
        algo->MultByMonomialInPlace(ctxtEncI, M / 4);
        ctxtDec = context->EvalAdd(ctxtEnc, ctxtEncI);
        
        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            algo->MultByIntegerInPlace(ctxtDec, scalar);
            std::cout << "scalar: " << scalar << std::endl;
            // ct_test = ctxtDec->Clone();
        }
        std::cout << "Approximate Mod Reduction Finish!" << std::endl;

    }
    else{
        //------------------------------------------------------------------------------
        // SPARSELY PACKED CASE
        //------------------------------------------------------------------------------
        std::cout << "Sparse Packed Case!" << std::endl;

        //------------------------------------------------------------------------------
        // Running PartialSum
        //------------------------------------------------------------------------------
        std::cout << "Start PartialSum..." << std::endl;
        partialsum(raised, N, slots);
        std::cout << "PartialSum Finish!" << std::endl;
        //------------------------------------------------------------------------------
        // Running CoeffsToSlots
        //------------------------------------------------------------------------------
        std::cout << "Start CoeffsToSlots..." << std::endl;
        algo->ModReduceInternalInPlace(raised, compositeDegree);

        auto ctxtEnc = (isLTBootstrap) ? LinearTransformCtS(raised) : CoeffsToSlots(raised);
        std::cout << "CoeffsToSlots Finish!" << std::endl;
        std::cout << "Start Process Conjugate..." << std::endl;
        auto conj = conjugate(ctxtEnc);
        context->EvalAddInPlace(ctxtEnc, conj);

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                context->ModReduceInPlace(ctxtEnc);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            }
        }
        std::cout << "Process Conjugate Finish!" << std::endl;

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------
        std::cout << "Start Approximate Mod Reduction..." << std::endl;
        // Evaluate Chebyshev series for the sine wave
        ctxtEnc = context->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);

        // Double-angle iterations
        if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
            (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
            if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            }
            uint32_t numIter;
            if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                numIter = shadow->R_UNIFORM;
            else
                numIter = shadow->R_SPARSE;
            ApplyDoubleAngleIterations(ctxtEnc, numIter);
        }
        ctxtDec = ctxtEnc->Clone();
        // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (cryptoParams->GetScalingTechnique() != COMPOSITESCALINGAUTO &&
            cryptoParams->GetScalingTechnique() != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            algo->MultByIntegerInPlace(ctxtDec, scalar);
        }
        std::cout << "Approximate Mod Reduction Finish!" << std::endl;
    }

#if NATIVEINT != 128
    // 64-bit only: scale back the message to its original scale.
    uint64_t corFactor = static_cast<uint64_t>(1) << std::llround(correction);
    algo->MultByIntegerInPlace(ctxtDec, corFactor);
    // ct_test = ctxtDec->Clone();
#endif
    std::cout << "corFactor: " << corFactor << std::endl;

    auto bootstrappingNumTowers = ctxtDec->GetElements()[0].GetNumOfElements();

    // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
    if (bootstrappingNumTowers <= initSizeQ) {
        return ciphertext->Clone();
    }

    // return ct_test;
    return ctxtDec;
}


// ========== Function Approximations ==========

Cipher HomoEncryptCompute::chebyshev(std::function<double(double)> func,
                     const Cipher& ciphertext, double a,
                     double b, uint32_t degree){
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }
    return context->EvalChebyshevFunction(func, ciphertext, a, b, degree);
}

Cipher HomoEncryptCompute::chebyshevseries(const Cipher& ciphertext, const std::vector<double>& coeffcient,
                     double a, double b){
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }
    return context->EvalChebyshevSeries(ciphertext, coeffcient, a, b);                    
}

Cipher HomoEncryptCompute::sigmoid(const Cipher &in, int n, int degree, int scaling) {
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }
    return context->EvalChebyshevFunction([scaling, n](double x) -> double {
        return 1/(n + n * pow(2.71828182846, -scaling*x));

    }, in, -1, 1, degree);
}

Cipher HomoEncryptCompute::sigmoid_tight(const Cipher &in, int n, int degree, int scaling) {
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }
    auto func = [scaling, n](double x) -> double {
        return 1.0 - 1/(n + n * pow(2.71828182846, -scaling*(x/90.0)));
    };
    double a = -0.15 * 90;
    double b = 1.05 * 90;
    return context->EvalChebyshevFunction(func, in, a, b, degree);
}

Cipher HomoEncryptCompute::sigmoid_tanh(const Cipher &in, int degree, int scaling, int scaling_tanh) {
    if (!context) {
        throw std::runtime_error("Context not initialized");
    }

    // 定义sigmoid函数
    auto sigmoid = [scaling](double x) -> double {
        return 1.0 / (1.0 + std::exp(-scaling * x));
    };

    // 定义tanh函数
    auto tanh_func = [scaling_tanh](double x) -> double {
        double exp_pos = std::exp(scaling_tanh * x);
        double exp_neg = std::exp(-scaling_tanh * x);
        return (exp_pos - exp_neg) / (exp_pos + exp_neg);
    };

    // 复合函数：sigmoid(tanh(x))
    return context->EvalChebyshevFunction([sigmoid, tanh_func](double x) -> double {
        return sigmoid(tanh_func(x));
    }, in, -1.0, 1.0, degree);
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

// ========== Plaintext Creation Utilities ==========

Plain HomoEncryptCompute::MakePlaintext(double val, uint32_t level) {
    std::vector<std::complex<double>> x1;
    for (size_t i = 0; i < context->GetRingDimension()/2 ;i++) {
        x1.push_back(std::complex<double>(val, 0.0));
    }
    return context->MakeCKKSPackedPlaintext(x1, 1,level);
}

Plain HomoEncryptCompute::MakePlaintext(std::vector<std::complex<double>> v, uint32_t level) {
    std::vector<std::complex<double>> x1;
    uint32_t numSlots=v.size();
    for (size_t i = 0; i < context->GetRingDimension()/2/numSlots ;i++) {
        for(size_t j = 0; j < numSlots; j++) {
            x1.push_back(v[j]);
        }
    }
    return context->MakeCKKSPackedPlaintext(x1, 1,level);
}

Plain HomoEncryptCompute::MakePlaintext(std::vector<double> v, uint32_t level) {
    std::vector<std::complex<double>> x1;
    uint32_t numSlots=v.size();
    for (size_t i = 0; i < context->GetRingDimension()/2/numSlots ;i++) {
        for(size_t j = 0; j < numSlots; j++) {
            x1.push_back(std::complex<double>(v[j],0));
        }
    }
    return context->MakeCKKSPackedPlaintext(x1, 1,level);
}

// ========== Decrypt Utilities ==========

std::vector<std::complex<double>> HomoEncryptCompute::DecryptCKKSPackedValue(const Cipher& ciphertext, uint32_t numSlots) {
    return DecryptCKKSPackedValue(ciphertext, key_pair.secretKey, numSlots);
}

std::vector<std::complex<double>> HomoEncryptCompute::DecryptCKKSPackedValue(
    const Cipher& ciphertext, 
    const PrivateKey<DCRTPoly>& privateKey, 
    uint32_t numSlots) {

    Plaintext result;
    ciphertext->GetCryptoContext()->Decrypt(privateKey, ciphertext, &result);
    result->SetLength(numSlots);
    return result->GetCKKSPackedValue();
}

// ========== Utility Functions ==========

std::vector<double> HomoEncryptCompute::genUniformReal(uint32_t n) {
    std::vector<double> vec(n);
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, 1.0);
    for (size_t i = 0; i < n; i++) {
        vec[i] = dis(gen);
    }
    return vec;
}

double HomoEncryptCompute::estimatePrecision(std::vector<std::complex<double>>& v1, std::vector<std::complex<double>>& v2) {
    double precVal = 0.;
    uint32_t n = v1.size();
    uint32_t validCount = 0;
    
    for(size_t i = 0; i < n; i++) {
        double denom = abs(v2[i].real());
        double numer = abs(v1[i].real() - v2[i].real());
        
        // 跳过分母为零或接近零的情况
        if (denom < 1e-10) {
            std::cout << "Warning: v2[" << i << "].real() is too small: " << denom << std::endl;
            continue;
        }
        
        // 跳过分子为零或接近零的情况（完美匹配）
        if (numer < 1e-15) {
            // 如果误差极小，认为是非常高的精度（比如50位）
            precVal += 50.0;
            validCount++;
            continue;
        }
        
        double ratio = numer / denom;
        
        // 确保 ratio 在合理范围内
        if (ratio >= 1.0) {
            std::cout << "Warning: relative error >= 1 at index " << i << std::endl;
            // 如果误差大于原值，精度为负或0
            precVal += 0.0;
        } else {
            double prec = -std::log2(ratio);
            precVal += prec;
        }
        validCount++;
    }
    
    // 如果没有有效的数据点，返回0
    if (validCount == 0) {
        std::cout << "Error: No valid data points for precision estimation" << std::endl;
        return 0.0;
    }
    
    return precVal / validCount;
}

// ========== Print Functions ==========

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

void HomoEncryptCompute::print_moduli_chain_detail(const DCRTPoly& poly){
    int num_primes = poly.GetNumOfElements();
    double total_bit_len = 0.0;
    for (int i = 0; i < num_primes; i++) {
        auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
        std::cout << "q_" << i << ": " 
                    << qi
                    << ",  log q_" << i <<": " << log(qi.ConvertToDouble()) / log(2)
                    << std::endl;
        total_bit_len += log(qi.ConvertToDouble()) / log(2);
    }   
    std::cout << "Total bit length: " << total_bit_len << std::endl;
}

void HomoEncryptCompute::printIntegerMod(const BigInteger &x,const BigInteger &q) {
    if(x>q/2)
        // std::cout << "-" << (q-x) << " " << q;
        std::cout << "-" << (q-x);
    else
        // std::cout << x << " " << q;
        std::cout << x ;
}

void HomoEncryptCompute::printBigVectorMod(const BigVector &x) {
    std::cout << "[";
    for (size_t i = 0; i < x.GetLength(); i++) {
        printIntegerMod(x[i],x.GetModulus());
        if (i<x.GetLength()-1) std::cout << " ";
    }
    std::cout << "]";
}

void HomoEncryptCompute::printPoly(const Poly &p) {
    std::cout << "[";
    for (size_t i = 0; i < p.GetLength(); i++) {
        printIntegerMod(p[i],p.GetModulus());
        if (i<p.GetLength()-1) std::cout << " ";
    }
    std::cout << "]";
}

// ========== Test Functions ==========

void HomoEncryptCompute::testSlots2Coeffs() {
    std::cout << "Test SlotsToCoeffs" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 4; // This is the number of slots in the plaintext
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 4);

    // parameters.SetExecutionMode(EXEC_NOISE_ESTIMATION); // This is to allow computation with complex numbers
    parameters.SetScalingTechnique(FIXEDAUTO);
    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    uint32_t levelBudget = 2;
    
    usint depth = 3;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();
    // This is the maximum number of slots that can be used for full packing.

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=2; // number of remaining levels after SlotsToCoeffs
    EvalSlotsToCoeffsSetup(levelBudget, numSlots,lDec);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1),
        std::complex<double>(2., 0.2),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.4)
    };

    std::cout << "Input value: " << vec << std::endl;

    // To test SlotsToCoeffs, we must encode with twice the number of slots.
    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lDec-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph->SetSlots(numSlots);
   
    ciph = SlotsToCoeffs(ciph);

    std::cout << "After S2C: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    Poly xx1=myDecrypt(ciph, key_pair.secretKey);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ciph->GetSlots() << std::endl;
    std::cout << "modulus size: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;
    for(size_t i=0; i<xx1.GetLength(); i++) {
        std::cout << i << " ";
        std::cout << " ";
        q = xx1.GetModulus();
         if(xx1[i]>q/2)
            x = q - xx1[i];
        else
            x= xx1[i];
        std::cout << bigIntegerToDouble(x)/ciph->GetScalingFactor();

        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph, &result);
    result->SetLength(numSlots);
    std::cout << "Output after StC \n\t" << result << std::endl;
}

void HomoEncryptCompute::testStC_CtS() {
    std::cout << "Test SlotsToCoeffs & CoeffsToSlots" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 4; // This is the number of slots in the plaintext
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 5);

    parameters.SetScalingTechnique(FIXEDAUTO);
    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    std::vector<uint32_t> levelBudget = {2, 2};
    
    usint depth = 6;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=4; // number of remaining levels after StC
    uint32_t lEnc=3; // number of remaining levels after CtS
    // EvalStC_CtSSetup(levelBudget, numSlots, lDec, lEnc);
    EvalStC_CtSSetup(levelBudget, numSlots * 2, lDec, lEnc);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots * 2);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1),
        std::complex<double>(2., 0.2),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.4),
        std::complex<double>(1., 0.1),
        std::complex<double>(2., 0.2),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.4)
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lDec-1,nullptr,numSlots * 2);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    // ciph->SetSlots(numSlots);
   
    ciph = SlotsToCoeffs(ciph);
    auto ciph2 = ciph->Clone();

    std::cout << "After S2C: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    Poly xx1=myDecrypt(ciph, key_pair.secretKey);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ciph->GetSlots() << std::endl;
    std::cout << "modulus size: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;
    for(size_t i=0; i<xx1.GetLength(); i++) {
        std::cout << i << " ";
        std::cout << " ";
        q = xx1.GetModulus();
         if(xx1[i]>q/2)
            x = q - xx1[i];
        else
            x= xx1[i];
        std::cout << bigIntegerToDouble(x)/ciph->GetScalingFactor();

        std::cout << std::endl;
    }
    std::cout << std::endl;

    ciph2 = CoeffsToSlots(ciph2);
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph2, &result);
    result->SetLength(numSlots);
    std::cout << "Output after StC and CtS \n\t" << result << std::endl;
}

void HomoEncryptCompute::testCtS_StC() {
    std::cout << "Test CoeffsToSlots & SlotsToCoeffs" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 4;
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 4);

    parameters.SetScalingTechnique(FIXEDAUTO);
    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    std::vector<uint32_t> levelBudget = {2, 2};
    
    usint depth = 6;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=3;
    uint32_t lEnc=4;
    EvalStC_CtSSetup(levelBudget, numSlots, lDec, lEnc);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1),
        std::complex<double>(2., 0.2),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.4)
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lEnc-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph->SetSlots(numSlots);
   
    ciph = CoeffsToSlots(ciph);
    std::cout << "After CtS: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph = SlotsToCoeffs(ciph);
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph, &result);
    result->SetLength(numSlots);
    std::cout << "Output after CtS and StC \n\t" << result << std::endl;
}

void HomoEncryptCompute::mytestSlots2Coeffs() {
    std::cout << "Test SlotsToCoeffs" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t ring_dim = 1 << 4;
    uint32_t batchSize = ring_dim / 2;
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 4);

    parameters.SetScalingTechnique(FIXEDAUTO);
    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    uint32_t levelBudget = 2;
    
    usint depth = 3;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << "CKKS scheme is using slots " << numSlots << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=2;
    EvalSlotsToCoeffsSetup(levelBudget, numSlots,lDec);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1), std::complex<double>(5., 0.5),
        std::complex<double>(2., 0.2), std::complex<double>(6., 0.6),
        std::complex<double>(3., 0.3), std::complex<double>(7., 0.7),
        std::complex<double>(4., 0.4), std::complex<double>(8., 0.8)
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lDec-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    std::cout << "Start SlotsToCoeffs" << std::endl;
    ciph = SlotsToCoeffs(ciph);

    std::cout << "After S2C: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    Poly xx1=myDecrypt(ciph, key_pair.secretKey);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ciph->GetSlots() << std::endl;
    std::cout << "modulus size: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;
    for(size_t i=0; i<xx1.GetLength(); i++) {
        std::cout << i << " ";
        std::cout << " ";
        q = xx1.GetModulus();
         if(xx1[i]>q/2)
            x = q - xx1[i];
        else
            x= xx1[i];
        std::cout << bigIntegerToDouble(x)/ciph->GetScalingFactor();

        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void HomoEncryptCompute::mytestStC_CtS() {
    std::cout << "Test SlotsToCoeffs & CoeffsToSlots" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 32;
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 6);

    int dcrtBits = 45;
    int firstMod = 60; // This is the first modulus size in bits, it
    parameters.SetScalingTechnique(FIXEDAUTO);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(firstMod);

    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    std::vector<uint32_t> levelBudget = {2, 2};
    
    usint depth = 6;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=4;
    uint32_t lEnc=3;
    EvalStC_CtSSetup(levelBudget, numSlots, lDec, lEnc);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1),
        std::complex<double>(7., 0.7),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.3)
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lDec-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;
    std::cout << "Initial q: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;

    ciph->SetSlots(numSlots);
   
    ciph = SlotsToCoeffs(ciph);
    auto ciph2 = ciph->Clone();

    std::cout << "After StC: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    Poly xx1=myDecrypt(ciph, key_pair.secretKey);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ciph->GetSlots() << std::endl;
    std::cout << "modulus size: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;
    for(size_t i=0; i<xx1.GetLength(); i++) {
        std::cout << i << " ";
        std::cout << " ";
        q = xx1.GetModulus();
         if(xx1[i]>q/2)
            x = q - xx1[i];
        else
            x= xx1[i];
        // std::cout << "q: " << log2(q.ConvertToDouble()) << " " << "q: " << q << " ";
        // std::cout << "SaclingFactor: " << log2(ciph->GetScalingFactor()) << " " << "x: " << ciph->GetScalingFactor() << " ";
        std::cout << bigIntegerToDouble(x)/ciph->GetScalingFactor();

        std::cout << std::endl;
    }
    std::cout << std::endl;

    ciph2 = CoeffsToSlots(ciph2);
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph2, &result);
    result->SetLength(numSlots);
    std::cout << "Output after StC and CtS \n\t" << result << std::endl;
}

void HomoEncryptCompute::mytestCtS_StC() {
    std::cout << "Test CoeffsToSlots & SlotsToCoeffs" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 16;
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 5);

    parameters.SetScalingTechnique(FIXEDAUTO);
    CKKSDataType ckksDataType = COMPLEX;
    parameters.SetCKKSDataType(ckksDataType);

    std::vector<uint32_t> levelBudget = {2, 2};
    
    usint depth = 6;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=3;
    uint32_t lEnc=4;
    EvalStC_CtSSetup(levelBudget, numSlots, lDec, lEnc);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<std::complex<double>> vec={
        std::complex<double>(1., 0.1),
        std::complex<double>(7., 0.7),
        std::complex<double>(3., 0.3),
        std::complex<double>(4., 0.3)
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lEnc-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph->SetSlots(numSlots);
   
    ciph = CoeffsToSlots(ciph);
    std::cout << "After CtS: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph = SlotsToCoeffs(ciph);
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph, &result);
    result->SetLength(numSlots);
    std::cout << "Output after CtS and StC \n\t" << result << std::endl;
}

void HomoEncryptCompute::testStC_CtSReal(){
    std::cout << "Test SlotsToCoeffs & CoeffsToSlots" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    uint32_t batchSize = 4;
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 4);

    parameters.SetScalingTechnique(FIXEDAUTO);

    std::vector<uint32_t> levelBudget = {2, 2};
    
    usint depth = 6;
    parameters.SetMultiplicativeDepth(depth);

    std::cout << "Depth: " << depth << std::endl;

    context = GenCryptoContext(parameters);

    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    usint ringDim = context->GetRingDimension();

    usint numSlots = batchSize;
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;

    std::vector<uint32_t> dim1 = {0, 0};
    uint32_t lDec=4;
    uint32_t lEnc=3;
    EvalStC_CtSSetup(levelBudget, numSlots, lDec, lEnc);

    key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);
    context->EvalBootstrapKeyGen(key_pair.secretKey, numSlots);

    std::vector<double> vec={
        1.,
        2.,
        3.,
        4.
    };

    std::cout << "Input value: " << vec << std::endl;

    Plaintext ptxt1 = context->MakeCKKSPackedPlaintext(vec, 1,depth-lDec-1,nullptr,numSlots);
    Ciphertext<DCRTPoly> ciph = context->Encrypt(key_pair.publicKey, ptxt1);

    std::cout << "After encrypt: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    ciph->SetSlots(numSlots);
   
    ciph = SlotsToCoeffs(ciph);
    std::cout << "After S2C: level:" << ciph->GetLevel() << " Encoding scale:" << log2(ciph->GetScalingFactor()) << std::endl;

    Poly xx1=myDecrypt(ciph, key_pair.secretKey);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ciph->GetSlots() << std::endl;
    std::cout << "modulus size: " << ciph->GetElements()[0].GetModulus().GetMSB() << std::endl;
    for(size_t i=0; i<xx1.GetLength(); i++) {
        std::cout << i << " ";
        std::cout << " ";
        q = xx1.GetModulus();
         if(xx1[i]>q/2)
            x = q - xx1[i];
        else
            x= xx1[i];
        std::cout << bigIntegerToDouble(x)/ciph->GetScalingFactor();

        std::cout << std::endl;
    }
    std::cout << std::endl;

    ciph = CoeffsToSlots(ciph);
    Plaintext result;
    context->Decrypt(key_pair.secretKey, ciph, &result);
    result->SetLength(numSlots);
    std::cout << "Output after StC and CtS \n\t" << result << std::endl;
}
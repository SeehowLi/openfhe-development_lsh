/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-03 21:55:24
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-07 16:04:48
 * @FilePath: \openfhe-development\src\pke\examples\discrete_ckks.cpp
 * @Description: 用于实现离散CKKS的示例代码
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-fhe.h"
using namespace lbcrypto;

void StC_CtS_example();
void discrete_ckks_example();

int main() {
    // discrete_ckks_example();
    StC_CtS_example();

    return 0;
}

void StC_CtS_example() {
    std::cout << "--------------------------------- STC CTS EXAMPLE ---------------------------------"
              << std::endl;
    // 参数容器的声明
    CCParams<CryptoContextCKKSRNS> parameters;
    // 设置密钥分布--注意：CKKS的密钥分布可以是SPARSE_TERNARY或UNIFORM_TERNARY。
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    // SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    // 设置安全等级，这里不设置是为了测试
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 12);
    // 设置密钥切换模式--混合切换方法，digit size默认是3，增加会减少复杂度但是增加存储密钥大小
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

    // 设置缩放因子参数
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits               = 59;
    usint firstMod               = 60;
#endif
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    // 设置levelBudget
    std::vector<uint32_t> levelBudget = {2, 2};
    // 设置bsgs的维度，值为0就是自动选
    std::vector<uint32_t> bsgsDim = {0, 0};
    // 设置层级
    uint32_t levelsAvailableAfterBootstrap = 8;
    usint depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);
    std::cout << "BTS depth is: " << FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist) << std::endl;
    std::cout << "Multiplicative Depth is:" << depth << std::endl;

    // 声明密文上下文容器,传递参数
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    // 功能启用
    cc ->Enable(PKE);
    cc ->Enable(KEYSWITCH);
    cc ->Enable(LEVELEDSHE);
    cc ->Enable(ADVANCEDSHE);
    cc ->Enable(FHE);
    cc ->Enable(DISCRETECKKS);

    // 设置槽数
    uint32_t numSlots = 1 << 3;
    std::cout << "Slot Number is:" << numSlots << std::endl;

    // 预计算参数
    cc->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots);
    std::cout << "Bootstrap Setup is done." << std::endl;

    // 生成密钥
    auto keypair = cc->KeyGen();
    cc->EvalMultKeyGen(keypair.secretKey);
    cc->EvalBootstrapKeyGen(keypair.secretKey, numSlots);

    // 显式预计算
    // cc->EvalBootstrapPrecompute(numSlots);

    // 生成加密向量
    std::vector<double> x = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};


    // 生成明文多项式--稀疏打包
    Plaintext pt = cc->MakeCKKSPackedPlaintext(x, 1, 1, nullptr, numSlots);
    pt->SetLength(numSlots);
    std::cout << "Input: " << pt;
    
    // 加密明文
    Ciphertext<DCRTPoly> ct = cc->Encrypt(keypair.publicKey, pt);
    // std::cout << "Ciphertext: " << ct << std::endl;
    std::cout << "Initial number of levels remaining: " << depth - ct->GetLevel() << std::endl;

    // // BTS
    // auto ciphertextAfter = cc->EvalBootstrap(ct);
    // std::cout << "Number of levels remaining after bootstrapping: " << depth - ciphertextAfter->GetLevel() << std::endl
    //           << std::endl;
    // StC
    auto ctAfterStC = cc->EvalStC(ct);
    std::cout << "SlotToCoeff SUCESS!" << std::endl;

    // CtS
    auto ctAfterCtS = cc->EvalCtS(ctAfterStC);
    std::cout << "CoeffToSlot SUCESS!" << std::endl;

    Plaintext result;
    cc->Decrypt(keypair.secretKey, ctAfterCtS, &result);
    result->SetLength(numSlots);
    std::cout << "Output after StC -> CtS \n\t" << result << std::endl;



}

void discrete_ckks_example() {
    std::cout << "--------------------------------- DISCRETE CKKS EXAMPLE ---------------------------------"
              << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 10);
    

}
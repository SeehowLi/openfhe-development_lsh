//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Simple examples for CKKS
 */

#define PROFILE

#include "openfhe.h"

using namespace lbcrypto;

int main() {
    // Step 1: Setup CryptoContext

    // A. 指定主要参数
    /* A1) 乘法深度:
     * 我们在此设置的 CKKS 方案将适用于任何乘法深度等于 'multDepth' 的计算。
     * 这是给定乘法的最大可能深度，但不是方案支持的总乘法次数。
     *
     * 例如，计算 f(x, y) = x^2 + x*y + y^2 + x + y 的乘法深度为 1，
     * 但需要总共 3 次乘法。另一方面，计算 g(x_i) = x1*x2*x3*x4 可以
     * 实现为乘法深度为 3 的计算，例如 g(x_i) = ((x1*x2)*x3)*x4，
     * 或实现为乘法深度为 2 的计算，例如 g(x_i) = (x1*x2)*(x3*x4)。
     *
     * 出于性能原因，通常优先选择以最小乘法深度执行操作。
     */
    uint32_t multDepth = 1;

    /* A2) 缩放因子的位长度。
     * CKKS 适用于实数，但这些数字被编码为整数。
     * 例如，实数 m=0.01 被编码为 m'=round(m*D)，其中 D 是一个称为缩放因子的方案参数。
     * 假设 D=1000，那么 m' 是 10（一个整数）。假设基于 m' 的计算结果是 130，
     * 那么在解密时，缩放因子被移除，因此用户看到的实数结果是 0.13。
     *
     * 参数 'scaleModSize' 决定了缩放因子 D 的位长度，而不是缩放因子本身。
     * 后者是实现特定的，并且在某些版本的 CKKS（例如 FLEXIBLEAUTO）中，
     * 它可能在不同的密文之间有所不同。
     *
     * 选择 'scaleModSize' 取决于计算所需的精度，以及其他参数如 multDepth 或安全标准。
     * 这是因为其他参数决定了计算过程中会引入多少噪声（记住 CKKS 是一个近似方案，
     * 每次操作都会引入少量噪声）。缩放因子应足够大，以同时容纳这些噪声并支持符合所需精度的结果。
     */
    uint32_t scaleModSize = 50;//50位的缩放因子

    /* A3) 明文槽的数量，用于密文中。
     * CKKS 在每个密文中打包多个明文值。
     * 槽的最大数量取决于一个称为环维度的安全参数。
     * 在此示例中，我们不直接指定环维度，而是让库根据我们选择的安全级别、
     * 我们希望支持的乘法深度以及缩放因子大小来选择它。
     *
     * 请使用方法 GetRingDimension() 来找出这些参数所使用的确切环维度。
     * 给定环维度 N，最大批处理大小为 N/2，这是由于 CKKS 的工作方式决定的。
     */
    uint32_t batchSize = 8;

    /* A4) 所需的安全级别基于 FHE 标准。
     * 此参数可以取四个值。三个可能的值分别对应于 128 位、192 位和 256 位安全性，
     * 第四个值对应于 "NotSet"，这意味着用户需要负责选择安全参数。
     * 自然地，"NotSet" 应仅在非生产环境中使用，或者由了解其选择安全性影响的专家使用。
     *
     * 如果选择了给定的安全级别，库将参考当前由 FHE 标准联盟定义的安全参数表
     * (https://homomorphicencryption.org/introduction/) 来自动选择安全参数。
     * 有关更多详细信息，请参阅以下参考文献中的 "推荐参数表"：
     * http://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
     */
    // CCParams<CryptoContextCKKSRNS> parameters是参数容器，包括各种配置参数
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetBatchSize(batchSize);

    // 支持复数运算
    // parameters.SetExecutionMode(EXEC_NOISE_ESTIMATION);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    // Enable the features that you wish to use
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

    // B. Step 2: Key Generation
    /* B1) Generate encryption keys.
   * These are used for encryption/decryption, as well as in generating
   * different kinds of keys.
   */
    auto keys = cc->KeyGen();

    /* B2) 生成数字大小
     * 在 CKKS 中，每当有人使用密钥 s 对两个密文进行乘法时，
     * 我们会得到一个结果，其中一些分量在密钥 s 下有效，
     * 而另一个附加分量在密钥 s^2 下有效。
     *
     * 在大多数情况下，我们希望对乘法结果进行重新线性化，
     * 即我们希望将密文的 s^2 分量转换为在原始密钥 s 下有效。
     * 为此，我们需要通过以下代码行创建一个重新线性化密钥。
     */
    cc->EvalMultKeyGen(keys.secretKey);

    /* B3) 生成旋转密钥
     * CKKS 支持对打包密文的内容进行旋转，但要实现这一点，
     * 我们需要创建一个称为旋转密钥的对象。这可以通过以下调用完成，
     * 它接受一个向量作为输入，其中的索引对应于我们希望支持的旋转偏移量。
     * 负索引对应右移，正索引对应左移。查看此演示的输出以了解示例。
     *
     * 请记住，旋转作用于批处理大小或整个环维度（如果未指定批处理大小）。
     * 这意味着，如果环维度为 8 且未指定批处理大小，
     * 那么输入 (1,2,3,4,0,0,0,0) 旋转 2 次将变为 (3,4,0,0,0,0,1,2)，
     * 而不是 (3,4,1,2,0,0,0,0)。
     * 如果环维度为 8 且批处理大小设置为 4，
     * 那么 (1,2,3,4) 的旋转 2 次将变为 (3,4,1,2)。
     * 此外，正如可以在此演示的输出中观察到的那样，
     * 由于 CKKS 是近似的，零并不是真正的零——它们只是非常小的数字。
     */
    cc->EvalRotateKeyGen(keys.secretKey, {1, -2});

    // Step 3: Encoding and encryption of inputs

    // Inputs--都是实数。。不是虚数-换成虚数试一下
    // std::vector<double> x1 = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    // std::vector<double> x2 = {5.0, 4.0, 3.0, 2.0, 1.0, 0.75, 0.5, 0.25};
    std::vector<std::complex<double>> x1 = {std::complex<double>(0.25, 0.1), std::complex<double>(0.5, 0.2),
                                            std::complex<double>(0.75, 0.3), std::complex<double>(1.0, 0.4),
                                            std::complex<double>(2.0, 0.5), std::complex<double>(3.0, 0.6), 
                                            std::complex<double>(4.0, 0.7), std::complex<double>(5.0, 0.8)};
    std::vector<std::complex<double>> x2 = {std::complex<double>(5.0, 0.1), std::complex<double>(4.0, 0.2),
                                            std::complex<double>(3.0, 0.3), std::complex<double>(2.0, 0.4),
                                            std::complex<double>(1.0, 0.5), std::complex<double>(0.75, 0.6), 
                                            std::complex<double>(0.5, 0.7), std::complex<double>(0.25, 0.8)};

    // Encoding as plaintexts--直接打包
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(x2);

    std::cout << "Input x1: " << ptxt1 << std::endl;
    std::cout << "Input x2: " << ptxt2 << std::endl;

    // Encrypt the encoded vectors
    auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto c2 = cc->Encrypt(keys.publicKey, ptxt2);

    // Step 4: Evaluation

    // Homomorphic addition
    auto cAdd = cc->EvalAdd(c1, c2);

    // Homomorphic subtraction
    auto cSub = cc->EvalSub(c1, c2);

    // Homomorphic scalar multiplication
    auto cScalar = cc->EvalMult(c1, 4.0);

    // Homomorphic multiplication--估计里面自动RS并且RLK了
    auto cMul = cc->EvalMult(c1, c2);

    // Homomorphic rotations
    auto cRot1 = cc->EvalRotate(c1, 1);
    auto cRot2 = cc->EvalRotate(c1, -2);

    // Step 5: Decryption and output
    Plaintext result;
    // We set the cout precision to 8 decimal digits for a nicer output.
    // If you want to see the error/noise introduced by CKKS, bump it up
    // to 15 and it should become visible.
    std::cout.precision(16);

    std::cout << std::endl << "Results of homomorphic computations: " << std::endl;

    cc->Decrypt(keys.secretKey, c1, &result);
    result->SetLength(batchSize);
    std::cout << "x1 = " << result;
    std::cout << "Estimated precision in bits: " << result->GetLogPrecision() << std::endl;

    // Decrypt the result of addition
    cc->Decrypt(keys.secretKey, cAdd, &result);
    result->SetLength(batchSize);
    std::cout << "x1 + x2 = " << result;
    std::cout << "Estimated precision in bits: " << result->GetLogPrecision() << std::endl;

    // Decrypt the result of subtraction
    cc->Decrypt(keys.secretKey, cSub, &result);
    result->SetLength(batchSize);
    std::cout << "x1 - x2 = " << result << std::endl;

    // Decrypt the result of scalar multiplication
    cc->Decrypt(keys.secretKey, cScalar, &result);
    result->SetLength(batchSize);
    std::cout << "4 * x1 = " << result << std::endl;

    // Decrypt the result of multiplication
    cc->Decrypt(keys.secretKey, cMul, &result);
    result->SetLength(batchSize);
    std::cout << "x1 * x2 = " << result << std::endl;

    // Decrypt the result of rotations

    cc->Decrypt(keys.secretKey, cRot1, &result);
    result->SetLength(batchSize);
    std::cout << std::endl << "In rotations, very small outputs (~10^-10 here) correspond to 0's:" << std::endl;
    std::cout << "x1 rotate by 1 = " << result << std::endl;

    cc->Decrypt(keys.secretKey, cRot2, &result);
    result->SetLength(batchSize);
    std::cout << "x1 rotate by -2 = " << result << std::endl;

    return 0;
}

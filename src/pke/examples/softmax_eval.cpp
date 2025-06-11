/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-06-02 23:27:17
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-06-11 09:53:25
 * @FilePath: \openfhe-development\src\pke\examples\softmax_eval.cpp
 * @Description: approximate softmax function
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

#include "openfhe.h"
#include "math/chebyshev.h"


using namespace lbcrypto;

void EvalSoftmaxExample();
void EvalSoftmaxExample_Cho();
std::vector<std::complex<double>> ComputeSoftmax(const std::vector<std::complex<double>>& input); 
void TestInverseSqrtApproximation_Chen();


int main(int argc, char* argv[]) {
    TestInverseSqrtApproximation_Chen();

    // EvalSoftmaxExample();
    
    // std::cout << "--------------------------------- EVAL SOFTMAX FUNCTION BY Cho---------------------------------"
    //           << std::endl;
    // EvalSoftmaxExample_Cho();
    return 0;
}

// In this example, we evaluate the softmax function
void EvalSoftmaxExample() {
    std::cout << "--------------------------------- EVAL SOFTMAX FUNCTION ---------------------------------"
              << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 10);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    usint scalingModSize = 50;
    usint firstModSize   = 60;
#endif
    parameters.SetScalingModSize(scalingModSize);
    parameters.SetFirstModSize(firstModSize);

    // Choosing a higher degree yields better precision, but a longer runtime.
    uint32_t polyDegree = 16;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 16;

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // We need to enable Advanced SHE to use the Chebyshev approximation.
    cc->Enable(ADVANCEDSHE);

    auto keyPair = cc->KeyGen();
    // We need to generate mult keys to run Chebyshev approximations.
    cc->EvalMultKeyGen(keyPair.secretKey);
 
    std::vector<std::complex<double>> input{-10.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 10.0};
    // std::vector<std::complex<double>> input{1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0};
    size_t encodedLength = input.size();
    // We need to generate rotate keys to evalute softmax function.
    std::vector<int32_t> rotationIndices;
    for (uint32_t i = 1; i <= encodedLength; i <<= 1) {
        rotationIndices.push_back(static_cast<int32_t>(i));
    }
    cc->EvalAtIndexKeyGen(keyPair.secretKey, rotationIndices);
    // sum key
    cc->EvalSumKeyGen(keyPair.secretKey);

    // full slots packing
    Plaintext plaintext  = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext      = cc->Encrypt(keyPair.publicKey, plaintext);

    double lowerBound = -11;
    double upperBound = 11;
    auto result       = cc->EvalSoftmaxUsingCheb(ciphertext, lowerBound, upperBound, polyDegree, encodedLength);

    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, result, &plaintextDec);
    plaintextDec->SetLength(encodedLength);

    std::vector<std::complex<double>> expectedOutput({
    0.000212, 0.000576, 0.001567, 0.004259, 0.011581, 0.031472, 0.085563, 0.232554, 0.632216});
    std::cout << "Expected output\n\t" << expectedOutput << std::endl;

    std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;
}

// In this example, we evaluate the softmax function using the method form paper:
// "Fast and Accurate Homomorphic Softmax Evaluation"--without BTS 
void EvalSoftmaxExample_Cho() {
    std::cout << "--------------------------------- EVAL SOFTMAX FUNCTION BY Cho---------------------------------"
              << std::endl;
    std::cout << "---------------------------------         Without BTS         ---------------------------------"
              << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_NotSet);
    //ring dimension is 1 << 10 = 1024
    parameters.SetRingDim(1 << 10);
#if NATIVEINT == 128
    usint scalingModSize = 78;
    usint firstModSize   = 89;
#else
    // ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint scalingModSize = 50;
    usint firstModSize   = 60;
#endif
    parameters.SetScalingModSize(scalingModSize);
    // parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstModSize);
    
    // Set the degree for polynomial approximation.
    uint32_t polyDegree = 16;
    // Set the depth properly to ensure the softmax function can be evaluated.
    uint32_t multDepth = 36;
    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    // cc->Enable(FHE);

    auto keyPair = cc->KeyGen();
    // We need to generate mult keys
    cc->EvalMultKeyGen(keyPair.secretKey);
    // Gen roate keys

    size_t logM_bound = 7;//cant be too small
    size_t m_bound = 1 << logM_bound;//128
    size_t loglen_softmax = 8;//4
    size_t len_softmax = 1 << loglen_softmax;//16
    double ln_len = std::log(static_cast<double>(len_softmax));    // ln(len_softmax)
    double loglnlen_softmax = std::log2(ln_len);  // log₂(ln(len_softmax))
    size_t k = static_cast<size_t>(std::round(logM_bound - loglnlen_softmax));//k≈logM-loglnn
    size_t softmax_num = 1;
    size_t slotnumInUse = softmax_num * len_softmax; // 16 slots in use

    std::cout << "Softmax num is :" << softmax_num << std::endl;
    std::cout << "Softmax dimension is :" << len_softmax << std::endl;
    std::cout << "k is :" << k << std::endl;

    std::vector<int32_t> rotationIndices;
    for (uint32_t i = 1; i <= len_softmax; i <<= 1) {
        rotationIndices.push_back(static_cast<int32_t>(i));
    }
    cc->EvalAtIndexKeyGen(keyPair.secretKey, rotationIndices);
    // sum key -- perhaps no use
    cc->EvalSumKeyGen(keyPair.secretKey);

    // Generate input vector
    std::vector<std::complex<double>> input;
    input.reserve(len_softmax);
    // random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    // Set norm paramters
    double mean = -static_cast<double>(m_bound) / 2.0;  // -64.0
    double stddev = static_cast<double>(m_bound) / 6.0; // 约21.33
    double input_lower_bound = -static_cast<double>(m_bound); // -128.0
    double input_upper_bound = 0.0;
    std::normal_distribution<double> normal_dist(mean, stddev);
    std::cout << "Generating " << len_softmax << " samples from normal distribution:" << std::endl;
    std::cout << "Mean: " << mean << ", Stddev: " << stddev << std::endl;
    std::cout << "Range: [" << input_lower_bound << ", " << input_upper_bound << "]" << std::endl;
    
    for (uint32_t i = 0; i < len_softmax; i++) {
        double value;
        do {
            value = normal_dist(gen);
        } while (value < input_lower_bound || value > input_upper_bound);  // 截断超出范围的值

        input.push_back(std::complex<double>(value, 0.0));
    }
    std::cout << "Generated input vector:" << std::endl;
    for (uint32_t i = 0; i < len_softmax; i++) {
        std::cout << "input[" << i << "] = " << input[i].real() << std::endl;
    }

    // full slots packing
    Plaintext plaintext  = cc->MakeCKKSPackedPlaintext(input);
    auto ciphertext      = cc->Encrypt(keyPair.publicKey, plaintext);

    auto result       = cc->EvalSoftmaxUsingCho_Algo2(ciphertext, polyDegree, 1, slotnumInUse, len_softmax);

    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, result, &plaintextDec);
    plaintextDec->SetLength(slotnumInUse);

    std::vector<std::complex<double>> expectedOutput = ComputeSoftmax(input);
    std::cout << "Expected output\n\t" << expectedOutput << std::endl;

    std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;


}

// Automatically compute softmax
std::vector<std::complex<double>> ComputeSoftmax(const std::vector<std::complex<double>>& input) {
    if (input.empty()) {
        return {};
    }
    
    size_t n = input.size();
    std::vector<double> real_parts(n);
    
    // 提取实部
    for (size_t i = 0; i < n; i++) {
        real_parts[i] = input[i].real();
    }
    
    // 数值稳定的softmax计算
    // 步骤1：找到最大值以避免数值溢出
    double max_val = *std::max_element(real_parts.begin(), real_parts.end());
    
    // 步骤2：计算 exp(x_i - max_val)
    std::vector<double> exp_values(n);
    double sum_exp = 0.0;
    
    for (size_t i = 0; i < n; i++) {
        exp_values[i] = std::exp(real_parts[i] - max_val);
        sum_exp += exp_values[i];
    }
    
    // 步骤3：计算softmax = exp(x_i - max_val) / sum(exp(x_j - max_val))
    std::vector<std::complex<double>> softmax_result(n);
    for (size_t i = 0; i < n; i++) {
        double softmax_val = exp_values[i] / sum_exp;
        softmax_result[i] = std::complex<double>(softmax_val, 0.0);
    }
    
    return softmax_result;
}

void TestInverseSqrtApproximation_Chen() {
    std::cout << "=== x^(-1/2) Chebyshev Approximation Analysis ===" << std::endl;
    
    // 定义 x^(-1/2) 函数
    auto inverseSqrt = [](double x) -> double { 
        if (x <= 0.0) {
            throw std::invalid_argument("x must be positive for x^(-1/2)");
        }
        return 1.0 / std::sqrt(x); 
    };
    
    // 设置逼近区间 [a, b] - 必须是正数区间
    double a = 0.125;   // 下界
    double b = 1.0;  // 上界
    uint32_t degree = 16;  // 多项式度数
    
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nFunction: x^(-1/2) on [" << a << ", " << b << "]" << std::endl;
    std::cout << "Polynomial degree: " << degree << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    // 计算逼近误差
    double maxError = EvalChebyshevApproximationError(inverseSqrt, a, b, degree);
    std::cout << "Maximum absolute error: " << maxError << std::endl;
    
    // 计算精度维数
    double precisionDigits = EvalChebyshevPrecisionDigits(inverseSqrt, a, b, degree);
    std::cout << "Precision (decimal digits): " << precisionDigits << std::endl;
    
    // 分析不同度数的精度变化
    std::cout << "\n=== Precision vs Degree Analysis ===" << std::endl;
    std::cout << "Degree\tMax Error\t\tPrecision (digits)" << std::endl;
    std::cout << "------\t---------\t\t------------------" << std::endl;
    
    for (uint32_t deg = 3; deg <= 10; ++deg) {
        double error = EvalChebyshevApproximationError(inverseSqrt, a, b, deg);
        double precision = EvalChebyshevPrecisionDigits(inverseSqrt, a, b, deg);
        std::cout << deg << "\t" << std::scientific << std::setprecision(6) << error 
                  << "\t\t" << std::fixed << std::setprecision(2) << precision << std::endl;
    }
    
    // 显示 Chebyshev 系数
    std::cout << "\n=== Chebyshev Coefficients for degree " << degree << " ===" << std::endl;
    auto coefficients = EvalChebyshevCoefficients(inverseSqrt, a, b, degree);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        std::cout << "c[" << i << "] = " << std::scientific << std::setprecision(8) 
                  << coefficients[i] << std::endl;
    }
}

// x^(-1/2) in [0.125, 1.0] degree 16 minimax coefficients -- using sollya
// index：c0, c1, c2, ..., c16
// std::vector<double> inverseSqrtCoefficients_16 = {
//     7.8651000577494472615545734433510411216485229572126,          // c0 
//     -118.334172206421257459319486176677528930798734959243,        // c1 (x)
//     1447.8283259073668051398481817155815702289220627149,          // c2 (x^2)
//     -12651.1998742763880719655984425890144783529292141994,        // c3 (x^3)
//     80878.850667914528669803087010966608333290920944268,          // c4 (x^4)
//     -3.87873659346168004903515404862896669730788474161e5,         // c5 (x^5)
//     1.4206871159263745662745587748872063714883992570628e6,        // c6 (x^6)
//     -4.0202586353765137660089811847775201169068288306554e6,       // c7 (x^7)
//     8.8417899176579085494832779559985687036802360002235e6,        // c8 (x^8)
//     -1.51249666967800698696156681944843054642021585090437e7,      // c9 (x^9)
//     2.0030712615711983811933412256791334158261327840513e7,        // c10 (x^10)
//     -2.03098683751543675710887104156972051550425627447135e7,      // c11 (x^11)
//     1.5454444464766109789359488104752212309438662008134e7,        // c12 (x^12)
//     -8.5333834474530857722730149460532410099773307495878e6,       // c13 (x^13)
//     3.22574925598832079715795327371008422202376639673525e6,        // c14 (x^14)
//     -7.4623003549269075371761697078648342227308357830233e5,       // c15 (x^15)
//     79633.46950683225426570444163916755408884981338578            // c16 (x^16)
// };

// x^(-1/2) in [0.125, 1.0] degree 8 minimax coefficients -- using sollya
// index：c0, c1, c2, ..., c8
// std::vector<double> inverseSqrtCoefficients_8 = {
//     5.7451840735225940750552800360033244445618602829813,         // c0
//     -42.927230735903967807200581394652935289743556751988,        // c1 (x)
//     234.76071979544915608532588925745475667362453475232,         // c2 (x^2)
//     -814.2398036540641779887674777253577059776689322037,         // c3 (x^3)
//     1800.6157499647011358128579077218671016395667559305,         // c4 (x^4)
//     -2526.0301792801780840024530925539305538332334599804,        // c5 (x^5)
//     2171.02078659047653933225993925651062898215426137,           // c6 (x^6)
//     -1040.90401129487980491019648887409102837959644338427,       // c7 (x^7)
//     212.959794033816828693563117004474060873978041917916         // c8 (x^8)
// };
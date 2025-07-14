/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-06-02 23:27:17
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-06-27 10:11:45
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
void ringPackingTest();


int main(int argc, char* argv[]) {
    TestInverseSqrtApproximation_Chen();

    // EvalSoftmaxExample();
    
    // EvalSoftmaxExample_Cho();

    // ringPackingTest();
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

//InvSqrt using chebyshev
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
    double a = 1.0/256.0;   // 下界
    double b = 1.0;  // 上界
    uint32_t degree = 32;  // 多项式度数
    
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nFunction: x^(-1/2) on [" << a << ", " << b << "]" << std::endl;
    std::cout << "Polynomial degree: " << degree << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // 计算系数
    auto coefficients = EvalChebyshevCoefficients(inverseSqrt, a, b, degree);
    
    // 计算逼近误差
    double maxError = EvalChebyshevApproximationError(inverseSqrt, a, b, degree, coefficients);
    std::cout << "Maximum absolute error: " << maxError << std::endl;
    
    // 计算精度维数
    double precisionDigits = EvalChebyshevPrecisionDigits(inverseSqrt, a, b, degree, coefficients);
    std::cout << "Precision (digits): " << precisionDigits << std::endl;
    
    // 显示 Chebyshev 系数
    // std::cout << "\n=== Chebyshev Coefficients for degree " << degree << " ===" << std::endl;
    // for (size_t i = 0; i < coefficients.size(); ++i) {
    //     std::cout << "c[" << i << "] = " << std::scientific << std::setprecision(8) 
    //               << coefficients[i] << std::endl;
    // }

    std::cout << "================== Remez-Newton ==================" << std::endl;
    // 加权Remez近似的结果
    std::vector<double> remez_poly8={13.8402265776084348688918600300353053395410080429087,
    -330.30231375541776259455014025942384818790992306028,
    3726.6482501724203606228326417338892533966200047035,
    -21060.443928621267844424940653903288363681320682549,
    65884.433122419486014043084306114819855206651576107,
    -1.1943620297249377029791605574530979171222686047278e5,
    1.24813897967780739332704467364084781037083649617862e5,
    -69734.080653062218310244918210064240972030685068457,
    16123.4224505302005410384246075221246062231994098755};

    std::vector<double> remez_poly16 = {19.178480355053953310447819770168738504950680450193,
    -1229.08034434872635058909893302850287765310972520668,
    45491.254938255845278300763402686409135153598260506,
    -9.451772284548826595953901516547732884612483654958e5,
    1.20845629273147707008970002647072889325379579626483e7,
    -1.01968108029664395559036491860350669118526962543327e8,
    5.9619714465561544328159236650328031958678154828313e8,
    -2.4973270253415289976168061263180512280581616679784e9,
    7.6595830949205634462466643691373512975801663280381e9,
    -1.7420528138084951910582373215817641727464918059521e10,
    2.950758214937015434923832238267296853990019212954e10,
    -3.706831538059775466347567180637783643494094984964e10,
    3.40377382079951728619555481613062647751480064715826e10,
    -2.21854230366203524752469282926333984468739681452306e10,
    9.7173054240518971201165482218051981232382743610754e9,
    -2.5639036019826462488286810021436661617468817652622e9,
    3.078756036729417720526553508507125892488442894146e8};


    double maxerror_poly = EvalPolynomialApproximationRelativeError(inverseSqrt, a, b, remez_poly16);
    std::cout << "Maximum relative error for Remez polynomial: " << maxerror_poly << std::endl;
    double precision_poly = EvalPolynomialPrecisionDigits(inverseSqrt, a, b, remez_poly8);
    std::cout << "Precision (digits) for Remez polynomial: " << precision_poly << std::endl;
    
    // ========== Newton迭代增强测试 ==========
    std::cout << "\n=== Newton Iteration Enhancement ===" << std::endl;
    
    // 单个值测试
    double test_x = 0.5;
    std::cout << "Testing at x = " << test_x << std::endl;
    
    // double remez_result = EvalPolynomial(remez_poly16, test_x);
    // double newton_result = EvalNewtonIteration(remez_poly16, test_x, 5);
    double remez_result = EvalPolynomial(remez_poly8, test_x);
    double newton_result = EvalNewtonIteration(remez_poly8, test_x, 5);
    double exact_value = inverseSqrt(test_x);
    
    std::cout << "Remez result:    " << std::fixed << std::setprecision(15) << remez_result << std::endl;
    std::cout << "Newton result:   " << std::fixed << std::setprecision(15) << newton_result << std::endl;
    std::cout << "Exact value:     " << std::fixed << std::setprecision(15) << exact_value << std::endl;
    std::cout << "Remez error:     " << std::scientific << std::setprecision(6) << std::abs(remez_result - exact_value) << std::endl;
    std::cout << "Newton error:    " << std::scientific << std::setprecision(6) << std::abs(newton_result - exact_value) << std::endl;
    
    // 区间精度分析
    // auto newton_analysis = EvalNewtonEnhancedApproximation(remez_poly16, inverseSqrt, a, b, 5, 1000);
    auto newton_analysis = EvalNewtonEnhancedApproximation(remez_poly8, inverseSqrt, a, b, 5, 1000);
    std::cout << "\nNewton-enhanced approximation on [" << a << ", " << b << "]:" << std::endl;
    std::cout << "Max error: " << std::scientific << std::setprecision(6) << newton_analysis.maxError << std::endl;
    std::cout << "Precision: " << std::fixed << std::setprecision(2) << newton_analysis.precisionDigits << " digits" << std::endl;
    std::cout << "Improvement: " << std::scientific << newton_analysis.improvementFactor << "x" << std::endl;
    
    // 迭代次数分析
    std::cout << "\n--- Newton Iterations Analysis ---" << std::endl;
    std::cout << "Iterations\tMax Error\t\tPrecision (digits)" << std::endl;
    std::cout << "----------\t---------\t\t------------------" << std::endl;
    
    // auto iter_analysis = EvalNewtonIterationAnalysis(remez_poly16, inverseSqrt, a, b, 6, 100);
    auto iter_analysis = EvalNewtonIterationAnalysis(remez_poly8, inverseSqrt, a, b, 6, 100);
    for (size_t i = 0; i < iter_analysis.size(); ++i) {
        std::cout << (i+1) << "\t\t" << std::scientific << std::setprecision(4) << iter_analysis[i].maxError 
                  << "\t\t" << std::fixed << std::setprecision(2) << iter_analysis[i].precisionDigits << std::endl;
    }
}

void ringPackingTest() {
    // Test the sparse ring packing
    std::cout << "=== Ring Packing Test ===" << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 4);

    uint32_t numSlots = 2; // Number of slots to pack
    parameters.SetBatchSize(numSlots);

    #if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        // Currently, only FIXEDMANUAL and FIXEDAUTO modes are supported for 128-bit CKKS bootstrapping.
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits               = 78;
        usint firstMod               = 89;
    #else
        // All modes are supported for 64-bit CKKS bootstrapping.
        ScalingTechnique rescaleTech = FLEXIBLEAUTO;
        usint dcrtBits               = 59;
        usint firstMod               = 60;
    #endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    uint32_t depth = 2; // Set a small depth for testing
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);

    uint32_t ringDim = cc->GetRingDimension();
    std::cout << "Ring dimension: " << ringDim << std::endl;
    std::cout << "Number of slots: " << numSlots << std::endl;

    // Generate random input
    std::vector<double> x;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (size_t i = 0; i < numSlots; i++) {
        x.push_back(dis(gen));
    }

    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x, 1, 0, nullptr, numSlots);
    ptxt->SetLength(numSlots);
    std::cout << "Input: " << ptxt << std::endl;

    //Get to plaintext polynomial coefficients
    auto poly_ptxt = ptxt->GetElement<DCRTPoly>();
    std::cout << "Number of towers: " << poly_ptxt.GetNumOfElements() << std::endl;

    auto tower0 = poly_ptxt.GetElementAtIndex(0);
    std::cout << "Tower 0 coefficients (first 16):" << std::endl;
    for (size_t i = 0; i < std::min(usint(16), tower0.GetLength()); ++i) {
        std::cout << "coeff[" << i << "] = " << tower0[i] << std::endl;
    }
    
}


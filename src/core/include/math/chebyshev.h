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
 * This code provides Chebyshev approximation utilities.
 */

#ifndef LBCRYPTO_INC_MATH_CHEBYSHEV_H
#define LBCRYPTO_INC_MATH_CHEBYSHEV_H

#include <cstdint>
#include <functional>
#include <vector>

/**
 * @namespace lbcrypto
 * The namespace of lbcrypto
 */
namespace lbcrypto {

// move to cryptocontext?

/**
 * Newton增强结果结构体
 */
struct NewtonEnhancedResult {
    double maxError;
    double avgError;
    double maxRelativeError;
    double avgRelativeError;  
    double precisionDigits;
    double improvementFactor;
};

/**
 * Method for calculating Chebyshev coefficients for an input function
 * over the range [a,b]. These coefficents can be input into
 * EValChebyshevSeries to evaluate the function.
 *
 * @param func is the function to be approximated
 * @param a - lower bound of argument for which the coefficients were found
 * @param b - upper bound of argument for which the coefficients were found
 * @param degree Desired degree of approximation
 * @return the coefficients of the Chebyshev approximation.
 */
std::vector<double> EvalChebyshevCoefficients(std::function<double(double)> func, double a, double b, uint32_t degree);

/**
 * Calculate the approximation error of Chebyshev polynomial approximation
 * for a given function over the range [a,b] with specified degree.
 *
 * @anchor SeehowLi
 * @param func is the function to be approximated
 * @param a - lower bound of the approximation range
 * @param b - upper bound of the approximation range  
 * @param degree - degree of the Chebyshev polynomial
 * @param coefficients - coefficients of the Chebyshev polynomial
 * @param numTestPoints - number of test points to evaluate error (default: 1000)
 * @return the maximum absolute error of the approximation
 */
double EvalChebyshevApproximationError(std::function<double(double)> func, double a, double b, 
                                      uint32_t degree, std::vector<double> coefficients, size_t numTestPoints = 1000);

/**
 * Calculate the precision (in decimal digits) of Chebyshev polynomial approximation
 * for a given function over the range [a,b] with specified degree.
 *
 * @author Seehow Li
 * @param func is the function to be approximated
 * @param a - lower bound of the approximation range
 * @param b - upper bound of the approximation range
 * @param degree - degree of the Chebyshev polynomial  
 * @param coefficients - coefficients of the Chebyshev polynomial
 * @param numTestPoints - number of test points to evaluate precision (default: 1000)
 * @return the precision in decimal digits (higher value means better precision)
 */
double EvalChebyshevPrecisionDigits(std::function<double(double)> func, double a, double b,
                                   uint32_t degree, std::vector<double> coefficients, size_t numTestPoints = 1000);

/**
 * 计算多项式逼近的最大误差
 * @param targetFunc 目标函数
 * @param a 区间下界
 * @param b 区间上界
 * @param coefficients 多项式系数数组，从0次项到最高次项
 * @param numTestPoints 测试点数量（默认1000）
 * @return 最大绝对误差
 */
double EvalPolynomialApproximationError(std::function<double(double)> targetFunc,
                                       double a, double b, const std::vector<double>& coefficients,
                                       size_t numTestPoints = 1000);

/**
 * 计算多项式逼近的最大误差
 * @param targetFunc 目标函数
 * @param a 区间下界
 * @param b 区间上界
 * @param coefficients 多项式系数数组，从0次项到最高次项
 * @param numTestPoints 测试点数量（默认1000）
 * @return 最大相对误差
 */
double EvalPolynomialApproximationRelativeError(std::function<double(double)> targetFunc,
                                       double a, double b, const std::vector<double>& coefficients,
                                       size_t numTestPoints = 1000);

/**
 * 计算多项式逼近的精度（二进制）
 * @param targetFunc 目标函数
 * @param a 区间下界
 * @param b 区间上界
 * @param coefficients 多项式系数数组，从0次项到最高次项
 * @param numTestPoints 测试点数量（默认1000）
 * @return 近似的二进制位数
 */
double EvalPolynomialPrecisionDigits(std::function<double(double)> func, double a, double b,
                                std::vector<double> coefficients, size_t numTestPoints = 1000);


/**
 * 使用多项式系数计算多项式值（Horner方法）
 * @param coeffs 多项式系数数组，从0次项到最高次项
 * @param x 计算点
 * @return 多项式在x处的值
 */
double EvalPolynomial(const std::vector<double>& coeffs, double x);

/**
 * 对单个值进行Newton迭代增强计算 1/sqrt(x)
 * @param remez_coeffs Remez多项式系数（作为初始逼近）
 * @param x 计算点
 * @param num_iterations Newton迭代次数（默认5次）
 * @return Newton增强后的结果
 */
double EvalNewtonIteration(const std::vector<double>& remez_coeffs, double x, int num_iterations = 5);

/**
 * 优化的Newton迭代（4次乘法版本）
 * @param remez_coeffs Remez多项式系数
 * @param x 计算点  
 * @param num_iterations Newton迭代次数（默认5次）
 * @return 优化的Newton增强结果
 */
double EvalNewtonIterationOptimized(const std::vector<double>& remez_coeffs, double x, int num_iterations = 5);

/**
 * 在指定区间上评估Newton增强的逼近效果
 * @param remez_coeffs Remez多项式系数
 * @param target_func 目标函数
 * @param a 区间下界
 * @param b 区间上界
 * @param num_iterations Newton迭代次数（默认5次）
 * @param num_test_points 测试点数量（默认1000）
 * @return Newton增强的详细结果
 */
NewtonEnhancedResult EvalNewtonEnhancedApproximation(
    const std::vector<double>& remez_coeffs,
    std::function<double(double)> target_func,
    double a, double b,
    int num_iterations = 5,
    size_t num_test_points = 1000);

/**
 * 批量计算Newton增强结果
 * @param remez_coeffs Remez多项式系数
 * @param x_values 输入值数组
 * @param num_iterations Newton迭代次数（默认5次）
 * @return Newton增强结果数组
 */
std::vector<double> EvalNewtonIterationBatch(
    const std::vector<double>& remez_coeffs,
    const std::vector<double>& x_values,
    int num_iterations = 5);

/**
 * 分析不同Newton迭代次数的精度效果
 * @param remez_coeffs Remez多项式系数
 * @param target_func 目标函数
 * @param a 区间下界
 * @param b 区间上界
 * @param max_iterations 最大迭代次数
 * @param num_test_points 测试点数量
 * @return 每个迭代次数对应的精度结果
 */
std::vector<NewtonEnhancedResult> EvalNewtonIterationAnalysis(
    const std::vector<double>& remez_coeffs,
    std::function<double(double)> target_func,
    double a, double b,
    int max_iterations = 6,
    size_t num_test_points = 100);

}  // namespace lbcrypto

#endif


//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2023, NJIT, Duality Technologies Inc. and other contributors
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
  This code provides Chebyshev approximation utilities
 */

#include "math/chebyshev.h"
#include "utils/exception.h"

#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

namespace lbcrypto {

std::vector<double> EvalChebyshevCoefficients(std::function<double(double)> func, double a, double b, uint32_t degree) {
    if (!degree) {
        OPENFHE_THROW("The degree of approximation can not be zero");
    }
    // the number of coefficients to be generated should be degree+1 as zero is also included
    size_t coeffTotal{degree + 1};
    double bMinusA = 0.5 * (b - a);
    double bPlusA  = 0.5 * (b + a);
    double PiByDeg = M_PI / static_cast<double>(coeffTotal);
    std::vector<double> functionPoints(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i)
        functionPoints[i] = func(std::cos(PiByDeg * (i + 0.5)) * bMinusA + bPlusA);

    double multFactor = 2.0 / static_cast<double>(coeffTotal);
    std::vector<double> coefficients(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i) {
        for (size_t j = 0; j < coeffTotal; ++j)
            coefficients[i] += functionPoints[j] * std::cos(PiByDeg * i * (j + 0.5));
        coefficients[i] *= multFactor;
    }
    return coefficients;
}


double EvalChebyshevApproximationError(std::function<double(double)> func, double a, double b, 
                                      uint32_t degree, std::vector<double> coefficients, size_t numTestPoints) {
    if (!degree) {
        OPENFHE_THROW("The degree of approximation can not be zero");
    }
    
    if (numTestPoints == 0) {
        OPENFHE_THROW("Number of test points must be greater than zero");
    }
    
    if (a >= b) {
        OPENFHE_THROW("Lower bound must be less than upper bound");
    }
    
    // 获取 Chebyshev 系数
    // std::vector<double> coefficients = EvalChebyshevCoefficients(func, a, b, degree);
    
    // 传递系数

    // 计算测试点上的最大误差
    double maxError = 0.0;
    double stepSize = (b - a) / (numTestPoints - 1);
    
    for (size_t i = 0; i < numTestPoints; ++i) {
        double x = a + i * stepSize;
        
        // 计算真实函数值
        double trueValue = func(x);
        
        // 计算 Chebyshev 逼近值
        // 将 x 映射到 [-1, 1] 区间
        double mappedX = (2.0 * x - a - b) / (b - a);
        
        // 使用 Chebyshev 多项式计算逼近值
        double approxValue = 0.0;
        double T0 = 1.0;              // T_0(x) = 1
        double T1 = mappedX;          // T_1(x) = x
        
        // 处理 T_0 项 (第一个系数需要除以2)
        approxValue += coefficients[0] / 2.0;
        
        if (degree >= 1) {
            approxValue += coefficients[1] * T1;
        }
        
        // 使用递推关系计算高阶 Chebyshev 多项式: T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
        for (uint32_t j = 2; j <= degree; ++j) {
            double T2 = 2.0 * mappedX * T1 - T0;
            approxValue += coefficients[j] * T2;
            T0 = T1;
            T1 = T2;
        }
        
        // 计算绝对误差
        double error = std::abs(trueValue - approxValue);
        if (error > maxError) {
            maxError = error;
        }
    }
    
    return maxError;
}

double EvalChebyshevPrecisionDigits(std::function<double(double)> func, double a, double b,
                                   uint32_t degree, std::vector<double> coefficients, size_t numTestPoints) {
    double maxError = EvalChebyshevApproximationError(func, a, b, degree, coefficients, numTestPoints);
    
    if (maxError <= 0.0) {
        return std::numeric_limits<double>::infinity();  // 完美逼近
    }
    
    // 精度位数 = -log2(最大误差)
    return -std::log2(maxError);
}

double EvalPolynomialApproximationError(std::function<double(double)> targetFunc,
                                        double a, double b, const std::vector<double>& coefficients,
                                       size_t numTestPoints) {
    if (coefficients.empty()) {
        OPENFHE_THROW("Coefficients vector cannot be empty");
    }
    
    if (numTestPoints == 0) {
        OPENFHE_THROW("Number of test points must be greater than zero");
    }
    
    if (a >= b) {
        OPENFHE_THROW("Lower bound must be less than upper bound");
    }
    
    double maxError = 0.0;
    double stepSize = (b - a) / (numTestPoints - 1);
    
    for (size_t i = 0; i < numTestPoints; ++i) {
        double x = a + i * stepSize;
        
        // 计算目标函数的真实值
        double trueValue = targetFunc(x);
        
        // 计算多项式逼近值: P(x) = c₀ + c₁x + c₂x² + ... + cₙxⁿ
        double polyValue = 0.0;
        double xPower = 1.0;  // x^0 = 1
        
        for (size_t j = 0; j < coefficients.size(); ++j) {
            polyValue += coefficients[j] * xPower;
            xPower *= x;  // 计算下一个幂次
        }
        
        // 计算绝对误差
        double error = std::abs(trueValue - polyValue);
        if (error > maxError) {
            maxError = error;
            std::cout << "Max error occurs at x= " << x << std::endl;
        }
    }
    
    return maxError;
}

double EvalPolynomialApproximationRelativeError(std::function<double(double)> targetFunc,
                                        double a, double b, const std::vector<double>& coefficients,
                                       size_t numTestPoints) {
    if (coefficients.empty()) {
        OPENFHE_THROW("Coefficients vector cannot be empty");
    }
    
    if (numTestPoints == 0) {
        OPENFHE_THROW("Number of test points must be greater than zero");
    }
    
    if (a >= b) {
        OPENFHE_THROW("Lower bound must be less than upper bound");
    }
    
    double maxRelativeError = 0.0;
    double stepSize = (b - a) / (numTestPoints - 1);
    
    for (size_t i = 0; i < numTestPoints; ++i) {
        double x = a + i * stepSize;
        
        // 计算目标函数的真实值
        double trueValue = targetFunc(x);
        
        // 计算多项式逼近值: P(x) = c₀ + c₁x + c₂x² + ... + cₙxⁿ
        double polyValue = 0.0;
        double xPower = 1.0;  // x^0 = 1
        
        for (size_t j = 0; j < coefficients.size(); ++j) {
            polyValue += coefficients[j] * xPower;
            xPower *= x;  // 计算下一个幂次
        }
        
        // 计算相对误差
        double relative_error = (std::abs(trueValue - polyValue)) / (std::abs(trueValue) + 1e-15);  // 防止除以零
        if (relative_error > maxRelativeError) {
            maxRelativeError = relative_error;
            std::cout << "Max relative error occurs at x= " << x << std::endl;
        }
    }
    
    return maxRelativeError;
}

double EvalPolynomialPrecisionDigits(std::function<double(double)> func, double a, double b,
                                     std::vector<double> coefficients, size_t numTestPoints) {
    double maxError = EvalPolynomialApproximationError(func, a, b, coefficients, numTestPoints);
    
    if (maxError <= 0.0) {
        return std::numeric_limits<double>::infinity();  // 完美逼近
    }
    
    // 精度位数 = -log2(最大误差)
    return -std::log2(maxError);
}

double EvalPolynomial(const std::vector<double>& coeffs, double x) {
    if (coeffs.empty()) return 0.0;
    
    // 使用Horner方法提高数值稳定性
    double result = coeffs.back();
    for (int i = static_cast<int>(coeffs.size()) - 2; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

double EvalNewtonIteration(const std::vector<double>& remez_coeffs, double x, int num_iterations) {
    if (remez_coeffs.empty()) {
        OPENFHE_THROW("Remez coefficients vector cannot be empty");
    }
    
    if (x <= 0.0) {
        OPENFHE_THROW("Input x must be positive for inverse square root calculation");
    }
    
    // 使用Remez多项式作为初始近似
    double y = EvalPolynomial(remez_coeffs, x);
    
    // Newton迭代: y_{n+1} = y_n * (3 - x * y_n^2) / 2
    // 这是计算 1/sqrt(x) 的Newton迭代公式
    for (int i = 0; i < num_iterations; ++i) {
        double y_squared = y * y;
        double x_y_squared = x * y_squared;
        y = y * (3.0 - x_y_squared) * 0.5;
    }
    
    return y;
}

double EvalNewtonIterationOptimized(const std::vector<double>& remez_coeffs, double x, int num_iterations) {
    if (remez_coeffs.empty()) {
        OPENFHE_THROW("Remez coefficients vector cannot be empty");
    }
    
    if (x <= 0.0) {
        OPENFHE_THROW("Input x must be positive for inverse square root calculation");
    }
    
    double y = EvalPolynomial(remez_coeffs, x);
    double x_half = x * 0.5;
    
    for (int i = 0; i < num_iterations; ++i) {
        // 论文中的优化算法：减少乘法次数
        double z1 = x_half * y;       // l(z1) = max(l(x/2), l(y)) + 1
        double z2 = y * y;            // l(z2) = l(y) + 1
        double z3 = 1.5 * y;          // l(z3) = l(y) + 1
        y = z3 - z1 * z2;             // level max(l(x/2), l(y)) + 2
    }
    
    return y;
}

NewtonEnhancedResult EvalNewtonEnhancedApproximation(
    const std::vector<double>& remez_coeffs,
    std::function<double(double)> target_func,
    double a, double b,
    int num_iterations,
    size_t num_test_points) {
    
    if (remez_coeffs.empty()) {
        OPENFHE_THROW("Remez coefficients vector cannot be empty");
    }
    
    if (num_test_points == 0) {
        OPENFHE_THROW("Number of test points must be greater than zero");
    }
    
    if (a >= b) {
        OPENFHE_THROW("Lower bound must be less than upper bound");
    }
    
    NewtonEnhancedResult result = {};
    
    double total_abs_error = 0.0;
    double total_rel_error = 0.0;
    double step_size = (b - a) / (num_test_points - 1);
    
    // 同时计算原始Remez误差用于比较
    double max_remez_error = 0.0;
    
    for (size_t i = 0; i < num_test_points; ++i) {
        double x = a + i * step_size;
        
        try {
            // 计算各种结果
            double exact = target_func(x);
            double remez_only = EvalPolynomial(remez_coeffs, x);
            double newton_enhanced = EvalNewtonIteration(remez_coeffs, x, num_iterations);
            
            // 计算误差
            double remez_error = std::abs(remez_only - exact);
            double newton_error = std::abs(newton_enhanced - exact);
            double rel_error = (std::abs(exact) > 1e-15) ? newton_error / std::abs(exact) : 0.0;
            
            // 更新统计信息
            max_remez_error = std::max(max_remez_error, remez_error);
            result.maxError = std::max(result.maxError, newton_error);
            result.maxRelativeError = std::max(result.maxRelativeError, rel_error);
            
            total_abs_error += newton_error;
            total_rel_error += rel_error;
        } catch (const std::exception& e) {
            // 跳过有问题的点
            continue;
        }
    }
    
    result.avgError = total_abs_error / num_test_points;
    result.avgRelativeError = total_rel_error / num_test_points;
    result.precisionDigits = (result.maxError > 0) ? -std::log2(result.maxError) : std::numeric_limits<double>::infinity();
    result.improvementFactor = (result.maxError > 0) ? max_remez_error / result.maxError : std::numeric_limits<double>::infinity();
    
    return result;
}

std::vector<double> EvalNewtonIterationBatch(
    const std::vector<double>& remez_coeffs,
    const std::vector<double>& x_values,
    int num_iterations) {
    
    std::vector<double> results;
    results.reserve(x_values.size());
    
    for (double x : x_values) {
        try {
            results.push_back(EvalNewtonIteration(remez_coeffs, x, num_iterations));
        } catch (const std::exception& e) {
            // 对于有问题的输入，返回NaN
            results.push_back(std::numeric_limits<double>::quiet_NaN());
        }
    }
    
    return results;
}

std::vector<NewtonEnhancedResult> EvalNewtonIterationAnalysis(
    const std::vector<double>& remez_coeffs,
    std::function<double(double)> target_func,
    double a, double b,
    int max_iterations,
    size_t num_test_points) {
    
    std::vector<NewtonEnhancedResult> results;
    results.reserve(max_iterations);
    
    for (int iter = 1; iter <= max_iterations; ++iter) {
        NewtonEnhancedResult result = EvalNewtonEnhancedApproximation(
            remez_coeffs, target_func, a, b, iter, num_test_points);
        results.push_back(result);
    }
    
    return results;
}

}  // namespace lbcrypto

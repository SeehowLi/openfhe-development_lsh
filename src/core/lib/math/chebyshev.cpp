
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
                                      uint32_t degree, size_t numTestPoints) {
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
    std::vector<double> coefficients = EvalChebyshevCoefficients(func, a, b, degree);
    
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
                                   uint32_t degree, size_t numTestPoints) {
    double maxError = EvalChebyshevApproximationError(func, a, b, degree, numTestPoints);
    
    if (maxError <= 0.0) {
        return std::numeric_limits<double>::infinity();  // 完美逼近
    }
    
    // 精度位数 = -log10(最大误差)
    return -std::log10(maxError);
}

}  // namespace lbcrypto

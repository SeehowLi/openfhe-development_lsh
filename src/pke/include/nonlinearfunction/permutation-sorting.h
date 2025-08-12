/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 21:20:40
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 12:18:27
 * @FilePath: \openfhe-development\src\pke\include\nonlinearfunction\permutation-sorting.h
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#ifndef SRC_NONLINEARFUNCTION_PERMUTATION_SORTING_H_
#define SRC_NONLINEARFUNCTION_PERMUTATION_SORTING_H_

#include "sorting-method.h"
#include <cmath>            // ✅ 添加这个（用于std::log2等数学函数）

namespace lbcrypto{

class PermutationSorting : public SortingMethod {
private:
    int m_sigmoidScaling;
    int m_sigmoidDegree;
    int m_sincDegree;
    bool m_enableTieoffset;
    int m_n;
    double m_delta;
    
public:
    explicit PermutationSorting(HomoEncryptCompute& hec, const SortingParams& params);
    
    // 实现基类接口 - 使用SortingInput
    Cipher Sort(const SortingInput& input, const SortingParams& params) override;

    std::string GetMethodName() const override { return "Permutation-based"; }
    int EstimateDepthRequired(const SortingParams& params) const override;
    
private:
    // 主要接口 - 使用两个密文
    Cipher DoPermutationSort(const Cipher& in_exp, const Cipher& in_rep);
    
    // 内部排序实现（使用扩展的输入）
    Cipher InternalSort(const Cipher& in_exp, const Cipher& in_rep);
    
    // 核心算法步骤
    Cipher ComputeIndexing(const Cipher& in_exp, const Cipher& in_rep);
    Cipher ComputeTieOffset(const Cipher& in_exp, const Cipher& in_rep);
    Cipher ComputeSorting(const Cipher& indexes, const Cipher& in_rep);
    
    // 参数自动配置
    void AutoConfigureParameters();
};

}

#endif  // SRC_NONLINEARFUNCTION_PERMUTATION_SORTING_H_
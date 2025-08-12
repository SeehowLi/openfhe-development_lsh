/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 21:21:02
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 11:10:34
 * @FilePath: \openfhe-development\src\pke\include\nonlinearfunction\network-sorting.h
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#ifndef SRC_NONLINEARFUNCTION_NETWORK_SORTING_H_
#define SRC_NONLINEARFUNCTION_NETWORK_SORTING_H_

#include "sorting-method.h"
#include "vector"
#include <cmath>

namespace lbcrypto {

class NetworkSorting : public SortingMethod {
private:
    int m_n;
    int m_reluDegree;
    double m_inputScale;
    
public:
    explicit NetworkSorting(HomoEncryptCompute& hec, const SortingParams& params);
    
    
    // 实现基类接口 - 使用SortingInput（只用standard部分）
    Cipher Sort(const SortingInput& input, const SortingParams& params) override;
    
    std::string GetMethodName() const override { 
        return "Network-based Sorting"; 
    }
    
    int EstimateDepthRequired(const SortingParams& params) const override;
    
private:
    // 主要实现 - 单个密文输入
    Cipher DoNetworkSort(const Cipher& input, const SortingParams& params);

    // 原有的私有方法
    Cipher Swap(const Cipher& in, int arrowsDelta, int round, int stage);
    std::vector<Plain> GenerateLayerMasks(int level, int slots, 
                                        int round, int stage,
                                        double maskValue = 1.0);
    
    void AutoConfigureReLUDegree(double precision);
};

} // namespace lbcrypto

#endif  // SRC_NONLINEARFUNCTION_NETWORK_SORTING_H_

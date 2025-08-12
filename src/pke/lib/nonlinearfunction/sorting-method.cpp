/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 09:01:33
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 12:37:41
 * @FilePath: \openfhe-development\src\pke\lib\nonlinearfunction\sorting-method.cpp
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#include "nonlinearfunction/sorting-method.h"
#include "nonlinearfunction/permutation-sorting.h"
#include "nonlinearfunction/network-sorting.h"
#include "nonlinearfunction/sorting-utils.h"

#include <iostream>
#include <stdexcept>

namespace lbcrypto {

void SortingMethod::PrintIfVerbose(const std::string& message) {
    if (m_verbose) {
        std::cout << "[" << GetMethodName() << "] " << message << std::endl;
    }
}

std::unique_ptr<SortingMethod> SortingMethod::Create(
    SortingParams::Method method, 
    HomoEncryptCompute& hec,
    const SortingParams& params) {
    
    switch (method) {
        case SortingParams::PERMUTATION:
            return std::make_unique<PermutationSorting>(hec, params);
            
        case SortingParams::NETWORK:
            return std::make_unique<NetworkSorting>(hec, params);
            
        default:
            throw std::runtime_error("Unknown sorting method");
    }
}

} // namespace lbcrypto
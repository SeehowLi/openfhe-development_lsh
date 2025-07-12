/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 23:39:57
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 16:47:17
 * @FilePath: \openfhe-development\src\pke\lib\nonlinearfunction\network-sorting.cpp
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#include "nonlinearfunction/network-sorting.h"
#include "nonlinearfunction/sorting-utils.h"
#include <iostream>
#include <chrono>

namespace lbcrypto {

NetworkSorting::NetworkSorting(HomoEncryptCompute& hec, const SortingParams& params)
    : SortingMethod(hec, params.verbose),
      m_n(params.vectorSize),
      m_reluDegree(params.reluDegree),
      m_inputScale(params.inputScale) {
    
    // 如果未设置ReLU度数，根据精度自动配置
    if (m_reluDegree == 0) {
        AutoConfigureReLUDegree(params.precision);
    }
    
    // 如果未设置输入缩放，使用默认值
    if (m_inputScale == 0) {
        m_inputScale = 0.95;
    }
}

// 使用SortingInput（只用standard部分）
Cipher NetworkSorting::Sort(const SortingInput& input, const SortingParams& params) {
    if (!input.standard) {
        throw std::runtime_error("NetworkSorting requires standard input");
    }
    return DoNetworkSort(*input.standard, params);
}

// 主要实现
Cipher NetworkSorting::DoNetworkSort(const Cipher& input, const SortingParams& params) {
    PrintIfVerbose("Starting network-based sorting");
    PrintIfVerbose("Vector size: " + std::to_string(m_n) + 
                   ", ReLU degree: " + std::to_string(m_reluDegree) +
                   ", Input scale: " + std::to_string(m_inputScale));
    
    // 计算总迭代次数
    int iterations = (std::log2(m_n) * (std::log2(m_n) + 1)) / 2;
    
    // 复制输入（避免修改原始数据）
    Cipher result = m_hec.add(input, 0.0);
    
    int currentIteration = 1;
    
    // Bitonic sorting network主循环
    for (int i = 0; i < std::log2(m_n); i++) {
        for (int j = 0; j < i + 1; j++) {
            int arrowsDelta = std::pow(2, i - j);
            int stage = i - j;
            int round = j;
            
            auto start = std::chrono::steady_clock::now();
            
            PrintIfVerbose("Layer " + std::to_string(currentIteration) + 
                          ": arrowsDelta=" + std::to_string(arrowsDelta) +
                          ", stage=" + std::to_string(stage) + 
                          ", round=" + std::to_string(round));
            
            // 执行比较-交换操作
            result = Swap(result, arrowsDelta, stage, round);
            
            if (m_verbose) {
                print_duration(start, "Swap operation");
            }
            
            // 如果不是最后一次迭代，执行bootstrapping
            if (currentIteration < iterations) {
                start = std::chrono::steady_clock::now();
                result = m_hec.bootstrap(result);
                
                if (m_verbose) {
                    print_duration(start, "Bootstrapping");
                    m_hec.print(result, m_n);
                    std::cout << "Layer " << currentIteration 
                             << " / " << iterations << " completed." << std::endl;
                }
            }
            
            currentIteration++;
        }
    }
    
    PrintIfVerbose("Network-based sorting completed");
    
    return result;
}

Cipher NetworkSorting::Swap(const Cipher& in, int arrowsDelta, int round, int stage) {
    // 旋转以获取比较对
    Cipher rot_pos = m_hec.rot(in, arrowsDelta);
    Cipher rot_neg = m_hec.rot(in, -arrowsDelta);
    
    if (m_verbose) {
        m_hec.print(in, m_n, "  Input: ");
        m_hec.print(rot_pos, m_n, "  Rot+: ");
    }
    
    // 使用ReLU近似计算min函数
    Cipher m1 = m_hec.relu(m_hec.sub(in, rot_pos), m_reluDegree, m_n);
    
    // 根据m1计算其他值
    Cipher m3 = m_hec.sub(m_hec.add(in, rot_pos), m1);
    Cipher m4 = m_hec.rot(m1, -arrowsDelta);
    Cipher m2 = m_hec.sub(m_hec.add(in, rot_neg), m4);

    if (m_verbose) {
        m_hec.print(m1, m_n, "  m1: ");
        m_hec.print(m2, m_n, "  m2: ");
        m_hec.print(m3, m_n, "  m3: ");
        m_hec.print(m4, m_n, "  m4: ");
    }

    // 生成掩码
    auto masks = GenerateLayerMasks(m1->GetLevel(), m1->GetSlots(), round, stage);
    
    // 应用掩码并组合结果
    Cipher result = m_hec.add_tree({
        m_hec.mult(m1, masks[0]),
        m_hec.mult(m2, masks[1]),
        m_hec.mult(m3, masks[2]),
        m_hec.mult(m4, masks[3])
    });
    
    if (m_verbose) {
        m_hec.print(result, m_n, "  Result: ");
    }
    
    return result;
}

// std::vector<Plain> NetworkSorting::GenerateLayerMasks(int encodingLevel, int numSlots, 
//                                                      int round, int stage, 
//                                                      double maskValue) {
//     std::vector<double> mask_1, mask_2, mask_3, mask_4;
    
//     // 原始实现的逻辑更简单
//     for (int i = 0; i < numSlots / (std::pow(2, round + 2)); i++) {
        
//         // 第一部分
//         for (int times = 0; times < std::pow(2, stage); times++) {
//             // mask_1 active
//             for (int j = 0; j < std::pow(2, round); j++) {
//                 mask_1.push_back(maskValue);
//                 mask_2.push_back(0);
//                 mask_3.push_back(0);
//                 mask_4.push_back(0);
//             }
            
//             // mask_2 active
//             for (int j = 0; j < std::pow(2, round); j++) {
//                 mask_1.push_back(0);
//                 mask_2.push_back(maskValue);
//                 mask_3.push_back(0);
//                 mask_4.push_back(0);
//             }
//         }
        
//         if ((i + 1) * std::pow(2, stage + round + 1) >= numSlots) {
//             break;
//         }
        
//         // 第二部分
//         for (int times = 0; times < std::pow(2, stage); times++) {
//             // mask_3 active
//             for (int j = 0; j < std::pow(2, round); j++) {
//                 mask_1.push_back(0);
//                 mask_2.push_back(0);
//                 mask_3.push_back(maskValue);
//                 mask_4.push_back(0);
//             }
            
//             // mask_4 active
//             for (int j = 0; j < std::pow(2, round); j++) {
//                 mask_1.push_back(0);
//                 mask_2.push_back(0);
//                 mask_3.push_back(0);
//                 mask_4.push_back(maskValue);
//             }
//         }
        
//         if ((i + 1) * std::pow(2, stage + round + 2) >= numSlots) {
//             break;
//         }
//     }
    
//     // 确保掩码大小正确
//     mask_1.resize(numSlots, 0);
//     mask_2.resize(numSlots, 0);
//     mask_3.resize(numSlots, 0);
//     mask_4.resize(numSlots, 0);
    
//     return {
//         m_hec.encode(mask_1, encodingLevel, numSlots),
//         m_hec.encode(mask_2, encodingLevel, numSlots),
//         m_hec.encode(mask_3, encodingLevel, numSlots),
//         m_hec.encode(mask_4, encodingLevel, numSlots)
//     };
// }

std::vector<Plain> NetworkSorting::GenerateLayerMasks(int encodingLevel, int numSlots, 
                                                     int round, int stage, 
                                                     double maskValue) {
    std::vector<double> mask_1, mask_2, mask_3, mask_4;
    
    // 使用原始实现的逻辑
    for (int i = 0; i < numSlots / (std::pow(2, round + 2)); i++) {
        
        for (int times = 0; times < std::pow(2, stage); times++) {
            for (int j = 0; j < std::pow(2, round); j++) {
                mask_1.push_back(maskValue);
                mask_2.push_back(0);
                mask_3.push_back(0);
                mask_4.push_back(0);
            }
            
            for (int j = 0; j < std::pow(2, round); j++) {
                mask_1.push_back(0);
                mask_2.push_back(maskValue);
                mask_3.push_back(0);
                mask_4.push_back(0);
            }
        }
        
        if ((i + 1) * std::pow(2, stage + round + 1) >= numSlots) {
            break;
        }
        
        for (int times = 0; times < std::pow(2, stage); times++) {
            for (int j = 0; j < std::pow(2, round); j++) {
                mask_1.push_back(0);
                mask_2.push_back(0);
                mask_3.push_back(maskValue);
                mask_4.push_back(0);
            }
            
            for (int j = 0; j < std::pow(2, round); j++) {
                mask_1.push_back(0);
                mask_2.push_back(0);
                mask_3.push_back(0);
                mask_4.push_back(maskValue);
            }
        }
        
        if ((i + 1) * std::pow(2, stage + round + 2) >= numSlots) {
            break;
        }
    }
    
    // 确保大小正确
    mask_1.resize(numSlots, 0);
    mask_2.resize(numSlots, 0);
    mask_3.resize(numSlots, 0);
    mask_4.resize(numSlots, 0);
    
    return {
        m_hec.encode(mask_1, encodingLevel, numSlots),
        m_hec.encode(mask_2, encodingLevel, numSlots),
        m_hec.encode(mask_3, encodingLevel, numSlots),
        m_hec.encode(mask_4, encodingLevel, numSlots)
    };
}


int NetworkSorting::EstimateDepthRequired(const SortingParams& params) const {
    // ReLU多项式的深度消耗
    int depth = poly_evaluation_cost(m_reluDegree);
    
    // 掩码操作的深度
    depth += 1;
    
    // 加上安全余量
    return depth + 1;
}

void NetworkSorting::AutoConfigureReLUDegree(double precision) {
    if (precision >= 0.1) {
        m_reluDegree = 119;
        PrintIfVerbose("Auto-configured ReLU degree: 119 (for precision >= 0.1)");
    } else if (precision >= 0.01) {
        m_reluDegree = 495;
        PrintIfVerbose("Auto-configured ReLU degree: 495 (for precision >= 0.01)");
    } else if (precision >= 0.001) {
        m_reluDegree = 495;
        PrintIfVerbose("Auto-configured ReLU degree: 495 (for precision >= 0.001)");
    } else {
        m_reluDegree = 495;
        PrintIfVerbose("Auto-configured ReLU degree: 495 (for precision < 0.001)");
    }
}

} // namespace lbcrypto
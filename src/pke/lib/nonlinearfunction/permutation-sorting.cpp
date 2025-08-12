/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 23:02:53
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-27 11:05:07
 * @FilePath: \openfhe-131\src\pke\lib\nonlinearfunction\permutation-sorting.cpp
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#include "nonlinearfunction/permutation-sorting.h"
#include "nonlinearfunction/sorting-utils.h"
#include <cmath>
#include <iostream>

namespace lbcrypto {

PermutationSorting::PermutationSorting(HomoEncryptCompute& hec, const SortingParams& params)
    : SortingMethod(hec, params.verbose),
      m_sigmoidScaling(params.sigmoidScaling),
      m_sigmoidDegree(params.sigmoidDegree),
      m_sincDegree(params.sincDegree),
      m_enableTieoffset(params.enableTieoffset),
      m_n(params.vectorSize),
      m_delta(params.precision) {
    
    // 如果参数未设置，自动配置
    if (m_sigmoidScaling == 0 || m_sigmoidDegree == 0 || m_sincDegree == 0) {
        AutoConfigureParameters();
    }
}

// 实现基类接口 - 使用SortingInput
Cipher PermutationSorting::Sort(const SortingInput& input, const SortingParams& params) {
    std::cout << "DEBUG: PermutationSorting::Sort called" << std::endl;
    
    if (!input.hasExpandedInputs()) {
        throw std::runtime_error("PermutationSorting requires expanded and repeated inputs");
    }
    
    std::cout << "DEBUG: About to call DoPermutationSort" << std::endl;
    Cipher result = DoPermutationSort(*input.expanded, *input.repeated);
    std::cout << "DEBUG: DoPermutationSort returned successfully" << std::endl;
    
    return result;
}

// Permutation方法的具体实现
// Cipher PermutationSorting::DoPermutationSort(const Cipher& in_exp, const Cipher& in_rep) {
//     PrintIfVerbose("Starting permutation-based sorting");
    
//     auto start = std::chrono::steady_clock::now();
    
//     // 步骤1：计算索引
//     PrintIfVerbose("Computing sorting indices...");
//     Cipher indexing = ComputeIndexing(in_exp, in_rep);
//     PrintIfVerbose("Indexing computation completed");
    
//     if (m_verbose) {
//         // print_duration(start, "Index computation");
//         // 打印一些索引信息
//         // m_hec.print(indexing, 16, "Indices (first 16 values): ");
//     }
    
//     // 步骤2：处理相等元素
//     if (m_enableTieoffset) {
//         start = std::chrono::steady_clock::now();
//         PrintIfVerbose("Computing tie-breaking offsets...");
        
//         Cipher offset = ComputeTieOffset(in_exp, in_rep);
//         indexing = m_hec.add(indexing, offset);
        
//         if (m_verbose) {
//             // print_duration(start, "Tie-breaking");
//             // m_hec.print(indexing, m_n * m_n, "Final indices: ");
//         }
//     } else {
//         PrintIfVerbose("Tie-offset disabled, skipping...");
//     }
    
//     // 步骤3：应用排列
//     start = std::chrono::steady_clock::now();
//     PrintIfVerbose("Applying permutation matrix...");
//     PrintIfVerbose("  Calling ComputeSorting with n=" + std::to_string(m_n));
    
//     Cipher sorted = ComputeSorting(indexing, in_rep);
//     PrintIfVerbose("ComputeSorting completed");
    
//     if (m_verbose) {
//         // print_duration(start, "Permutation application");
//         // m_hec.print(sorted, m_n, "Sorted result: ");
//     }
    
//     return sorted;
// }

Cipher PermutationSorting::DoPermutationSort(const Cipher& in_exp, const Cipher& in_rep) {
    std::cout << "DEBUG 1: Entering DoPermutationSort" << std::endl;
    
    // 步骤1：计算索引
    Cipher indexing = ComputeIndexing(in_exp, in_rep);
    std::cout << "DEBUG 2: ComputeIndexing returned" << std::endl;
    
    // 直接进入步骤3
    std::cout << "DEBUG 3: About to call ComputeSorting" << std::endl;
    Cipher sorted = ComputeSorting(indexing, in_rep);
    std::cout << "DEBUG 4: ComputeSorting returned" << std::endl;
    
    return sorted;
}

Cipher PermutationSorting::ComputeIndexing(const Cipher& in_exp, const Cipher& in_rep) {
    PrintIfVerbose("ComputeIndexing: Starting");
    PrintIfVerbose("  n = " + std::to_string(m_n));
    PrintIfVerbose("  sigmoid scaling = " + std::to_string(m_sigmoidScaling));
    PrintIfVerbose("  sigmoid degree = " + std::to_string(m_sigmoidDegree));
    // 计算差值
    Cipher difference = m_hec.sub(in_exp, in_rep);
    PrintIfVerbose("ComputeIndexing: Computed difference");
    
    // 应用sigmoid近似进行比较
    PrintIfVerbose("ComputeIndexing: Applying sigmoid...");
    Cipher cmp = m_hec.sigmoid(difference, m_n, m_sigmoidDegree, -m_sigmoidScaling);
    PrintIfVerbose("ComputeIndexing: Sigmoid completed");
    
    // 旋转求和计算每个元素的位置
    PrintIfVerbose("ComputeIndexing: Computing rotsum...");
    Cipher indexes = m_hec.rotsum(cmp, m_n);
    PrintIfVerbose("ComputeIndexing: Rotsum completed");

    // 减去偏移量
    PrintIfVerbose("ComputeIndexing: Encoding offset value...");
    double offset_value = 0.5 / m_n;
    PrintIfVerbose("  Offset value = " + std::to_string(offset_value));
    PrintIfVerbose("  Level = " + std::to_string(indexes->GetLevel()));
    PrintIfVerbose("  Slots = " + std::to_string(m_n * m_n));
    
    Plain offset = m_hec.encode(offset_value, indexes->GetLevel(), m_n * m_n);
    PrintIfVerbose("ComputeIndexing: Encoded offset");
    
    Cipher result = m_hec.sub(indexes, offset);
    PrintIfVerbose("ComputeIndexing: Subtracted offset");
    
    return result;
}

Cipher PermutationSorting::ComputeTieOffset(const Cipher& in_exp, const Cipher& in_rep) {
    // 根据精度选择合适的度数
    int d_tie = 247;
    if (m_delta >= 0.1) {
        d_tie = 247;
    } else if (m_delta >= 0.01) {
        d_tie = 495;
    } else if (m_delta >= 0.001) {
        d_tie = 1007;
    } else if (m_delta >= 0.0001) {
        d_tie = 4007;
    }
    
    // 使用sinc函数检测相等元素
    Cipher eq = m_hec.sinc(m_hec.sub(in_exp, in_rep), d_tie, 1.0 / m_delta);
    
    // 复制用于后续计算
    Plain pt_eqclone = m_hec.encode(0, eq->GetLevel(), m_n*m_n);
    Cipher eqclone = m_hec.add(eq, pt_eqclone);

    // 旋转累加计算ri
    for (int i = 0; i < std::log2(m_n); i++) {
        int rotindex = std::pow(2, i);
        eqclone = m_hec.add(eqclone, m_hec.rot(eqclone, m_n * rotindex));
    }
    
    Cipher sx = m_hec.mult(eqclone, 0.5 / m_n);
    
    // 创建三角矩阵
    std::vector<double> triangular_matrix;
    for (int rows = 0; rows < m_n; rows++) {
        for (int cols = m_n - rows; cols < m_n; cols++) {
            triangular_matrix.push_back(0);
        }
        for (int cols = 0; cols < m_n - rows; cols++) {
            triangular_matrix.push_back(1.0 / (double) m_n);
        }
    }
    
    Plain triang = m_hec.encode(triangular_matrix, eq->GetLevel(), m_n * m_n);
    Cipher dx = m_hec.mult(eq, triang);
    
    // 旋转累加计算ki
    for (int i = 0; i < std::log2(m_n); i++) {
        int rotindex = std::pow(2, i);
        dx = m_hec.add(dx, m_hec.rot(dx, m_n * rotindex));
    }
    
    Cipher offset = m_hec.sub(sx, dx);
    return m_hec.add(offset, 0.5 / m_n);
}

Cipher PermutationSorting::ComputeSorting(const Cipher& indexes, const Cipher& in_rep) {
    PrintIfVerbose("ComputeSorting: Starting");
    PrintIfVerbose("  n = " + std::to_string(m_n));
    PrintIfVerbose("  sinc degree = " + std::to_string(m_sincDegree));
    
    // 创建位置编码
    PrintIfVerbose("ComputeSorting: Creating position encoding...");
    std::vector<double> positions;
    positions.reserve(m_n * m_n);  // 预分配
    
    for (int i = 0; i < m_n; i++) {
        for (int j = 0; j < m_n; j++) {
            positions.push_back(i / (double)m_n);
        }
    }
    PrintIfVerbose("  Position vector size = " + std::to_string(positions.size()));
    
    // 编码位置向量
    PrintIfVerbose("ComputeSorting: Encoding positions...");
    Plain positions_plain = m_hec.encode(positions, indexes->GetLevel(), m_n * m_n);
    PrintIfVerbose("ComputeSorting: Positions encoded");
    
    // 计算排列增量
    PrintIfVerbose("ComputeSorting: Computing permutation delta...");
    Cipher permutation_delta = m_hec.sub(indexes, positions_plain);
    PrintIfVerbose("ComputeSorting: Delta computed");
    
    // 使用sinc函数创建排列矩阵
    PrintIfVerbose("ComputeSorting: Applying sinc function...");
    Cipher permutation_matrix = m_hec.sinc(permutation_delta, m_sincDegree, m_n);
    PrintIfVerbose("ComputeSorting: Sinc completed");
    
    // 应用排列矩阵
    PrintIfVerbose("ComputeSorting: Multiplying with input...");
    Cipher sorted = m_hec.mult(in_rep, permutation_matrix);
    PrintIfVerbose("ComputeSorting: Multiplication completed");
    
    // 旋转累加得到最终结果
    PrintIfVerbose("ComputeSorting: Final rotation and accumulation...");
    for (int i = 0; i < std::log2(m_n); i++) {
        int rotindex = std::pow(2, i);
        PrintIfVerbose("  Rotating by " + std::to_string(rotindex));
        sorted = m_hec.add(sorted, m_hec.rot(sorted, rotindex));
    }
    PrintIfVerbose("ComputeSorting: Completed");
    
    return sorted;

}

int PermutationSorting::EstimateDepthRequired(const SortingParams& params) const {
    // 与原始实现的set_permutation_parameters保持一致
    int partial_depth = 0;
    
    if (m_delta >= 0.1) {
        partial_depth = 10;
    } else if (m_delta >= 0.01) {
        partial_depth = 10;
    } else if (m_delta >= 0.001) {
        partial_depth = 14;
    } else if (m_delta >= 0.0001) {
        partial_depth = 15;
    }
    
    // 根据n添加sinc的深度
    if (m_n <= 8) {
        partial_depth += 6;
    } else if (m_n == 16) {
        partial_depth += 6;
    } else if (m_n == 32) {
        partial_depth += 7;
    } else if (m_n == 64) {
        partial_depth += 8;
    } else if (m_n == 128) {
        partial_depth += 9;
    }
    
    // 原始实现：+1 for the last matrix mult
    return partial_depth + 1;
}

void PermutationSorting::AutoConfigureParameters() {
    // 根据精度要求设置sigmoid参数
    if (m_delta >= 0.1) {
        m_sigmoidScaling = 650;
        m_sigmoidDegree = 1006;
    } else if (m_delta >= 0.01) {
        m_sigmoidScaling = 650;
        m_sigmoidDegree = 1006;
    } else if (m_delta >= 0.001) {
        m_sigmoidScaling = 9170;
        m_sigmoidDegree = 16000;
    } else if (m_delta >= 0.0001) {
        m_sigmoidScaling = 16000;
        m_sigmoidDegree = 32000;
    } else {
        throw std::runtime_error("Precision " + std::to_string(m_delta) + " is too small");
    }
    
    // 根据向量大小设置sinc度数
    if (m_n <= 8) {
        m_sincDegree = 59;
    } else if (m_n <= 16) {
        m_sincDegree = 59;
    } else if (m_n <= 32) {
        m_sincDegree = 119;
    } else if (m_n <= 64) {
        m_sincDegree = 247;
    } else if (m_n <= 128) {
        m_sincDegree = 495;
    } else {
        m_sincDegree = 1007;
    }
    
    PrintIfVerbose("Auto-configured parameters:");
    PrintIfVerbose("  Sigmoid scaling: " + std::to_string(m_sigmoidScaling));
    PrintIfVerbose("  Sigmoid degree: " + std::to_string(m_sigmoidDegree));
    PrintIfVerbose("  Sinc degree: " + std::to_string(m_sincDegree));
}

} // namespace lbcrypto
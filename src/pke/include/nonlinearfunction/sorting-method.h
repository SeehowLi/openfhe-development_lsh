/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 09:04:27
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-12 12:28:53
 * @FilePath: \openfhe-development\src\pke\include\nonlinearfunction\sorting-method.h
 * @Description: 
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

// 文件保护，为了防止重复包含
#ifndef SRC_NONLINEARFUNCTION_SORTING_METHOD_H_
#define SRC_NONLINEARFUNCTION_SORTING_METHOD_H_

#include "homoencrypt-compute.h"
#include <memory>           // ✅ 添加这个（用于std::shared_ptr）
#include <string>           // ✅ 添加这个（用于std::string）
#include <stdexcept>        // ✅ 添加这个（用于异常）

using namespace std;
using namespace std::chrono;

namespace lbcrypto{

// 参数
struct SortingParams{
    enum Method {
        PERMUTATION,
        NETWORK
    };
    
    // 通用参数
    Method method = PERMUTATION;
    int vectorSize = 0;
    double precision = 0.001;
    bool verbose = false;
    
    // Permutation特有参数
    int sigmoidScaling = 0;
    int sigmoidDegree = 0;
    int sincDegree = 0;
    bool enableTieoffset = false;
    
    // Network特有参数
    int reluDegree = 0;
    double inputScale = 1.0;

};

// 排序输入数据 - 支持不同算法的输入需求
struct SortingInput {
    std::shared_ptr<Cipher> standard;      // 标准输入（所有算法都使用）
    std::shared_ptr<Cipher> expanded;      // 扩展编码（Permutation使用）
    std::shared_ptr<Cipher> repeated;      // 重复编码（Permutation使用）
    
    // 单输入构造（Network使用）
    explicit SortingInput(const Cipher& input) 
        : standard(std::make_shared<Cipher>(input)), 
          expanded(nullptr), 
          repeated(nullptr) {}

    // 三输入构造（Permutation使用）
    SortingInput(const Cipher& std, const Cipher& exp, const Cipher& rep)
        : standard(std::make_shared<Cipher>(std)), 
          expanded(std::make_shared<Cipher>(exp)), 
          repeated(std::make_shared<Cipher>(rep)) {}

    // 检查是否有有效的扩展输入
    bool hasStandardInput() const { return standard != nullptr; }
    bool hasExpandedInputs() const { return expanded != nullptr && repeated != nullptr; }
};

// 排序方法基类
class SortingMethod {
protected:
    HomoEncryptCompute& m_hec;
    bool m_verbose;

public:
    explicit SortingMethod(HomoEncryptCompute& hec, bool verbose = false)
        : m_hec(hec), m_verbose(verbose) {}
    
    virtual ~SortingMethod() = default;

    // 通用接口 - 使用SortingInput支持不同的输入需求
    virtual Cipher Sort(const SortingInput& input, const SortingParams& params) = 0;
    
    // 便捷接口 - 自动处理输入准备
    Cipher Sort(const Cipher& input, const SortingParams& params) {
        return Sort(SortingInput(input), params);
    }

    // 为 PermutationSorting 提供的便利接口
    Cipher Sort(const Cipher& standard, const Cipher& expanded, const Cipher& repeated, const SortingParams& params) {
        return Sort(SortingInput(standard, expanded, repeated), params);
    }

    // 核心接口
    virtual std::string GetMethodName() const = 0;
    virtual int EstimateDepthRequired(const SortingParams& params) const = 0;

    // 工具方法
    void PrintIfVerbose(const std::string& message);

    // 工厂方法
    static std::unique_ptr<SortingMethod> Create(
        SortingParams::Method method, 
        HomoEncryptCompute& hec,
        const SortingParams& params);

};

}  // namespace lbcrypto

#endif  // SRC_NONLINEARFUNCTION_SORTING_METHOD_H_
// 该文件包含排序相关的非线性函数定义和实现
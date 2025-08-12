/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-11 19:59:21
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-08-02 13:10:12
 * @FilePath: \openfhe-development\src\pke\include\homoencrypt-compute.h
 * @Description: Unified homomorphic encryption computation library
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */
#ifndef SRC_LIB_HOMOENCRYPT_COMPUTE_H_
#define SRC_LIB_HOMOENCRYPT_COMPUTE_H_

#include "openfhe.h"
#include "ciphertext-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include <random>
#include <typeinfo>

using namespace lbcrypto;
using namespace std;
using namespace std::chrono;

// 自定义类型简化书写
using Plain = Plaintext;
using Cipher = Ciphertext<DCRTPoly>;

/**
 * @class HomoEncryptCompute
 * @brief 统一的同态加密计算库，封装了OpenFHE的CKKS方案功能
 * @author SeehowLi
 * 
 * 该类提供了完整的同态加密操作接口，包括：
 * - 加密上下文生成
 * - 密钥管理
 * - 基本同态运算
 * - 高级功能（如SlotsToCoeffs转换）
 * - 函数近似计算
 * - 各种工具函数
 * 
 * @note 本实现借鉴了以下团队的工作：
 * - OpenFHE开发团队的CKKS实现和示例代码
 * - [Jean-Sebastien Coron and Robin Kostler]关于SlotsToCoeffs提取的工作
 * - [Lorenzo]等人在函数近似方面的算法
 * 
 */
class HomoEncryptCompute {
private:
    CryptoContext<DCRTPoly> context;  // 加密上下文
    KeyPair<DCRTPoly> key_pair;       // 密钥对（公钥和私钥）
    
    /**
     * @brief 打印模数链信息（简化版）
     * @param poly DCRT多项式
     * 输出格式：log(QP): xxx
     */
    void print_moduli_chain(const DCRTPoly& poly);
    
    /**
     * @brief 获取FHE算法实例
     * @param cryptoContext 加密上下文
     * @return FHECKKSRNS算法的共享指针
     * @throw runtime_error 如果无法获取FHE算法
     */
    std::shared_ptr<FHECKKSRNS> getFHEAlgorithm(const CryptoContext<DCRTPoly>& cryptoContext);
    
    /**
     * @brief 将大整数转换为double类型
     * @tparam BigIntType 大整数类型
     * @param bigInt 输入的大整数
     * @return 转换后的double值
     */
    template<typename BigIntType>
    double bigIntegerToDouble(const BigIntType& bigInt) {
        std::stringstream ss;
        ss << bigInt;
        return std::stod(ss.str());
    }

public:
    /**
     * @brief 默认构造函数
     */
    HomoEncryptCompute(){};

    // ========== Context Generation Functions ==========
    
    /**
     * @brief 为KNN算法生成加密上下文
     * @param num_slots 密文中的槽数（数据打包数）
     * @param levels_required 所需的电路深度（乘法深度）
     * @param ring_dim 环维度（必须是2的幂）
     * @param toy_parameters true使用测试参数（不安全但快速），false使用128位安全参数
     * 
     * 该函数会自动生成所需的旋转密钥，特别优化了KNN算法的需求
     */
    void generate_context_knn(int num_slots, int levels_required, uint32_t ring_dim, bool toy_parameters);

    /**
     * @brief 生成用于测试的加密上下文（非安全参数）
     * @param num_slots 密文中的槽数
     * @param levels_required 所需的电路深度
     * @param ring_dim 环维度
     * 
     * 警告：仅用于测试和开发，不要在生产环境使用
     */
    void generate_context_toy(int num_slots = 1 << 3, int levels_required = 6, uint32_t ring_dim = 1 << 4);

    /**
     * @brief 生成用于测试的加密上下文（非安全参数）--CtS first BTS
     * @param num_slots 密文中的槽数
     * @param levels_required 所需的电路深度
     * @param ring_dim 环维度
     * 
     * 警告：仅用于测试和开发，不要在生产环境使用
     */
    void generate_context_ctsfirstbts_toy(int num_slots = 1 << 11, int levels_extra = 10, uint32_t ring_dim = 1 << 12);

    /**
     * @brief 生成用于测试的加密上下文（非安全参数）--StC first BTS
     * @param num_slots 密文中的槽数
     * @param levels_required 所需的电路深度
     * @param ring_dim 环维度
     * 
     * 警告：仅用于测试和开发，不要在生产环境使用
     */
    void generate_context_stcfirstbts_toy(int num_slots = 1 << 11, int levels_extra = 10, uint32_t ring_dim = 1 << 12);

    /**
     * @brief 生成具有128位安全强度的加密上下文
     * @param num_slots 密文中的槽数
     * @param levels_required 所需的电路深度
     * @param ring_dim 环维度（通常需要>=2^16以满足安全要求）
     * 
     * 用于生产环境的安全参数设置
     */
    void generate_context_128bit(int num_slots = 1 << 3, int levels_required = 6, uint32_t ring_dim = 1 << 16,
                                 int dcrtBits = 59, int firstMod = 60);
    
    /**
     * @brief 生成具有bts的加密上下文，用于测试
     * @param num_slots 密文中的槽数
     * @param levels_required 所需的电路深度
     * @param ring_dim 环维度（通常需要>=2^16以满足安全要求）
     * 
     * 用于生产环境的安全参数设置
     */
    void generate_context_bts(std::vector<uint32_t> levelBudget, int num_slots = 1 << 3, int levels_required = 6, uint32_t ring_dim = 1 << 16,\
                                int lStC = 1, int lCtS = 1,
                                 int dcrtBits = 59, int firstMod = 60);

    // ========== Key Generation Functions ==========
    
    /**
     * @brief 生成网络排序所需的旋转密钥
     * @param num_slots 槽数
     * 
     * 生成2的幂次的正负旋转密钥，用于实现高效的排序网络
     */
    void generate_rotation_keys_network(int num_slots);

    /**
     * @brief 生成单个旋转密钥
     * @param index 旋转索引（正数表示左旋，负数表示右旋）
     */
    void generate_rotation_key(int index);

    /**
     * @brief 批量生成多个旋转密钥
     * @param rotations 旋转索引向量
     * 
     * 一次性生成多个旋转密钥，比多次调用单个版本更高效
     */
    void generate_rotation_key(vector<int> rotations);

    // ========== Basic FHE Operations ==========
    
    /**
     * @brief 生成自举密钥
     * @param secretkey 私钥
     * @param num_slots 使用的槽数
     */
    void genbootkey(PrivateKey<DCRTPoly> secretkey, int num_slots);
    
     /**
     * @brief 获取私钥
     * @param secretkey 私钥
     */
    PrivateKey<DCRTPoly> getsecretkey();
    
    // --- Encode operations ---
    
    /**
     * @brief 编码double向量为明文
     * @param vec 要编码的数据向量
     * @param level 编码层级（影响精度）
     * @param num_slots 使用的槽数
     * @return 编码后的明文对象
     */
    Plain encode(const vector<double>& vec, int level, int num_slots);
    
    /**
     * @brief 编码double向量为明文（指定缩放度）
     * @param vec 要编码的数据向量
     * @param scale_deg 缩放度（2^scale_deg）
     * @param level 编码层级
     * @param num_slots 使用的槽数
     * @return 编码后的明文对象
     */
    Plain encode(const vector<double>& vec, size_t scale_deg, int level, int num_slots);
    
    /**
     * @brief 编码单个值为明文（复制到所有槽）
     * @param value 要编码的值
     * @param level 编码层级
     * @param num_slots 槽数
     * @return 编码后的明文对象，所有槽包含相同的值
     */
    Plain encode(double value, int level, int num_slots);

    /**
     * @brief 编码单个值为明文（全槽打包）
     * @param vec 要编码的值向量
     * @param level 编码层级
     * @return 编码后的明文对象，所有槽包含相同的值
     */
    Plain encode(const vector<double>& vec, int level);

    // --- Encrypt operations ---
    
    /**
     * @brief 加密明文
     * @param p 明文对象
     * @return 密文对象
     */
    Cipher encrypt(const Plain& p);
    
    /**
     * @brief 直接加密double向量
     * @param vec 要加密的数据向量
     * @param level 编码层级（默认0）
     * @param plaintext_num_slots 槽数（默认0表示自动）
     * @return 密文对象
     */
    Cipher encrypt(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0);
    
    /**
     * @brief 扩展编码加密（每个元素重复多次）
     * @param vec 原始数据向量
     * @param level 编码层级
     * @param plaintext_num_slots 总槽数
     * @param repetitions 每个元素的重复次数
     * @return 密文对象
     * 
     * 例如：[1,2] 重复2次 -> [1,1,2,2]
     */
    Cipher encrypt_expanded(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0, int repetitions = 1);
    
    /**
     * @brief 重复编码加密（整个向量重复多次）
     * @param vec 原始数据向量
     * @param level 编码层级
     * @param plaintext_num_slots 总槽数
     * @param repetitions 向量的重复次数
     * @return 密文对象
     * 
     * 例如：[1,2] 重复2次 -> [1,2,1,2]
     */
    Cipher encrypt_repeated(const vector<double>& vec, int level = 0, int plaintext_num_slots = 0, int repetitions = 1);

    // --- Decrypt and decode operations ---
    
    /**
     * @brief 解码明文为double向量
     * @param p 明文对象
     * @return 解码后的double向量
     */
    vector<double> decode(const Plain& p);
    
    /**
     * @brief 解密密文
     * @param c 密文对象
     * @return 解密后的明文对象
     */
    Plain decrypt(const Cipher& c);

    // --- Arithmetic operations ---
    
    /**
     * @brief 密文加法
     * @param c1 第一个密文
     * @param c2 第二个密文
     * @return 结果密文 (c1 + c2)
     */
    Cipher add(const Cipher& c1, const Cipher& c2);
    
    /**
     * @brief 密文与明文加法
     * @param c 密文
     * @param p 明文
     * @return 结果密文 (c + p)
     */
    Cipher add(const Cipher& c, Plain& p);
    
    /**
     * @brief 密文与常数加法
     * @param c 密文
     * @param d 常数（会复制到所有槽）
     * @return 结果密文 (c + d)
     */
    Cipher add(const Cipher& c, double d);
    
    /**
     * @brief 原地密文加法
     * @param c1 第一个密文（会被修改）
     * @param c2 第二个密文
     * 
     * 执行后 c1 = c1 + c2
     */
    void add_inplace(Cipher& c1, const Cipher& c2);
    
    /**
     * @brief 多个密文的树形加法
     * @param v 密文向量
     * @return 所有密文的和
     * 
     * 使用树形结构优化多个密文相加的性能
     */
    Cipher add_tree(vector<Cipher> v);
    
    /**
     * @brief 密文减法
     * @param c1 被减数密文
     * @param c2 减数密文
     * @return 结果密文 (c1 - c2)
     */
    Cipher sub(const Cipher& c1, const Cipher& c2);
    
    /**
     * @brief 密文与明文减法
     * @param c 被减数密文
     * @param p 减数明文
     * @return 结果密文 (c - p)
     */
    Cipher sub(const Cipher& c, Plain& p);
    
    /**
     * @brief 密文与明文乘法
     * @param c 密文
     * @param p 明文
     * @return 结果密文 (c * p)
     */
    Cipher mult(const Cipher& c, const Plain& p);
    
    /**
     * @brief 密文与常数乘法
     * @param c 密文
     * @param v 常数（会复制到所有槽）
     * @return 结果密文 (c * v)
     */
    Cipher mult(const Cipher& c, double v);
    
    /**
     * @brief 密文乘法
     * @param c1 第一个密文
     * @param c2 第二个密文
     * @return 结果密文 (c1 * c2)
     * 
     * 注意：会消耗一层乘法深度，增加噪声
     */
    Cipher mult(const Cipher& c1, const Cipher& c2);
    
    /**
     * @brief 密文平方
     * @param c1 密文
     * @return 结果密文 (c1^2)
     * 
     * 比自乘更高效
     */
    Cipher square(const Cipher& c1);

    // --- Rotation ---
    
    /**
     * @brief 旋转密文槽
     * @param c 密文
     * @param index 旋转索引（正数左旋，负数右旋）
     * @return 旋转后的密文
     * 
     * 需要预先生成对应的旋转密钥
     */
    Cipher rot(const Cipher& c, int index);
    
    /**
     * @brief 旋转求和
     * @param in 输入密文
     * @param n 求和距离
     * @return 结果密文
     * 
     * 将相距n的元素累加，用于实现高效的归约操作
     */
    Cipher rotsum(const Cipher& in, int n);

    /**
     * @brief 模约简/重缩放
     * @param ciphertext 输入密文
     * @param degree 重缩放的度数 
     * @return 结果密文
     * 
     * 将密文重缩放degree个level
     */
    Cipher modreduceinternal(Cipher& in, int degree);

    /**
     * @brief 模约简/重缩放
     * @param ciphertext 输入密文
     * @param degree 重缩放的度数 
     * @return 内部进行操作
     * 
     * 将密文重缩放degree个level
     */
    void modreduceinternalinplace(Cipher& in, int degree);

    // ========== Polynomial Operations ==========
    
    /**
     * @brief 将DCRT多项式转换为大模数多项式
     * @param poly DCRT多项式
     * @return 大模数多项式
     * 
     * 使用CRT插值重构单一大模数表示
     */
    Poly PolyFromDCRTPoly(const DCRTPoly& poly);
    
    /**
     * @brief 多项式乘以X^shift（原生多项式版本）
     * @param poly 输入多项式
     * @param shift 移位量
     * @return 移位后的多项式
     * 
     * 在模X^N+1的环中执行
     */
    NativePoly ShiftRight(const NativePoly& poly, uint32_t shift);
    
    /**
     * @brief 多项式乘以X^shift（DCRT多项式版本）
     * @param poly 输入DCRT多项式
     * @param shift 移位量
     * @return 移位后的DCRT多项式
     */
    DCRTPoly ShiftRight(const DCRTPoly& poly, uint32_t shift);
    
    // ========== Ciphertext Utilities ==========
    
    /**
     * @brief 密文乘以X^shift
     * @param ciphertext 输入密文
     * @param shift 移位量
     * @return 移位后的密文
     * 
     * 对密文的所有多项式分量执行移位
     */
    Cipher ShiftRight(const Cipher& ciphertext, uint32_t shift);
    
    /**
     * @brief 特殊解密函数（返回多项式形式）
     * @param c 密文
     * @param privateKey 私钥
     * @return 解密后的多项式
     * 
     * 用于调试和分析，不执行解码
     */
    Poly myDecrypt(const Cipher& c, const PrivateKey<DCRTPoly> privateKey);
    
    /**
     * @brief 特殊解密函数（使用内部私钥）
     * @param c 密文
     * @return 解密后的多项式
     */
    Poly myDecrypt(const Cipher& c);
    
    /**
     * @brief 从密文提取LWE形式
     * @param ciphertext 输入密文
     * @return LWE向量表示
     * 
     * 用于密文的第0个系数
     */
    BigVector LWEfromCiph(const Cipher& ciphertext);
    
    /**
     * @brief 解密LWE密文（单个值）
     * @param ciphertext 密文
     * @param privateKey 私钥
     * @return 解密后的大整数
     */
    BigInteger DecryptLWE(const Cipher& ciphertext, const PrivateKey<DCRTPoly>& privateKey);
    
    /**
     * @brief 解密LWE密文（使用内部私钥）
     * @param ciphertext 密文
     * @return 解密后的大整数
     */
    BigInteger DecryptLWE(const Cipher& ciphertext);
    
    /**
     * @brief 解密LWE密文（多个值）
     * @param ciphertext 密文
     * @param privateKey 私钥
     * @param n 要解密的值的数量
     * @return 解密后的大整数向量
     */
    BigVector DecryptLWE(const Cipher& ciphertext, const PrivateKey<DCRTPoly>& privateKey, int n);
    
    /**
     * @brief 解密LWE密文（多个值，使用内部私钥）
     * @param ciphertext 密文
     * @param n 要解密的值的数量
     * @return 解密后的大整数向量
     */
    BigVector DecryptLWE(const Cipher& ciphertext, int n);
    
    /**
     * @brief 计算密文的共轭
     * @param ciphertext 输入密文
     * @return 共轭后的密文
     * 
     * 对复数槽执行共轭操作
     */
    Cipher CiphertextConjugate(const Cipher& ciphertext);


    // ========== Bootstrapping Operations ==========
    // 将大模数转换为double类型
    double GetBigModulus(const std::shared_ptr<lbcrypto::CryptoParametersCKKSRNS> cryptoParams);
    // 二倍角公式迭代
    void ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ciphertext, uint32_t numIter);
    // 通过调整密文的缩放因子来降低噪声--待定  
    void AdjustCiphertext(Ciphertext<DCRTPoly>& ciphertext, double correction);
    
    void ExtendCiphertext(std::vector<DCRTPoly>& ctxtDCRT, const CryptoContextImpl<DCRTPoly>& cc,
                          const std::shared_ptr<DCRTPoly::Params> elementParamsRaisedPtr);
    Cipher conjugate(const Cipher& ciphertext);
                          
    /**
     * @brief 自举操作（刷新密文），原始的自举操作--CtS first
     * @param c 密文
     * @return 刷新后的密文
     * 
     * 降低噪声，恢复可用的乘法深度
     */
    Cipher bootstrap(Cipher& c);

    /**
     * @brief 自举前除了预计算之外的一些其他设置
     * @param no
     * @todo 后续可能删除或者与其他的函数合并
     * @return 设置的context
     * 
     * context的一些设置
     */
    void bootothersetup();

    /**
     * @brief 部分和，用于稀疏打包，有效槽位数据变为N/n倍
     * @param c 密文
     * @return 部分和后的密文
     * 
     * 稀疏打包的对齐操作
     */
    Cipher modraise(Cipher& c, const uint32_t& correction, const uint32_t& compositeDegree = 1);

    /**
     * @brief 部分和，用于稀疏打包，有效槽位数据变为N/n倍
     * @param c 密文
     * @return 部分和后的密文
     * 
     * 稀疏打包的对齐操作
     */
    Cipher partialsum(Cipher& c, const size_t& N, const size_t& n);

    /**
     * @brief 自举操作（刷新密文），原始的自举操作--StC first
     * @param c 密文
     * @return 刷新后的密文
     * 
     * 降低噪声，恢复可用的乘法深度
     */
    Cipher bootstrapCtSfirst(Cipher&ciphertext, uint32_t numIterations = 1 , uint32_t precision = 0);

    /**
     * @brief 自举操作（刷新密文），原始的自举操作--StC first
     * @param c 密文
     * @return 刷新后的密文
     * 
     * 降低噪声，恢复可用的乘法深度
     */
    Cipher bootstrapStCfirst(Cipher&ciphertext, uint32_t numIterations = 1, uint32_t precision = 0);


    // ========== SlotsToCoeffs Operations ==========
    
    /**
     * @brief 设置槽到系数转换
     * @param levelBudget 层级预算（转换消耗的层数）
     * @param numSlots 槽数
     * @param lDec 转换后剩余的层数
     * 
     * 必须在执行SlotsToCoeffs前调用
     */
    void EvalSlotsToCoeffsSetup(uint32_t levelBudget, uint32_t numSlots, uint32_t lDec);

    /**
     * @brief 设置双向转换（StC和CtS）
     * @param levelBudget 层级预算
     * @param numSlots 槽数
     * @param lStC StC后剩余的层数
     * @param lCtS CtS后剩余的层数
     * 
     * 同时设置两个方向的转换
     */
    void EvalStC_CtSSetup(std::vector<uint32_t> levelBudget, uint32_t numSlots, uint32_t lStC, uint32_t lCtS);

    /**
     * @brief 设置双向转换（StC和CtS）--这个是一步到位的，不是全定制
     * @param levelBudget 层级预算
     * @param dim1 维度（默认{0,0}表示自动）
     * @param numSlots 槽数
     * @param correctionFactor 修正因子（默认0）
     * @param isStCFirst 是否StC-first BTS（默认false表示CtS-first BTS）
     * @param lDec StC后剩余的层数
     * @param lEnc CtS后剩余的层数
     * 
     * 同时设置两个方向的转换
     */
    void EvalBTSSetup(vector<uint32_t> levelBudget, std::vector<uint32_t> dim1 = {0,0}, 
                      uint32_t numSlots = 0, uint32_t correctionFactor = 0, bool isStCFirst = false);

    /**
     * @brief 槽到系数转换
     * @param ciph 输入密文（槽表示）
     * @return 输出密文（系数表示）
     * 
     * 将CKKS的槽表示转换为多项式系数表示
     */
    Cipher SlotsToCoeffs(const Cipher& ciph);

    /**
     * @brief 系数到槽转换
     * @param ciph 输入密文（系数表示）
     * @return 输出密文（槽表示）
     * 
     * 将多项式系数表示转换回CKKS的槽表示
     */
    Cipher CoeffsToSlots(const Cipher& ciph);

    /**
     * @brief 线性变换--CtS
     * @param ciph 输入密文（系数表示）
     * @return 输出线性变换之后的密文
     * 
     * 矩阵乘法
     */
    Cipher LinearTransformCtS(const Cipher& ciph);

    /**
     * @brief 线性变换--StC
     * @param ciph 输入密文（系数表示）
     * @return 输出线性变换之后的密文
     * 
     * 矩阵乘法
     */
    Cipher LinearTransformStC(const Cipher& ciph);

    // ========== Function Approximations ==========
    
    /**
     * @brief 切比雪夫多项式近似
     * @param ciphertext 输入密文
     * @param coeffcient 切比雪夫多项式系数
     * @param a 区间下界
     * @param b 区间上界
     * @return 函数值的密文
     * 
     * 使用切比雪夫多项式在[a,b]区间近似func
     */
    Cipher chebyshevseries(const Cipher& ciphertext, const std::vector<double>& coeffcient,
                     double a, double b);

    /**
     * @brief 切比雪夫多项式近似
     * @param func 要近似的函数
     * @param ciphertext 输入密文
     * @param a 区间下界
     * @param b 区间上界
     * @param degree 多项式次数
     * @return 函数值的密文
     * 
     * 使用切比雪夫多项式在[a,b]区间近似func
     */
    Cipher chebyshev(std::function<double(double)> func,
                     const Cipher& ciphertext, double a,
                     double b, uint32_t degree);

    /**
     * @brief Sigmoid函数近似
     * @param in 输入密文
     * @param n 缩放参数（影响输出范围）
     * @param degree 多项式次数
     * @param scaling 输入缩放因子
     * @return sigmoid(scaling*x)/n 的密文
     */
    Cipher sigmoid(const Cipher& in, int n, int degree, int scaling);
    
    /**
     * @brief 紧凑Sigmoid函数近似
     * @param in 输入密文
     * @param n 缩放参数
     * @param degree 多项式次数
     * @param scaling 输入缩放因子
     * @return 1 - sigmoid(scaling*x)/n 的密文
     * 
     * 在更紧凑的区间[-0.15, 1.05]内近似
     */
    Cipher sigmoid_tight(const Cipher& in, int n, int degree, int scaling);
    
    /**
     * @brief Sigmoid和Tanh的复合函数近似
     * @param in 输入密文
     * @param degree 多项式次数
     * @param scaling Sigmoid的缩放因子
     * @param scaling_tanh Tanh的缩放因子
     * @return sigmoid(tanh(x)) 的密文
     * 
     * 用于实现更陡峭的激活函数
     */
    Cipher sigmoid_tanh(const Cipher& in, int degree, int scaling, int scaling_tanh);

    /**
     * @brief Sinc函数近似
     * @param in 输入密文
     * @param degree 多项式次数
     * @param n 频率参数
     * @return sin(πnx)/(πnx) 的密文
     */
    Cipher sinc(const Cipher& in, int degree, int n);
    
    /**
     * @brief ReLU函数近似
     * @param in 输入密文
     * @param degree 多项式次数
     * @param n 未使用参数（保留接口一致性）
     * @return max(0, x) 的密文
     */
    Cipher relu(const Cipher& in, int degree, int n);

    // ========== Plaintext Creation Utilities ==========
    
    /**
     * @brief 创建常数明文（所有槽相同值）
     * @param val 常数值
     * @param level 编码层级
     * @return 明文对象
     */
    Plain MakePlaintext(double val, uint32_t level);
    
    /**
     * @brief 从复数向量创建明文
     * @param v 复数向量
     * @param level 编码层级
     * @return 明文对象
     * 
     * 向量会被重复以填满所有槽
     */
    Plain MakePlaintext(std::vector<std::complex<double>> v, uint32_t level);
    
    /**
     * @brief 从实数向量创建明文
     * @param v 实数向量
     * @param level 编码层级
     * @return 明文对象
     * 
     * 实数会被转换为复数（虚部为0）
     */
    Plain MakePlaintext(std::vector<double> v, uint32_t level);

    // ========== Decrypt Utilities ==========
    
    /**
     * @brief 解密并获取CKKS打包值（使用内部私钥）
     * @param ciphertext 密文
     * @param numSlots 要提取的槽数
     * @return 复数向量
     */
    std::vector<std::complex<double>> DecryptCKKSPackedValue(const Cipher& ciphertext, uint32_t numSlots);
    
    /**
     * @brief 解密并获取CKKS打包值（指定私钥）
     * @param ciphertext 密文
     * @param privateKey 私钥
     * @param numSlots 要提取的槽数
     * @return 复数向量
     */
    std::vector<std::complex<double>> DecryptCKKSPackedValue(const Cipher& ciphertext, 
                                                             const PrivateKey<DCRTPoly>& privateKey, 
                                                             uint32_t numSlots);

    // ========== Utility Functions ==========
    
    /**
     * @brief 生成均匀分布的随机实数
     * @param n 向量长度
     * @return 随机数向量，范围[0.5, 1.0]
     */
    std::vector<double> genUniformReal(uint32_t n);
    
    /**
     * @brief 估计两个向量之间的精度
     * @param v1 第一个向量（通常是计算结果）
     * @param v2 第二个向量（通常是期望值）
     * @return 平均精度（以比特为单位）
     * 
     * 计算相对误差并转换为精度比特数
     */
    double estimatePrecision(std::vector<std::complex<double>>& v1, 
                           std::vector<std::complex<double>>& v2);

    // --- Print functions ---
    
    /**
     * @brief 打印密文内容
     * @param c 密文
     * @param slots 要打印的槽数（0表示全部）
     * @param prefix 打印前缀
     * 
     * 解密并格式化输出密文内容
     */
    void print(const Cipher& c, int slots = 0, string prefix = "");
    
    /**
     * @brief 以模数形式打印大整数
     * @param x 大整数
     * @param q 模数
     * 
     * 如果x > q/2，打印为负数形式
     */
    void printIntegerMod(const BigInteger& x, const BigInteger& q);
    
    /**
     * @brief 打印大整数向量（模数形式）
     * @param x 大整数向量
     * 
     * 格式：[x1 x2 ... xn]
     */
    void printBigVectorMod(const BigVector& x);
    
    /**
     * @brief 打印多项式（模数形式）
     * @param p 多项式
     * 
     * 打印所有系数
     */
    void printPoly(const Poly& p);
    
    /**
     * @brief 打印详细的模数链信息
     * @param poly DCRT多项式
     * 
     * 输出每个模数及其比特长度
     */
    void print_moduli_chain_detail(const DCRTPoly& poly);

    // ========== Test Functions ==========
    
    /**
     * @brief 测试SlotsToCoeffs功能（稀疏打包）
     * 
     * 使用小参数测试基本的槽到系数转换
     */
    void testSlots2Coeffs();
    
    /**
     * @brief 测试SlotsToCoeffs功能（全打包）
     * 
     * 使用全槽打包测试转换
     */
    void mytestSlots2Coeffs();
    
    /**
     * @brief 测试StC后接CtS（稀疏打包）
     * 
     * 验证双向转换的正确性
     */
    void testStC_CtS();
    
    /**
     * @brief 测试CtS后接StC（稀疏打包）
     * 
     * 验证反向转换的正确性
     */
    void testCtS_StC();
    
    /**
     * @brief 测试StC后接CtS（全打包）
     */
    void mytestStC_CtS();
    
    /**
     * @brief 测试CtS后接StC（全打包）
     */
    void mytestCtS_StC();
    
    /**
     * @brief 测试实数情况下的StC和CtS
     * 
     * 使用实数而非复数进行测试
     */
    void testStC_CtSReal();
    
    // ========== Getters ==========
    
    /**
     * @brief 获取加密上下文
     * @return 当前的加密上下文
     */
    CryptoContext<DCRTPoly> getContext() const { return context; }
    
    /**
     * @brief 获取密钥对
     * @return 当前的密钥对（公钥和私钥）
     */
    KeyPair<DCRTPoly> getKeyPair() const { return key_pair; }

};

#endif  // SRC_LIB_HOMOENCRYPT_COMPUTE_H_
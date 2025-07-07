# OpenFHE 修改记录

## 概述
本文档记录所有对OpenFHE库进行的自定义修改、功能添加和实验性特性。

**作者**: [李世豪] 
**贡献者：**[李世豪 etc.] 
**开始日期**: 2025-01-07  
**基础版本**: OpenFHE [v1.2.4]  
**开发环境**: MinGW64 / Windows 11/Ubuntu 22.04

---

## 修改汇总表

| 日期 | 功能/模块 | 类型 | 状态 | 相关文件 | 备注 |
|------|----------|------|------|----------|------|
| 2025-01-07 | 自定义日志系统 | 功能增强 | 进行中 | `src/core/lib/utils/logger.cpp` | 用于调试 |
| 2025-06-07 | softmax的“归一化-平方” | 论文复现 | 进行中 | `src/pke/examples/softmax_eval.cpp`  | baseline |
| 2025-07-07 | 离散CKKS的模约简       | 论文复现 | 进行中 | `src/pke/examples/discrete_ckks.cpp` | baseline |
|  |  |  |  |  |  |

---

## 详细修改记录

### 2025-01-07: 自定义日志系统
**类型**: 功能增强  
**状态**: 进行中  
**优先级**: 中  
**开发耗时**: 3小时  

#### 功能描述
添加了一个自定义日志系统，用于追踪库的操作过程和调试信息，特别是在进行同态加密操作时能够记录详细的性能数据。

#### 修改的文件
- `src/core/include/utils/logger.h` - 新增日志类定义
- `src/core/lib/utils/logger.cpp` - 日志功能实现
- `CMakeLists.txt` - 添加日志编译选项

#### 核心代码实现
```cpp
// 添加的主要类和函数：
class CustomLogger {
private:
    std::string moduleName;
    std::ofstream logFile;
    
public:
    // 记录操作日志
    void logOperation(const std::string& operation, LogLevel level);
    // 记录性能数据
    void logPerformance(const std::string& operation, double time_ms);
    // 记录加密参数
    void logCryptoParams(const CryptoParameters& params);
};
```

#### 使用示例
```cpp
// 在BGV加密中使用
CustomLogger logger("BGV加密模块");
auto start = std::chrono::high_resolution_clock::now();

// 执行加密操作
ciphertext = crypto_context->Encrypt(publicKey, plaintext);

auto end = std::chrono::high_resolution_clock::now();
double duration = std::chrono::duration<double, std::milli>(end - start).count();
logger.logPerformance("BGV加密", duration);
```

#### 测试情况
- [x] 单元测试：`tests/unit/test_logger.cpp`
- [x] BGV方案集成测试
- [ ] BFV方案集成测试
- [ ] CKKS方案集成测试
- [x] 性能影响测试：额外开销 <0.1%

#### 已知问题
1. **多线程安全性**：需要在并发环境下验证互斥锁的有效性
2. **日志文件大小**：未实现日志轮转功能，长时间运行可能产生巨大日志文件

#### 后续改进计划
- [ ] 实现日志文件自动轮转（超过100MB时创建新文件）
- [ ] 添加日志级别过滤功能
- [ ] 支持输出到不同目标（文件/控制台/网络）
- [ ] 添加日志格式自定义功能

#### 相关设计决策
- 选择单例模式确保全局只有一个日志实例
- 使用C++17的filesystem库处理文件路径
- 日志格式采用：`[时间戳] [级别] [模块] 消息`

---

### 2025-06-06: softmax “平方-归一化”实现
**类型**: 性能优化  
**状态**: 设计阶段  
**优先级**: 高  

#### 功能描述
开发一个自动参数选择器，根据安全级别和性能需求自动选择最优的BGV参数...

---

## 编译配置修改

### 新增的编译选项
```cmake
# 添加到 CMakeLists.txt
option(ENABLE_CUSTOM_LOGGING "启用自定义日志功能" ON)
option(ENABLE_PERFORMANCE_TRACKING "启用性能追踪" ON)

if(ENABLE_CUSTOM_LOGGING)
    add_definitions(-DCUSTOM_LOGGING_ENABLED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOG_LEVEL=2")
endif()
```

### 编译命令
```bash
# 启用所有自定义功能
cmake .. -DENABLE_CUSTOM_LOGGING=ON -DENABLE_PERFORMANCE_TRACKING=ON

# 发布版本（禁用日志）
cmake .. -DENABLE_CUSTOM_LOGGING=OFF -DCMAKE_BUILD_TYPE=Release
```

---

## 性能影响评估

| 功能 | 内存开销 | CPU开销 | 说明 |
|------|---------|---------|------|
| 日志系统 | ~1MB | <0.1% | 仅在启用时产生开销 |
| 参数优化器 | ~5MB | 初始化时2-3秒 | 缓存优化结果 |

---

## 与上游版本的兼容性

### 合并策略
1. 所有修改保持在独立分支：`feature/custom-modifications`
2. 使用条件编译确保可以完全禁用自定义功能
3. 保持API向后兼容

### 版本兼容性测试
- [x] OpenFHE v1.2.x - 完全兼容
- [ ] OpenFHE v1.4.0 - 待测试

---

## 参考资料
- [OpenFHE官方文档](https://openfhe.org/docs)
- [同态加密性能优化论文](链接)
- 内部设计文档：`docs/开发文档/架构设计.md`
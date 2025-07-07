#!/bin/bash
###
 # @Author: SeehowLi lsh0126@nudt.edu.cn
 # @Date: 2025-07-07 16:43:57
 # @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 # @LastEditTime: 2025-07-07 16:44:23
 # @FilePath: \openfhe-development\lsh\smart_collect.sh
 # @Description: 
 # 
 # Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
### 

# 创建输出目录
mkdir -p code_chunks

echo "Starting code collection for OpenFHE project..."

# 1. 项目概览和摘要
echo "Creating project summary..."
> code_chunks/0_summary.txt
echo "====== OPENFHE PROJECT SUMMARY ======" >> code_chunks/0_summary.txt
echo -e "\nGenerated on: $(date)" >> code_chunks/0_summary.txt
echo -e "\n=== Directory Structure ===" >> code_chunks/0_summary.txt
find . -type d -name ".git" -prune -o -type d -name "build" -prune -o -type d -print | head -50 | sed 's|[^/]*/|- |g' >> code_chunks/0_summary.txt

echo -e "\n\n=== File Statistics ===" >> code_chunks/0_summary.txt
echo "Header files (.h, .hpp):" >> code_chunks/0_summary.txt
find . -name "*.h" -o -name "*.hpp" | grep -v "build\|\.git" | wc -l >> code_chunks/0_summary.txt
echo "Source files (.cpp):" >> code_chunks/0_summary.txt
find . -name "*.cpp" | grep -v "build\|\.git" | wc -l >> code_chunks/0_summary.txt

echo -e "\n\n=== Core Header Files ===" >> code_chunks/0_summary.txt
find . \( -name "*.h" -o -name "*.hpp" \) -path "*/src/*" -o -path "*/include/*" | grep -v "build\|\.git" | sort >> code_chunks/0_summary.txt

# 2. 核心头文件（类定义和接口）
echo "Collecting headers and interfaces..."
> code_chunks/1_headers_and_interfaces.txt
echo "====== CORE HEADERS AND INTERFACES ======" >> code_chunks/1_headers_and_interfaces.txt
echo "Generated on: $(date)" >> code_chunks/1_headers_and_interfaces.txt

find . \( -name "*.h" -o -name "*.hpp" \) | grep -v "build\|\.git\|test" | sort | while read -r file; do
    echo -e "\n\n//========================================" >> code_chunks/1_headers_and_interfaces.txt
    echo "//===== File: $file =====" >> code_chunks/1_headers_and_interfaces.txt
    echo "//========================================" >> code_chunks/1_headers_and_interfaces.txt
    
    # 提取关键内容：类定义、重要函数声明等
    cat "$file" | grep -E -A3 -B1 "^[[:space:]]*(template|class|struct|enum|typedef|namespace|public:|protected:|private:|virtual|static|explicit)" | grep -v "^--$" >> code_chunks/1_headers_and_interfaces.txt
done

# 3. 实现文件摘要
echo "Collecting implementation signatures..."
> code_chunks/2_implementations.txt
echo "====== KEY IMPLEMENTATIONS ======" >> code_chunks/2_implementations.txt
echo "Generated on: $(date)" >> code_chunks/2_implementations.txt

find . -name "*.cpp" | grep -v "build\|\.git\|test\|example" | sort | while read -r file; do
    # 检查文件是否有实质内容
    if [ -s "$file" ]; then
        echo -e "\n\n//===== File: $file =====" >> code_chunks/2_implementations.txt
        # 提取函数实现签名和重要的类方法
        grep -E "^[[:space:]]*[a-zA-Z_][a-zA-Z0-9_:<>]*::[a-zA-Z_][a-zA-Z0-9_]*\s*\(" "$file" | head -20 >> code_chunks/2_implementations.txt
    fi
done

# 4. CMake和构建配置
echo "Collecting build configurations..."
> code_chunks/3_build_config.txt
echo "====== BUILD CONFIGURATION ======" >> code_chunks/3_build_config.txt
echo "Generated on: $(date)" >> code_chunks/3_build_config.txt

# 主CMakeLists.txt
if [ -f "CMakeLists.txt" ]; then
    echo -e "\n//===== Main CMakeLists.txt =====" >> code_chunks/3_build_config.txt
    cat CMakeLists.txt >> code_chunks/3_build_config.txt
fi

# 其他CMake文件
find . -name "*.cmake" -o -name "CMakeLists.txt" | grep -v "build\|\.git" | grep -v "^./CMakeLists.txt$" | while read -r file; do
    echo -e "\n\n//===== $file =====" >> code_chunks/3_build_config.txt
    cat "$file" >> code_chunks/3_build_config.txt
done

# 5. 创建关键算法文件列表（OpenFHE特定）
echo "Identifying key OpenFHE components..."
> code_chunks/4_key_algorithms.txt
echo "====== KEY OPENFHE ALGORITHMS ======" >> code_chunks/4_key_algorithms.txt

# 查找关键的加密算法实现
for keyword in "BGV" "BFV" "CKKS" "TFHE" "scheme" "crypto" "lattice" "encrypt" "decrypt" "key"; do
    echo -e "\n=== Files containing '$keyword' ===" >> code_chunks/4_key_algorithms.txt
    find . \( -name "*.h" -o -name "*.hpp" -o -name "*.cpp" \) | xargs grep -l "$keyword" 2>/dev/null | grep -v "build\|\.git" | head -10 >> code_chunks/4_key_algorithms.txt
done

# 生成文件大小报告
echo -e "\n\n=== Generated file sizes ===" 
ls -lh code_chunks/

echo -e "\nCode collection complete! Files are in ./code_chunks/"
echo "You can now upload these files one by one, starting with 0_summary.txt"
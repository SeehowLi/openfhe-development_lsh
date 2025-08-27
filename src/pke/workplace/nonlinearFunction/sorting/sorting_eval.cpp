#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "nonlinearfunction/sorting-method.h"
#include "nonlinearfunction/sorting-utils.h"

using namespace lbcrypto;
using namespace std;
using namespace std::chrono;

int main() {
    // 这里是主函数的入口
    // 由于代码片段中没有提供完整的实现细节，这里仅作为占位符
    cout << "This is a placeholder for the main function." << endl;
    return 0;
}

// // 全局变量
// HomoEncryptCompute hec;
// vector<double> input_values;
// int n;
// double delta = 0.001;
// int precision_digits = 3;
// bool toy = false;
// SortingType sortingType = NONE;

// // 函数声明
// void read_arguments(int argc, char *argv[]);
// void evaluate_sorting_accuracy(const Cipher& result, int circuit_depth);
// void print_usage();

// int main(int argc, char *argv[]) {
//     // 读取命令行参数
//     read_arguments(argc, argv);
    
//     if (sortingType == NONE) {
//         cerr << "You must pick a sorting method. Add either --permutation or --network" << endl;
//         print_usage();
//         return 1;
//     }
    
//     cout << "Selected sorting type: " << to_string(sortingType) << endl;
    
//     // 设置精度显示位数
//     if (delta >= 0.1) precision_digits = 1;
//     else if (delta >= 0.01) precision_digits = 2;
//     else if (delta >= 0.001) precision_digits = 3;
//     else if (delta >= 0.0001) precision_digits = 4;
    
//     cout << setprecision(precision_digits) << fixed;
    
//     // 调试输出
//     cout << "\n=== Input Parameters ===" << endl;
//     cout << "Vector size (n): " << n << endl;
//     cout << "Precision (delta): " << delta << endl;
//     cout << "Toy parameters: " << (toy ? "Yes" : "No") << endl;
//     cout << "Input values: ";
//     for (const auto& val : input_values) {
//         cout << val << " ";
//     }
//     cout << endl;
    
//     auto start_time = steady_clock::now();
    
//     Cipher result;
//     int circuit_depth = 0;
    
//     // 创建排序参数
//     SortingParams params;
//     params.vectorSize = n;
//     params.precision = delta;
//     params.verbose = toy;  // 使用toy参数控制verbose
    
//     if (sortingType == PERMUTATION) {
//         cout << "\n=== Permutation-based Sorting ===" << endl;

//         // 先设置参数
//         int sigmoid_scaling = 0;
//         int degree_sigmoid = 0;
//         int degree_sinc = 0;
//         int partial_depth = 0;
        
//         if (delta >= 0.1) {
//             precision_digits = 1;
//             sigmoid_scaling = 650;
//             degree_sigmoid = 1006;
//             partial_depth = 10;
//         } else if (delta >= 0.01) {
//             precision_digits = 2;
//             sigmoid_scaling = 650;
//             degree_sigmoid = 1006;
//             partial_depth = 10;
//         } else if (delta >= 0.001) {
//             precision_digits = 3;
//             sigmoid_scaling = 9170;
//             degree_sigmoid = 16000;
//             partial_depth = 14;
//         } else if (delta >= 0.0001) {
//             precision_digits = 4;
//             sigmoid_scaling = 16000;
//             degree_sigmoid = 32000;
//             partial_depth = 15;
//             degree_sinc = 495;
//         } else {
//             cerr << "The required min distance '" << delta << "' is too small!" << endl;
//             return 1;
//         }
        
//         if (n <= 8) {
//             degree_sinc = 59;
//             partial_depth += 6;
//         } else if (n == 16) {
//             degree_sinc = 59;
//             partial_depth += 6;
//         } else if (n == 32) {
//             degree_sinc = 119;
//             partial_depth += 7;
//         } else if (n == 64) {
//             degree_sinc = 247;
//             partial_depth += 8;
//         } else if (n == 128) {
//             degree_sinc = 495;
//             partial_depth += 9;
//         }
        
//         circuit_depth = partial_depth + 1;
        
//         cout << "Circuit depth: " << circuit_depth << endl;
//         cout << setprecision(precision_digits) << fixed;
//         cout << endl << "Ciphertext: " << endl << input_values << endl << endl << "δ: " << delta << ", ";
        
//         // 生成加密上下文
//         hec.generate_context_permutation(n * n, circuit_depth, toy);
        
//         // 生成旋转密钥
//         for (int i = 0; i < log2(n); i++) {
//             hec.generate_rotation_key(pow(2, i) * n);
//             hec.generate_rotation_key(pow(2, i));
//         }
        
//         // 加密数据
//         Cipher in_exp = hec.encrypt_expanded(input_values, 0, n*n, n);
//         Cipher in_rep = hec.encrypt_repeated(input_values, 0, n*n, n);
        
//         // 创建排序参数
//         SortingParams params;
//         params.method = SortingParams::PERMUTATION;
//         params.vectorSize = n;
//         params.precision = delta;
//         params.verbose = toy;
//         params.sigmoidScaling = sigmoid_scaling;
//         params.sigmoidDegree = degree_sigmoid;
//         params.sincDegree = degree_sinc;
//         params.enableTieoffset = false;
        
//         // 检查--tieoffset参数
//         for (int i = 1; i < argc; i++) {
//             if (string(argv[i]) == "--tieoffset") {
//                 params.enableTieoffset = true;
//                 break;
//             }
//         }
//         // 创建排序器并执行排序
//         auto sorter = SortingMethod::Create(params.method, hec, params);
//         std::cout << "DEBUG: Sorter created" << std::endl;

//         std::cout << "DEBUG: About to call Sort" << std::endl;
//         result = sorter->Sort(SortingInput(hec.encrypt(input_values), in_exp, in_rep), params);
//         std::cout << "DEBUG: Sort returned successfully" << std::endl;
//     } else if (sortingType == NETWORK) {
//         cout << "\n=== Network-based Sorting ===" << endl;
        
//         params.method = SortingParams::NETWORK;
//         params.inputScale = 0.95;
//         params.reluDegree = 0;  // 让NetworkSorting自动配置
        
//         // 创建临时排序器来估算深度
//         auto temp_sorter = SortingMethod::Create(params.method, hec, params);
//         int levels_consumption = temp_sorter->EstimateDepthRequired(params);
        
//         cout << "Levels consumption estimated: " << levels_consumption << endl;
        
//         // 生成加密上下文
//         circuit_depth = hec.generate_context_network(n, levels_consumption, toy);
//         cout << "Total circuit depth (with bootstrapping): " << circuit_depth << endl;
        
//         // 生成旋转密钥
//         cout << "Generating rotation keys..." << endl;
//         hec.generate_rotation_keys_network(n);
        
//         // 缩放输入值
//         cout << "Scaling input values by " << params.inputScale << "..." << endl;
//         for (auto& val : input_values) {
//             val *= params.inputScale;
//         }
        
//         // 加密数据 - 注意level设置
//         int initial_level = circuit_depth - levels_consumption - 3;
//         cout << "Encrypting data at level " << initial_level << "..." << endl;
//         Cipher in = hec.encrypt(input_values, initial_level, n);
        
//         // 创建排序器并执行排序
//         cout << "Starting sorting..." << endl;
//         auto sorter = SortingMethod::Create(params.method, hec, params);
//         result = sorter->Sort(in, params);
//     }
    
//     cout << "\n=== Sorting Completed ===" << endl;
//     std::cout << "DEBUG: About to print duration" << std::endl;
//     print_duration(start_time, "The sorting took:");
//     std::cout << "DEBUG: Duration printed" << std::endl;

//     // 评估结果
//     std::cout << "DEBUG: About to call evaluate_sorting_accuracy" << std::endl;
//     evaluate_sorting_accuracy(result, circuit_depth);
//     std::cout << "DEBUG: evaluate_sorting_accuracy completed" << std::endl;

//     return 0;
// }

// void evaluate_sorting_accuracy(const Cipher& result, int circuit_depth) {
//     std::cout << "DEBUG: Entering evaluate_sorting_accuracy" << std::endl;
    
//     cout << "Level final: " << result->GetLevel() << "/" << circuit_depth << endl << endl;
    
//     std::cout << "DEBUG: About to decrypt result" << std::endl;
//     vector<double> sorted_fhe = hec.decode(hec.decrypt(result));
//     std::cout << "DEBUG: Decryption successful, size = " << sorted_fhe.size() << std::endl;
    
//     vector<double> results_fhe;
    
//     if (sortingType == PERMUTATION) {
//         // 原始实现：从n*n向量中提取结果
//         for (int i = 0; i < n * n; i += n) {
//             results_fhe.push_back(sorted_fhe[i]);  // 注意：原始实现没有除以input_scale
//         }
//     } else if (sortingType == NETWORK) {
//         // 原始实现：直接取前n个元素
//         for (int i = 0; i < n; i += 1) {
//             results_fhe.push_back(sorted_fhe[i]);  // 结果已经是缩放过的
//         }
//     }
    
//     // 对原始输入排序（注意：对于Network，input_values已经被缩放了）
//     sort(input_values.begin(), input_values.end());
    
//     cout << endl << "Expected:  " << input_values << endl;
//     cout << endl << "Obtained:  " << results_fhe << endl << endl;
    
//     int corrects = 0;
    
//     for (int i = 0; i < n; i++) {
//         if (abs(input_values[i] - results_fhe[i]) < delta) corrects++;
//     }
//     cout << "Corrects (up to " << delta << "): " << GREEN_TEXT << corrects << RESET_COLOR 
//          << "/" << GREEN_TEXT << n << RESET_COLOR << endl;
    
//     cout << "Precision bits: " << GREEN_TEXT << precision_bits(input_values, results_fhe) << RESET_COLOR << endl;
// }

// void read_arguments(int argc, char *argv[]) {
//     // 默认值
//     n = 64;
//     bool random_elements = true;
    
//     // 解析参数
//     for (int i = 1; i < argc; i++) {
//         string arg = argv[i];
        
//         if (arg == "--permutation") {
//             sortingType = PERMUTATION;
//         }
//         else if (arg == "--network") {
//             sortingType = NETWORK;
//         }
//         else if (arg == "--toy") {
//             toy = true;
//         }
//         else if (arg == "--delta" && i + 1 < argc) {
//             delta = stod(argv[++i]);
//         }
//         else if (arg == "--random" && i + 1 < argc) {
//             n = stoi(argv[++i]);
//             if ((n & (n - 1)) != 0) {
//                 cerr << "The number of values must be a power of two" << endl;
//                 exit(1);
//             }
//         }
//         else if (arg == "--file" && i + 1 < argc) {
//             ifstream f(argv[++i]);
//             if (f) {
//                 stringstream ss;
//                 ss << f.rdbuf();
//                 input_values = parse_input_vector("[ " + ss.str() + " ]");
//                 n = input_values.size();
//                 random_elements = false;
                
//                 // 计算最小差值
//                 vector<double> sorted_vals = input_values;
//                 sort(sorted_vals.begin(), sorted_vals.end());
//                 double min_diff = 1.0;
//                 for (size_t j = 1; j < sorted_vals.size(); j++) {
//                     double diff = sorted_vals[j] - sorted_vals[j-1];
//                     if (diff > 0 && diff < min_diff) {
//                         min_diff = diff;
//                     }
//                 }
//                 delta = min_diff;
//             } else {
//                 cout << "Could not find file" << endl;
//                 exit(1);
//             }
//         }
//         else if (arg == "--vector" && i + 1 < argc) {
//             input_values = parse_input_vector(argv[++i]);
//             n = input_values.size();
//             random_elements = false;
//         }
//         else if (arg.front() == '[' && arg.back() == ']') {
//             input_values = parse_input_vector(arg);
//             n = input_values.size();
//             random_elements = false;
//         }
//     }
    
//     // 生成随机数据（如果需要）
//     if (random_elements || input_values.empty()) {
//         input_values = generate_close_randoms(n, delta);
//     }
// }

// void print_usage() {
//     cout << "\nUsage: sorting_eval [options]" << endl;
//     cout << "\nRequired (choose one):" << endl;
//     cout << "  --permutation        Use permutation-based sorting" << endl;
//     cout << "  --network           Use network-based sorting" << endl;
//     cout << "\nOptional:" << endl;
//     cout << "  --random <n>        Generate n random values (must be power of 2)" << endl;
//     cout << "  --vector <v>        Input vector, e.g., \"[3,1,4,1,5,9,2,6]\"" << endl;
//     cout << "  --file <f>          Read input from file" << endl;
//     cout << "  --delta <d>         Precision (default: 0.001)" << endl;
//     cout << "  --toy               Use toy security parameters" << endl;
//     cout << "\nExamples:" << endl;
//     cout << "  ./sorting_eval --network --random 32" << endl;
//     cout << "  ./sorting_eval --permutation --vector \"[5,3,8,1,9,2,7,4]\"" << endl;
// }
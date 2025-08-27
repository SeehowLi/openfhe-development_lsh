#include "openfhe.h"
#include "homoencrypt-compute.h"
#include "math/chebyshev.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>

// int main(int argc, char* argv[]) {
//     // 默认值
//     int k = 1;
//     int degree = 45;
    
//     // 解析命令行参数
//     if (argc >= 2) {
//         k = std::atoi(argv[1]);
//     }
//     if (argc >= 3) {
//         degree = std::atoi(argv[2]);
//     }
    
//     // 打印参数信息
//     std::cout << "Approximating function: sin(2π·" << k << "·x)/(2π)" << std::endl;
//     std::cout << "Using Chebyshev polynomial of degree: " << degree << std::endl;
//     std::cout << "Interval: [-1, 1]" << std::endl;
//     std::cout << std::string(60, '-') << std::endl;
    
//     // 计算切比雪夫系数
//     // 注意：lambda表达式需要捕获k的值
//     auto coefficients = EvalChebyshevCoefficients(
//         [k](double x) -> double {
//             if (k == 0) {
//                 // 特殊处理 k=0 的情况：sin(0)/2π = 0
//                 return 0.0;
//             }
//             return sin(2 * M_PI * k * x) / (2 * M_PI);
//         }, 
//         -1.0, 1.0, degree
//     );
    
//     // 输出系数
//     std::cout << "Chebyshev coefficients:" << std::endl;
//     std::cout << std::scientific << std::setprecision(15);
    
//     for (size_t i = 0; i < coefficients.size(); i++) {
//         std::cout << "coefficients[" << std::setw(2) << i << "] = " 
//                   << std::setw(22) << coefficients[i] << std::endl;
//     }
    
//     // 输出一些统计信息
//     std::cout << std::string(60, '-') << std::endl;
//     std::cout << "Total number of coefficients: " << coefficients.size() << std::endl;
    
//     // 找出最大系数
//     double max_coeff = 0.0;
//     size_t max_index = 0;
//     for (size_t i = 0; i < coefficients.size(); i++) {
//         if (std::abs(coefficients[i]) > std::abs(max_coeff)) {
//             max_coeff = coefficients[i];
//             max_index = i;
//         }
//     }
//     std::cout << "Largest coefficient (by absolute value): coefficients[" 
//               << max_index << "] = " << max_coeff << std::endl;
    
//     // 验证一些特殊点的值（可选）
//     std::cout << std::string(60, '-') << std::endl;
//     std::cout << "Verification at special points:" << std::endl;
    
//     // 在x=0处的值（使用切比雪夫多项式的性质）
//     double value_at_0 = 0.0;
//     for (size_t n = 0; n < coefficients.size(); n++) {
//         if (n % 2 == 0) {
//             value_at_0 += coefficients[n] * ((n/2) % 2 == 0 ? 1.0 : -1.0);
//         }
//     }
//     std::cout << "Approximation at x=0: " << value_at_0 << std::endl;
//     std::cout << "Exact value at x=0: " << 0.0 << std::endl;
    
//     return 0;
// }

int main() {
    int scaling = 1260;
    // int scaling_tanh = 10;
    int degree = 450;
    // int n = 1;

    // 定义sigmoid函数
    auto sigmoid = [scaling](double x) -> double {
        return 1.0 / (1.0 + std::exp(-scaling * (x / 90.0)));
    };

    // // 定义tanh函数
    // auto tanh_func = [scaling_tanh](double x) -> double {
    //     double exp_pos = std::exp(scaling_tanh * x);
    //     double exp_neg = std::exp(-scaling_tanh * x);
    //     return (exp_pos - exp_neg) / (exp_pos + exp_neg);
    // };

    // auto func = [sigmoid, tanh_func](double x) -> double {
    //     return sigmoid(tanh_func(x));
    // };
    // double a = -1.0;
    // double b =  1.0;
    double a = -0.15 * 90;
    double b =  1.05 * 90;
    
    // 计算切比雪夫系数
    auto coefficients = EvalChebyshevCoefficients(sigmoid, a, b, degree);
    
    // 输出系数为可直接复制的格式
    std::cout << "coefficients3{" << std::endl;
    std::cout << std::scientific << std::setprecision(17);
    
    // 每行输出4个系数
    for (size_t i = 0; i < coefficients.size(); i++) {
        if (i % 4 == 0 && i > 0) {
            std::cout << std::endl;  // 每4个系数换行
        }
        
        // 输出系数，保持对齐
        std::cout << std::setw(25) << coefficients[i];
        
        // 添加逗号（最后一个元素除外）
        if (i < coefficients.size() - 1) {
            std::cout << ",";
        }
        
        // 每行最后一个元素后添加适当的空格
        if (i % 4 < 3 && i < coefficients.size() - 1) {
            std::cout << " ";
        }
    }
    
    std::cout << "};" << std::endl;

    return 0;
}

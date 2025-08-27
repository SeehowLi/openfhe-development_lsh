#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <omp.h>
#include <chrono>

#include "homoencrypt-compute.h"

using namespace lbcrypto;

// 全局变量
HomoEncryptCompute hec;
vector<double> input_values;
bool toy = false;

// 全局变量来存储数据
vector<double> query_data;
vector<vector<double>> training_data;
// 平展成一维
vector<double> query_data_exp32;  //拓展的query_data
vector<double> training_data_1d;  // 前32个

// 时间测量辅助函数
class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::string task_name;
    
public:
    Timer(const std::string& name) : task_name(name) {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    ~Timer() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        cout << "[TIMING] " << task_name << ": " << duration.count() << " ms" << std::endl;
    }
};

// 辅助函数：提取方括号内的数字
vector<double> extract_numbers(const string& str) {
    vector<double> numbers;
    stringstream ss(str);
    string token;
    
    while (getline(ss, token, ',')) {
        // 移除空格和方括号
        token.erase(remove_if(token.begin(), token.end(), 
                    [](char c) { return c == '[' || c == ']' || c == ' '; }), 
                    token.end());
        
        if (!token.empty()) {
            try {
                numbers.push_back(stod(token));
            } catch (...) {
                // 忽略无法转换的值
            }
        }
    }
    
    return numbers;
}

// 归一化函数
void normalize_query_and_training() {
    // 现在假设数据分布的最大值是10.0
    double max_val = 10.0;
    double max_distance = sqrt(10.0) * 2 * max_val; // 计算最大距离，用于归一化
    // OpenMP并行归一化
    for (double& v : query_data) {
        v /= max_distance;
    }
    // 归一化training_data
    for (size_t i = 0; i < training_data.size(); i++) {
        for (double& v : training_data[i]) {
            v /= max_distance;
        }
    }
}

// 读取jsonl数据
void read_data_json(string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        // 查找 "query" 字段
        size_t query_pos = line.find("\"query\":");
        if (query_pos != string::npos) {
            size_t start = line.find('[', query_pos);
            size_t end = line.find(']', start);
            if (start != string::npos && end != string::npos) {
                string query_str = line.substr(start, end - start + 1);
                query_data = extract_numbers(query_str);
            }
        }
        
        // 查找 "data" 字段
        size_t data_pos = line.find("\"data\":");
        if (data_pos != string::npos) {
            training_data.clear();
            size_t start = line.find('[', data_pos);
            size_t pos = start + 1;
            
            while (pos < line.length()) {
                size_t row_start = line.find('[', pos);
                if (row_start == string::npos) break;
                
                size_t row_end = line.find(']', row_start);
                if (row_end == string::npos) break;
                
                string row_str = line.substr(row_start, row_end - row_start + 1);
                vector<double> row_data = extract_numbers(row_str);
                
                if (!row_data.empty()) {
                    training_data.push_back(row_data);
                }
                
                pos = row_end + 1;
            }
        }
    }

    file.close();
    
}

// ========================================================================================== //
// 将data拓展成128维度的向量
void expand_groups_dimension_inplace() {
    for (auto& vec : training_data) {
        vec.resize(128, 0.0);
    }
}

// 拓展query成128维度的向量
void expand_query_dimension_inplace() {
    query_data.resize(128, 0.0);  // resize会在后面补0
}

// 将group数据转换为一维向量的函数
vector<double> convert_group_to_1d(const vector<vector<double>>& group_data, int elements_per_row = 128) {
    vector<double> result;
    result.reserve(group_data.size() * elements_per_row);
    
    for (const auto& row : group_data) {
        // 取每行的前elements_per_row个元素（或所有元素如果不足）
        int count = min(elements_per_row, (int)row.size());
        for (int i = 0; i < count; i++) {
            result.push_back(row[i]);
        }
        
        // 如果不足elements_per_row个元素，用0填充
        for (int i = count; i < elements_per_row; i++) {
            result.push_back(0.0);
        }
    }
    
    return result;
}

// 处理所有group数据的函数
void process_groups_to_1d() {
    training_data_1d = convert_group_to_1d(training_data, 128);
}

// 将query_data扩展重复128次的函数
vector<double> expand_query_data(int repeat_times = 128) {
    vector<double> expanded_query;
    expanded_query.reserve(query_data.size() * repeat_times);
    
    for (int i = 0; i < repeat_times; i++) {
        for (double value : query_data) {
            expanded_query.push_back(value);
        }
    }
    
    return expanded_query;
}

// 创建距离掩码函数
void create_distance_mask(vector<double>& distance_mask, int num_groups = 128, int elements_per_group = 128) {
    int total_length = num_groups * elements_per_group;
    distance_mask.clear();
    distance_mask.resize(total_length, 0.0);
    
    // 每32个数字的第一个设置为1
    for (int i = 0; i < num_groups; i++) {
        distance_mask[i * elements_per_group] = 1.0;
    }
    
}

// 内存检测
void printPeakMemoryUsage() {
    std::ifstream file("/proc/self/status");
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.find("VmPeak:") == 0) {
            std::cout << "\n[MEMORY] Peak Virtual Memory: " << line.substr(7) << std::endl;
        } else if (line.find("VmHWM:") == 0) {
            std::cout << "[MEMORY] Peak Physical Memory (RSS): " << line.substr(6) << std::endl;
        }
    }
    file.close();
}

int main(int argc, char* argv[]) {
    auto total_start = chrono::high_resolution_clock::now();
    // 解析命令行参数
    string dataset_file = "/public1/home/m8s001097/openfhe-development/lsh/train.jsonl";
    // string dataset_file = "./data/data.jsonl";
    string predictions_file = "predictions.jsonl";
    
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--dataset" && i + 1 < argc) {
            dataset_file = argv[++i];
        } else if (arg == "--predictions" && i + 1 < argc) {
            predictions_file = argv[++i];
        }
    }

    // 设置OpenMP线程数
    int num_threads = omp_get_max_threads();
    num_threads = 64;
    omp_set_num_threads(num_threads);
    
    read_data_json(dataset_file);
    
    /**
     * @brief 数据预处理
     */
    {
        // 扩展数据维度
        expand_groups_dimension_inplace();
        // 扩展query维度
        expand_query_dimension_inplace();
        // 归一化查询和训练数据
        normalize_query_and_training();
        // 处理输入向量成一维->128*128
        process_groups_to_1d();
        
        // 并行重复拓展query_data
        query_data_exp32 = expand_query_data(128);
    }
        

    // 加密参数设置
    int num_slot_dis1 = 128;// 前三组group中有效距离的个数
    int num_slot = num_slot_dis1 * num_slot_dis1; // 最终的槽数16384
    // =============== 重要参数 =============== //
    int levels_required = 21;
    uint32_t ring_dim = 1 << 16; // 65536
    int sigmod_degree1 = 200; // Sigmoid的近似多项式的阶数
    int sigmod_scale1 = 900; // Sigmoid的近似多项式的缩放因子
    int sigmod_degree2 = 247; // Sigmoid的近似多项式的阶数
    int sigmod_scale2 = 900; // Sigmoid的近似多项式的缩放因子
    vector<int> rotations_distance;
    // 距离转换的掩码
    vector<double> distance_mask1;
    vector<double> distance_mask2;

    // 创建距离掩码--每128个数字的第一个为1
    create_distance_mask(distance_mask1, num_slot_dis1, num_slot_dis1);

    // ============================ 针对128槽的加密上下文 ============================ //
    {
        hec.generate_context_knn(num_slot, levels_required, ring_dim, toy);
    }
    
    /**
    *@brief 初始数据编码成明文向量
    */
    Plain pt_data ,ptxt_query32;
    pt_data = hec.encode(training_data_1d, 1, 0, num_slot);
    ptxt_query32 =  hec.encode(query_data_exp32, 1, 0, num_slot);


    /**
    *@brief 初始数据加密
    */
    Cipher ctdata, ctqry32;
    ctdata = hec.encrypt(pt_data);
    ctqry32 = hec.encrypt(ptxt_query32);


    /**
    *@brief 欧式距离计算
    */
    // 减法做差
    Cipher ctminus;
    ctminus = hec.sub(ctdata, ctqry32);

    // 平方
    Cipher ctminus_square;
    ctminus_square = hec.square(ctminus);

    // 新的处理方式
    auto ctdistance_square_before = ctminus_square->Clone();
    // auto ctdistance_square_after = ctminus_square->Clone();
    // TODO 可以并行
    for(int i = 0; i < log2(num_slot_dis1); i++) {
        hec.add_inplace(ctdistance_square_before, hec.rot(ctdistance_square_before, pow(2, i)));
        // hec.add_inplace(ctdistance_square_after,  hec.rot(ctdistance_square_after, -1*pow(2, i)));
    }

    vector<double> distance_square_bemask(num_slot, 0.0);
    // vector<double> distance_square_afmask(num_slot, 0.0);
    for (int i = 0; i < num_slot_dis1; i++)
    {
        // 每一组的第一个为1.0
        distance_square_bemask[i * num_slot_dis1] = 1.0; 
        // 每一组的最后一个为1.0
        // distance_square_afmask[i * num_slot_dis1 + num_slot_dis1 - 1] = 1.0;

    }
    auto pt_distance_square_bemask = hec.encode(distance_square_bemask, 1, ctdistance_square_before->GetLevel(), num_slot);
    // auto pt_distance_square_afmask = hec.encode(distance_square_afmask, 1, ctdistance_square_after->GetLevel(), num_slot);
    auto ctdistance_square_bemask = hec.mult(ctdistance_square_before, pt_distance_square_bemask);
    // auto ctdistance_square_afmask = hec.mult(ctdistance_square_after, pt_distance_square_afmask);
    for(int i = 0; i < log2(num_slot_dis1); i++) {
        hec.add_inplace(ctdistance_square_bemask, hec.rot(ctdistance_square_bemask, -1 * pow(2, i)));
        // hec.add_inplace(ctdistance_square_afmask,  hec.rot(ctdistance_square_afmask, pow(2, i)));
    }
    // auto ctdistance_rep = hec.add(ctdistance_square_bemask, ctdistance_square_afmask);
    auto ctdistance_rep = ctdistance_square_bemask->Clone();

    vector<double> rep2exp_mask(num_slot, 0.0);
    // 掩码的生成，每一组128个数，第x组的第x个数是1，其余是0，只看前100组有效数字
    for (int i = 0; i < 100; i++) {
        rep2exp_mask[i * num_slot_dis1 + i] = 1.0;
    }

    // rep的mask，后28组全部是0
    vector<double> rep_mask(num_slot_dis1 * num_slot_dis1, 0.0);

    for (int i = 0; i < num_slot_dis1 * num_slot_dis1; i++) {
        if (i / num_slot_dis1 < 100) {  // 判断是第几组
            rep_mask[i] = 1.0;
        }
    }

    /**
     * @brief 计算repeat和expand
     * @todo 可以进行并行优化
     */
    Cipher ctdistance_exp;

    auto pt_rep2exp_mask = hec.encode(rep2exp_mask, 1, ctdistance_rep->GetLevel(), num_slot);
    ctdistance_exp = hec.mult(ctdistance_rep, pt_rep2exp_mask)->Clone();

    auto pt_rep_mask = hec.encode(rep_mask, 1, ctdistance_rep->GetLevel(), num_slot);
    ctdistance_rep = hec.mult(ctdistance_rep, pt_rep_mask)->Clone();
                                   
    /**
     * @brief 计算expand
     */
    // 利用repeat生成expand--不需要算第四组的
    // 旋转求和
    for (int i = 0; i < log2(num_slot_dis1); i++) {
        hec.add_inplace(ctdistance_exp, hec.rot(ctdistance_exp, num_slot_dis1 * pow(2, i)));
    }

    double distance_bound = 0.5;
    vector<double> const06_exp(num_slot_dis1 * num_slot_dis1, distance_bound);
    vector<double> const06_rep(num_slot_dis1 * num_slot_dis1, distance_bound);
    // 将每组的前100个数据设为0
    for (int group = 0; group < num_slot_dis1; group++) {
        for (int i = 0; i < 100; i++) {
            const06_exp[group * num_slot_dis1 + i] = 0.0;
        }
    }
    // const06_rep: 前100组设为0，后28组保持0.6
    for (int group = 0; group < 100; group++) {
        for (int i = 0; i < num_slot_dis1; i++) {
            const06_rep[group * num_slot_dis1 + i] = 0.0;
        }
    }
    auto pt_const06_exp = hec.encode(const06_exp, 1, ctdistance_exp->GetLevel(), num_slot);
    auto pt_const06_rep = hec.encode(const06_rep, 1, ctdistance_rep->GetLevel(), num_slot);
    ctdistance_exp = hec.add(ctdistance_exp, pt_const06_exp);
    ctdistance_rep = hec.add(ctdistance_rep, pt_const06_rep);

    // rep和exp做差得到delta,利用sigmod将delta计算成mask,以此来讲distance转换成index
    auto ctdistance_del = hec.sub(ctdistance_exp, ctdistance_rep);
    
    // 合并计算sigmod(tanh(10*x))
    auto ct_distance_mask = hec.sigmoid_tanh(ctdistance_del, sigmod_degree1, sigmod_scale1, 10);

    auto ctindex = ct_distance_mask->Clone();
    
    // 旋转求和，得到整体数据的index
    for(int i = 0; i < log2(num_slot_dis1); i++) {
        hec.add_inplace(ctindex, hec.rot(ctindex, num_slot_dis1 * pow(2, i)));
    }

    auto ctdistance2index_mask = hec.add(ctindex, -11.0);

    vector<double> backminus30(num_slot_dis1 * num_slot_dis1, 0.0); // 先全部初始化为0
    // 填充每组的后28个元素为-30
    for (int i = 0; i < num_slot_dis1; i++) {
        for (int j = 100; j < 128; j++) {
            backminus30[i * 128 + j] = -30.0;
        }
    }
    auto pt_backminus30 = hec.encode(backminus30, ctdistance2index_mask->GetLevel(), num_slot);

    ctdistance2index_mask = hec.add(ctdistance2index_mask, pt_backminus30);
    // 变成0和1
    // TODO 可能会有0.5，这里先不考虑
    auto ctindex_exp_mask = hec.sigmoid_tight(ctdistance2index_mask, 1, sigmod_degree2, sigmod_scale2);

    // 原始数据的排列索引，从1~128的重复
    vector<double> ctoriginindex(num_slot_dis1 * num_slot_dis1, 0.0);
    for (int group = 0; group < num_slot_dis1; ++group) {
        for (int i = 0; i < num_slot_dis1; ++i) {
            ctoriginindex[group * num_slot_dis1 + i] = static_cast<double>(i + 1);
        }
    }
    auto pt_originindex = hec.encode(ctoriginindex, 1, ctindex_exp_mask->GetLevel(), num_slot);

    auto ctindex_exp = hec.mult(ctindex_exp_mask, pt_originindex);

    // 修改输出部分，写入文件而不是标准输出
    ofstream output_file(predictions_file);
    if (!output_file.is_open()) {
        cerr << "Error: Cannot open output file: " << predictions_file << endl;
        return 1;
    }

    // 解密结果
    Plain result;
    result = hec.decrypt(ctindex_exp);

    double index_threshold = 0.98; // 阈值,小于此值的判断为噪声
    vector<int> result_index;
    int j = 0;
    
    for(size_t i=0; i < result->GetLength(); i++) {
        if(result->GetCKKSPackedValue()[i].real() > index_threshold){
            result_index.push_back(static_cast<int>(round(result->GetCKKSPackedValue()[i].real())));
        }
    }
    // 后处理：检查并清理result_index--清理异常值
    if(result_index.size() > 10) {
        // 使用迭代器删除非递增元素
        auto it = result_index.begin();
        while(it != result_index.end() && result_index.size() > 10) {
            // 检查当前元素是否比前一个元素小或相等（非严格递增）
            if(it != result_index.begin() && *it <= *(it - 1)) {
                // 发现异常值，删除当前元素
                it = result_index.erase(it);
            } else {
                // 当前元素正常，继续下一个
                ++it;
            }
        }
        result_index.resize(10); // 确保结果最多10个
    }


    output_file << "{ \"answer\": [ ";
    for(size_t i=0; i < result_index.size(); i++) {
        output_file << result_index[i];
        j++;
        if(j < 10){
            output_file << ", ";
        }
        if(j >= 10) break;  // 只输出前10个结果
    }
    output_file << " ] }" << endl;
    output_file.close();

    auto total_end = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::milliseconds>(total_end - total_start);
    std::cout << "\n[TIMING] Total execution time: " << total_duration.count()/1000.0 << " s" << std::endl;

    // 添加这一行来打印峰值内存使用
    printPeakMemoryUsage();
    return 0;
}
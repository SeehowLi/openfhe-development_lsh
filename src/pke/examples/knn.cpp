#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <string>

#include "nonlinearfunction/sorting-method.h"
#include "nonlinearfunction/sorting-utils.h"

using namespace lbcrypto;
using namespace std;
using namespace std::chrono;

// 全局变量
HomoEncryptCompute hec;
vector<double> input_values;
int n = 32;
double delta = 0.001;
int precision_digits = 3;
bool toy = false;
SortingType sortingType = NONE;

// 全局变量来存储数据
vector<double> query_data;
vector<vector<double>> training_data;
vector<vector<double>> group1_data;  // 前32个
vector<vector<double>> group2_data;  // 中间32个
vector<vector<double>> group3_data;  // 后32个
vector<vector<double>> group4_data;  // 最后4个

vector<double> query_data_exp32;  //拓展的query_data
vector<double> query_data_exp4;  //拓展的query_data
vector<double> group1_1d;  // 前32个
vector<double> group2_1d;  // 中间32个
vector<double> group3_1d;  // 后32个
vector<double> group4_1d;  // 最后4个

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
    std::cout << "Max distance for normalization: " << max_distance << std::endl;

    // // 找到query_data中的最大值
    // for (double v : query_data) {
    //     if (abs(v) > max_val) max_val = abs(v);
    // }
    // // 找到training_data中的最大值
    // for (const auto& row : training_data) {
    //     for (double v : row) {
    //         if (abs(v) > max_training_val) max_training_val = abs(v);
    //     }
    // }
    // if (max_val == 0.0 && max_training_val == 0.0) return; // 防止除零

    // 归一化query_data
    for (double& v : query_data) {
        v /= max_distance;
    }
    // 归一化training_data
    for (auto& row : training_data) {
        for (double& v : row) {
            v /= max_distance;
        }
    }
    // 归一化group1_data
    for (auto& row : group1_data) {
        for (double& v : row) {
            v /= max_distance;
        }
    }
    // 归一化group2_data
    for (auto& row : group2_data) {
        for (double& v : row) {
            v /= max_distance;
        }
    }
    // 归一化group3_data
    for (auto& row : group3_data) {
        for (double& v : row) {
            v /= max_distance;
        }
    }
    // 归一化group4_data
    for (auto& row : group4_data) {
        for (double& v : row) {
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
    
    // 将数据分成三组
    if (training_data.size() >= 100) {
        // 清空之前的分组数据
        group1_data.clear();
        group2_data.clear();
        group3_data.clear();

        // 第一组：前32个 (索引 0-31)
        for (int i = 0; i < 32; i++) {
            group1_data.push_back(training_data[i]);
        }

        // 第二组：中间32个 (索引 32-63)
        for (int i = 32; i < 64; i++) {
            group2_data.push_back(training_data[i]);
        }

        // 第三组：后36个 (索引 64-99)
        for (int i = 64; i < 100; i++) {
            group3_data.push_back(training_data[i]);
        }
    }
    
    // 输出读取结果的统计信息
    cout << "Successfully read data from " << filename << endl;
    cout << "Query size: " << query_data.size() << endl;
    cout << "Total data size: " << training_data.size() << " x " 
         << (training_data.empty() ? 0 : training_data[0].size()) << endl;
    cout << "Group 1 size: " << group1_data.size() << endl;
    cout << "Group 2 size: " << group2_data.size() << endl;
    cout << "Group 3 size: " << group3_data.size() << endl;
}

// ================== 聚类算法，使得数据更有分布特征 ================== //

// 计算两个向量之间的欧氏距离
double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// 计算向量加法
std::vector<double> vector_add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// 计算向量除以标量
std::vector<double> vector_divide(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] / scalar;
    }
    return result;
}

// 简化的K-means实现
void kmeans_clustering(const std::vector<std::vector<double>>& data, 
                      int k, 
                      std::vector<int>& assignments,
                      std::vector<std::vector<double>>& centers) {
    int n = data.size();
    int dim = data[0].size();
    
    // 初始化随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, n - 1);
    
    // 随机初始化聚类中心
    centers.clear();
    centers.resize(k);
    std::vector<bool> selected(n, false);
    
    // K-means++初始化
    int first_center = dis(gen);
    centers[0] = data[first_center];
    selected[first_center] = true;
    
    for (int i = 1; i < k; ++i) {
        std::vector<double> min_distances(n, std::numeric_limits<double>::max());
        
        // 计算每个点到最近中心的距离
        for (int j = 0; j < n; ++j) {
            if (!selected[j]) {
                for (int c = 0; c < i; ++c) {
                    double dist = euclidean_distance(data[j], centers[c]);
                    min_distances[j] = std::min(min_distances[j], dist);
                }
            }
        }
        
        // 根据距离概率选择下一个中心
        double sum_dist = 0.0;
        for (double d : min_distances) {
            if (d < std::numeric_limits<double>::max()) {
                sum_dist += d * d;
            }
        }
        
        std::uniform_real_distribution<> dis_real(0.0, sum_dist);
        double threshold = dis_real(gen);
        double cumulative = 0.0;
        
        for (int j = 0; j < n; ++j) {
            if (!selected[j] && min_distances[j] < std::numeric_limits<double>::max()) {
                cumulative += min_distances[j] * min_distances[j];
                if (cumulative >= threshold) {
                    centers[i] = data[j];
                    selected[j] = true;
                    break;
                }
            }
        }
    }
    
    // K-means迭代
    assignments.resize(n);
    int max_iterations = 100;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        bool changed = false;
        
        // 分配步骤
        for (int i = 0; i < n; ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int best_cluster = 0;
            
            for (int j = 0; j < k; ++j) {
                double dist = euclidean_distance(data[i], centers[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_cluster = j;
                }
            }
            
            if (assignments[i] != best_cluster) {
                changed = true;
                assignments[i] = best_cluster;
            }
        }
        
        if (!changed) break;
        
        // 更新步骤
        for (int j = 0; j < k; ++j) {
            std::vector<double> sum(dim, 0.0);
            int count = 0;
            
            for (int i = 0; i < n; ++i) {
                if (assignments[i] == j) {
                    sum = vector_add(sum, data[i]);
                    count++;
                }
            }
            
            if (count > 0) {
                centers[j] = vector_divide(sum, count);
            }
        }
    }
}

// 球形打包分组函数
void sphere_packing_grouping(const std::vector<std::vector<double>>& training_data) {
    int n = training_data.size();  // 100
    int dim = training_data[0].size();  // 10
    
    // 清空全局变量
    group1_data.clear();
    group2_data.clear();
    group3_data.clear();
    group4_data.clear();
    
    // 执行K-means聚类，分成32个簇
    std::vector<int> cluster_assignments;
    std::vector<std::vector<double>> cluster_centers;
    kmeans_clustering(training_data, 32, cluster_assignments, cluster_centers);
    
    // 计算聚类中心的主轴（用于排序）
    // 简化版PCA：使用所有中心的平均值和第一个变化最大的维度
    std::vector<double> center_mean(dim, 0.0);
    for (const auto& center : cluster_centers) {
        center_mean = vector_add(center_mean, center);
    }
    center_mean = vector_divide(center_mean, cluster_centers.size());
    
    // 计算每个维度的方差
    std::vector<double> variances(dim, 0.0);
    for (const auto& center : cluster_centers) {
        for (int i = 0; i < dim; ++i) {
            double diff = center[i] - center_mean[i];
            variances[i] += diff * diff;
        }
    }
    
    // 找到方差最大的维度
    int max_var_dim = 0;
    double max_var = variances[0];
    for (int i = 1; i < dim; ++i) {
        if (variances[i] > max_var) {
            max_var = variances[i];
            max_var_dim = i;
        }
    }
    
    // 按最大方差维度对聚类中心排序
    std::vector<std::pair<double, int>> center_positions;
    for (int i = 0; i < static_cast<int>(cluster_centers.size()); ++i) {
        center_positions.push_back({cluster_centers[i][max_var_dim], i});
    }
    std::sort(center_positions.begin(), center_positions.end());
    
    // 创建数据点索引和它们的聚类分配
    std::vector<std::pair<int, int>> point_cluster_pairs;
    for (int i = 0; i < n; ++i) {
        point_cluster_pairs.push_back({i, cluster_assignments[i]});
    }
    
    // 按聚类中心的排序顺序分配数据点
    std::vector<std::vector<double>>* groups[4] = {&group1_data, &group2_data, &group3_data, &group4_data};
    int group_sizes[4] = {32, 32, 32, 4};
    int current_group = 0;
    int points_in_current_group = 0;
    
    // 按排序后的聚类顺序处理
    for (const auto& pos : center_positions) {
        int cluster_id = pos.second;
        
        // 找到属于当前聚类的所有点
        std::vector<int> points_in_cluster;
        for (const auto& pc : point_cluster_pairs) {
            if (pc.second == cluster_id) {
                points_in_cluster.push_back(pc.first);
            }
        }
        
        // 将这些点分配到当前组
        for (int point_idx : points_in_cluster) {
            // 检查当前组是否已满
            while (current_group < 4 && points_in_current_group >= group_sizes[current_group]) {
                current_group++;
                points_in_current_group = 0;
            }
            
            // 分配点到组
            if (current_group < 4) {
                groups[current_group]->push_back(training_data[point_idx]);
                points_in_current_group++;
            }
        }
    }
    
    // 处理可能的剩余点（由于聚类可能不完美均匀）
    // 确保每个组都有正确的大小
    int total_assigned = group1_data.size() + group2_data.size() + 
                        group3_data.size() + group4_data.size();
    
    if (total_assigned < n) {
        // 收集未分配的点
        std::vector<bool> assigned(n, false);
        for (const auto& vec : group1_data) {
            for (int i = 0; i < n; ++i) {
                if (training_data[i] == vec) {
                    assigned[i] = true;
                    break;
                }
            }
        }
        for (const auto& vec : group2_data) {
            for (int i = 0; i < n; ++i) {
                if (!assigned[i] && training_data[i] == vec) {
                    assigned[i] = true;
                    break;
                }
            }
        }
        for (const auto& vec : group3_data) {
            for (int i = 0; i < n; ++i) {
                if (!assigned[i] && training_data[i] == vec) {
                    assigned[i] = true;
                    break;
                }
            }
        }
        for (const auto& vec : group4_data) {
            for (int i = 0; i < n; ++i) {
                if (!assigned[i] && training_data[i] == vec) {
                    assigned[i] = true;
                    break;
                }
            }
        }
        
        // 将未分配的点添加到未满的组
        for (int i = 0; i < n; ++i) {
            if (!assigned[i]) {
                if (group1_data.size() < 32) {
                    group1_data.push_back(training_data[i]);
                } else if (group2_data.size() < 32) {
                    group2_data.push_back(training_data[i]);
                } else if (group3_data.size() < 32) {
                    group3_data.push_back(training_data[i]);
                } else if (group4_data.size() < 4) {
                    group4_data.push_back(training_data[i]);
                }
            }
        }
    }
}

// ========================================================================================== //

// 将data拓展成16维度的向量
void expand_groups_dimension_inplace() {
    // 扩展group1_data
    for (auto& vec : group1_data) {
        vec.resize(16, 0.0);  // resize会在后面补0
    }
    
    // 扩展group2_data
    for (auto& vec : group2_data) {
        vec.resize(16, 0.0);
    }
    
    // 扩展group3_data
    for (auto& vec : group3_data) {
        vec.resize(16, 0.0);
    }
    
    // 扩展group4_data
    for (auto& vec : group4_data) {
        vec.resize(16, 0.0);
    }
}

// 拓展query成16维度的向量
void expand_query_dimension_inplace() {
    query_data.resize(16, 0.0);  // resize会在后面补0
}

// 将group数据转换为一维向量的函数
vector<double> convert_group_to_1d(const vector<vector<double>>& group_data, int elements_per_row = 16) {
    vector<double> result;
    
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
    // 转换各组数据为一维向量
    group1_1d = convert_group_to_1d(group1_data);
    group2_1d = convert_group_to_1d(group2_data);
    group3_1d = convert_group_to_1d(group3_data);
    group4_1d = convert_group_to_1d(group4_data);

    // 打印转换结果的统计信息
    cout << "\n=== 1D Conversion Results ===" << endl;
    cout << "Group 1: " << group1_data.size() << " rows -> " << group1_1d.size() << " elements" << endl;
    cout << "Group 2: " << group2_data.size() << " rows -> " << group2_1d.size() << " elements" << endl;
    cout << "Group 3: " << group3_data.size() << " rows -> " << group3_1d.size() << " elements" << endl;
    cout << "Group 4: " << group4_data.size() << " rows -> " << group4_1d.size() << " elements" << endl;

    // 验证前几个元素
    cout << "\nGroup 1 (first 20 elements): ";
    for (int i = 0; i < min(32, (int)group1_1d.size()); i++) {
        cout << fixed << setprecision(3) << group1_1d[i] << " ";
    }
    cout << endl;
}

// 将query_data扩展重复32次的函数
vector<double> expand_query_data(int repeat_times = 32) {
    vector<double> expanded_query;
    
    for (int i = 0; i < repeat_times; i++) {
        for (double value : query_data) {
            expanded_query.push_back(value);
        }
    }
    
    return expanded_query;
}

// 可选：添加一个打印函数来验证分组
void print_group_samples() {
    cout << "\n=== Group 1 (first 2 rows) ===" << endl;
    for (int i = 0; i < min(2, (int)group1_data.size()); i++) {
        cout << "Row " << i << ": ";
        for (int j = 0; j < min(5, (int)group1_data[i].size()); j++) {
            cout << fixed << setprecision(2) << group1_data[i][j] << " ";
        }
        cout << "..." << endl;
    }
    
    cout << "\n=== Group 2 (first 2 rows) ===" << endl;
    for (int i = 0; i < min(2, (int)group2_data.size()); i++) {
        cout << "Row " << i << ": ";
        for (int j = 0; j < min(5, (int)group2_data[i].size()); j++) {
            cout << fixed << setprecision(2) << group2_data[i][j] << " ";
        }
        cout << "..." << endl;
    }
    
    cout << "\n=== Group 3 (first 2 rows) ===" << endl;
    for (int i = 0; i < min(2, (int)group3_data.size()); i++) {
        cout << "Row " << i << ": ";
        for (int j = 0; j < min(5, (int)group3_data[i].size()); j++) {
            cout << fixed << setprecision(2) << group3_data[i][j] << " ";
        }
        cout << "..." << endl;
    }
}

int main() {
    read_data_json("../lsh/train.jsonl");
    /**
     * @brief 数据预处理
     */
    // 分组
    sphere_packing_grouping(training_data);
    // 扩展数据维度
    expand_groups_dimension_inplace();
    // 扩展query维度
    expand_query_dimension_inplace();
    // 归一化查询和训练数据
    normalize_query_and_training();
    // 处理输入向量成一维->group1_1d, group2_1d, group3_1d->32*16/group4_1d->4*16
    process_groups_to_1d();
    // 重复拓展query_data->32*16
    query_data_exp32 = expand_query_data(32);
    // 重复拓展query_data->4*16
    query_data_exp4 = expand_query_data(4);

    // 加密参数设置
    int num_slots1 = 32 * 16;
    int num_slots2 =  4 * 16; 
    int levels_required = 45;
    uint32_t ring_dim = 1 << 15; // 32768
    int gap1 = static_cast<int>(ring_dim)/ ( 2 * num_slots1); // 512
    int gap2 = static_cast<int>(ring_dim)/ ( 2 * num_slots2); // 4096
    vector<int> rotations_distance;
    // ============================ 针对32槽的加密上下文 ============================ //
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);

    int dcrtBits = 45;
    int firstMod = 48;

    if (toy) {
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);

        if (num_slots1 <= 1 << 14) parameters.SetRingDim(1 << 15);
        if (num_slots1 <= 1 << 13) parameters.SetRingDim(1 << 14);
        if (num_slots1 <= 1 << 12) parameters.SetRingDim(1 << 13);
        if (num_slots1 <= 1 << 11) parameters.SetRingDim(1 << 12);

        cout << "n: " << num_slots1 << endl;
    } else {
        // parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(1 << 15);
    }

    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slots1);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(2);

    parameters.SetMultiplicativeDepth(levels_required);

    // 生成加密上下文
    auto context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    auto key_pair = context->KeyGen();

    context->EvalMultKeyGen(key_pair.secretKey);

    // 生成旋转密钥--for permutation
    for (int i = 0; i < log2(n); i++) {
        vector<int> rotations;
        rotations.push_back(pow(2, i) * n);
        context->EvalRotateKeyGen(key_pair.secretKey, rotations);
        vector<int> rotations2;
        rotations2.push_back(pow(2, i));
        context->EvalRotateKeyGen(key_pair.secretKey, rotations2);
    }
    // 生成求欧式距离的旋转密钥--16
    for (int i = 0; i < log2(16); i++) {
        
        rotations_distance.push_back(pow(2, i) * gap1);
        rotations_distance.push_back(pow(2, i) * gap2);
    }
    context->EvalRotateKeyGen(key_pair.secretKey, rotations_distance);

    /**
    *@brief 初始数据编码成明文向量
    */
    // 先测试一下group1
    Plain ptxt1 = context->MakeCKKSPackedPlaintext(group1_1d, 1, 1, nullptr, static_cast<size_t>(num_slots1));
    ptxt1->SetLength(num_slots1);
    cout << "Group1 encode sucess!" << endl;
    // group2
    Plain ptxt2 = context->MakeCKKSPackedPlaintext(group2_1d, 1, 1, nullptr, static_cast<size_t>(num_slots1));
    ptxt2->SetLength(num_slots1);
    cout << "Group2 encode sucess!" << endl;
    // group3
    Plain ptxt3 = context->MakeCKKSPackedPlaintext(group3_1d, 1, 1, nullptr, static_cast<size_t>(num_slots1));
    ptxt3->SetLength(num_slots1);
    cout << "Group3 encode sucess!" << endl;
    // group4
    Plain ptxt4 = context->MakeCKKSPackedPlaintext(group4_1d, 1, 1, nullptr, static_cast<size_t>(num_slots2));
    ptxt4->SetLength(num_slots2);
    cout << "Group4 encode sucess!" << endl;
    // query
    Plain ptxt_query32 = context->MakeCKKSPackedPlaintext(query_data_exp32, 1, 1, nullptr, static_cast<size_t>(num_slots1));
    ptxt_query32->SetLength(num_slots1);
    Plain ptxt_query4 = context->MakeCKKSPackedPlaintext(query_data_exp4, 1, 1, nullptr, static_cast<size_t>(num_slots2));
    ptxt_query4->SetLength(num_slots2);
    cout << "Query encode sucess!" << endl;
    // std::cout << "Input: " << ptxt1;

    /**
    *@brief 初始数据加密
    */
    // group1加密
    Cipher ctgr1 = context->Encrypt(key_pair.publicKey, ptxt1);
    uint32_t current_level = ctgr1->GetLevel();
    uint32_t current_depth = levels_required - current_level;
    cout << "Group1 encrypt sucess!" <<endl;
    // group2加密
    Cipher ctgr2 = context->Encrypt(key_pair.publicKey, ptxt2);
    cout << "Group2 encrypt sucess!" <<endl;
    // group3加密
    Cipher ctgr3 = context->Encrypt(key_pair.publicKey, ptxt3);
    cout << "Group3 encrypt sucess!" <<endl;
    // group4加密
    Cipher ctgr4 = context->Encrypt(key_pair.publicKey, ptxt4);
    cout << "Group4 encrypt sucess!" <<endl;
    // query加密
    Cipher ctqry32 = context->Encrypt(key_pair.publicKey, ptxt_query32);
    Cipher ctqry4 = context->Encrypt(key_pair.publicKey, ptxt_query4);
    cout << "Query encrypt sucess!" <<endl;
    std::cout << "Current depth is: " << current_depth << std::endl;

    /**
    *@brief 欧式距离计算
    */
    // 减法做差
    auto ctminus1 = context->EvalSub(ctgr1, ctqry32);
    auto ctminus2 = context->EvalSub(ctgr2, ctqry32);
    auto ctminus3 = context->EvalSub(ctgr3, ctqry32);
    auto ctminus4 = context->EvalSub(ctgr4, ctqry4);
    cout << "Subtraction done!" << endl;

    // 平方
    auto ctminus1_square = context->EvalSquare(ctminus1);
    auto ctminus2_square = context->EvalSquare(ctminus2);
    auto ctminus3_square = context->EvalSquare(ctminus3);
    auto ctminus4_square = context->EvalSquare(ctminus4);
    cout << "Square done!" << endl;
    
    auto ctdistance1_square = ctminus1_square->Clone();
    auto ctdistance2_square = ctminus2_square->Clone();
    auto ctdistance3_square = ctminus3_square->Clone();
    auto ctdistance4_square = ctminus4_square->Clone();
    // 旋转相加--计算欧式距离
    for (int i = 0; i < log2(16); i++) {
        ctdistance1_square = context->EvalAdd(ctdistance1_square, context->EvalRotate(ctdistance1_square, gap1 * pow(2, i)));
    }

    

    // Cipher ctqry = hec.encrypt(query_data_exp, 0, num_slots1);// query拓展的加密

    // auto ct_gr1_minus_qry = hec.sub(ctgr1, ctqry);// group1-query

    // auto ct_gr1_minus_qry_square = hec.mult(ct_gr1_minus_qry, ct_gr1_minus_qry);// 平方

    // auto ct_gr1_distance1 = context->EvalSum(ct_gr1_minus_qry_square, 10);

    // auto ct_gr1_distance2 = context->EvalSum(ct_gr1_minus_qry_square, num_slots1);

    // auto ct_gr1_distance3 = context->EvalSum(ct_gr1_minus_qry_square, 10 * num_slots1);

    Plain result1;
    Plain result4;
    context->Decrypt(key_pair.secretKey, ctdistance1_square, &result1);
    context->Decrypt(key_pair.secretKey, ctgr4, &result4);

    // auto plain2 = hec.decrypt(ct_gr1_distance2);

    // auto plain3 = hec.decrypt(ct_gr1_distance3);

    result1->SetLength(static_cast<size_t>(num_slots1));
    result4->SetLength(static_cast<size_t>(num_slots2));

    cout << "result1:" << result1 << endl;
    // cout << "result4:" << result4 << endl;


    // 打印分组样本以验证
    // print_group_samples();
    
    return 0;
}

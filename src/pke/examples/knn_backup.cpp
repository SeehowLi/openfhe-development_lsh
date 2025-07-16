#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <omp.h>
#include <chrono>

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
vector<vector<double>> group3_data;  // 后36个
vector<vector<double>> group4_data;  // 最后4个

vector<double> query_data_exp32;  //拓展的query_data
vector<double> query_data_exp4;  //拓展的query_data
vector<double> group1_1d;  // 前32个
vector<double> group2_1d;  // 中间32个
vector<double> group3_1d;  // 后32个
vector<double> group4_1d;  // 最后4个

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
        std::cout << "[TIMING] " << task_name << ": " << duration.count() << " ms" << std::endl;
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
    Timer timer("Normalization");
    
    // 现在假设数据分布的最大值是10.0
    double max_val = 10.0;
    double max_distance = sqrt(10.0) * 2 * max_val; // 计算最大距离，用于归一化
    std::cout << "Max distance for normalization: " << max_distance << std::endl;

    // OpenMP并行归一化
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // 归一化query_data
            for (double& v : query_data) {
                v /= max_distance;
            }
        }
        
        #pragma omp section
        {
            // 归一化training_data
            #pragma omp parallel for
            for (size_t i = 0; i < training_data.size(); i++) {
                for (double& v : training_data[i]) {
                    v /= max_distance;
                }
            }
        }
        
        #pragma omp section
        {
            // 归一化group1_data
            #pragma omp parallel for
            for (size_t i = 0; i < group1_data.size(); i++) {
                for (double& v : group1_data[i]) {
                    v /= max_distance;
                }
            }
        }
        
        #pragma omp section
        {
            // 归一化group2_data
            #pragma omp parallel for
            for (size_t i = 0; i < group2_data.size(); i++) {
                for (double& v : group2_data[i]) {
                    v /= max_distance;
                }
            }
        }
        
        #pragma omp section
        {
            // 归一化group3_data
            #pragma omp parallel for
            for (size_t i = 0; i < group3_data.size(); i++) {
                for (double& v : group3_data[i]) {
                    v /= max_distance;
                }
            }
        }
        
        #pragma omp section
        {
            // 归一化group4_data
            #pragma omp parallel for
            for (size_t i = 0; i < group4_data.size(); i++) {
                for (double& v : group4_data[i]) {
                    v /= max_distance;
                }
            }
        }
    }
}

// 读取jsonl数据
void read_data_json(string filename) {
    Timer timer("Reading JSON data");
    
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
        group4_data.clear();

        // 第一组：前32个 (索引 0-31)
        for (int i = 0; i < 32; i++) {
            group1_data.push_back(training_data[i]);
        }

        // 第二组：中间32个 (索引 32-63)
        for (int i = 32; i < 64; i++) {
            group2_data.push_back(training_data[i]);
        }

        // 第三组：后36个 (索引 64-95)
        for (int i = 64; i < 96; i++) {
            group3_data.push_back(training_data[i]);
        }

        // 第四组：后4个 (索引 96-99)
        for (int i = 96; i < 100; i++) {
            group4_data.push_back(training_data[i]);
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
    cout << "Group 4 size: " << group4_data.size() << endl;
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
    #pragma omp parallel for
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// 计算向量除以标量
std::vector<double> vector_divide(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    #pragma omp parallel for
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
    Timer timer("K-means clustering");
    
    int n = data.size();
    int dim = data[0].size();
    
    // 使用固定的种子值，确保每次运行结果相同,便于调试
    std::mt19937 gen(26);  // 固定种子为26
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
        
        // 并行计算每个点到最近中心的距离
        #pragma omp parallel for
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
        
        // 并行分配步骤
        #pragma omp parallel for reduction(||:changed)
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
        
        // 并行更新步骤
        #pragma omp parallel for
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
    Timer timer("Sphere packing grouping");
    
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

// 将data拓展成32维度的向量
void expand_groups_dimension_inplace() {
    Timer timer("Expanding groups dimension");
    
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // 扩展group1_data
            for (auto& vec : group1_data) {
                vec.resize(32, 0.0);  // resize会在后面补0
            }
        }
        
        #pragma omp section
        {
            // 扩展group2_data
            for (auto& vec : group2_data) {
                vec.resize(32, 0.0);
            }
        }
        
        #pragma omp section
        {
            // 扩展group3_data
            for (auto& vec : group3_data) {
                vec.resize(32, 0.0);
            }
        }
        
        #pragma omp section
        {
            // 扩展group4_data
            for (auto& vec : group4_data) {
                vec.resize(32, 0.0);
            }
        }
    }
}

// 拓展query成16维度的向量
void expand_query_dimension_inplace() {
    Timer timer("Expanding query dimension");
    query_data.resize(32, 0.0);  // resize会在后面补0
}

// 将group数据转换为一维向量的函数
vector<double> convert_group_to_1d(const vector<vector<double>>& group_data, int elements_per_row = 32) {
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
    Timer timer("Processing groups to 1D");
    
    // 并行转换各组数据为一维向量
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            group1_1d = convert_group_to_1d(group1_data, 32);
        }
        
        #pragma omp section
        {
            group2_1d = convert_group_to_1d(group2_data, 32);
        }
        
        #pragma omp section
        {
            group3_1d = convert_group_to_1d(group3_data, 32);
        }
        
        #pragma omp section
        {
            group4_1d = convert_group_to_1d(group4_data, 32);
        }
    }

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
    expanded_query.reserve(query_data.size() * repeat_times);
    
    for (int i = 0; i < repeat_times; i++) {
        for (double value : query_data) {
            expanded_query.push_back(value);
        }
    }
    
    return expanded_query;
}

// 创建距离掩码函数
void create_distance_mask(vector<double>& distance_mask, int num_groups = 32, int elements_per_group = 32) {
    Timer timer("Creating distance mask");
    
    int total_length = num_groups * elements_per_group;
    distance_mask.clear();
    distance_mask.resize(total_length, 0.0);
    
    // 每32个数字的第一个设置为1
    for (int i = 0; i < num_groups; i++) {
        distance_mask[i * elements_per_group] = 1.0;
    }
    
    // // 验证前几个值
    // cout << "Distance mask (first 48 elements): ";
    // for (int i = 0; i < min(48, total_length); i++) {
    //     cout << distance_mask[i] << " ";
    //     if ((i + 1) % 16 == 0) cout << "| ";  // 每16个数字用|分隔
    // }
    // cout << endl;
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
    
    // 设置OpenMP线程数
    int num_threads = omp_get_max_threads();
    cout << "Using " << num_threads << " OpenMP threads" << endl;
    omp_set_num_threads(num_threads);
    
    auto total_start = chrono::high_resolution_clock::now();
    
    read_data_json("../lsh/train.jsonl");
    
    /**
     * @brief 数据预处理
     */
    {
        Timer timer("Total data preprocessing");
        
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
        group4_1d.resize(32 * 32, 0.0); // 最后4个group扩展成1024维度
        
        // 并行重复拓展query_data
        query_data_exp32 = expand_query_data(32);
    }
        

    // 加密参数设置
    int num_slot_dis1 = 32;// 前三组group中有效距离的个数
    int num_slot_dis2 = 4; // 最后一组group中有效距离的个数
    int num_dim = 16; // 每个group的维度
    int num_slot = num_slot_dis1 * num_slot_dis1; // 最终的槽数1024
    int num_slots1 = num_slot_dis1 * num_dim;
    // int num_slots2 = num_slot_dis2 * num_dim;
    // =============== 重要参数 =============== //
    int levels_required = 45;
    uint32_t ring_dim = 1 << 15; // 32768
    double lowbound_dis = 0.0;
    double upbound_dis = 0.9;
    uint32_t sqrt_cheb_degree = 31; // 五层
    vector<int> rotations_distance;
    // 距离转换的掩码
    vector<double> distance_mask1;
    vector<double> distance_mask2;

    // 创建距离掩码--每32个数字的第一个为1
    create_distance_mask(distance_mask1, num_slot_dis1, num_slot_dis1);
    create_distance_mask(distance_mask2, num_slot_dis2, num_slot_dis1);
    distance_mask2.resize(num_slot, 0.0);

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
        // 安全强度设置
        // parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(ring_dim);
    }

    cout << "N: " << parameters.GetRingDim() << endl;

    parameters.SetBatchSize(num_slot);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    //This keeps memory small, at the cost of increasing the modulus
    parameters.SetNumLargeDigits(2);

    parameters.SetMultiplicativeDepth(levels_required);

    // 生成加密上下文
    {
        Timer timer("Generating crypto context");
        auto context = GenCryptoContext(parameters);
        context->Enable(PKE);
        context->Enable(KEYSWITCH);
        context->Enable(LEVELEDSHE);
        context->Enable(ADVANCEDSHE);
        context->Enable(FHE);

        auto key_pair = context->KeyGen();

        context->EvalMultKeyGen(key_pair.secretKey);

        // 生成旋转密钥--for permutation,计算欧氏距离,获得紧凑距离
        // TODO:后面需要调整,有一些不需要或者生成的密钥多了
        {
            Timer timer("Generating rotation keys");
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    vector<int> rotations;
                    for (int i = 0; i < log2(num_slot_dis1); i++) {
                        rotations.push_back(pow(2, i) * num_slot_dis1);
                    }
                    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
                }
                
                #pragma omp section
                {
                    vector<int> rotations2;
                    for (int i = 0; i < log2(num_slot); i++) {
                        rotations2.push_back(pow(2, i));
                    }
                    context->EvalRotateKeyGen(key_pair.secretKey, rotations2);
                }
                // #pragma omp section // 好像不需要了
                // {
                //     vector<int> rotations3;
                //     for (int i = 0; i < log2(num_slot) + 1; i++) {
                //         rotations3.push_back(pow(2, i) * (num_dim - 1));
                //         if(i == log2(num_slot_dis1)){
                //             rotations3.push_back((num_dim - 1) * (num_dim - 1));
                //         }
                //     }
                //     context->EvalRotateKeyGen(key_pair.secretKey, rotations3);
                // }
                #pragma omp section
                {
                    context->EvalRotateKeyGen(key_pair.secretKey, {-32, 32});
                }
            }
        }

        /**
        *@brief 初始数据编码成明文向量
        */
        Plain ptxt1, ptxt2, ptxt3, ptxt4, ptxt_query32, ptxt_query4;
        
        {
            Timer timer("Encoding plaintexts");
            
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    ptxt1 = context->MakeCKKSPackedPlaintext(group1_1d, 1, 1, nullptr, static_cast<size_t>(num_slot));
                    ptxt1->SetLength(num_slot);
                    #pragma omp critical
                    cout << "Group1 encode success!" << endl;
                }
                
                #pragma omp section
                {
                    ptxt2 = context->MakeCKKSPackedPlaintext(group2_1d, 1, 1, nullptr, static_cast<size_t>(num_slot));
                    ptxt2->SetLength(num_slot);
                    #pragma omp critical
                    cout << "Group2 encode success!" << endl;
                }
                
                #pragma omp section
                {
                    ptxt3 = context->MakeCKKSPackedPlaintext(group3_1d, 1, 1, nullptr, static_cast<size_t>(num_slot));
                    ptxt3->SetLength(num_slot);
                    #pragma omp critical
                    cout << "Group3 encode success!" << endl;
                }
                
                #pragma omp section
                {
                    ptxt4 = context->MakeCKKSPackedPlaintext(group4_1d, 1, 1, nullptr, static_cast<size_t>(num_slot));
                    ptxt4->SetLength(num_slot);
                    #pragma omp critical
                    cout << "Group4 encode success!" << endl;
                }
                
                #pragma omp section
                {
                    ptxt_query32 = context->MakeCKKSPackedPlaintext(query_data_exp32, 1, 1, nullptr, static_cast<size_t>(num_slot));
                    ptxt_query32->SetLength(num_slot);
                    #pragma omp critical
                    cout << "Query32 encode success!" << endl;
                }
            }
        }
        // cout << "ptxt1:" << ptxt1 << endl;
        // cout << "ptxt_query32:" << ptxt_query32 << endl;

        /**
        *@brief 初始数据加密
        */
        Cipher ctgr1, ctgr2, ctgr3, ctgr4, ctqry32, ctqry4;
        
        {
            Timer timer("Encrypting data");
            
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    ctgr1 = context->Encrypt(key_pair.publicKey, ptxt1);
                    #pragma omp critical
                    cout << "Group1 encrypt success!" << endl;
                }
                
                #pragma omp section
                {
                    ctgr2 = context->Encrypt(key_pair.publicKey, ptxt2);
                    #pragma omp critical
                    cout << "Group2 encrypt success!" << endl;
                }
                
                #pragma omp section
                {
                    ctgr3 = context->Encrypt(key_pair.publicKey, ptxt3);
                    #pragma omp critical
                    cout << "Group3 encrypt success!" << endl;
                }
                
                #pragma omp section
                {
                    ctgr4 = context->Encrypt(key_pair.publicKey, ptxt4);
                    #pragma omp critical
                    cout << "Group4 encrypt success!" << endl;
                }
                
                #pragma omp section
                {
                    ctqry32 = context->Encrypt(key_pair.publicKey, ptxt_query32);
                    #pragma omp critical
                    cout << "Query32 encrypt success!" << endl;
                }
            }
        }

        cout << "Current depth is: " << levels_required - ctgr1->GetLevel() << endl;

        /**
        *@brief 欧式距离计算
        */
        // 减法做差
        Cipher ctminus1, ctminus2, ctminus3, ctminus4;
        
        {
            Timer timer("Subtraction");
            
            #pragma omp parallel sections
            {
                #pragma omp section
                ctminus1 = context->EvalSub(ctgr1, ctqry32);
                
                #pragma omp section
                ctminus2 = context->EvalSub(ctgr2, ctqry32);
                
                #pragma omp section
                ctminus3 = context->EvalSub(ctgr3, ctqry32);
                
                #pragma omp section
                ctminus4 = context->EvalSub(ctgr4, ctqry32);
            }
            cout << "Subtraction done!" << endl;
        }

        // 平方
        Cipher ctminus1_square, ctminus2_square, ctminus3_square, ctminus4_square;
        
        {
            Timer timer("Squaring");
            
            #pragma omp parallel sections
            {
                #pragma omp section
                ctminus1_square = context->EvalSquare(ctminus1);
                
                #pragma omp section
                ctminus2_square = context->EvalSquare(ctminus2);
                
                #pragma omp section
                ctminus3_square = context->EvalSquare(ctminus3);
                
                #pragma omp section
                ctminus4_square = context->EvalSquare(ctminus4);
            }
            cout << "Square done!" << endl;
        }
        cout << "Current depth is: " << levels_required - ctminus1_square->GetLevel() << endl;
        
        // 旋转预处理--折半处理
        // TODO:需要新的旋转密钥--放到前面的并行处理
        vector<double> bemask(1024),afmask(1024);
        for (int i = 0; i < 1024; i++) {
            bemask[i] = ((i / 32) % 2 == 0) ? 1.0 : 0.0;  // 先32个1，再32个0
            afmask[i] = ((i / 32) % 2 == 1) ? 1.0 : 0.0;  // 先32个0，再32个1
        }
        auto pt_bemask = context->MakeCKKSPackedPlaintext(bemask, 1, ctminus1_square->GetLevel(), nullptr, num_slot);
        auto pt_afmask = context->MakeCKKSPackedPlaintext(afmask, 1, ctminus1_square->GetLevel(), nullptr, num_slot);
        
        vector<double> rep2exp_mask(num_slot, 0.0);
        // 掩码的生成，每一组32个数，第x组的第x个数是1，其余是0
        for (int i = 0; i < 32; i++) {
            rep2exp_mask[i * 32 + i] = 1.0;
        }

        /**
         * @brief 计算repeat和expand
         */
        Cipher ct_bemask1, ct_afmask1; Cipher ct_bemask2, ct_afmask2;
        Cipher ct_bemask3, ct_afmask3; Cipher ct_bemask4, ct_afmask4;
        Cipher ctdistance1_rep,ctdistance2_rep,ctdistance3_rep,ctdistance4_rep;
        Cipher ctdistance1_exp,ctdistance2_exp,ctdistance3_exp,ctdistance4_exp;
        {
            Timer timer("Expand Ciphertext");
            // #1先填满
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    ct_bemask1 = context->EvalMult(ctminus1_square, pt_bemask);
                    ct_afmask1 = context->EvalMult(ctminus1_square, pt_afmask);
                    
                    context->EvalAddInPlace(ct_bemask1, context->EvalRotate(ct_bemask1, -32));
                    context->EvalAddInPlace(ct_afmask1, context->EvalRotate(ct_afmask1,  32));
                }
                #pragma omp section
                {
                    ct_bemask2 = context->EvalMult(ctminus2_square, pt_bemask);
                    ct_afmask2 = context->EvalMult(ctminus2_square, pt_afmask);
                    
                    context->EvalAddInPlace(ct_bemask2, context->EvalRotate(ct_bemask2, -32));
                    context->EvalAddInPlace(ct_afmask2, context->EvalRotate(ct_afmask2,  32));
                }
                #pragma omp section
                {
                    ct_bemask3 = context->EvalMult(ctminus3_square, pt_bemask);
                    ct_afmask3 = context->EvalMult(ctminus3_square, pt_afmask);
                    
                    context->EvalAddInPlace(ct_bemask3, context->EvalRotate(ct_bemask3, -32));
                    context->EvalAddInPlace(ct_afmask3, context->EvalRotate(ct_afmask3,  32));
                }
                #pragma omp section
                {
                    ct_bemask4 = context->EvalMult(ctminus4_square, pt_bemask);
                    ct_afmask4 = context->EvalMult(ctminus4_square, pt_afmask);
                    
                    context->EvalAddInPlace(ct_bemask4, context->EvalRotate(ct_bemask4, -32));
                    context->EvalAddInPlace(ct_afmask4, context->EvalRotate(ct_afmask4,  32));
                }
            }
            // #2旋转求和
            for (int i = 0; i < log2(num_slot_dis1); i++) {
                #pragma omp parallel sections
                {
                    #pragma omp section
                    {
                        context->EvalAddInPlace(ct_bemask1, context->EvalRotate(ct_bemask1, pow(2, i)));
                        context->EvalAddInPlace(ct_afmask1, context->EvalRotate(ct_afmask1, pow(2, i)));
                    }
                    #pragma omp section
                    {
                        context->EvalAddInPlace(ct_bemask2, context->EvalRotate(ct_bemask2, pow(2, i)));
                        context->EvalAddInPlace(ct_afmask2, context->EvalRotate(ct_afmask2, pow(2, i)));
                    }
                    #pragma omp section
                    {
                        context->EvalAddInPlace(ct_bemask3, context->EvalRotate(ct_bemask3, pow(2, i)));
                        context->EvalAddInPlace(ct_afmask3, context->EvalRotate(ct_afmask3, pow(2, i)));
                    }
                    #pragma omp section
                    {
                        context->EvalAddInPlace(ct_bemask4, context->EvalRotate(ct_bemask4, pow(2, i)));
                        context->EvalAddInPlace(ct_afmask4, context->EvalRotate(ct_afmask4, pow(2, i)));
                    }
                }
            }
            // #3后处理
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    ct_bemask1 = context->EvalMult(ct_bemask1, pt_bemask);
                    ct_afmask1 = context->EvalMult(ct_afmask1, pt_bemask);
                    ct_afmask1 = context->EvalRotate(ct_afmask1, -32);

                    ctdistance1_rep = context->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); }, 
                                                context->EvalAdd(ct_bemask1, ct_afmask1), lowbound_dis, upbound_dis, sqrt_cheb_degree);
                    auto pt_rep2exp_mask = context->MakeCKKSPackedPlaintext(rep2exp_mask, 1, ctdistance1_rep->GetLevel(), nullptr, num_slot);
                    ctdistance1_exp = context->EvalMult(ctdistance1_rep, pt_rep2exp_mask);
                }
                #pragma omp section
                {
                    ct_bemask2 = context->EvalMult(ct_bemask2, pt_bemask);
                    ct_afmask2 = context->EvalMult(ct_afmask2, pt_bemask);
                    ct_afmask2 = context->EvalRotate(ct_afmask2, -32);

                    ctdistance2_rep = context->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); }, 
                                                context->EvalAdd(ct_bemask2, ct_afmask2), lowbound_dis, upbound_dis, sqrt_cheb_degree);
                    auto pt_rep2exp_mask = context->MakeCKKSPackedPlaintext(rep2exp_mask, 1, ctdistance2_rep->GetLevel(), nullptr, num_slot);
                    ctdistance2_exp = context->EvalMult(ctdistance2_rep, pt_rep2exp_mask);
                }
                #pragma omp section
                {
                    ct_bemask3 = context->EvalMult(ct_bemask3, pt_bemask);
                    ct_afmask3 = context->EvalMult(ct_afmask3, pt_bemask);
                    ct_afmask3 = context->EvalRotate(ct_afmask3, -32);

                    ctdistance3_rep = context->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); }, 
                                                context->EvalAdd(ct_bemask3, ct_afmask3), lowbound_dis, upbound_dis, sqrt_cheb_degree);
                    auto pt_rep2exp_mask = context->MakeCKKSPackedPlaintext(rep2exp_mask, 1, ctdistance3_rep->GetLevel(), nullptr, num_slot);
                    ctdistance3_exp = context->EvalMult(ctdistance3_rep, pt_rep2exp_mask);
                }
                #pragma omp section
                {
                    ct_bemask4 = context->EvalMult(ct_bemask4, pt_bemask);
                    ct_afmask4 = context->EvalMult(ct_afmask4, pt_bemask);
                    ct_afmask4 = context->EvalRotate(ct_afmask4, -32);

                    ctdistance4_rep = context->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); }, 
                                                context->EvalAdd(ct_bemask4, ct_afmask4), lowbound_dis, upbound_dis, sqrt_cheb_degree);
                }
            }
        }
        cout << "Comsumed depth is: " << ctdistance1_rep->GetLevel() << endl;
        cout << "Current depth is: " << levels_required - ctdistance1_rep->GetLevel() << endl;

        /**
         * @brief 计算expand
         */
        // 利用repeat生成expand--不需要算第四组的
        {
            Timer timer("Expand Ciphertext");
            // 旋转求和
            for (int i = 0; i < log2(num_slot_dis1); i++) {
                    #pragma omp parallel sections
                    {
                        #pragma omp section
                        {
                            context->EvalAddInPlace(ctdistance1_exp, context->EvalRotate(ctdistance1_exp, num_slot_dis1 * pow(2, i)));
                        }
                        #pragma omp section
                        {
                            context->EvalAddInPlace(ctdistance2_exp, context->EvalRotate(ctdistance2_exp, num_slot_dis1 * pow(2, i)));
                        }
                        #pragma omp section
                        {
                            context->EvalAddInPlace(ctdistance3_exp, context->EvalRotate(ctdistance3_exp, num_slot_dis1 * pow(2, i)));
                        }
                    }
                }
        }    
        cout << "Comsumed depth is: " << ctdistance1_exp->GetLevel() << endl;
        cout << "Current depth is: " << levels_required - ctdistance1_exp->GetLevel() << endl;


        /**
         * @brief 分别计算前三个group的sorting
         */
        // 创建排序参数
        SortingParams sorting_params;
        sorting_params.method = SortingParams::PERMUTATION;
        sorting_params.vectorSize = 32;  // 每个group有32个距离值
        sorting_params.precision = 0.001;  // 使用您设置的delta
        sorting_params.verbose = false;    // 或者设为true查看详细信息
        
        // 解密结果
        Plain result1, result2, result3, result4;
        
        {
            Timer timer("Decryption");
            
            context->Decrypt(key_pair.secretKey, ctdistance1_exp, &result1);
            context->Decrypt(key_pair.secretKey, ctdistance2_exp, &result2);
            context->Decrypt(key_pair.secretKey, ctdistance3_exp, &result3);
            context->Decrypt(key_pair.secretKey, ctdistance1_rep, &result4);
        }

        // result1->SetLength(static_cast<size_t>(num_slot_dis1));
        // result2->SetLength(static_cast<size_t>(num_slots1));
        // result4->SetLength(static_cast<size_t>(num_slots2));


        cout << "result1:" << result1 << endl;
        cout << "result2:" << result2 << endl;
        cout << "result3:" << result3 << endl;
        cout << "result4:" << result4 << endl;
    }

    
    auto total_end = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::milliseconds>(total_end - total_start);
    cout << "\n[TIMING] Total execution time: " << total_duration.count() << " ms" << endl;
    
    return 0;
}
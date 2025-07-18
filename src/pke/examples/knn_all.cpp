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

// 全局变量
HomoEncryptCompute hec;
vector<double> input_values;
int n = 32;
double delta = 0.001;
int precision_digits = 3;
bool toy = true;
SortingType sortingType = NONE;

// 全局变量来存储数据
vector<double> query_data;
vector<vector<double>> training_data;
vector<vector<double>> group1_data;  // 前32个
vector<vector<double>> group2_data;  // 中间32个
vector<vector<double>> group3_data;  // 后36个
vector<vector<double>> group4_data;  // 最后4个

vector<double> query_data_exp32;  //拓展的query_data
vector<double> training_data_1d;  // 前32个
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
    Timer timer("Normalization");
    
    // 现在假设数据分布的最大值是10.0
    double max_val = 10.0;
    double max_distance = sqrt(10.0) * 2 * max_val; // 计算最大距离，用于归一化
    cout << "Max distance for normalization: " << max_distance << std::endl;

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

// 将data拓展成128维度的向量
void expand_groups_dimension_inplace() {
    Timer timer("Expanding groups dimension");
    
    for (auto& vec : training_data) {
        vec.resize(128, 0.0);
    }
}

// 拓展query成128维度的向量
void expand_query_dimension_inplace() {
    Timer timer("Expanding query dimension");
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
    Timer timer("Processing groups to 1D");

    training_data_1d = convert_group_to_1d(training_data, 128);

    // 打印转换结果的统计信息
    // cout << "\n=== 1D Conversion Results ===" << endl;
    // cout << "Training Data: " << training_data.size() << " rows -> " << training_data_1d.size() << " elements" << endl;

    // 验证前几个元素
    // cout << "\nTraining Data (first 20 elements): ";
    // for (int i = 0; i < min(128, (int)training_data_1d.size()); i++) {
    //     cout << fixed << setprecision(3) << training_data_1d[i] << " ";
    // }
    // cout << endl;
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
    int levels_required = 50;
    uint32_t ring_dim = 1 << 16; // 65536
    double lowbound_dis = -0.6;
    double upbound_dis = 0.6;
    uint32_t sqrt_cheb_degree = 31; // 7层
    vector<int> rotations_distance;
    // 距离转换的掩码
    vector<double> distance_mask1;
    vector<double> distance_mask2;

    // 创建距离掩码--每128个数字的第一个为1
    create_distance_mask(distance_mask1, num_slot_dis1, num_slot_dis1);

    // ============================ 针对128槽的加密上下文 ============================ //
    {
        Timer timer("Generating crypto context and rotation keys");
        hec.generate_context_knn(num_slot, levels_required, ring_dim, toy);
    }
    
    /**
    *@brief 初始数据编码成明文向量
    */
    Plain pt_data ,ptxt_query32;
    {
        Timer timer("Encoding plaintexts");
        
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                pt_data = hec.encode(training_data_1d, 1, 1, num_slot);
                #pragma omp critical
                cout << "Training Data encode success!" << endl;
            }
            
            #pragma omp section
            {
                ptxt_query32 =  hec.encode(query_data_exp32, 1, 1, num_slot);
                #pragma omp critical
                cout << "Query32 encode success!" << endl;
            }
        }
    }

    /**
    *@brief 初始数据加密
    */
    Cipher ctdata, ctqry32;
    
    {
        Timer timer("Encrypting data");
        
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                ctdata = hec.encrypt(pt_data);
                #pragma omp critical
                cout << "Training Data encrypt success!" << endl;
            }
            
            #pragma omp section
            {
                ctqry32 = hec.encrypt(ptxt_query32);
                #pragma omp critical
                cout << "Query32 encrypt success!" << endl;
            }
        }
    }
    cout << "Current depth is: " << levels_required - ctdata->GetLevel() << endl;

    /**
    *@brief 欧式距离计算
    */
    // 减法做差
    Cipher ctminus;
    
    {
        Timer timer("Subtraction");

        ctminus = hec.sub(ctdata, ctqry32)->Clone();
        cout << "Subtraction done!" << endl;
    }

    // 平方
    Cipher ctminus_square;
    
    {
        Timer timer("Squaring");

        ctminus_square = hec.square(ctminus)->Clone();
        cout << "Square done!" << endl;
    }
    cout << "Current depth is: " << levels_required - ctminus_square->GetLevel() << endl;

    // 旋转预处理--折半处理
    // TODO:需要新的旋转密钥--放到前面的并行处理
    vector<double> bemask(num_slot_dis1 * num_slot_dis1);
    vector<double> afmask(num_slot_dis1 * num_slot_dis1);
    vector<double> bemaskafter(num_slot_dis1 * num_slot_dis1);

    for (int i = 0; i < num_slot_dis1 * num_slot_dis1; i++) {
        int row = i / num_slot_dis1;  // 当前是第几行（0-127）
        int col = i % num_slot_dis1;  // 当前行内的位置（0-127）
        
        // bemask: 偶数行原本是1，但只有前100个保持为1
        bemask[i] = (row % 2 == 0 && col < 100) ? 1.0 : 0.0;
        // afmask: 奇数行原本是1，但只有前100个保持为1
        afmask[i] = (row % 2 == 1 && col < 100) ? 1.0 : 0.0;

        // bemask: 偶数行1
        bemaskafter[i] = (row % 2 == 0 && col < 128) ? 1.0 : 0.0;
        // afmask: 奇数行1
        // afmaskafter[i] = (row % 2 == 1 && col < 128) ? 1.0 : 0.0;

    }
    auto pt_bemask = hec.encode(bemask, 1, ctminus_square->GetLevel(), num_slot);
    auto pt_afmask = hec.encode(afmask, 1, ctminus_square->GetLevel(), num_slot);
    auto pt_bemaskafter = hec.encode(bemaskafter, 1, ctminus_square->GetLevel(), num_slot);
    // auto pt_afmaskafter = hec.encode(afmaskafter, 1, ctminus_square->GetLevel(), num_slot);

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
    Cipher ct_bemask, ct_afmask;
    Cipher ctdistance_rep;
    Cipher ctdistance_exp;
    {
        Timer timer("Expand Ciphertext");
        // #1先填满
        ct_bemask= hec.mult(ctminus_square, pt_bemask)->Clone();
        ct_afmask = hec.mult(ctminus_square, pt_afmask)->Clone();

        hec.add_inplace(ct_bemask, hec.rot(ct_bemask, -1 * num_slot_dis1));
        hec.add_inplace(ct_afmask, hec.rot(ct_afmask,      num_slot_dis1));
        // #2旋转求和
        for (int i = 0; i < log2(num_slot_dis1); i++) {
            hec.add_inplace(ct_bemask, hec.rot(ct_bemask, pow(2, i)));
            hec.add_inplace(ct_afmask, hec.rot(ct_afmask, pow(2, i)));
        }
        // #3后处理
        ct_bemask = hec.mult(ct_bemask, pt_bemaskafter)->Clone();
        ct_afmask = hec.mult(ct_afmask, pt_bemaskafter)->Clone();
        ct_afmask = hec.rot(ct_afmask, -1 * num_slot_dis1)->Clone();
        // 开四次方跟，拉大不同数值之间的距离
        // ctdistance_rep = hec.chebyshev([](double x) -> double { return std::sqrt(x); }, 
        //                             hec.add(ct_bemask, ct_afmask), lowbound_dis, upbound_dis, sqrt_cheb_degree);
        ctdistance_rep = hec.add(ct_bemask, ct_afmask)->Clone();

        auto pt_rep2exp_mask = hec.encode(rep2exp_mask, 1, ctdistance_rep->GetLevel(), num_slot);
        ctdistance_exp = hec.mult(ctdistance_rep, pt_rep2exp_mask)->Clone();

        auto pt_rep_mask = hec.encode(rep_mask, 1, ctdistance_rep->GetLevel(), num_slot);
        ctdistance_rep = hec.mult(ctdistance_rep, pt_rep_mask)->Clone();
    }
    cout << "Comsumed depth is: " << ctdistance_rep->GetLevel() << endl;
    cout << "Current depth is: " << levels_required - ctdistance_rep->GetLevel() << endl;
    
                                   
    /**
     * @brief 计算expand
     */
    // 利用repeat生成expand--不需要算第四组的
    {
        Timer timer("Expand Ciphertext");
        // 旋转求和
        for (int i = 0; i < log2(num_slot_dis1); i++) {
            hec.add_inplace(ctdistance_exp, hec.rot(ctdistance_exp, num_slot_dis1 * pow(2, i)));
        }
    }
    
    cout << "Comsumed depth is: " << ctdistance_exp->GetLevel() << endl;
    cout << "Current depth is: " << levels_required - ctdistance_exp->GetLevel() << endl;

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
    
    auto ctdistance_del_pro = hec.chebyshev([](double x) -> double { return (std::exp(10 * x) - std::exp(-10 * x))/(std::exp(10 * x) + std::exp(-10 * x)); }, 
                                    ctdistance_del, lowbound_dis, upbound_dis, sqrt_cheb_degree);

    int distance_sig_degree = 2000; // Sigmoid的degree
    int distance_sig_scaling = 9170; // Sigmoid的scaling
    auto ct_distance_mask = hec.sigmoid_tight(ctdistance_del_pro, 1, distance_sig_degree, distance_sig_scaling);
    
    cout << "Comsumed depth is: " << ct_distance_mask->GetLevel() << endl;
    cout << "Current depth is: " << levels_required - ct_distance_mask->GetLevel() << endl;

    auto ctindex = ct_distance_mask->Clone();
    // 旋转求和，得到整体数据的index
    for(int i = 0; i < log(num_slot_dis1); i++) {
        hec.add_inplace(ctindex, hec.rot(ctindex, num_slot_dis1 * pow(2, i)));
    }

    // vector<double> const05(num_slot, 0.5);
    // auto pt_const05 = hec.encode(const05, 1, ctindex->GetLevel(), num_slot);
    // ctindex = hec.add(ctindex, pt_const05);


    /**
     * @brief 分别计算前三个group的sorting
     */
    // 创建排序参数
    SortingParams sorting_params;
    sorting_params.method = SortingParams::PERMUTATION;
    sorting_params.vectorSize = 32;  // 每个group有32个距离值
    sorting_params.precision = 0.001;  // 使用您设置的delta
    sorting_params.verbose = false;    // 或者设为true查看详细信息
    
    // 根据precision自动设置参数
    if (sorting_params.precision >= 0.001) {
        sorting_params.sigmoidScaling = 10000;
        sorting_params.sigmoidDegree = 16000;
        sorting_params.sincDegree = 119;  // 对于n=32
    }

    // 3. 验证提取的标准距离（可选，用于调试）
    // if (sorting_params.verbose) {
    //     Plain test_standard1 = hec.decrypt(ctdistance_standard);
    //     vector<double> test_values = hec.decode(test_standard1);
    //     cout << "Standard distances for Group 1 (first 10): ";
    //     for (int i = 0; i < min(10, 32); i++) {
    //         cout << fixed << setprecision(4) << test_values[i] << " ";
    //     }
    //     cout << endl;
    // }

    // 4. 创建排序器
    cout << "\n=== Creating Permutation Sorter ===" << endl;
    auto sorter = SortingMethod::Create(sorting_params.method, hec, sorting_params);

    // 解密结果
    Plain result1, result2, result3, result4, result5;
    
    {
        Timer timer("Decryption");
        
        result1 = hec.decrypt(ctdistance_exp);
        result2 = hec.decrypt(ctdistance_rep);
        result3 = hec.decrypt(ctdistance_del_pro);
        result4 = hec.decrypt(ct_distance_mask);
        result5 = hec.decrypt(ctindex);
    }
    
    result1->SetLength(static_cast<size_t>(num_slot_dis1));
    result2->SetLength(static_cast<size_t>(num_slot_dis1));
    result3->SetLength(static_cast<size_t>(num_slot_dis1));
    // result4->SetLength(static_cast<size_t>(num_slot_dis1));
    // result5->SetLength(static_cast<size_t>(num_slot_dis1));


    // cout << "result1:" << result1 << endl;
    // cout << "result2:" << result2 << endl;
    // cout << "result3:" << result3 << endl;
    // cout << "result4:" << result4 << endl;
    // cout << "result5:" << result5 << endl;
    for(size_t i=127; i < result4->GetLength(); i+= 128) {
            cout << "result4[" << i << "]: " << result4->GetCKKSPackedValue()[i] << " ";
        }
    cout << endl;

    for(size_t i=100; i < result4->GetLength(); i+= 128) {
            cout << "result4[" << i << "]: " << result4->GetCKKSPackedValue()[i] << " ";
        }
    cout << endl;

    auto total_end = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::milliseconds>(total_end - total_start);
    cout << "\n[TIMING] Total execution time: " << total_duration.count() << " ms" << endl;
    
    return 0;
}
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cmath>
#include <queue>
#include <functional>
#include <ctime>
#include <set>
#include <fstream>
using namespace std;

// 平衡准确率和覆盖率的参数
const int MIN_SEED_LEN = 12;              // 适当提高种子长度要求
const int MAX_CHAIN_GAP = 30000;          // 合理链间隙阈值
const int MAX_LOCAL_GAP = 10000;          // 局部比对间隙
const int BAND_WIDTH = 300;               // 带宽
const double MIN_SEED_SIM = 0.75;         // 提高最小种子相似度
const int KMER_SIZE = 11;                 // 折中的k-mer尺寸
const int MIN_CHAIN_SCORE = 20;           // 最小链得分阈值
const int MIN_ANCHOR_LEN = 20;            // 最小锚点长度
const int MAX_SEEDS = 500000;             // 合理种子数量
const int MAX_KMER_MISMATCH = 1;          // 只允许1个错配
const int MAX_EXTEND_LEN = 300;           // 最大种子扩展长度
const int MIN_MATCHES_FOR_ANCHOR = 15;    // 锚点最小匹配数

struct Anchor {
    int q_start, q_end;
    int r_start, r_end;
    bool is_reverse;
    double score;
    int matches;  // 实际匹配数
    
    Anchor(int qs, int qe, int rs, int re, bool rev = false, double s = 0.0, int m = 0)
        : q_start(qs), q_end(qe), r_start(rs), r_end(re), is_reverse(rev), score(s), matches(m) {}
    
    int length() const { return q_end - q_start; }
    int diag() const { return is_reverse ? (r_end - q_start) : (q_start - r_start); }
    int q_center() const { return (q_start + q_end) / 2; }
    int r_center() const { return (r_start + r_end) / 2; }
};

// 高效生成种子（更严格的质量控制）
vector<Anchor> generate_seeds(const string& ref, const string& query) {
    unordered_map<string, vector<int>> ref_index;
    vector<Anchor> seeds;
    
    // 构建参考索引（允许少量错配）
    for (int i = 0; i <= (int)ref.size() - KMER_SIZE; i++) {
        string kmer = ref.substr(i, KMER_SIZE);
        ref_index[kmer].push_back(i);
        
        // 生成1个错配的变体
        for (int j = 0; j < KMER_SIZE; j++) {
            string variant = kmer;
            for (char c : {'A', 'T', 'C', 'G'}) {
                if (c == kmer[j]) continue;
                variant[j] = c;
                ref_index[variant].push_back(i);
            }
        }
    }
    
    // 生成正向种子
    for (int q_pos = 0; q_pos <= (int)query.size() - KMER_SIZE; q_pos++) {
        if (seeds.size() > MAX_SEEDS) break;
        
        string kmer = query.substr(q_pos, KMER_SIZE);
        set<int> processed_positions;
        
        // 处理精确匹配和1个错配
        for (int mismatch = 0; mismatch <= MAX_KMER_MISMATCH; mismatch++) {
            vector<string> variants;
            if (mismatch == 0) {
                variants.push_back(kmer);
            } else {
                for (int j = 0; j < KMER_SIZE; j++) {
                    string variant = kmer;
                    for (char c : {'A', 'T', 'C', 'G'}) {
                        if (c == kmer[j]) continue;
                        variant[j] = c;
                        variants.push_back(variant);
                    }
                }
            }
            
            for (const auto& var : variants) {
                if (ref_index.find(var) == ref_index.end()) continue;
                
                for (int r_pos : ref_index[var]) {
                    if (processed_positions.find(r_pos) != processed_positions.end()) continue;
                    processed_positions.insert(r_pos);
                    
                    // 向左扩展（要求连续匹配）
                    int left_ext = 0;
                    int left_matches = 0;
                    while (q_pos - left_ext > 0 && r_pos - left_ext > 0 && left_ext < MAX_EXTEND_LEN) {
                        if (query[q_pos - left_ext - 1] != ref[r_pos - left_ext - 1]) {
                            // 允许少量不匹配但不超过1个
                            if (left_matches > 0) break;
                            left_matches++;
                        }
                        left_ext++;
                    }
                    
                    // 向右扩展（要求连续匹配）
                    int right_ext = 0;
                    int right_matches = 0;
                    while (q_pos + KMER_SIZE + right_ext < query.size() && 
                           r_pos + KMER_SIZE + right_ext < ref.size() && 
                           right_ext < MAX_EXTEND_LEN) {
                        if (query[q_pos + KMER_SIZE + right_ext] != ref[r_pos + KMER_SIZE + right_ext]) {
                            if (right_matches > 0) break;
                            right_matches++;
                        }
                        right_ext++;
                    }
                    
                    int total_len = KMER_SIZE + left_ext + right_ext;
                    int total_matches = KMER_SIZE - mismatch + left_ext + right_ext - left_matches - right_matches;
                    double similarity = (double)total_matches / total_len;
                    
                    if (total_len >= MIN_SEED_LEN && similarity >= MIN_SEED_SIM) {
                        seeds.push_back({
                            q_pos - left_ext, 
                            q_pos + KMER_SIZE + right_ext,
                            r_pos - left_ext, 
                            r_pos + KMER_SIZE + right_ext,
                            false,
                            (double)total_matches,
                            total_matches
                        });
                    }
                    
                    if (seeds.size() > MAX_SEEDS) break;
                }
                if (seeds.size() > MAX_SEEDS) break;
            }
            if (seeds.size() > MAX_SEEDS) break;
        }
        
        // 生成反向种子（更严格的质量控制）
        string rc_kmer = "";
        for (int i = KMER_SIZE-1; i >= 0; --i) {
            switch(kmer[i]) {
                case 'A': rc_kmer += 'T'; break;
                case 'T': rc_kmer += 'A'; break;
                case 'C': rc_kmer += 'G'; break;
                case 'G': rc_kmer += 'C'; break;
                default: rc_kmer += 'N';
            }
        }
        
        processed_positions.clear();
        for (int mismatch = 0; mismatch <= MAX_KMER_MISMATCH; mismatch++) {
            vector<string> variants;
            if (mismatch == 0) {
                variants.push_back(rc_kmer);
            } else {
                for (int j = 0; j < KMER_SIZE; j++) {
                    string variant = rc_kmer;
                    for (char c : {'A', 'T', 'C', 'G'}) {
                        if (c == rc_kmer[j]) continue;
                        variant[j] = c;
                        variants.push_back(variant);
                    }
                }
            }
            
            for (const auto& var : variants) {
                if (ref_index.find(var) == ref_index.end()) continue;
                
                for (int r_pos : ref_index[var]) {
                    if (processed_positions.find(r_pos) != processed_positions.end()) continue;
                    processed_positions.insert(r_pos);
                    
                    // 改进反向种子扩展逻辑
                    int left_ext = 0;
                    int right_ext = 0;
                    int left_matches = 0;
                    int right_matches = 0;
                    
                    // 同时向两侧扩展
                    while (left_ext < MAX_EXTEND_LEN || right_ext < MAX_EXTEND_LEN) {
                        bool extended = false;
                        
                        // 向左扩展 (query向左, ref向右)
                        if (q_pos - left_ext - 1 >= 0 && 
                            r_pos + KMER_SIZE + left_ext < ref.size() &&
                            left_ext < MAX_EXTEND_LEN) {
                            char q_char = query[q_pos - left_ext - 1];
                            char r_char = ref[r_pos + KMER_SIZE + left_ext];
                            char rc_r_char = 
                                r_char == 'A' ? 'T' :
                                r_char == 'T' ? 'A' :
                                r_char == 'C' ? 'G' : 'C';
                            
                            if (q_char == rc_r_char) {
                                left_ext++;
                                extended = true;
                            } else if (left_matches == 0) {
                                left_matches++;
                                left_ext++;
                                extended = true;
                            }
                        }
                        
                        // 向右扩展 (query向右, ref向左)
                        if (q_pos + KMER_SIZE + right_ext < query.size() && 
                            r_pos - right_ext - 1 >= 0 &&
                            right_ext < MAX_EXTEND_LEN) {
                            char q_char = query[q_pos + KMER_SIZE + right_ext];
                            char r_char = ref[r_pos - right_ext - 1];
                            char rc_r_char = 
                                r_char == 'A' ? 'T' :
                                r_char == 'T' ? 'A' :
                                r_char == 'C' ? 'G' : 'C';
                            
                            if (q_char == rc_r_char) {
                                right_ext++;
                                extended = true;
                            } else if (right_matches == 0) {
                                right_matches++;
                                right_ext++;
                                extended = true;
                            }
                        }
                        
                        if (!extended) break;
                    }
                    
                    int total_len = KMER_SIZE + left_ext + right_ext;
                    int total_matches = KMER_SIZE - mismatch + left_ext + right_ext - left_matches - right_matches;
                    double similarity = (double)total_matches / total_len;
                    
                    if (total_len >= MIN_SEED_LEN && similarity >= MIN_SEED_SIM) {
                        seeds.push_back({
                            q_pos - left_ext, 
                            q_pos + KMER_SIZE + right_ext,
                            r_pos - right_ext, 
                            r_pos + KMER_SIZE + left_ext,
                            true,
                            (double)total_matches,
                            total_matches
                        });
                    }
                    
                    if (seeds.size() > MAX_SEEDS) break;
                }
                if (seeds.size() > MAX_SEEDS) break;
            }
            if (seeds.size() > MAX_SEEDS) break;
        }
    }
    
    // 移除重叠种子（保留高质量种子）
    if (seeds.empty()) return seeds;
    
    sort(seeds.begin(), seeds.end(), [](const Anchor& a, const Anchor& b) {
        if (a.q_start != b.q_start) return a.q_start < b.q_start;
        return a.score > b.score; // 优先保留高质量种子
    });
    
    vector<Anchor> unique_seeds;
    unique_seeds.push_back(seeds[0]);
    
    for (int i = 1; i < seeds.size(); i++) {
        if (seeds[i].q_start >= unique_seeds.back().q_end) {
            unique_seeds.push_back(seeds[i]);
        } else if (seeds[i].score > unique_seeds.back().score * 1.2) {
            // 如果质量显著更高，则替换
            unique_seeds.back() = seeds[i];
        }
    }
    
    return unique_seeds;
}

// 链式对齐（增加质量检查）
vector<vector<Anchor>> chain_anchors(vector<Anchor>& seeds) {
    if (seeds.empty()) return {};
    
    // 按query位置排序
    sort(seeds.begin(), seeds.end(), [](const Anchor& a, const Anchor& b) {
        return a.q_start < b.q_start;
    });
    
    // DP数组和路径跟踪
    vector<double> dp(seeds.size(), 0.0);
    vector<int> prev(seeds.size(), -1);
    
    for (int i = 0; i < seeds.size(); i++) {
        dp[i] = seeds[i].score;
    }
    
    // 动态规划 - 只连接高质量锚点
    for (int i = 0; i < seeds.size(); i++) {
        // 只考虑高质量种子
        if (seeds[i].matches < MIN_MATCHES_FOR_ANCHOR) continue;
        
        for (int j = 0; j < i; j++) {
            if (seeds[j].matches < MIN_MATCHES_FOR_ANCHOR) continue;
            if (seeds[i].is_reverse != seeds[j].is_reverse) continue;
            
            int q_gap = seeds[i].q_start - seeds[j].q_end;
            int r_gap;
            if (!seeds[i].is_reverse) {
                r_gap = seeds[i].r_start - seeds[j].r_end;
            } else {
                r_gap = seeds[j].r_start - seeds[i].r_end;
            }
            
            if (q_gap < 0 || r_gap < 0) continue;
            if (q_gap > MAX_CHAIN_GAP || r_gap > MAX_CHAIN_GAP) continue;
            
            // 间隙惩罚（更严格的线性惩罚）
            double gap_penalty = 0.05 * (q_gap + r_gap) + 0.1 * abs(q_gap - r_gap);
            double new_score = dp[j] + seeds[i].score - gap_penalty;
            
            if (new_score > dp[i]) {
                dp[i] = new_score;
                prev[i] = j;
            }
        }
    }
    
    // 构建链（只包含高质量链）
    vector<vector<Anchor>> chains;
    vector<bool> used(seeds.size(), false);
    
    while (true) {
        int best_idx = -1;
        double best_score = MIN_CHAIN_SCORE;
        
        for (int i = 0; i < seeds.size(); i++) {
            if (!used[i] && seeds[i].matches >= MIN_MATCHES_FOR_ANCHOR && dp[i] > best_score) {
                best_score = dp[i];
                best_idx = i;
            }
        }
        
        if (best_idx == -1) break;
        
        vector<Anchor> chain;
        int cur = best_idx;
        while (cur != -1 && !used[cur]) {
            chain.push_back(seeds[cur]);
            used[cur] = true;
            cur = prev[cur];
        }
        
        // 只保留足够长的链
        if (chain.size() >= 2 || (chain.size() == 1 && chain[0].matches >= MIN_MATCHES_FOR_ANCHOR * 2)) {
            reverse(chain.begin(), chain.end());
            chains.push_back(chain);
        }
    }
    
    // 按链得分排序
    sort(chains.begin(), chains.end(), [](const vector<Anchor>& a, const vector<Anchor>& b) {
        double sum_a = 0, sum_b = 0;
        for (const auto& anchor : a) sum_a += anchor.score;
        for (const auto& anchor : b) sum_b += anchor.score;
        return sum_a > sum_b;
    });
    
    return chains;
}

// 局部比对（更严格的质量控制）
vector<Anchor> banded_local_align(const string& ref, const string& query, 
                                 int q_start, int q_end, 
                                 int r_start, int r_end,
                                 bool is_reverse) {
    vector<Anchor> anchors;
    int len_q = q_end - q_start;
    int len_r = r_end - r_start;
    
    if (len_q <= 0 || len_r <= 0) return anchors;
    
    // 动态调整带宽
    int dynamic_band = min(max(len_q, len_r) / 2, BAND_WIDTH);
    dynamic_band = max(dynamic_band, 50);
    
    // 初始化DP矩阵
    vector<vector<int>> dp(len_q + 1, vector<int>(len_r + 1, 0));
    vector<vector<char>> path(len_q + 1, vector<char>(len_r + 1, ' '));
    
    // 填充DP矩阵（增加匹配得分）
    for (int i = 1; i <= len_q; i++) {
        int min_j = max(1, i - dynamic_band);
        int max_j = min(len_r, i + dynamic_band);
        
        for (int j = min_j; j <= max_j; j++) {
            char q_char = query[q_start + i - 1];
            char r_char;
            if (is_reverse) {
                r_char = ref[r_end - j];
            } else {
                r_char = ref[r_start + j - 1];
            }
            
            int match_score = (q_char == r_char) ? 4 : -3;  // 增加匹配/不匹配惩罚
            int diag = dp[i-1][j-1] + match_score;
            int up = dp[i-1][j] - 4;    // 增加空位惩罚
            int left = dp[i][j-1] - 4;  // 增加空位惩罚
            
            if (diag >= up && diag >= left && diag > 0) {
                dp[i][j] = diag;
                path[i][j] = 'D';
            } else if (up >= left && up > 0) {
                dp[i][j] = up;
                path[i][j] = 'U';
            } else if (left > 0) {
                dp[i][j] = left;
                path[i][j] = 'L';
            }
        }
    }
    
    // 回溯寻找最佳路径（要求最小匹配数）
    int best_i = 0, best_j = 0, best_score = 0;
    for (int i = 1; i <= len_q; i++) {
        for (int j = 1; j <= len_r; j++) {
            if (dp[i][j] > best_score) {
                best_score = dp[i][j];
                best_i = i;
                best_j = j;
            }
        }
    }
    
    // 计算最小所需分数
    int min_score = MIN_ANCHOR_LEN * 2;
    if (best_score < min_score) return anchors;
    
    // 回溯
    int i = best_i, j = best_j;
    vector<pair<int, int>> path_points;
    int total_matches = 0;
    
    while (i > 0 && j > 0 && dp[i][j] > 0) {
        path_points.push_back({i, j});
        if (path[i][j] == 'D') {
            char q_char = query[q_start + i - 1];
            char r_char;
            if (is_reverse) {
                r_char = ref[r_end - j];
            } else {
                r_char = ref[r_start + j - 1];
            }
            if (q_char == r_char) total_matches++;
            i--;
            j--;
        } else if (path[i][j] == 'U') {
            i--;
        } else if (path[i][j] == 'L') {
            j--;
        } else {
            break;
        }
    }
    
    reverse(path_points.begin(), path_points.end());
    
    if (path_points.empty()) return anchors;
    
    // 创建锚点（要求最小匹配数）
    if (total_matches < MIN_MATCHES_FOR_ANCHOR) return anchors;
    
    int start_i = path_points.front().first;
    int start_j = path_points.front().second;
    int end_i = path_points.back().first;
    int end_j = path_points.back().second;
    
    double similarity = (double)total_matches / path_points.size();
    if (similarity < MIN_SEED_SIM) return anchors;
    
    anchors.push_back({
        q_start + start_i - 1, 
        q_start + end_i,
        is_reverse ? (r_end - end_j) : (r_start + start_j - 1),
        is_reverse ? (r_end - start_j + 1) : (r_start + end_j),
        is_reverse,
        (double)path_points.size() * similarity,
        total_matches
    });
    
    return anchors;
}

// 间隙填补（增加质量检查）
vector<Anchor> global_gap_filling(const string& ref, const string& query, 
                                 const vector<vector<Anchor>>& chains) {
    vector<Anchor> all_anchors;
    
    for (const auto& chain : chains) {
        for (const auto& anchor : chain) {
            all_anchors.push_back(anchor);
        }
    }
    
    sort(all_anchors.begin(), all_anchors.end(), [](const Anchor& a, const Anchor& b) {
        return a.q_start < b.q_start;
    });
    
    vector<Anchor> filled_anchors;
    if (all_anchors.empty()) return filled_anchors;
    
    filled_anchors.push_back(all_anchors[0]);
    
    for (int i = 1; i < all_anchors.size(); i++) {
        const Anchor& prev = filled_anchors.back();
        const Anchor& curr = all_anchors[i];
        
        if (prev.is_reverse == curr.is_reverse) {
            int q_gap = curr.q_start - prev.q_end;
            int r_gap;
            if (!prev.is_reverse) {
                r_gap = curr.r_start - prev.r_end;
            } else {
                r_gap = prev.r_start - curr.r_end;
            }
            
            // 只填补合理大小的间隙
            if (q_gap > 0 && q_gap < MAX_LOCAL_GAP && abs(r_gap) < MAX_LOCAL_GAP) {
                vector<Anchor> new_anchors = banded_local_align(
                    ref, query,
                    prev.q_end, curr.q_start,
                    prev.is_reverse ? curr.r_end : prev.r_end,
                    prev.is_reverse ? prev.r_start : curr.r_start,
                    prev.is_reverse
                );
                
                // 只添加高质量锚点
                for (const auto& na : new_anchors) {
                    if (na.matches >= MIN_MATCHES_FOR_ANCHOR) {
                        filled_anchors.push_back(na);
                    }
                }
            }
        }
        
        // 只添加高质量锚点
        if (curr.matches >= MIN_MATCHES_FOR_ANCHOR) {
            filled_anchors.push_back(curr);
        }
    }
    
    return filled_anchors;
}

// 覆盖未比对区域（更严格的质量控制）
vector<Anchor> cover_remaining(const string& ref, const string& query, 
                              const vector<Anchor>& anchors) {
    vector<Anchor> final_anchors;
    int last_q_end = 0;
    
    vector<Anchor> sorted_anchors = anchors;
    sort(sorted_anchors.begin(), sorted_anchors.end(), [](const Anchor& a, const Anchor& b) {
        return a.q_start < b.q_start;
    });
    
    // 起始区域
    if (!sorted_anchors.empty() && sorted_anchors.front().q_start > 0) {
        vector<Anchor> new_anchors = banded_local_align(
            ref, query,
            0, sorted_anchors.front().q_start,
            0, ref.size(),
            false
        );
        
        for (const auto& na : new_anchors) {
            if (na.matches >= MIN_MATCHES_FOR_ANCHOR) {
                final_anchors.push_back(na);
                last_q_end = na.q_end;
            }
        }
    }
    
    for (const auto& anchor : sorted_anchors) {
        if (anchor.q_start > last_q_end) {
            vector<Anchor> new_anchors = banded_local_align(
                ref, query,
                last_q_end, anchor.q_start,
                0, ref.size(),
                false
            );
            
            for (const auto& na : new_anchors) {
                if (na.matches >= MIN_MATCHES_FOR_ANCHOR) {
                    final_anchors.push_back(na);
                    last_q_end = na.q_end;
                }
            }
        }
        
        if (anchor.matches >= MIN_MATCHES_FOR_ANCHOR) {
            final_anchors.push_back(anchor);
            last_q_end = anchor.q_end;
        }
    }
    
    // 结束区域
    if (last_q_end < query.size()) {
        vector<Anchor> new_anchors = banded_local_align(
            ref, query,
            last_q_end, query.size(),
            0, ref.size(),
            false
        );
        
        for (const auto& na : new_anchors) {
            if (na.matches >= MIN_MATCHES_FOR_ANCHOR) {
                final_anchors.push_back(na);
            }
        }
    }
    
    // 最终合并：连接相邻的高质量锚点
    vector<Anchor> merged_anchors;
    for (const auto& anchor : final_anchors) {
        if (merged_anchors.empty()) {
            merged_anchors.push_back(anchor);
            continue;
        }
        
        Anchor& last = merged_anchors.back();
        if (anchor.is_reverse == last.is_reverse &&
            anchor.q_start <= last.q_end + 5) {
            
            // 检查参考位置连续性
            bool ref_continuous = false;
            if (!last.is_reverse) {
                ref_continuous = abs(anchor.r_start - last.r_end) < 10;
            } else {
                ref_continuous = abs(anchor.r_end - last.r_start) < 10;
            }
            
            if (ref_continuous) {
                // 合并锚点
                last.q_end = anchor.q_end;
                if (!last.is_reverse) {
                    last.r_end = anchor.r_end;
                } else {
                    last.r_start = anchor.r_start;
                }
                last.score += anchor.score;
                last.matches += anchor.matches;
                continue;
            }
        }
        merged_anchors.push_back(anchor);
    }
    
    return merged_anchors;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <reference_file> <query_file>" << endl;
        return 1;
    }

    ifstream ref_file(argv[1]);
    ifstream query_file(argv[2]);
    
    if (!ref_file.is_open() || !query_file.is_open()) {
        cerr << "Error opening files!" << endl;
        return 1;
    }

    string ref, query;
    string line;
    
    // 读取参考序列（支持多行）
    while (getline(ref_file, line)) {
        // 忽略FASTA格式的头部
        if (line[0] != '>') {
            ref += line;
        }
    }
    
    // 读取查询序列（支持多行）
    while (getline(query_file, line)) {
        if (line[0] != '>') {
            query += line;
        }
    }

    clock_t start = clock();
    
    // 步骤1: 生成高质量种子
    cout << "Generating seeds..." << endl;
    vector<Anchor> seeds = generate_seeds(ref, query);
    cout << "Generated " << seeds.size() << " seeds | Time: " 
         << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    if (!seeds.empty()) {
        cout << "First seed: q(" << seeds[0].q_start << "-" << seeds[0].q_end 
             << ") r(" << seeds[0].r_start << "-" << seeds[0].r_end 
             << ") matches: " << seeds[0].matches << endl;
    }

    // 步骤2: 链式对齐
    cout << "Performing chaining..." << endl;
    vector<vector<Anchor>> chains = chain_anchors(seeds);
    cout << "Formed " << chains.size() << " chains | Time: " 
         << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    if (!chains.empty()) {
        cout << "Largest chain: " << chains[0].size() << " anchors, total matches: " 
             << chains[0][0].matches << endl;
    }

    // 步骤3: 全局间隙填补
    cout << "Global gap filling..." << endl;
    vector<Anchor> filled_anchors = global_gap_filling(ref, query, chains);
    cout << "After gap filling: " << filled_anchors.size() << " anchors | Time: " 
         << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    
    // 步骤4: 覆盖剩余区域
    cout << "Covering remaining regions..." << endl;
    vector<Anchor> final_anchors = cover_remaining(ref, query, filled_anchors);
    cout << "Final anchor count: " << final_anchors.size() << " | Time: " 
         << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    
    // 输出结果
    cout << "[";
    for (int i = 0; i < final_anchors.size(); i++) {
        const auto& a = final_anchors[i];
        cout << "(" << a.q_start << ", " << a.q_end
             << ", " << a.r_start << ", " << a.r_end << ")";
        if (i < final_anchors.size() - 1) 
            cout << ", ";
    }
    cout << "]" << endl;

    return 0;
}

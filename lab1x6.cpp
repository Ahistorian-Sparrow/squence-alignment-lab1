#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

struct RepeatInfo {
    int position; 
    int size;     
    int count;    
    bool is_reverse;
};

string reverse_complement(const string &s) {
    string rc;
    for (int i = s.size()-1; i >= 0; --i) {
        switch(s[i]) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            default: rc += 'X';
        }
    }
    return rc;
}

long long compute_hash(const string &s) {
    const int base = 4;
    const long long mod = 1e18+3;
    long long hash = 0;
    for(char c : s) {
        int val = 0;
        switch(c) {
            case 'A': val=0; break;
            case 'T': val=1; break;
            case 'C': val=2; break;
            case 'G': val=3; break;
        }
        hash = (hash * base + val) % mod;
    }
    return hash;
}

vector<RepeatInfo> find_repeats(const string &ref, const string &query) {
    vector<RepeatInfo> result;
    
    // 寻找首次差异位置
    int split_pos = 0;
    while (split_pos < min(ref.size(), query.size()) && 
           ref[split_pos] == query[split_pos]) split_pos++;
    
    // 预处理哈希表
    unordered_map<int, unordered_map<long long, vector<pair<int, bool>>>> hash_map;
    for(int L=1; L<=split_pos; L++) {
        for(int i=0; i<=split_pos-L; i++) {
            string s = ref.substr(i, L);
            long long h = compute_hash(s);
            hash_map[L][h].push_back({i, false});
            
            string rc = reverse_complement(s);
            long long h_rc = compute_hash(rc);
            hash_map[L][h_rc].push_back({i, true});
        }
    }

    int i = split_pos;
    while(i < query.size()) {
        int max_L = min(split_pos, (int)(query.size()-i));
        bool found = false;

        // 从最大长度开始尝试
        for(int L=max_L; L>=10; L--) {
            if(i+L > query.size()) continue;

            string q_sub = query.substr(i, L);
            string q_rc = reverse_complement(q_sub);
            long long h_q = compute_hash(q_sub);
            long long h_rc = compute_hash(q_rc);

            if(!hash_map.count(L)) continue;

            for(auto h : {h_q, h_rc}) {
                if(!hash_map[L].count(h)) continue;

                for(auto &[ref_start, is_rev] : hash_map[L][h]) {
                    string ref_sub = ref.substr(ref_start, L);
                    string expected = is_rev ? reverse_complement(ref_sub) : ref_sub;
                    string actual = (h == h_q) ? q_sub : q_rc;

                    if(expected != actual) continue;

                    // 计算基础重复次数
                    int count = 1;
                    int j = i + L;
                    while(j+L <= query.size()) {
                        if(query.substr(j, L) == expected) {
                            count++;
                            j += L;
                        } else break;
                    }

                    // 尾差异处理
                    int remaining = query.size() - j;
                    if(remaining > L) {
                      if(!is_rev){
                        for(int Ltail=0; Ltail<L; Ltail++) {
                            string overlap = query.substr(j, L - Ltail);
                            string tail = query.substr(j - Ltail, Ltail);
                            string ref_overlap;
                            string ref_tail;
                            
                            ref_overlap = ref.substr(ref_start, L - Ltail);
                            ref_tail = ref.substr(ref_start - Ltail, Ltail);

                            if(overlap == ref_overlap) {
                                if(tail == ref_tail){
                                    count++;
                                    j += L - Ltail;
                                    ref_start -= Ltail;
                                }
                                break;
                            }
                        }
                      }else{
                        for(int Ltail=0; Ltail<L; Ltail++) {
                            string overlap = query.substr(j, L - 2 * Ltail);
                            string tail = query.substr(j - Ltail, Ltail);
                            string q_overlap;
                            string q_tail;
                            
                            q_overlap = reverse_complement(ref.substr(ref_start + Ltail, L - 2 * Ltail));
                            q_tail = reverse_complement(ref.substr(ref_start, Ltail));

                            if(overlap == q_overlap) {
                                if(tail == q_tail){
                                    count++;
                                    L -= Ltail;
                                    j += L - Ltail;
                                    ref_start += Ltail;
                                }
                                break;
                            }
                        }
                      }
                    }

                    result.push_back({
                        ref_start + L,  // 调整后的位置
                        L,
                        count,      // 额外次数
                        is_rev
                    });
                    i = j;
                    found = true;
                    goto next_segment;
                }
            }
        }
        next_segment:
        if(!found) break;
    }

    return result;
}

int main() {
    // 测试数据
    string ref = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"; 
    string query = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"; 

    auto repeats = find_repeats(ref, query);
    
    for(auto &r : repeats) {
        cout << "position: " << r.position
             << " size: " << r.size
             << " count: " << r.count 
             << " reverse: " << (r.is_reverse ? "yes" : "no") << endl;
    }
}
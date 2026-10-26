/**
 * tools/fast_census.cpp — High-Performance Parallel C++ Census Engine for Isotypic Composition Transversals
 * 
 * Generates all set partitions of Z_n via Knuth's Algorithm S (Restricted Growth Strings)
 * using register bitmasks and parallel OpenMP threads to exhaustively verify transversal
 * counts against the closed-form formula T(n).
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <map>
#include <string>
#include <sstream>
#include <omp.h>
#include <getopt.h>

using namespace std;

// Closed-form formula calculation in C++
uint64_t factorial_u64(int n) {
    uint64_t r = 1;
    for (int i = 2; i <= n; ++i) r *= i;
    return r;
}

__int128 factorial_128(int n) {
    __int128 r = 1;
    for (int i = 2; i <= n; ++i) r *= i;
    return r;
}

__int128 multinomial_128(int n, const vector<int>& parts) {
    __int128 num = factorial_128(n);
    __int128 den = 1;
    for (int p : parts) den *= factorial_128(p);
    return num / den;
}

void get_partitions(int n, int max_p, vector<int>& cur, vector<vector<int>>& res) {
    if (n == 0) {
        res.push_back(cur);
        return;
    }
    for (int p = min(n, max_p); p >= 1; --p) {
        cur.push_back(p);
        get_partitions(n - p, p, cur, res);
        cur.pop_back();
    }
}

__int128 formula_shape_128(int n, const vector<int>& shape) {
    int k = (int)shape.size();
    map<int, int> counts;
    for (int p : shape) counts[p]++;
    
    vector<int> mults;
    for (auto const& [val, count] : counts) mults.push_back(count);
    __int128 c_lambda = multinomial_128(k, mults);

    vector<int> caps;
    for (int p : shape) caps.push_back(p - 1);
    __int128 cap_mult = multinomial_128(n - k, caps);

    return (__int128)n * c_lambda * c_lambda * cap_mult;
}

__int128 formula_total_128(int n) {
    vector<vector<int>> shapes;
    vector<int> cur;
    get_partitions(n, n, cur, shapes);
    __int128 tot = 0;
    for (auto const& sh : shapes) {
        tot += formula_shape_128(n, sh);
    }
    return tot;
}

string to_string_128(__int128 v) {
    if (v == 0) return "0";
    string s;
    bool neg = false;
    if (v < 0) { neg = true; v = -v; }
    while (v > 0) {
        s += (char)('0' + (v % 10));
        v /= 10;
    }
    if (neg) s += '-';
    reverse(s.begin(), s.end());
    return s;
}

// Generate unique compositions for a given shape
void generate_compositions(vector<int> shape, vector<vector<int>>& comps) {
    sort(shape.begin(), shape.end());
    do {
        comps.push_back(shape);
    } while (next_permutation(shape.begin(), shape.end()));
}

// Precomputed transversal bitmasks for a shape
struct ShapeMasks {
    vector<int> shape;
    vector<uint32_t> masks;
};

// Shape index mapping
struct ShapeRegistry {
    int n;
    vector<vector<int>> shapes;
    vector<ShapeMasks> shape_data;
    // Map from compact shape encoding to shape index
    map<vector<int>, int> shape_to_idx;

    ShapeRegistry(int n_) : n(n_) {
        vector<int> cur;
        get_partitions(n, n, cur, shapes);
        shape_data.resize(shapes.size());
        for (size_t i = 0; i < shapes.size(); ++i) {
            shape_to_idx[shapes[i]] = (int)i;
            shape_data[i].shape = shapes[i];

            // Generate compositions
            vector<vector<int>> comps;
            generate_compositions(shapes[i], comps);

            // Generate masks
            int k = (int)shapes[i].size();
            for (auto const& gamma : comps) {
                for (int o = 0; o < n; ++o) {
                    uint32_t m = 0;
                    int cur_v = o;
                    for (int r = 0; r < k; ++r) {
                        m |= (1U << cur_v);
                        cur_v = (cur_v + gamma[r]) % n;
                    }
                    shape_data[i].masks.push_back(m);
                }
            }
        }
    }
};

struct CensusResult {
    uint64_t total_partitions = 0;
    __int128 total_transversals = 0;
    uint64_t zero_solution_partitions = 0;
    vector<uint64_t> partitions_by_shape;
    vector<__int128> transversals_by_shape;
    vector<uint64_t> zeros_by_shape;

    void init(size_t num_shapes) {
        partitions_by_shape.assign(num_shapes, 0);
        transversals_by_shape.assign(num_shapes, 0);
        zeros_by_shape.assign(num_shapes, 0);
    }

    void merge(const CensusResult& o) {
        total_partitions += o.total_partitions;
        total_transversals += o.total_transversals;
        zero_solution_partitions += o.zero_solution_partitions;
        for (size_t i = 0; i < partitions_by_shape.size(); ++i) {
            partitions_by_shape[i] += o.partitions_by_shape[i];
            transversals_by_shape[i] += o.transversals_by_shape[i];
            zeros_by_shape[i] += o.zeros_by_shape[i];
        }
    }
};

// Worker function to evaluate a sub-tree of Algorithm S
void explore_rg_subtree(int n, int prefix_len, const vector<int>& prefix,
                        const ShapeRegistry& reg, CensusResult& res) {
    int a[32];
    int m[32]; // m[j] = max(a[0], ..., a[j])
    for (int i = 0; i < prefix_len; ++i) {
        a[i] = prefix[i];
        m[i] = (i == 0) ? a[0] : max(m[i - 1], a[i]);
    }

    int j = prefix_len;
    if (j == n) {
        // Evaluate single prefix
        int k = m[n - 1] + 1;
        uint32_t b_mask[32] = {0};
        int b_size[32] = {0};
        for (int i = 0; i < n; ++i) {
            b_mask[a[i]] |= (1U << i);
            b_size[a[i]]++;
        }
        vector<int> sh(b_size, b_size + k);
        sort(sh.rbegin(), sh.rend());
        int s_idx = reg.shape_to_idx.at(sh);

        const auto& masks = reg.shape_data[s_idx].masks;
        uint64_t sols = 0;
        for (uint32_t g : masks) {
            bool ok = true;
            for (int b = 0; b < k; ++b) {
                if ((g & b_mask[b]) == 0) {
                    ok = false;
                    break;
                }
            }
            if (ok) sols++;
        }

        res.total_partitions++;
        res.total_transversals += sols;
        res.partitions_by_shape[s_idx]++;
        res.transversals_by_shape[s_idx] += sols;
        if (sols == 0) {
            res.zero_solution_partitions++;
            res.zeros_by_shape[s_idx]++;
        }
        return;
    }

    // Iterative Algorithm S for elements j to n-1
    a[j] = 0;
    while (j >= prefix_len) {
        if (j == n - 1) {
            // Leaf level: iterate over all valid values for a[n-1]
            int limit = m[n - 2] + 1;
            for (int val = 0; val <= limit; ++val) {
                a[n - 1] = val;
                int k = max(m[n - 2], val) + 1;

                uint32_t b_mask[32] = {0};
                int b_size[32] = {0};
                for (int i = 0; i < n; ++i) {
                    b_mask[a[i]] |= (1U << i);
                    b_size[a[i]]++;
                }
                vector<int> sh(b_size, b_size + k);
                sort(sh.rbegin(), sh.rend());
                int s_idx = reg.shape_to_idx.at(sh);

                const auto& masks = reg.shape_data[s_idx].masks;
                uint64_t sols = 0;
                for (uint32_t g : masks) {
                    bool ok = true;
                    for (int b = 0; b < k; ++b) {
                        if ((g & b_mask[b]) == 0) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) sols++;
                }

                res.total_partitions++;
                res.total_transversals += sols;
                res.partitions_by_shape[s_idx]++;
                res.transversals_by_shape[s_idx] += sols;
                if (sols == 0) {
                    res.zero_solution_partitions++;
                    res.zeros_by_shape[s_idx]++;
                }
            }
            // Backtrack
            --j;
            while (j >= prefix_len && a[j] >= m[j - 1] + 1) {
                --j;
            }
            if (j >= prefix_len) {
                a[j]++;
                m[j] = max(m[j - 1], a[j]);
                ++j;
                a[j] = 0;
            }
        } else {
            m[j] = max(m[j - 1], a[j]);
            ++j;
            a[j] = 0;
        }
    }
}

// Collect prefixes of depth L
void generate_prefixes(int n, int target_depth, vector<int>& cur, vector<vector<int>>& prefixes) {
    if ((int)cur.size() == target_depth || (int)cur.size() == n) {
        prefixes.push_back(cur);
        return;
    }
    int j = (int)cur.size();
    int m_val = 0;
    for (int x : cur) m_val = max(m_val, x);
    int limit = m_val + 1;
    for (int val = 0; val <= limit; ++val) {
        cur.push_back(val);
        generate_prefixes(n, target_depth, cur, prefixes);
        cur.pop_back();
    }
}

int main(int argc, char* argv[]) {
    int n = 12;
    int num_threads = omp_get_max_threads();
    int prefix_depth = 4;

    static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"threads", required_argument, 0, 't'},
        {"depth", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "n:t:d:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'n': n = atoi(optarg); break;
            case 't': num_threads = atoi(optarg); break;
            case 'd': prefix_depth = atoi(optarg); break;
            case 'h':
            default:
                cout << "Usage: " << argv[0] << " [-n N] [-t THREADS] [-d PREFIX_DEPTH]\n";
                return 0;
        }
    }

    omp_set_num_threads(num_threads);

    cout << "=======================================================================\n";
    cout << " FAST TRANSVERSAL CENSUS ENGINE (C++ / OpenMP)\n";
    cout << " Lattice: Z_" << n << ", Threads: " << num_threads << ", Prefix Depth: " << prefix_depth << "\n";
    cout << "=======================================================================\n";

    auto t0 = chrono::high_resolution_clock::now();

    cout << "Initialising shape registry and precomputing polygon bitmasks...\n";
    ShapeRegistry reg(n);
    cout << "  Total shapes: " << reg.shapes.size() << "\n";

    // Prefixes for parallel distribution
    vector<vector<int>> prefixes;
    vector<int> initial_cur = {0};
    generate_prefixes(n, min(prefix_depth, n), initial_cur, prefixes);
    cout << "  Generated " << prefixes.size() << " work chunks for " << num_threads << " threads.\n";

    CensusResult global_res;
    global_res.init(reg.shapes.size());

    cout << "Executing parallel partition search...\n";
    auto t_search_start = chrono::high_resolution_clock::now();

    #pragma omp parallel
    {
        CensusResult local_res;
        local_res.init(reg.shapes.size());

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < prefixes.size(); ++i) {
            explore_rg_subtree(n, (int)prefixes[i].size(), prefixes[i], reg, local_res);
        }

        #pragma omp critical
        {
            global_res.merge(local_res);
        }
    }

    auto t1 = chrono::high_resolution_clock::now();
    double total_sec = chrono::duration<double>(t1 - t0).count();
    double search_sec = chrono::duration<double>(t1 - t_search_start).count();

    __int128 formula_val = formula_total_128(n);

    cout << "\n=======================================================================\n";
    cout << " CENSUS RESULTS FOR Z_" << n << "\n";
    cout << "=======================================================================\n";
    cout << "Total Partitions (B_" << n << "):      " << global_res.total_partitions << "\n";
    cout << "Total Rooted Transversals:  " << to_string_128(global_res.total_transversals) << "\n";
    cout << "Formula Prediction T(" << n << "):    " << to_string_128(formula_val) << "\n";
    cout << "Formula Exact Match?        " << (global_res.total_transversals == formula_val ? "EXACT MATCH (VERIFIED)" : "MISMATCH!") << "\n";
    cout << "Zero-Solution Partitions:   " << global_res.zero_solution_partitions 
         << " (" << fixed << setprecision(2) << (double)global_res.zero_solution_partitions / global_res.total_partitions * 100.0 << "%)\n";
    cout << "Search Time:                " << fixed << setprecision(3) << search_sec << " s\n";
    cout << "Total Elapsed Time:         " << fixed << setprecision(3) << total_sec << " s\n";
    cout << "Throughput:                 " << fixed << setprecision(1) << (double)global_res.total_partitions / search_sec << " partitions/sec\n";
    cout << "=======================================================================\n";

    return 0;
}

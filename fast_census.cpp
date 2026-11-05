/**
 * tools/fast_census.cpp — High-Performance Parallel C++ Census Engine for Isotypic Composition Transversals
 * 
 * Generates all set partitions of Z_n via Knuth's Algorithm S (Restricted Growth Strings)
 * using register bitmasks and parallel OpenMP threads to exhaustively verify transversal
 * counts against the closed-form formula T(n).
 *
 * Features:
 *  - Real-time periodic progress flushing with elapsed time, ETA, and running metrics.
 *  - Atomic checkpointing and automatic resumption across system interruptions.
 */

#include <iostream>
#include <fstream>
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
#include <unistd.h>

using namespace std;

const uint32_t CHECKPOINT_MAGIC = 0x43454E53; // "CENS"

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

struct ShapeData {
    vector<int> shape;
    vector<uint32_t> masks;
};

struct ShapeRegistry {
    int n;
    vector<vector<int>> shapes;
    map<vector<int>, int> shape_to_idx;
    vector<ShapeData> shape_data;

    ShapeRegistry(int n_) : n(n_) {
        vector<int> cur;
        get_partitions(n, n, cur, shapes);
        shape_data.resize(shapes.size());
        for (size_t i = 0; i < shapes.size(); ++i) {
            shape_to_idx[shapes[i]] = (int)i;
            shape_data[i].shape = shapes[i];

            // Compositions
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
        total_partitions = 0;
        total_transversals = 0;
        zero_solution_partitions = 0;
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

bool save_checkpoint(const string& filename, int n, int prefix_depth,
                     double elapsed_sec, const CensusResult& res,
                     const vector<uint8_t>& chunk_done) {
    string tmp_file = filename + ".tmp";
    ofstream out(tmp_file, ios::binary);
    if (!out.is_open()) return false;

    uint32_t magic = CHECKPOINT_MAGIC;
    out.write((const char*)&magic, sizeof(magic));
    out.write((const char*)&n, sizeof(n));
    out.write((const char*)&prefix_depth, sizeof(prefix_depth));
    out.write((const char*)&elapsed_sec, sizeof(elapsed_sec));

    out.write((const char*)&res.total_partitions, sizeof(res.total_partitions));
    out.write((const char*)&res.total_transversals, sizeof(res.total_transversals));
    out.write((const char*)&res.zero_solution_partitions, sizeof(res.zero_solution_partitions));

    uint64_t num_shapes = res.partitions_by_shape.size();
    out.write((const char*)&num_shapes, sizeof(num_shapes));
    out.write((const char*)res.partitions_by_shape.data(), num_shapes * sizeof(uint64_t));
    out.write((const char*)res.transversals_by_shape.data(), num_shapes * sizeof(__int128));
    out.write((const char*)res.zeros_by_shape.data(), num_shapes * sizeof(uint64_t));

    uint64_t num_chunks = chunk_done.size();
    out.write((const char*)&num_chunks, sizeof(num_chunks));
    out.write((const char*)chunk_done.data(), num_chunks * sizeof(uint8_t));

    out.close();
    return (rename(tmp_file.c_str(), filename.c_str()) == 0);
}

bool load_checkpoint(const string& filename, int expected_n, int expected_depth,
                     size_t expected_shapes, size_t expected_chunks,
                     double& elapsed_sec, CensusResult& res,
                     vector<uint8_t>& chunk_done) {
    ifstream in(filename, ios::binary);
    if (!in.is_open()) return false;

    uint32_t magic = 0;
    in.read((char*)&magic, sizeof(magic));
    if (magic != CHECKPOINT_MAGIC) return false;

    int file_n = 0, file_depth = 0;
    in.read((char*)&file_n, sizeof(file_n));
    in.read((char*)&file_depth, sizeof(file_depth));
    if (file_n != expected_n || file_depth != expected_depth) return false;

    in.read((char*)&elapsed_sec, sizeof(elapsed_sec));

    in.read((char*)&res.total_partitions, sizeof(res.total_partitions));
    in.read((char*)&res.total_transversals, sizeof(res.total_transversals));
    in.read((char*)&res.zero_solution_partitions, sizeof(res.zero_solution_partitions));

    uint64_t num_shapes = 0;
    in.read((char*)&num_shapes, sizeof(num_shapes));
    if (num_shapes != expected_shapes) return false;

    res.partitions_by_shape.resize(num_shapes);
    res.transversals_by_shape.resize(num_shapes);
    res.zeros_by_shape.resize(num_shapes);

    in.read((char*)res.partitions_by_shape.data(), num_shapes * sizeof(uint64_t));
    in.read((char*)res.transversals_by_shape.data(), num_shapes * sizeof(__int128));
    in.read((char*)res.zeros_by_shape.data(), num_shapes * sizeof(uint64_t));

    uint64_t num_chunks = 0;
    in.read((char*)&num_chunks, sizeof(num_chunks));
    if (num_chunks != expected_chunks) return false;

    chunk_done.resize(num_chunks);
    in.read((char*)chunk_done.data(), num_chunks * sizeof(uint8_t));

    return in.good();
}

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

    a[j] = 0;
    while (j >= prefix_len) {
        if (j == n - 1) {
            // Leaf level
            int limit = m[j - 1] + 1;
            for (int val = 0; val <= limit; ++val) {
                a[j] = val;
                int k = max(m[j - 1], val) + 1;

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
    int prefix_depth = -1; // default set based on n
    string checkpoint_path = "";
    bool fresh_start = false;
    bool enable_checkpoint = true;
    int progress_interval = -1;
    string json_out_path = "";

    static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"threads", required_argument, 0, 't'},
        {"depth", required_argument, 0, 'd'},
        {"checkpoint", required_argument, 0, 'c'},
        {"interval", required_argument, 0, 'i'},
        {"fresh", no_argument, 0, 'f'},
        {"json-out", required_argument, 0, 'j'},
        {"no-checkpoint", no_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, "n:t:d:c:i:fj:h", long_options, &opt_index)) != -1) {
        switch (opt) {
            case 'n': n = atoi(optarg); break;
            case 't': num_threads = atoi(optarg); break;
            case 'd': prefix_depth = atoi(optarg); break;
            case 'c': checkpoint_path = optarg; break;
            case 'i': progress_interval = atoi(optarg); break;
            case 'f': fresh_start = true; break;
            case 'j': json_out_path = optarg; break;
            case 0:
                if (string(long_options[opt_index].name) == "no-checkpoint") {
                    enable_checkpoint = false;
                }
                break;
            case 'h':
            default:
                cout << "Usage: " << argv[0] << " [options]\n"
                     << "  -n N                  Lattice size Zn (default: 12)\n"
                     << "  -t THREADS            OpenMP worker threads (default: max)\n"
                     << "  -d DEPTH              Algorithm S prefix depth (default: 4 for n<=12, 9 for n>=16)\n"
                     << "  -c, --checkpoint PATH Checkpoint file path\n"
                     << "  -i, --interval N      Progress log & checkpoint interval in chunks\n"
                     << "  -f, --fresh           Ignore any existing checkpoint and start fresh\n"
                     << "  -j, --json-out PATH   Export per-shape census to JSON\n"
                     << "      --no-checkpoint   Disable checkpoint file generation\n"
                     << "  -h, --help            Show help\n";
                return 0;
        }
    }

    if (prefix_depth < 0) {
        prefix_depth = (n >= 15) ? 9 : 4;
    }

    if (checkpoint_path.empty()) {
        checkpoint_path = "census_z" + to_string(n) + "_d" + to_string(prefix_depth) + "_checkpoint.dat";
    }

    omp_set_num_threads(num_threads);

    cout << "=======================================================================\n";
    cout << " FAST TRANSVERSAL CENSUS ENGINE (C++ / OpenMP)\n";
    cout << " Lattice: Z_" << n << ", Threads: " << num_threads << ", Prefix Depth: " << prefix_depth << "\n";
    if (enable_checkpoint) {
        cout << " Checkpoint File: " << checkpoint_path << "\n";
    }
    cout << "=======================================================================\n";

    auto t0 = chrono::high_resolution_clock::now();

    cout << "Initialising shape registry and precomputing polygon bitmasks...\n" << flush;
    ShapeRegistry reg(n);
    cout << "  Total shapes: " << reg.shapes.size() << "\n" << flush;

    // Prefixes for parallel distribution
    vector<vector<int>> prefixes;
    vector<int> initial_cur = {0};
    generate_prefixes(n, min(prefix_depth, n), initial_cur, prefixes);
    size_t total_chunks = prefixes.size();
    cout << "  Generated " << total_chunks << " work chunks for " << num_threads << " threads.\n" << flush;

    if (progress_interval < 0) {
        progress_interval = max(1, (int)(total_chunks / 100)); // every ~1%
    }

    CensusResult global_res;
    global_res.init(reg.shapes.size());
    vector<uint8_t> chunk_done(total_chunks, 0);
    double prior_elapsed_sec = 0.0;
    size_t resumed_chunks = 0;

    if (enable_checkpoint && !fresh_start) {
        if (load_checkpoint(checkpoint_path, n, prefix_depth, reg.shapes.size(),
                            total_chunks, prior_elapsed_sec, global_res, chunk_done)) {
            for (uint8_t d : chunk_done) if (d) resumed_chunks++;
            cout << ">> Found valid checkpoint! Resuming from " << resumed_chunks << " / " << total_chunks
                 << " chunks (" << fixed << setprecision(1) << (double)resumed_chunks / total_chunks * 100.0 << "% completed).\n"
                 << "   Prior elapsed: " << fixed << setprecision(1) << prior_elapsed_sec / 60.0 << " min, "
                 << "Partitions so far: " << global_res.total_partitions << "\n" << flush;
        }
    }

    cout << "Executing parallel partition search...\n" << flush;
    auto t_search_start = chrono::high_resolution_clock::now();
    uint64_t completed_chunks = resumed_chunks;

    #pragma omp parallel
    {
        CensusResult local_res;
        local_res.init(reg.shapes.size());

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < total_chunks; ++i) {
            if (chunk_done[i]) continue;

            local_res.init(reg.shapes.size());
            explore_rg_subtree(n, (int)prefixes[i].size(), prefixes[i], reg, local_res);

            #pragma omp critical
            {
                global_res.merge(local_res);
                chunk_done[i] = 1;
                completed_chunks++;

                if (completed_chunks % progress_interval == 0 || completed_chunks == total_chunks) {
                    auto now = chrono::high_resolution_clock::now();
                    double this_elapsed = chrono::duration<double>(now - t_search_start).count();
                    double total_run_sec = prior_elapsed_sec + this_elapsed;
                    double speed_chunks = (this_elapsed > 0) ? (double)(completed_chunks - resumed_chunks) / this_elapsed : 0;
                    double eta_sec = (speed_chunks > 0) ? (double)(total_chunks - completed_chunks) / speed_chunks : 0;
                    double pct = (double)completed_chunks / total_chunks * 100.0;
                    double zero_pct = (global_res.total_partitions > 0)
                                      ? (double)global_res.zero_solution_partitions / global_res.total_partitions * 100.0
                                      : 0.0;

                    cout << "[" << fixed << setprecision(1) << setw(5) << pct << "%] "
                         << "Chunks: " << setw(5) << completed_chunks << "/" << total_chunks
                         << " | Elapsed: " << setprecision(1) << setw(5) << total_run_sec / 60.0 << "m"
                         << " | ETA: " << setprecision(1) << setw(5) << eta_sec / 60.0 << "m"
                         << " | Partitions: " << global_res.total_partitions
                         << " | Zeros: " << setprecision(2) << zero_pct << "%"
                         << endl << flush;

                    if (enable_checkpoint && completed_chunks < total_chunks) {
                        save_checkpoint(checkpoint_path, n, prefix_depth, total_run_sec, global_res, chunk_done);
                    }
                }
            }
        }
    }

    auto t1 = chrono::high_resolution_clock::now();
    double search_sec = prior_elapsed_sec + chrono::duration<double>(t1 - t_search_start).count();
    double total_sec = chrono::duration<double>(t1 - t0).count() + prior_elapsed_sec;

    // Remove checkpoint on successful completion
    if (enable_checkpoint) {
        unlink(checkpoint_path.c_str());
    }

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
    cout << "=======================================================================\n" << flush;

    if (!json_out_path.empty()) {
        ofstream jout(json_out_path);
        jout << "{\n";
        jout << "  \"n\": " << n << ",\n";
        jout << "  \"bell_number\": " << global_res.total_partitions << ",\n";
        jout << "  \"total_transversals\": \"" << to_string_128(global_res.total_transversals) << "\",\n";
        jout << "  \"krell_number\": \"" << to_string_128(global_res.total_transversals / n) << "\",\n";
        jout << "  \"zero_partitions\": " << global_res.zero_solution_partitions << ",\n";
        jout << "  \"shapes\": [\n";
        for (size_t i = 0; i < reg.shapes.size(); ++i) {
            const auto& sh = reg.shapes[i];
            int k = (int)sh.size();
            uint64_t p_lam = global_res.partitions_by_shape[i];
            __int128 t_lam = global_res.transversals_by_shape[i];
            __int128 s_star = t_lam / n;
            uint64_t z_lam = global_res.zeros_by_shape[i];
            uint64_t solv = p_lam - z_lam;
            double mean_k = (p_lam > 0) ? (double)s_star / p_lam : 0.0;

            map<int, int> mults;
            for (int x : sh) mults[x]++;
            __int128 c_lam = factorial_128(k);
            for (auto const& kv : mults) c_lam /= factorial_128(kv.second);

            vector<int> caps;
            for (int x : sh) caps.push_back(x - 1);
            __int128 cap_mult = multinomial_128(n - k, caps);

            jout << "    {\n";
            jout << "      \"shape\": [";
            for (size_t j = 0; j < sh.size(); ++j) {
                jout << sh[j] << (j + 1 < sh.size() ? ", " : "");
            }
            jout << "],\n";
            jout << "      \"k\": " << k << ",\n";
            jout << "      \"p_lambda\": " << p_lam << ",\n";
            jout << "      \"c_lambda\": \"" << to_string_128(c_lam) << "\",\n";
            jout << "      \"capacity\": \"" << to_string_128(cap_mult) << "\",\n";
            jout << "      \"striding\": \"" << to_string_128(s_star) << "\",\n";
            jout << "      \"mean\": " << fixed << setprecision(4) << mean_k << ",\n";
            jout << "      \"zero_count\": " << z_lam << ",\n";
            jout << "      \"solvable\": " << solv << "\n";
            jout << "    }" << (i + 1 < reg.shapes.size() ? "," : "") << "\n";
        }
        jout << "  ]\n";
        jout << "}\n";
        jout.close();
        cout << "Census shape data written to: " << json_out_path << "\n";
    }

    return 0;
}

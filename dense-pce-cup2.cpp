#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <filesystem>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cassert>
#include <deque>
#include "EBBkC/src/edge_oriented.h"

// Globals required by EBBkC in library mode
int K = 0; // target clique size
int L = 1;   // early-termination t-plex parameter (strictest setting)
unsigned long long N = 0ULL; // number of vertices in the graph

// Ablation gating via command-line flags
// Defaults: full FPCE path (order+edge+Turan)
static bool g_enable_order_bound = true;
static bool g_enable_edge_bound  = true;
static bool g_enable_turan       = true;

enum class CupMode {
    Off,
    Stats,
    Compare,
    Assist,
    StrictAssist
};

static CupMode g_cup_mode = CupMode::Assist;
static bool g_enable_cup_auto_disable = true;

static constexpr int CUP_MIN_WINDOW_SIZE = 4;
static constexpr int CUP_MAX_WINDOW_SIZE = 256;
static constexpr int CUP_GREEN_MAX_CHANGED_FOOTPRINT = 32;
static constexpr int CUP_YELLOW_MAX_CHANGED_FOOTPRINT = 96;
static constexpr int CUP_HARD_MAX_CHANGED_FOOTPRINT = 128;
static constexpr int CUP_HIGH_CHANGED_FOOTPRINT = 64;
static constexpr int CUP_RED_SMALL_WINDOW_MAX = 4;
static constexpr int CUP_CHANGED_TO_WINDOW_RATIO_LIMIT = 8;
static constexpr int CUP_CHANGED_TO_SHORTLIST_RATIO_LIMIT = 12;
static constexpr int CUP_YELLOW_MAX_RELEVANT_UNCLES = 2;
static constexpr int CUP_ROLLING_WINDOW_SIZE = 64;
static constexpr int CUP_ROLLING_MIN_SAMPLES = 16;
static constexpr int CUP_MAX_FRAME_CANDIDATES = 256;
static constexpr int CUP_MAX_SHORTLIST_PER_UNCLE = 64;

static const char* cup_mode_name(CupMode mode) {
    switch (mode) {
        case CupMode::Stats:   return "stats";
        case CupMode::Compare: return "compare";
        case CupMode::Assist:  return "assist";
        case CupMode::StrictAssist: return "strict-assist";
        default:               return "off";
    }
}

// counting sort with tracklist
class Graph {
public:
    Graph() = default;
    // Temporary adjacency for building; cleared after CSR finalize
    std::vector<std::vector<int>> adj_map;
    // CSR storage
    std::vector<uint32_t> csr_offsets;   // size n+1 (total nodes + 1)
    std::vector<uint32_t> csr_neighbors; // size m (total edges)

    // CSR view for EBBkC zero-copy entry (from PKT/BSR binaries)
    std::vector<uint32_t> csr_node_off;  // size n+1
    std::vector<int>      csr_edge_dst;  // size m
    uint32_t csr_n = 0; // number of nodes in the induced CSR
    uint32_t csr_m = 0; // number of edges in the induced CSR
    // Optional reverse map when CSR is core-shrunk: new-id -> original id
    std::vector<int> csr_rev_map;

    int total_nodes = 0 ;
    // --- Order-bound state ---
    int degeneracy = 0;                 // ξ_G  (max core)
    std::vector<int> core_numbers;      // per-vertex core number (optional for later bounds)

    void add_edge(int u, int v) {
        if (u >= (int)adj_map.size() || v >= (int)adj_map.size()) {
            adj_map.resize(std::max(u, v) + 1);
        }
        adj_map[u].push_back(v); // add edge u-v to the adjacency list
        adj_map[v].push_back(u); // add edge v-u to the adjacency list
    }

    void finalize_adjacency() { // converts adj_map into CSR format
        // Sort and unique temporary adjacency
        for (auto& nbrs : adj_map) { // [[], [1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
            std::sort(nbrs.begin(), nbrs.end()); // sort each vertex's neighbor list
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end()); // remove duplicate edges
        }
        // Build CSR offset array
        int n = static_cast<int>(adj_map.size());
        total_nodes = std::max(total_nodes, n);
        csr_offsets.assign(static_cast<size_t>(total_nodes + 1), 0u); 
        for (int u = 0; u < total_nodes; ++u) { 
            csr_offsets[u + 1] = csr_offsets[u] + static_cast<uint32_t>(u < n ? adj_map[u].size() : 0u); // cumulative vertex count
        }
        // Populate CSR neighbor array
        uint32_t m = csr_offsets[total_nodes]; // number of edges in the graph
        csr_neighbors.clear(); // clear the neighbor array
        csr_neighbors.resize(m); // resize the neighbor array to the number of edges in the graph
        for (int u = 0; u < total_nodes; ++u) {
            uint32_t off = csr_offsets[u]; // offset of the neighbor list for node u
            if (u < n) {
                const auto& nbrs = adj_map[u]; // neighbor list for node u
                for (size_t k = 0; k < nbrs.size(); ++k) csr_neighbors[off + static_cast<uint32_t>(k)] = static_cast<uint32_t>(nbrs[k]); // add the neighbor to the neighbor array
            }
        }
        // Free temporary adjacency to save memory
        adj_map.clear(); // clear the adjacency list
        adj_map.shrink_to_fit();
    }

    // Matula & Beck k-core in O(n+m) over CSR to get core_numbers and ξ_G
    void compute_core_numbers_from_csr() {
        const int n = total_nodes;
        core_numbers.assign(n, 0);

        if (n == 0) { degeneracy = 0; return; }

        std::vector<int> deg(n);
        int maxdeg = 0;
        for (int u = 0; u < n; ++u) {
            int du = static_cast<int>(csr_offsets[u + 1] - csr_offsets[u]);
            deg[u] = du;
            if (du > maxdeg) maxdeg = du;
        }

        std::vector<int> bin(maxdeg + 1, 0);
        for (int u = 0; u < n; ++u) bin[deg[u]]++;

        int start = 0;
        for (int d = 0; d <= maxdeg; ++d) {
            int num = bin[d];
            bin[d] = start;
            start += num;
        }

        std::vector<int> pos(n), vert(n);
        for (int u = 0; u < n; ++u) {
            pos[u] = bin[deg[u]];
            vert[pos[u]] = u;
            bin[deg[u]]++;
        }
        for (int d = maxdeg; d > 0; --d) bin[d] = bin[d - 1];
        bin[0] = 0;

        int degen = 0;
        for (int i = 0; i < n; ++i) {
            int v = vert[i];
            degen = std::max(degen, deg[v]);
            core_numbers[v] = deg[v];
            // decrease-neighbor step
            for (int64_t it = csr_offsets[v]; it < csr_offsets[v + 1]; ++it) {
                int u = csr_neighbors[it];
                if (deg[u] > deg[v]) {
                    int du = deg[u];
                    int pu = pos[u];
                    int pw = bin[du];
                    int w  = vert[pw];
                    if (u != w) {
                        pos[u] = pw;        vert[pu] = w;
                        pos[w] = pu;        vert[pw] = u;
                    }
                    bin[du]++;             deg[u]--;
                }
            }
        }
        degeneracy = degen;
    }

    // μ(θ, ξG): min of the two corollaries (exactly what we used in mod-edge-order)
    int compute_order_bound(double theta, int xi) const {
        if (theta <= 0.0) return total_nodes;
        const long double th = static_cast<long double>(theta);
        const long double x  = static_cast<long double>(xi);
        const long double eps = 1e-12L;

        // Corollary 1: |S| ≤ ⌊ 2 ξG / θ ⌋
        int b1 = static_cast<int>(std::floor((2.0L * x) / th + eps));

        // Corollary 2: if θ > ξG / (ξG + 1), |S| ≤ ⌊ 1 / (1 − ξG / ((ξG + 1) θ)) ⌋
        int b2 = std::numeric_limits<int>::max();
        long double thresh = x / (x + 1.0L);
        if (th > thresh + eps) {
            long double denom = 1.0L - x / ((x + 1.0L) * th);
            if (denom <= 0.0L) denom = eps;
            b2 = static_cast<int>(std::floor(1.0L / denom + eps));
        }
        return std::min(b1, b2);
    }

    bool has_edge(int u, int v) const {
        if (u < 0 || v < 0 || u >= total_nodes) return false; // u and v are out of bounds
        if (csr_offsets.empty()) return false;
        uint32_t begin = csr_offsets[static_cast<size_t>(u)]; // begin of the neighbor list for node u
        uint32_t end = csr_offsets[static_cast<size_t>(u + 1)]; // end of the neighbor list for node u
        uint32_t vv = static_cast<uint32_t>(v); // type-casting
        return std::binary_search(csr_neighbors.begin() + begin, csr_neighbors.begin() + end, vv); // binary search is O(log degree) instead of O(degree).
    }

    void read_graph_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        int current_vertex = 0;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line == "[EOF]") break; // end of file
            total_nodes++;
            if (line.empty()) {
                current_vertex++; // skip empty lines
                continue;
            }

            std::stringstream ss(line);
            int neighbor;
            while (ss >> neighbor) {
                add_edge(current_vertex, neighbor); // add edge current_vertex-neighbor to the adjacency list
            }
            current_vertex++; // increment the current vertex
        }
        // Ensure isolated trailing vertices exist as empty adjacency lists
        if (adj_map.size() < static_cast<size_t>(total_nodes)) {
            adj_map.resize(static_cast<size_t>(total_nodes)); // resize the adjacency list to the number of nodes in the graph
        }
        finalize_adjacency(); // convert the adjacency list to the CSR format
    }

    // Load graph from BBkC binaries and build only our own sorted CSR.
    void read_graph_from_bbkcbinaries(const std::string& dir, const std::string& original_graph_path) {
        // file setup
        std::filesystem::path deg_path = std::filesystem::path(dir) / "b_degree.bin"; // path to the degree file
        std::filesystem::path adj_path = std::filesystem::path(dir) / "b_adj.bin"; // path to the adjacency file
        std::ifstream deg_file(deg_path, std::ios::binary); // open the degree file
        std::ifstream adj_file(adj_path, std::ios::binary); // open the adjacency file
        if (!deg_file.is_open() || !adj_file.is_open()) {
            std::cerr << "Error: Could not open BBkC binaries in directory " << dir << std::endl;
            return;
        }

        // Read degress file header
        // b_degree.bin: [ui tt][ui n][ui m][ui degree[0..n-1]]
        uint32_t tt=0, n=0, m=0; // tt (type size), n (vertices), m (edges)
        deg_file.read(reinterpret_cast<char*>(&tt), sizeof(uint32_t)); // read the type size
        deg_file.read(reinterpret_cast<char*>(&n), sizeof(uint32_t)); // read the number of vertices
        deg_file.read(reinterpret_cast<char*>(&m), sizeof(uint32_t)); // read the number of edges
        if (!deg_file.good()) { std::cerr << "Error: Failed reading degree header\n"; return; } // check if the degree file is good
        if (tt != sizeof(uint32_t)) std::cerr << "Warning: ui size header " << tt << "\n"; // check if the type size is correct

        // Read degree array
        std::vector<uint32_t> degrees(n); // degree array
        if (n>0) {
            deg_file.read(reinterpret_cast<char*>(degrees.data()), static_cast<std::streamsize>(n*sizeof(uint32_t))); // read the degree array
            if (!deg_file.good()) { std::cerr << "Error: Failed reading degree array\n"; return; } // check if the degree array is good
        }

        // Build only our enumerator CSR (sorted rows for binary_search in has_edge())
        total_nodes = static_cast<int>(n); // set the number of nodes in the graph
        csr_offsets.assign(n+1, 0); // initialize the CSR offsets
        for (uint32_t i=0;i<n;++i) csr_offsets[i+1] = csr_offsets[i] + degrees[i];
        csr_neighbors.resize(csr_offsets[n]); // resize the CSR neighbors to the number of edges in the graph

        // Stream neighbors and sort each row slice
        for (uint32_t u=0; u<n; ++u) {
            const uint32_t du  = degrees[u]; // degree of node u
            const uint32_t off = csr_offsets[u]; // offset of the neighbor list for node u
            const uint32_t end = off + du; // end of the neighbor list for node u
            for (uint32_t k=0; k<du; ++k) {
                uint32_t v; // neighbor of node u
                adj_file.read(reinterpret_cast<char*>(&v), sizeof(uint32_t));
                if (!adj_file.good()) { std::cerr << "Error: Unexpected EOF in b_adj.bin (u="<<u<<")\n"; return; } // check if the adjacency file is good
                csr_neighbors[off + k] = v; // add the neighbor to the CSR neighbors
            }
            std::sort(csr_neighbors.begin() + off, csr_neighbors.begin() + end); // sort the neighbor list for node u
        }
    }
};

class PseudoCliqueEnumerator {
private:
    Graph& graph;
    float theta;
    int min_size;
    int max_size;
    float theta_P;
    int total_nodes_in_P;
    int total_edges_in_P;
    int iter_count = 0;

    // --- Edge-bound state ---
    int lhs_required_edges = 0;             // ceil(theta * l*(l-1)/2)   [fixed for a run]
    std::vector<int> core_hist;             // histogram of core numbers among vertices in P
    int cur_min_core = INT_MAX;             // c(P) = min core within P
    int max_core_seen = 0;                  // histogram capacity

    struct Track {
        uint8_t  inP;     // 0/1 flag
        uint16_t degP;    // ∈ [0, μ] 
        uint16_t degNP;   // ∈ [0, μ] 
    };

    // Intrusive bucketed lists for degrees 0..max_size
    std::vector<int32_t> headP, headNP;   // heads per degree (-1 if empty)
    std::vector<int32_t> prevP, nextP;    // per-vertex links in P buckets
    std::vector<int32_t> prevNP, nextNP;  // per-vertex links in NP buckets
    std::vector<int32_t> countP, countNP; // counts per degree

    int maxNPdeg = 0; // maximum non-empty degree bucket in NP (outside-P only)
    std::vector<int> pseudo_cliques_count;
    std::vector<Track> tracks; //
    int order_ub = -1;   // μ (order bound); -1 means "disabled"

    enum class ChildBucketKind : uint8_t {
        None = 0,
        Low = 1,
        EqDelta = 2,
        EqDeltaPlus1 = 3
    };

    struct CupGenerated {
        std::vector<int> accepted_children;
        std::vector<int> fallback_candidates;
        int relevant_uncle_count = 0;
        int resolved_window_count = 0;
        int direct_accept_count = 0;
        int fallback_count = 0;
        int estimated_shortlist_size = 0;
        int handled_low_count = 0;
        int handled_eq_delta_count = 0;
        int handled_eq_delta_plus1_count = 0;
        bool covers_low_bucket = false;
        bool covers_eq_delta_bucket = false;
        bool covers_eq_delta_plus1_bucket = false;

        void reset() {
            accepted_children.clear();
            fallback_candidates.clear();
            relevant_uncle_count = 0;
            resolved_window_count = 0;
            direct_accept_count = 0;
            fallback_count = 0;
            estimated_shortlist_size = 0;
            handled_low_count = 0;
            handled_eq_delta_count = 0;
            handled_eq_delta_plus1_count = 0;
            covers_low_bucket = false;
            covers_eq_delta_bucket = false;
            covers_eq_delta_plus1_bucket = false;
        }
    };

    struct NodeScratch {
        std::vector<int> children;
        std::vector<int> baseline_children;
        std::vector<int> compare_children;
        CupGenerated cup_generated;

        void reset() {
            children.clear();
            baseline_children.clear();
            compare_children.clear();
            cup_generated.reset();
        }
    };

    struct CupFrame {
        std::vector<int> low_bucket_candidates;
        std::vector<int> eq_delta_candidates;
        std::vector<int> eq_delta_plus1_candidates;
        int source_vertex = -1;
        int previous_uncle_count = 0;
        int contributing_uncle_count = 0;
        int changed_footprint = 0;
        int guard_summary = INT_MAX;

        void reset() {
            low_bucket_candidates.clear();
            eq_delta_candidates.clear();
            eq_delta_plus1_candidates.clear();
            source_vertex = -1;
            previous_uncle_count = 0;
            contributing_uncle_count = 0;
            changed_footprint = 0;
            guard_summary = INT_MAX;
        }
    };

    struct RelevantUncles {
        const CupFrame* frame = nullptr;
        int relevant_uncle_count = 0;
        int estimated_candidate_count = 0;

        bool empty() const {
            return frame == nullptr || relevant_uncle_count == 0 || estimated_candidate_count == 0;
        }
    };

    struct CupNodeContext {
        bool attempted = false;
        bool active = false;
        bool auto_disabled = false;
        int exact_window_size = 0;
        int relevant_uncle_count = 0;
        int changed_footprint = 0;
        int recovered_from_uncles = 0;
        int uncle_candidates_needing_fallback = 0;
        int direct_accept_count = 0;
        int resolved_window_count = 0;
        bool covered_low_bucket = false;
        bool covered_eq_delta_bucket = false;
        bool covered_eq_delta_plus1_bucket = false;
    };

    struct RollingCupMetrics {
        std::deque<int> changed_footprints;
        std::deque<int> shortlist_sizes;
        std::deque<int> recovered_children;
        std::deque<int> fallback_counts;
        std::deque<int> direct_accept_counts;
        std::deque<int> covered_low;
        std::deque<int> covered_eq_delta;
        std::deque<int> covered_any;

        long long sum_changed_footprints = 0;
        long long sum_shortlist_sizes = 0;
        long long sum_recovered_children = 0;
        long long sum_fallback_counts = 0;
        long long sum_direct_accept_counts = 0;
        long long sum_covered_low = 0;
        long long sum_covered_eq_delta = 0;
        long long sum_covered_any = 0;

        int size() const {
            return static_cast<int>(changed_footprints.size());
        }
    };

    struct CupStats {
        unsigned long long nodes_with_previous_uncles = 0ULL;
        unsigned long long nodes_without_previous_uncles = 0ULL;
        unsigned long long assist_attempts = 0ULL;
        unsigned long long assist_successes = 0ULL;
        unsigned long long fallback_validations = 0ULL;
        unsigned long long direct_accepts = 0ULL;
        unsigned long long total_shortlist_size = 0ULL;
        unsigned long long total_child_window_size = 0ULL;
        unsigned long long total_relevant_uncle_count = 0ULL;
        unsigned long long total_changed_footprint = 0ULL;
        unsigned long long baseline_children_recovered = 0ULL;
        unsigned long long auto_disables = 0ULL;
        unsigned long long covered_low_buckets = 0ULL;
        unsigned long long covered_eq_delta_buckets = 0ULL;
        unsigned long long covered_eq_delta_plus1_buckets = 0ULL;
        unsigned long long compare_parity_checks = 0ULL;
    };

    std::vector<NodeScratch> scratch_stack;
    std::vector<CupFrame> cup_frame_stack;
    std::vector<CupFrame*> parent_cup_frame_stack;
    std::vector<uint32_t> cup_generated_mark;
    std::vector<uint32_t> cup_handled_mark;
    std::vector<uint32_t> baseline_mark;
    uint32_t cup_generated_token = 1;
    uint32_t cup_handled_token = 1;
    uint32_t baseline_mark_token = 1;
    bool cup_runtime_disabled = false;
    CupStats cup_stats;
    RollingCupMetrics rolling_cup_metrics;

    inline static double ratio_or_zero(long long num, long long den) {
        if (den <= 0) return 0.0;
        return static_cast<double>(num) / static_cast<double>(den);
    }

    void push_rolling_cup_metrics(const CupNodeContext& cup_ctx, const CupGenerated& generated) {
        auto pop_front_if_full = [&]() {
            if (rolling_cup_metrics.size() < CUP_ROLLING_WINDOW_SIZE) return;

            rolling_cup_metrics.sum_changed_footprints -= rolling_cup_metrics.changed_footprints.front();
            rolling_cup_metrics.changed_footprints.pop_front();

            rolling_cup_metrics.sum_shortlist_sizes -= rolling_cup_metrics.shortlist_sizes.front();
            rolling_cup_metrics.shortlist_sizes.pop_front();

            rolling_cup_metrics.sum_recovered_children -= rolling_cup_metrics.recovered_children.front();
            rolling_cup_metrics.recovered_children.pop_front();

            rolling_cup_metrics.sum_fallback_counts -= rolling_cup_metrics.fallback_counts.front();
            rolling_cup_metrics.fallback_counts.pop_front();

            rolling_cup_metrics.sum_direct_accept_counts -= rolling_cup_metrics.direct_accept_counts.front();
            rolling_cup_metrics.direct_accept_counts.pop_front();

            rolling_cup_metrics.sum_covered_low -= rolling_cup_metrics.covered_low.front();
            rolling_cup_metrics.covered_low.pop_front();

            rolling_cup_metrics.sum_covered_eq_delta -= rolling_cup_metrics.covered_eq_delta.front();
            rolling_cup_metrics.covered_eq_delta.pop_front();

            rolling_cup_metrics.sum_covered_any -= rolling_cup_metrics.covered_any.front();
            rolling_cup_metrics.covered_any.pop_front();
        };

        pop_front_if_full();

        const int shortlist = static_cast<int>(generated.accepted_children.size() + generated.fallback_candidates.size());
        const int covered_low = generated.covers_low_bucket ? 1 : 0;
        const int covered_eq_delta = generated.covers_eq_delta_bucket ? 1 : 0;
        const int covered_any = (generated.covers_low_bucket || generated.covers_eq_delta_bucket || generated.covers_eq_delta_plus1_bucket) ? 1 : 0;

        rolling_cup_metrics.changed_footprints.push_back(cup_ctx.changed_footprint);
        rolling_cup_metrics.sum_changed_footprints += cup_ctx.changed_footprint;

        rolling_cup_metrics.shortlist_sizes.push_back(shortlist);
        rolling_cup_metrics.sum_shortlist_sizes += shortlist;

        rolling_cup_metrics.recovered_children.push_back(cup_ctx.recovered_from_uncles);
        rolling_cup_metrics.sum_recovered_children += cup_ctx.recovered_from_uncles;

        rolling_cup_metrics.fallback_counts.push_back(generated.fallback_count);
        rolling_cup_metrics.sum_fallback_counts += generated.fallback_count;

        rolling_cup_metrics.direct_accept_counts.push_back(generated.direct_accept_count);
        rolling_cup_metrics.sum_direct_accept_counts += generated.direct_accept_count;

        rolling_cup_metrics.covered_low.push_back(covered_low);
        rolling_cup_metrics.sum_covered_low += covered_low;

        rolling_cup_metrics.covered_eq_delta.push_back(covered_eq_delta);
        rolling_cup_metrics.sum_covered_eq_delta += covered_eq_delta;

        rolling_cup_metrics.covered_any.push_back(covered_any);
        rolling_cup_metrics.sum_covered_any += covered_any;
    }

    double rolling_avg_changed_footprint() const {
        const int n = rolling_cup_metrics.size();
        if (n == 0) return 0.0;
        return static_cast<double>(rolling_cup_metrics.sum_changed_footprints) / static_cast<double>(n);
    }

    double rolling_recovery_rate() const {
        return ratio_or_zero(rolling_cup_metrics.sum_recovered_children, rolling_cup_metrics.sum_shortlist_sizes);
    }

    double rolling_bucket_coverage_rate() const {
        const int n = rolling_cup_metrics.size();
        if (n == 0) return 0.0;
        return static_cast<double>(rolling_cup_metrics.sum_covered_any) / static_cast<double>(n);
    }

    double rolling_low_bucket_coverage_rate() const {
        const int n = rolling_cup_metrics.size();
        if (n == 0) return 0.0;
        return static_cast<double>(rolling_cup_metrics.sum_covered_low) / static_cast<double>(n);
    }

    double rolling_eq_delta_coverage_rate() const {
        const int n = rolling_cup_metrics.size();
        if (n == 0) return 0.0;
        return static_cast<double>(rolling_cup_metrics.sum_covered_eq_delta) / static_cast<double>(n);
    }

    double rolling_fallback_ratio() const {
        return ratio_or_zero(rolling_cup_metrics.sum_fallback_counts, rolling_cup_metrics.sum_shortlist_sizes);
    }

    double rolling_direct_accept_rate() const {
        return ratio_or_zero(rolling_cup_metrics.sum_direct_accept_counts, rolling_cup_metrics.sum_shortlist_sizes);
    }

    int adaptive_changed_footprint_limit() const {
        int limit = CUP_HARD_MAX_CHANGED_FOOTPRINT;
        if (rolling_cup_metrics.size() < CUP_ROLLING_MIN_SAMPLES) return limit;

        const double avg_changed = rolling_avg_changed_footprint();
        const double recovery = rolling_recovery_rate();
        const double coverage = rolling_bucket_coverage_rate();
        const double fallback = rolling_fallback_ratio();

        if (avg_changed > 120.0 && recovery < 0.35 && coverage < 0.10) {
            return 64;
        }
        if (avg_changed > 96.0 && recovery < 0.45 && coverage < 0.15) {
            return 96;
        }
        if (avg_changed > 80.0 && fallback > 0.70 && recovery < 0.50) {
            return 96;
        }
        return limit;
    }

    inline void bump_token(std::vector<uint32_t>& mark, uint32_t& token) {
        if (++token == 0) {
            std::fill(mark.begin(), mark.end(), 0);
            token = 1;
        }
    }

    inline NodeScratch& scratch_for_depth(int depth) {
        if (depth >= static_cast<int>(scratch_stack.size())) {
            scratch_stack.resize(static_cast<size_t>(depth + 1));
        }
        scratch_stack[static_cast<size_t>(depth)].reset();
        return scratch_stack[static_cast<size_t>(depth)];
    }

    inline CupFrame& frame_for_depth(int depth) {
        if (depth >= static_cast<int>(cup_frame_stack.size())) {
            cup_frame_stack.resize(static_cast<size_t>(depth + 1));
        }
        cup_frame_stack[static_cast<size_t>(depth)].reset();
        return cup_frame_stack[static_cast<size_t>(depth)];
    }

    inline CupFrame* parent_cup_frame() const {
        return parent_cup_frame_stack.empty() ? nullptr : parent_cup_frame_stack.back();
    }

    bool is_current_maximal() const {
        // NP stores OUTSIDE-P vertices only.
        // Current P is maximal iff no outside vertex has degNP >= ceil(theta_P).
        int min_add_deg = static_cast<int>(std::ceil(theta_P)); // minimum degree of a vertex that can be  added to the current pseudo-clique P
        if (min_add_deg < 0) min_add_deg = 0; // if the minimum degree is less than 0, set it to 0
        // If min_add_deg==0, any outside vertex can extend; maximal only when P covers all vertices.
        if (min_add_deg == 0) return (total_nodes_in_P == graph.total_nodes);
        return maxNPdeg < min_add_deg;
    }

    // Optional accounting (like fpce.c/mod-edge-order)
    unsigned long long numcalls_saved_by_edge_bound = 0ULL;

    bool is_vertex_in_exact_child_window(int u, int min_add_deg, int deltaP) const {
        if (u < 0 || u >= graph.total_nodes || tracks[u].inP) return false;
        const int deg = static_cast<int>(tracks[u].degNP);
        if (deg < min_add_deg) return false;
        return deg <= (deltaP + 1);
    }

    ChildBucketKind classify_child_bucket(int u, int min_add_deg, int deltaP) const {
        if (!is_vertex_in_exact_child_window(u, min_add_deg, deltaP)) return ChildBucketKind::None;
        const int deg = static_cast<int>(tracks[u].degNP);
        if (deg < deltaP) return ChildBucketKind::Low;
        if (deg == deltaP) return ChildBucketKind::EqDelta;
        if (deg == deltaP + 1) return ChildBucketKind::EqDeltaPlus1;
        return ChildBucketKind::None;
    }

    bool can_direct_accept_from_cup(int u, int current_v, int min_add_deg, int deltaP, bool strict_mode) const {
        const ChildBucketKind bucket = classify_child_bucket(u, min_add_deg, deltaP);
        if (bucket == ChildBucketKind::Low) return true;
        if (bucket == ChildBucketKind::EqDelta && u < current_v) return true;
        if (!strict_mode) return false;
        return false;
    }

    bool validate_child_exact(int u, int current_v, int min_add_deg, int deltaP) const {
        if (!is_vertex_in_exact_child_window(u, min_add_deg, deltaP)) return false;

        const int deg = static_cast<int>(tracks[u].degNP);
        if (deg < deltaP) return true;

        if (deg == deltaP) {
            if (u < current_v) return true;
            if (!graph.has_edge(current_v, u)) return false;
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (x < u && !graph.has_edge(u, x)) return false;
            }
            return true;
        }

        if (deg == deltaP + 1) {
            if (!(current_v > u) || !graph.has_edge(current_v, u)) return false;
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (!graph.has_edge(u, x)) return false;
            }
            const int targetDeg = deltaP + 1;
            if (targetDeg < static_cast<int>(headP.size())) {
                for (int x = headP[targetDeg]; x != -1; x = nextP[x]) {
                    if (x < u && !graph.has_edge(u, x)) return false;
                }
            }
            return true;
        }

        return false;
    }

    int compute_exact_child_window_size(int min_add_deg, int deltaP) const {
        if (maxNPdeg < min_add_deg) return 0;
        const int lower = std::max(0, min_add_deg);
        const int upper = std::min(maxNPdeg, deltaP + 1);
        if (upper < lower) return 0;

        int count = 0;
        for (int deg = lower; deg <= upper && deg < static_cast<int>(countNP.size()); ++deg) {
            count += countNP[deg];
        }
        return count;
    }

    int count_low_bucket_vertices(int min_add_deg, int deltaP) const {
        if (deltaP <= min_add_deg) return 0;
        int count = 0;
        for (int deg = std::max(0, min_add_deg); deg < deltaP && deg < static_cast<int>(countNP.size()); ++deg) {
            count += countNP[deg];
        }
        return count;
    }

    int count_eq_delta_bucket_vertices(int min_add_deg, int deltaP) const {
        if (deltaP < min_add_deg || deltaP < 0 || deltaP >= static_cast<int>(countNP.size())) return 0;
        return countNP[deltaP];
    }

    int count_eq_delta_plus1_bucket_vertices(int min_add_deg, int deltaP) const {
        const int targetDeg = deltaP + 1;
        if (targetDeg < min_add_deg || targetDeg < 0 || targetDeg >= static_cast<int>(countNP.size())) return 0;
        return countNP[targetDeg];
    }

    int compute_changed_footprint(int v) const {
        if (v < 0 || v + 1 >= static_cast<int>(graph.csr_offsets.size())) return 0;
        int changed = 0;
        const uint32_t begin = graph.csr_offsets[static_cast<size_t>(v)];
        const uint32_t end = graph.csr_offsets[static_cast<size_t>(v) + 1];
        for (uint32_t ei = begin; ei < end; ++ei) {
            const int u = static_cast<int>(graph.csr_neighbors[ei]);
            if (!tracks[u].inP) ++changed;
        }
        return changed;
    }

    RelevantUncles collect_relevant_previous_uncles() const {
        RelevantUncles relevant;
        const CupFrame* previous_uncles = parent_cup_frame();
        if (previous_uncles == nullptr) return relevant;
        const int estimated = static_cast<int>(previous_uncles->low_bucket_candidates.size() +
                                               previous_uncles->eq_delta_candidates.size() +
                                               previous_uncles->eq_delta_plus1_candidates.size());
        if (previous_uncles->contributing_uncle_count == 0 || estimated == 0) return relevant;
        relevant.frame = previous_uncles;
        relevant.relevant_uncle_count = previous_uncles->contributing_uncle_count;
        relevant.estimated_candidate_count = estimated;
        return relevant;
    }

    bool should_activate_cup(const RelevantUncles& relevant_uncles, int exact_window_size, int changed_footprint) {
        if (g_cup_mode == CupMode::Off || g_cup_mode == CupMode::Stats) return false;
        if (cup_runtime_disabled && (g_cup_mode == CupMode::Assist || g_cup_mode == CupMode::StrictAssist)) return false;
        if (relevant_uncles.empty()) return false;
        if (relevant_uncles.relevant_uncle_count > CUP_MAX_SHORTLIST_PER_UNCLE) return false;
        if (exact_window_size < CUP_MIN_WINDOW_SIZE || exact_window_size > CUP_MAX_WINDOW_SIZE) return false;

        const int estimated_count = std::max(1, relevant_uncles.estimated_candidate_count);
        const int adaptive_limit = adaptive_changed_footprint_limit();

        // Fast practical rejects: hard cap and size-aware early exits.
        if (changed_footprint > adaptive_limit) return false;
        if (changed_footprint > CUP_YELLOW_MAX_CHANGED_FOOTPRINT &&
            exact_window_size <= CUP_RED_SMALL_WINDOW_MAX) {
            return false;
        }
        if (changed_footprint > CUP_HIGH_CHANGED_FOOTPRINT &&
            relevant_uncles.estimated_candidate_count >= exact_window_size) {
            return false;
        }
        if (changed_footprint > CUP_HIGH_CHANGED_FOOTPRINT &&
            relevant_uncles.relevant_uncle_count > CUP_YELLOW_MAX_RELEVANT_UNCLES) {
            return false;
        }

        // Relative-to-work guards: block when update cost dwarfs possible savings.
        if (changed_footprint > CUP_CHANGED_TO_WINDOW_RATIO_LIMIT * std::max(1, exact_window_size)) return false;
        if (changed_footprint > CUP_CHANGED_TO_SHORTLIST_RATIO_LIMIT * estimated_count) return false;

        if (relevant_uncles.estimated_candidate_count >= exact_window_size) return false;
        if (relevant_uncles.frame->guard_summary != INT_MAX &&
            relevant_uncles.frame->guard_summary > exact_window_size + 1) {
            return false;
        }

        const bool is_green_zone = changed_footprint <= CUP_GREEN_MAX_CHANGED_FOOTPRINT;
        const bool is_yellow_zone =
            changed_footprint > CUP_GREEN_MAX_CHANGED_FOOTPRINT &&
            changed_footprint <= CUP_YELLOW_MAX_CHANGED_FOOTPRINT;
        const bool is_red_zone = changed_footprint > CUP_YELLOW_MAX_CHANGED_FOOTPRINT;

        if (is_yellow_zone) {
            if (relevant_uncles.relevant_uncle_count > CUP_YELLOW_MAX_RELEVANT_UNCLES) return false;
            if (relevant_uncles.estimated_candidate_count * 4 >= exact_window_size * 3) return false;
        }

        const double prior_bucket_coverage_signal = rolling_bucket_coverage_rate();
        const double prior_direct_accept_signal = rolling_direct_accept_rate();
        const double fallback_risk = rolling_fallback_ratio();

        const int estimated_candidate_reduction =
            std::max(0, exact_window_size - relevant_uncles.estimated_candidate_count);
        const double alpha = is_green_zone ? 0.04 : (is_yellow_zone ? 0.08 : 0.12);
        const double score =
            3.0 * static_cast<double>(estimated_candidate_reduction) +
            2.0 * prior_bucket_coverage_signal +
            2.0 * prior_direct_accept_signal -
            static_cast<double>(relevant_uncles.relevant_uncle_count) -
            fallback_risk -
            alpha * static_cast<double>(changed_footprint);

        if (is_red_zone) {
            const bool has_meaningful_bucket_override =
                rolling_cup_metrics.size() >= CUP_ROLLING_MIN_SAMPLES &&
                (rolling_low_bucket_coverage_rate() >= 0.20 || rolling_eq_delta_coverage_rate() >= 0.20) &&
                rolling_recovery_rate() >= 0.45;
            if (!has_meaningful_bucket_override) return false;
            return score > 1.0;
        }

        if (is_yellow_zone) return score > 0.0;

        return score > -1.0;
    }

    void collect_children_bucket_lt_delta(int min_add_deg, int deltaP,
                                          const std::vector<uint32_t>* skip_marks,
                                          uint32_t skip_token,
                                          std::vector<int>& out) const {
        for (int deg = std::max(0, min_add_deg); deg < deltaP; ++deg) {
            if (deg < 0 || deg >= static_cast<int>(headNP.size())) continue;
            for (int u = headNP[deg]; u != -1; u = nextNP[u]) {
                if (skip_marks != nullptr &&
                    static_cast<size_t>(u) < skip_marks->size() &&
                    (*skip_marks)[static_cast<size_t>(u)] == skip_token) {
                    continue;
                }
                out.push_back(u);
            }
        }
    }

    void collect_children_bucket_eq_delta(int current_v, int min_add_deg, int deltaP,
                                          const std::vector<uint32_t>* skip_marks,
                                          uint32_t skip_token,
                                          std::vector<int>& out) const {
        if (deltaP < 0 || deltaP >= static_cast<int>(headNP.size()) || deltaP < min_add_deg) return;
        for (int u = headNP[deltaP]; u != -1; u = nextNP[u]) {
            if (skip_marks != nullptr &&
                static_cast<size_t>(u) < skip_marks->size() &&
                (*skip_marks)[static_cast<size_t>(u)] == skip_token) {
                continue;
            }
            if (validate_child_exact(u, current_v, min_add_deg, deltaP)) {
                out.push_back(u);
            }
        }
    }

    void collect_children_bucket_eq_delta_plus1(int current_v, int min_add_deg, int deltaP,
                                                const std::vector<uint32_t>* skip_marks,
                                                uint32_t skip_token,
                                                std::vector<int>& out) const {
        const int targetDeg = deltaP + 1;
        if (targetDeg < 0 || targetDeg >= static_cast<int>(headNP.size()) || targetDeg < min_add_deg) return;
        for (int u = headNP[targetDeg]; u != -1; u = nextNP[u]) {
            if (skip_marks != nullptr &&
                static_cast<size_t>(u) < skip_marks->size() &&
                (*skip_marks)[static_cast<size_t>(u)] == skip_token) {
                continue;
            }
            if (validate_child_exact(u, current_v, min_add_deg, deltaP)) {
                out.push_back(u);
            }
        }
    }

    void collect_children_baseline_full(int current_v, int min_add_deg, int deltaP,
                                        std::vector<int>& out) const {
        out.clear();
        collect_children_bucket_lt_delta(min_add_deg, deltaP, nullptr, 0, out);
        collect_children_bucket_eq_delta(current_v, min_add_deg, deltaP, nullptr, 0, out);
        collect_children_bucket_eq_delta_plus1(current_v, min_add_deg, deltaP, nullptr, 0, out);
    }

    void cup_generate_candidates_from_uncles(int current_v, int min_add_deg, int deltaP,
                                             const RelevantUncles& relevant_uncles,
                                             bool strict_direct_accept,
                                             CupGenerated& out) {
        out.reset();
        out.relevant_uncle_count = relevant_uncles.relevant_uncle_count;
        out.estimated_shortlist_size = relevant_uncles.estimated_candidate_count;
        if (relevant_uncles.empty()) return;

        const int low_total = count_low_bucket_vertices(min_add_deg, deltaP);
        const int eq_total = count_eq_delta_bucket_vertices(min_add_deg, deltaP);
        const int plus1_total = count_eq_delta_plus1_bucket_vertices(min_add_deg, deltaP);

        bump_token(cup_generated_mark, cup_generated_token);
        bump_token(cup_handled_mark, cup_handled_token);

        auto process_source = [&](const std::vector<int>& candidates) {
            for (int u : candidates) {
                if (u < 0 || u >= graph.total_nodes) continue;
                if (cup_generated_mark[static_cast<size_t>(u)] == cup_generated_token) continue;
                const ChildBucketKind bucket = classify_child_bucket(u, min_add_deg, deltaP);
                if (bucket == ChildBucketKind::None) continue;

                cup_generated_mark[static_cast<size_t>(u)] = cup_generated_token;
                cup_handled_mark[static_cast<size_t>(u)] = cup_handled_token;

                if (bucket == ChildBucketKind::Low) ++out.handled_low_count;
                else if (bucket == ChildBucketKind::EqDelta) ++out.handled_eq_delta_count;
                else if (bucket == ChildBucketKind::EqDeltaPlus1) ++out.handled_eq_delta_plus1_count;

                if (can_direct_accept_from_cup(u, current_v, min_add_deg, deltaP, strict_direct_accept)) {
                    out.accepted_children.push_back(u);
                    ++out.direct_accept_count;
                } else {
                    out.fallback_candidates.push_back(u);
                    ++out.fallback_count;
                }
            }
        };

        process_source(relevant_uncles.frame->low_bucket_candidates);
        process_source(relevant_uncles.frame->eq_delta_candidates);
        process_source(relevant_uncles.frame->eq_delta_plus1_candidates);

        out.resolved_window_count =
            out.handled_low_count + out.handled_eq_delta_count + out.handled_eq_delta_plus1_count;
        out.covers_low_bucket = (low_total > 0 && out.handled_low_count >= low_total);
        out.covers_eq_delta_bucket = (eq_total > 0 && out.handled_eq_delta_count >= eq_total);
        out.covers_eq_delta_plus1_bucket = (plus1_total > 0 && out.handled_eq_delta_plus1_count >= plus1_total);
    }

    void validate_cup_fallback_candidates(int current_v, int min_add_deg, int deltaP,
                                          const CupGenerated& generated,
                                          std::vector<int>& out) const {
        for (int u : generated.fallback_candidates) {
            if (validate_child_exact(u, current_v, min_add_deg, deltaP)) {
                out.push_back(u);
            }
        }
    }

    void baseline_complete_residual_low_bucket(int min_add_deg, int deltaP,
                                               const CupGenerated& generated,
                                               std::vector<int>& out) const {
        if (generated.covers_low_bucket) return;
        collect_children_bucket_lt_delta(min_add_deg, deltaP, &cup_handled_mark, cup_handled_token, out);
    }

    void baseline_complete_residual_eq_delta(int current_v, int min_add_deg, int deltaP,
                                             const CupGenerated& generated,
                                             std::vector<int>& out) const {
        if (generated.covers_eq_delta_bucket) return;
        collect_children_bucket_eq_delta(current_v, min_add_deg, deltaP, &cup_handled_mark, cup_handled_token, out);
    }

    void baseline_complete_residual_eq_delta_plus1(int current_v, int min_add_deg, int deltaP,
                                                   const CupGenerated& generated,
                                                   std::vector<int>& out) const {
        if (generated.covers_eq_delta_plus1_bucket) return;
        collect_children_bucket_eq_delta_plus1(current_v, min_add_deg, deltaP, &cup_handled_mark, cup_handled_token, out);
    }

    void publish_current_node_to_parent_frame(int current_v, const std::vector<int>& exact_children,
                                              int min_add_deg, int deltaP, int changed_footprint) {
        CupFrame* parent_frame_ptr = parent_cup_frame();
        if (parent_frame_ptr == nullptr || g_cup_mode == CupMode::Off) return;

        parent_frame_ptr->previous_uncle_count++;
        parent_frame_ptr->changed_footprint =
            std::max(parent_frame_ptr->changed_footprint, changed_footprint);
        parent_frame_ptr->guard_summary =
            std::min(parent_frame_ptr->guard_summary, compute_exact_child_window_size(min_add_deg, deltaP));
        parent_frame_ptr->source_vertex = current_v;

        if (exact_children.empty()) return;

        parent_frame_ptr->contributing_uncle_count++;
        size_t published = 0;
        for (int u : exact_children) {
            const ChildBucketKind bucket = classify_child_bucket(u, min_add_deg, deltaP);
            if (bucket == ChildBucketKind::Low) {
                parent_frame_ptr->low_bucket_candidates.push_back(u);
            } else if (bucket == ChildBucketKind::EqDelta) {
                parent_frame_ptr->eq_delta_candidates.push_back(u);
            } else if (bucket == ChildBucketKind::EqDeltaPlus1) {
                parent_frame_ptr->eq_delta_plus1_candidates.push_back(u);
            } else {
                continue;
            }

            ++published;
            if (published >= static_cast<size_t>(CUP_MAX_SHORTLIST_PER_UNCLE)) break;
            const size_t total_size = parent_frame_ptr->low_bucket_candidates.size() +
                                      parent_frame_ptr->eq_delta_candidates.size() +
                                      parent_frame_ptr->eq_delta_plus1_candidates.size();
            if (total_size >= static_cast<size_t>(CUP_MAX_FRAME_CANDIDATES)) break;
        }
    }

    void finalize_cup_stats(const std::vector<int>& baseline_children,
                            const CupGenerated& generated,
                            CupNodeContext& cup_ctx) {
        if (!cup_ctx.active) return;

        bump_token(baseline_mark, baseline_mark_token);
        for (int u : baseline_children) {
            baseline_mark[static_cast<size_t>(u)] = baseline_mark_token;
        }

        int recovered = 0;
        for (int u : generated.accepted_children) {
            if (baseline_mark[static_cast<size_t>(u)] == baseline_mark_token) ++recovered;
        }
        for (int u : generated.fallback_candidates) {
            if (baseline_mark[static_cast<size_t>(u)] == baseline_mark_token) ++recovered;
        }
        cup_ctx.recovered_from_uncles = recovered;
        cup_ctx.direct_accept_count = generated.direct_accept_count;
        cup_ctx.uncle_candidates_needing_fallback = generated.fallback_count;
        cup_ctx.resolved_window_count = generated.resolved_window_count;
        cup_ctx.covered_low_bucket = generated.covers_low_bucket;
        cup_ctx.covered_eq_delta_bucket = generated.covers_eq_delta_bucket;
        cup_ctx.covered_eq_delta_plus1_bucket = generated.covers_eq_delta_plus1_bucket;

        cup_stats.assist_attempts++;
        cup_stats.total_shortlist_size += static_cast<unsigned long long>(
            generated.accepted_children.size() + generated.fallback_candidates.size());
        cup_stats.total_relevant_uncle_count += static_cast<unsigned long long>(cup_ctx.relevant_uncle_count);
        cup_stats.total_changed_footprint += static_cast<unsigned long long>(cup_ctx.changed_footprint);
        cup_stats.fallback_validations += static_cast<unsigned long long>(generated.fallback_count);
        cup_stats.direct_accepts += static_cast<unsigned long long>(generated.direct_accept_count);
        cup_stats.baseline_children_recovered += static_cast<unsigned long long>(recovered);
        if (generated.covers_low_bucket) ++cup_stats.covered_low_buckets;
        if (generated.covers_eq_delta_bucket) ++cup_stats.covered_eq_delta_buckets;
        if (generated.covers_eq_delta_plus1_bucket) ++cup_stats.covered_eq_delta_plus1_buckets;
        if (recovered > 0) ++cup_stats.assist_successes;

        // Keep adaptive tuning tied to assist modes only; compare mode remains parity-only.
        if (g_cup_mode == CupMode::Assist || g_cup_mode == CupMode::StrictAssist) {
            push_rolling_cup_metrics(cup_ctx, generated);
        }

        if ((g_cup_mode == CupMode::Assist || g_cup_mode == CupMode::StrictAssist) &&
            g_enable_cup_auto_disable && !cup_runtime_disabled && cup_stats.assist_attempts >= 128ULL) {
            const bool low_recovery =
                (cup_stats.baseline_children_recovered * 3ULL) < cup_stats.total_shortlist_size;
            const bool high_fallback =
                cup_stats.fallback_validations > (cup_stats.direct_accepts + 1ULL) * 4ULL;
            const bool poor_bucket_coverage =
                (cup_stats.covered_low_buckets + cup_stats.covered_eq_delta_buckets +
                 cup_stats.covered_eq_delta_plus1_buckets) * 3ULL < cup_stats.assist_attempts;
            if (low_recovery && high_fallback && poor_bucket_coverage) {
                cup_runtime_disabled = true;
                cup_stats.auto_disables++;
            }
        }
    }

    // Update report_if_new to take size directly
    inline void report_if_new(int size) {
        if (size < min_size) return; // size of P is less than l still, so return
        if (size < (int)pseudo_cliques_count.size())
            pseudo_cliques_count[size] += 1;
    }
    
    // Small helpers that work with intrusive lists
    // returns the minimum degree in P
    inline int min_degree_in_P_fast() const {
        if (total_nodes_in_P == 0) return 0;
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            if (d >= 0 && d < (int)headP.size() && headP[d] != -1) return d; // if the degree is not -1 and less than the size of the P bucket and the head of the P bucket for degree d is not -1, return the degree
        }
        return 0;
    }
    inline int current_min_core_in_P_fast() const { // used by prune_by_edge_bound method
        return (cur_min_core == INT_MAX) ? 0 : cur_min_core;
    }

    bool prune_by_edge_bound() const {
        const int S = min_size;                 // ℓ
        const int P = total_nodes_in_P;         // |P|
        const int t = S - P;                    // vertices still needed
    
        // Only prune when 1 ≤ |P| < ℓ  (same gating as mod-edge-order / fpce.c)
        if (P < 1 || t <= 0) return false;
    
        const int delta   = min_degree_in_P_fast();         // δ(P)
        const int cP      = current_min_core_in_P_fast();   // c(P)
        const double sumE = static_cast<double>(total_edges_in_P);
        const int LHS     = lhs_required_edges;
    
        // Lemma 9 shape (fpce.c): g = c(P) - δ(P) - 1
        const int g = cP - delta - 1;
    
        if (g >= 0) {
            // RHS1 = |E[P]| + min(c(P), δ(P)+t)*(t + δ + 0.5) - 0.5*(δ+1)*(2δ+1)
            const int max_value = std::min(cP, delta + t);
            const double rhs1 = std::floor(
                sumE + max_value * (t + delta + 0.5)
                    - 0.5 * (delta + 1) * (2 * delta + 1)
            );
            if (LHS > rhs1) return true;   // EDGE_BOUND1
        } else {
            // RHS2 = |E[P]| + t*δ
            const double rhs2 = std::floor(sumE + t * static_cast<double>(delta));
            if (LHS > rhs2) return true;   // EDGE_BOUND2
        }
    
        // Lemma 8 fallback: RHS3 = |E[P]| + t*(δ + (t+1)/2)
        const double rhs3 = std::floor(sumE + t * (delta + (t + 1) * 0.5));
        if (LHS > rhs3) return true;       // EDGE_BOUND3
    
        return false;
    }

public:
    PseudoCliqueEnumerator(Graph& graph, float theta = 1.0, int min_size = 1, int max_size = -1, int order_ub = -1) : graph(graph), theta(theta), min_size(min_size) {        
        this->order_ub = order_ub; // μ (order bound); -1 means "disabled"
        std::cout << "total nodes: " << graph.total_nodes << " "   << "\n";

        this->max_size = max_size != -1 ? max_size : graph.total_nodes;

        //pseudo_cliques.resize(max_size + 1);
        pseudo_cliques_count.resize(max_size + 1, 0); // count of pseudo-cliques of size i


        // Initialize intrusive buckets, assign(count, value)
        headP.assign(max_size + 1, -1); // headP[k]: Points to first vertex inside P having degree k
        headNP.assign(max_size + 2, -1); // +1 for targetDeg access, headNP[k]: Points to first candidate vertex outside P (NP) connected to k vertices inside P
        countP.assign(max_size + 1, 0); // countP[k]: Stores the # of vertices inside P having degree k
        countNP.assign(max_size + 2, 0); // countNP[k]: Stores the # of vertices outside P (NP) connected to k vertices inside P
        prevP.assign(graph.total_nodes, -1); // prevP[v]: Points to the previous vertex in the P bucket for vertex v
        nextP.assign(graph.total_nodes, -1); // nextP[v]: Points to the next vertex in the P bucket for vertex v
        prevNP.assign(graph.total_nodes, -1); // prevNP[v]: Points to the previous vertex in the NP bucket for vertex v
        nextNP.assign(graph.total_nodes, -1); // nextNP[v]: Points to the next vertex in the NP bucket for vertex v
        
        tracks.clear();
        tracks.reserve(graph.total_nodes); // reserve space for all vertices
        tracks.resize(graph.total_nodes); // initialize all tracks to 0
        for (int i = 0; i < graph.total_nodes; ++i) {
            tracks[i].inP = 0; // 0: vertex is not in P
            tracks[i].degP = -1; // -1: vertex has no degree in P
            tracks[i].degNP = 0; // 0: vertex has no degree in P
            insertNP_front(i, 0); // insert vertex into NP bucket with degree 0
        }

        theta_P = 0; // theta_P: Theta * |P| * (|P| + 1) / 2 - |E_P|
        total_nodes_in_P = 0;
        total_edges_in_P = 0;

        // LHS is fixed for the (ℓ, θ) target
        // Use a small epsilon to prevent 12.0000001 -> 13
        lhs_required_edges = static_cast<int>(
            std::ceil(theta * (min_size * (min_size - 1) / 2.0) - 1e-7)
        );// minimum # of edges P must have to satisfy theta

        // Prepare core histogram only if edge-bound can use it
        max_core_seen = 0;
        if (g_enable_edge_bound) {
            for (int c : graph.core_numbers) max_core_seen = std::max(max_core_seen, c); // max core number in the graph
            core_hist.assign(max_core_seen + 1, 0); // initialize core histogram with max core number + 1
            cur_min_core = INT_MAX; // initialize current minimum core to INT_MAX
        } else {
            core_hist.clear();
            cur_min_core = INT_MAX;
        }

        scratch_stack.reserve(static_cast<size_t>(std::max(1, this->max_size + 1)));
        cup_frame_stack.reserve(static_cast<size_t>(std::max(1, this->max_size + 1)));
        parent_cup_frame_stack.reserve(static_cast<size_t>(std::max(1, this->max_size + 1)));
        cup_generated_mark.assign(static_cast<size_t>(graph.total_nodes), 0u);
        cup_handled_mark.assign(static_cast<size_t>(graph.total_nodes), 0u);
        baseline_mark.assign(static_cast<size_t>(graph.total_nodes), 0u);
    }

    // calculates the minimum connection required for any new node to join the current P
    void set_theta_P() { //θ(K) = θ * clq(|K| + 1) - |E[K]|
        theta_P = (theta * total_nodes_in_P * (total_nodes_in_P + 1) / 2.0) - total_edges_in_P;
    }

    // --- Intrusive list helpers ---
    // removes vertex v from the doubly linked list corresponding to degree deg inside P
    inline void unlinkP(int v, int deg) {
        int pv = prevP[v]; // previous vertex in the P bucket for vertex v
        int nv = nextP[v]; // next vertex in the P bucket for vertex v
        if (pv != -1) nextP[pv] = nv; else headP[deg] = nv; // if the previous vertex is not -1, set the next vertex of the previous vertex to the next vertex
        if (nv != -1) prevP[nv] = pv; // if the next vertex is not -1, set the previous vertex of the next vertex to the previous vertex
        prevP[v] = nextP[v] = -1; // set the previous and next vertex of vertex v to -1
        if (deg >= 0 && deg < (int)countP.size()) countP[deg]--; // decrement the count of vertices in the P bucket for degree deg
    }
    // inserts vertex v at the beginning of the doubly linked list for the degree bucket
    inline void insertP_front(int v, int deg) {
        prevP[v] = -1; // set the previous vertex of vertex v to -1
        nextP[v] = headP[deg]; // set the next vertex of vertex v to the head of the P bucket for degree deg
        if (headP[deg] != -1) prevP[headP[deg]] = v; // if the head of the P bucket for degree deg is not -1, set the previous vertex of the head of the P bucket for degree deg to vertex v
        headP[deg] = v; // set the head of the P bucket for degree deg to vertex v
        if (deg >= 0 && deg < (int)countP.size()) countP[deg]++; // increment the count of vertices in the P bucket for degree deg
    }
    inline void unlinkNP(int v, int deg) {
        int pv = prevNP[v]; // previous vertex in the NP bucket for vertex v
        int nv = nextNP[v]; // next vertex in the NP bucket for vertex v
        if (pv != -1) nextNP[pv] = nv; else headNP[deg] = nv; // if the previous vertex is not -1, set the next vertex of the previous vertex to the next vertex
        if (nv != -1) prevNP[nv] = pv; // if the next vertex is not -1, set the previous vertex of the next vertex to the previous vertex
        prevNP[v] = nextNP[v] = -1; // set the previous and next vertex of vertex v to -1
        if (deg >= 0 && deg < (int)countNP.size()) {
            countNP[deg]--; // decrement the count of vertices in the NP bucket for degree deg
            // Maintain O(1) maximality via maxNPdeg
            if (deg == maxNPdeg && countNP[deg] == 0) {
                while (maxNPdeg > 0 && countNP[maxNPdeg] == 0) --maxNPdeg;
            }
        }
    }
    inline void insertNP_front(int v, int deg) {
        prevNP[v] = -1; // set the previous vertex of vertex v to -1
        nextNP[v] = headNP[deg]; // set the next vertex of vertex v to the head of the NP bucket for degree deg
        if (headNP[deg] != -1) prevNP[headNP[deg]] = v; // if the head of the NP bucket for degree deg is not -1, set the previous vertex of the head of the NP bucket for degree deg to vertex v
        headNP[deg] = v; // set the head of the NP bucket for degree deg to vertex v
        if (deg >= 0 && deg < (int)countNP.size()) {
            countNP[deg]++; // increment the count of vertices in the NP bucket for degree deg
            // Maintain O(1) maximality via maxNPdeg
            if (deg > maxNPdeg) maxNPdeg = deg;
        }
    }


    void add_to_inside_P(int v, CupFrame* parent_frame = nullptr);
    void remove_from_inside_P(int v);
    void iter(int v);
    // Internal: insert v into P WITHOUT recursing (used by Turán seed build)
    void add_vertex_internal(int v);

    std::vector<int> get_pseudo_cliques_count(){
        return pseudo_cliques_count;
    };

    int get_iter_count(){
        return iter_count;
    };

    unsigned long long get_numcalls_saved_by_edge_bound() const {
        return numcalls_saved_by_edge_bound;
    }

    const CupStats& get_cup_stats() const {
        return cup_stats;
    }

    bool is_cup_runtime_disabled() const {
        return cup_runtime_disabled;
    }

    // Enumerate starting from a Turán seed (R-clique) provided by EBBkC, takes a "seed" clique (found by the external library) and "teleports" the algorithm to that state, effectively skipping the early levels of the recursion tree to prune the search space.
    void enumerate_with_turan(const std::vector<int>& turan_seed) {
        if (turan_seed.empty()) return;

        // A) Deterministic order: sort by (deg, id) to match mod-edge-order seed canonicalization
        auto sorted_seed = turan_seed;
        std::sort(sorted_seed.begin(), sorted_seed.end(),
                [&](int a, int b) {
                    // degree proxy from CSR slice lengths, sort by degree first, then by id to ensure deterministic ordering
                    const uint32_t da = graph.csr_offsets[a+1] - graph.csr_offsets[a]; // cumulative count of next node - current node
                    const uint32_t db = graph.csr_offsets[b+1] - graph.csr_offsets[b];
                    return std::make_pair(da, a) < std::make_pair(db, b);
                });

        // B) Build the whole seed into P with NO recursion
        for (int v : sorted_seed) add_vertex_internal(v);

        // C) Choose v* from the current δ(P) bucket (minimum-id in that bucket)
        int deltaP = 0;
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            if (d >= 0 && d < (int)headP.size() && headP[d] != -1) { deltaP = d; break; } // starting scanning from 0, first non-empty bucket corrsponds to min-degree
        }
        int v_star = -1;
        if (deltaP >= 0 && deltaP < (int)headP.size()) { // if the minimum degree in P is not -1 and less than the size of the P bucket
            int best = -1;
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) { // iterate through the P bucket for degree deltaP
                if (best == -1 || x < best) best = x; // find the vertex with the minimum id
            }
            v_star = (best != -1 ? best : sorted_seed.front()); // if the best vertex is not -1, set v_star to the best vertex, otherwise set v_star to the first vertex in the sorted seed
        } else {
            v_star = sorted_seed.front(); // if the minimum degree in P is -1 or greater than the size of the P bucket, set v_star to the first vertex in the sorted seed
        }

        // D) Recurse from v*
        iter(v_star);

        // E) Roll back the seed in reverse order
        for (auto it = sorted_seed.rbegin(); it != sorted_seed.rend(); ++it) {
            remove_from_inside_P(*it);
        }
    }


    
}; 

// add a vertex to the P bucket without recursing, used by Turán seed build
// add a vertex to the P bucket without recursing, used by Turán seed build
void PseudoCliqueEnumerator::add_vertex_internal(int v) {
    if (total_nodes_in_P > max_size) return;

    // NP stores OUTSIDE-P vertices only: remove v from NP before moving it into P
    {
        const int dnp_v = static_cast<int>(tracks[v].degNP);
        unlinkNP(v, dnp_v);
    }

    // move v into P
    tracks[v].inP  = 1; // set vertex v to be in P
    tracks[v].degP = static_cast<uint16_t>(tracks[v].degNP); // deg into current P
    insertP_front(v, tracks[v].degP); // insert vertex v into P with degree tracks[v].degP

    // core histogram (for edge-bound)
    int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0; // get the core number of vertex v
    if ((int)core_hist.size() <= cv) core_hist.resize(cv + 1, 0); // resize the core histogram to the core number of vertex v
    core_hist[cv] += 1; // increments the count of vertices in P that have core cv
    if (cv < cur_min_core) cur_min_core = cv; // update the minimum core number in P

    // update neighbors' P/NP degrees and |E[P]|
    uint32_t begin = (v >= 0 && v + 1 < (int)graph.csr_offsets.size())
                     ? graph.csr_offsets[(size_t)v] : 0u; // neighbor list for vertex v starts
    uint32_t end   = (v >= 0 && v + 1 < (int)graph.csr_offsets.size())
                     ? graph.csr_offsets[(size_t)v + 1] : begin; // neighbor list for vertex v ends
    for (uint32_t ei = begin; ei < end; ++ei) { // iterate through the neighbors of vertex v
        int u = (int)graph.csr_neighbors[ei]; // get the neighbor of vertex v
        // Update neighbors inside P
        if (tracks[u].inP) {
            total_edges_in_P += 1; // increment the number of edges in P
            int d = static_cast<int>(tracks[u].degP);
            unlinkP(u, d); // remove the neighbor from P
            insertP_front(u, d + 1); // insert the neighbor into P with degree d + 1
            tracks[u].degP = static_cast<uint16_t>(d + 1); // update the degree of the neighbor in P
        } else {
            // Update neighbors outside P (NP)
            int dnp = static_cast<int>(tracks[u].degNP); // get the degree of the neighbor in NP
            unlinkNP(u, dnp); // remove the neighbor from NP
            insertNP_front(u, dnp + 1); // insert the neighbor into NP with degree dnp + 1
            tracks[u].degNP = static_cast<uint16_t>(dnp + 1); // update the degree of the neighbor in NP
        }
    }

    total_nodes_in_P += 1; // increment the number of nodes in P
    set_theta_P(); // update the theta_P
}

void PseudoCliqueEnumerator::add_to_inside_P(int v, CupFrame* parent_frame) {
    parent_cup_frame_stack.push_back(parent_frame);
    add_vertex_internal(v);   // state updates only
    iter(v);                  // now recurse once (normal growth path)
    parent_cup_frame_stack.pop_back();
}

void PseudoCliqueEnumerator::remove_from_inside_P(int v) {
    const int v_deg_in_P = static_cast<int>(tracks[v].degP); // deg(v, P\{v}) after removal

    tracks[v].inP = 0; // set vertex v to be not in P
    int degree = static_cast<int>(tracks[v].degP);
    if (degree >= 0) unlinkP(v, degree); // remove vertex v from P

    // v becomes outside: re-insert into NP with its degree into the new P
    tracks[v].degNP = static_cast<uint16_t>(v_deg_in_P);
    insertNP_front(v, v_deg_in_P);

    // Iterate Neighbors of v
    {
        uint32_t begin = (v >= 0 && v + 1 < (int)graph.csr_offsets.size()) ? graph.csr_offsets[static_cast<size_t>(v)] : 0u; // neighbor list for vertex v starts
        uint32_t end = (v >= 0 && v + 1 < (int)graph.csr_offsets.size()) ? graph.csr_offsets[static_cast<size_t>(v + 1)] : begin; // neighbor list for vertex v ends
        for (uint32_t ei = begin; ei < end; ++ei) {
            int adj_u = static_cast<int>(graph.csr_neighbors[ei]); // neighbor of v
            if (tracks[adj_u].inP) {
                // Downgrade neighbors inside P
                total_edges_in_P--; // update |E[P]|
                int d = static_cast<int>(tracks[adj_u].degP); // degree of neighbor
                unlinkP(adj_u, d); // remove neighbor from P
                if (d > 0) {
                    insertP_front(adj_u, d - 1); // insert the neighbor into P with degree (degree - 1)
                    tracks[adj_u].degP = static_cast<uint16_t>(d - 1); // update the degree of the neighbor in P
                } else {
                    // Sanity: degree should not be negative; clamp to 0
                    tracks[adj_u].degP = 0; // set the degree of the neighbor in P to 0
                }
            } else {
                // Downgrade neighbors outside P (NP)
                int degree_np = static_cast<int>(tracks[adj_u].degNP);
                unlinkNP(adj_u, degree_np); // remove neighbor from NP
                if (degree_np > 0) {
                    insertNP_front(adj_u, degree_np - 1); // insert the neighbor into NP with degree (degree_np - 1)
                    tracks[adj_u].degNP = static_cast<uint16_t>(degree_np - 1); // update the degree of the neighbor in NP
                } else {
                    tracks[adj_u].degNP = 0; // set the degree of the neighbor in NP to 0
                }
            }
        }
    }

    // Update Edge Bound Statistics and c(P) on removal
    int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0; // get the core number of vertex v
    if (cv >= 0 && cv < (int)core_hist.size()) {
        if (core_hist[cv] > 0) core_hist[cv] -= 1; // decrement the core number of vertex v in the core histogram
        if (core_hist[cv] == 0 && cv == cur_min_core) {
            // advance to next non-empty bucket or reset
            int nextc = cv;
            while (nextc < (int)core_hist.size() && core_hist[nextc] == 0) ++nextc;
            cur_min_core = (nextc < (int)core_hist.size()) ? nextc : INT_MAX;
        }
    }

    total_nodes_in_P--;
    set_theta_P();
}

void PseudoCliqueEnumerator::iter(int v) {
    iter_count++;

    // EDGE bound: apply only when 1 ≤ |P| < ℓ (parity-correct ob-eb)
    if (g_enable_edge_bound && total_nodes_in_P >= 1 && total_nodes_in_P < min_size) {
        if (prune_by_edge_bound()) {
            numcalls_saved_by_edge_bound++;
            return;
        }
    }

    // Integer add-degree threshold must be ceil(theta_P)
    const int min_add_deg = static_cast<int>(std::ceil(theta_P));

    // Always partition with respect to δ(P)
    const int deltaP = min_degree_in_P_fast(); // min-deg currently inside P

    // --- maximality before child enumeration ---
    if (is_current_maximal()) {
        if (total_nodes_in_P >= min_size) report_if_new(total_nodes_in_P);
        return;
    }

    // Local order-bound stop at μ: do not expand past μ
    if (g_enable_order_bound && order_ub > 0 && total_nodes_in_P == order_ub) return;

    const int depth = total_nodes_in_P;
    NodeScratch& scratch = scratch_for_depth(depth);
    CupFrame& local_cup_frame = frame_for_depth(depth);
    CupNodeContext cup_ctx;
    cup_ctx.exact_window_size = compute_exact_child_window_size(min_add_deg, deltaP);
    const RelevantUncles relevant_uncles = collect_relevant_previous_uncles();

    if (!relevant_uncles.empty()) {
        cup_stats.nodes_with_previous_uncles++;
    } else {
        cup_stats.nodes_without_previous_uncles++;
    }
    cup_stats.total_child_window_size += static_cast<unsigned long long>(cup_ctx.exact_window_size);
    cup_ctx.changed_footprint = compute_changed_footprint(v);
    cup_ctx.attempted = !relevant_uncles.empty();
    cup_ctx.auto_disabled = cup_runtime_disabled &&
        (g_cup_mode == CupMode::Assist || g_cup_mode == CupMode::StrictAssist);
    cup_ctx.relevant_uncle_count = relevant_uncles.relevant_uncle_count;

    if (g_cup_mode == CupMode::Off) {
        collect_children_baseline_full(v, min_add_deg, deltaP, scratch.children);
    } else if (g_cup_mode == CupMode::Stats) {
        cup_stats.total_relevant_uncle_count += static_cast<unsigned long long>(relevant_uncles.relevant_uncle_count);
        cup_stats.total_changed_footprint += static_cast<unsigned long long>(cup_ctx.changed_footprint);
        cup_stats.total_shortlist_size += static_cast<unsigned long long>(relevant_uncles.estimated_candidate_count);
        collect_children_baseline_full(v, min_add_deg, deltaP, scratch.children);
    } else {
        cup_ctx.active = should_activate_cup(relevant_uncles, cup_ctx.exact_window_size, cup_ctx.changed_footprint);
        if (!cup_ctx.active) {
            collect_children_baseline_full(v, min_add_deg, deltaP, scratch.children);
        } else {
            const bool strict_direct_accept = (g_cup_mode == CupMode::StrictAssist);
            cup_generate_candidates_from_uncles(
                v,
                min_add_deg,
                deltaP,
                relevant_uncles,
                strict_direct_accept,
                scratch.cup_generated
            );

            scratch.children.insert(
                scratch.children.end(),
                scratch.cup_generated.accepted_children.begin(),
                scratch.cup_generated.accepted_children.end()
            );
            validate_cup_fallback_candidates(v, min_add_deg, deltaP, scratch.cup_generated, scratch.children);
            baseline_complete_residual_low_bucket(min_add_deg, deltaP, scratch.cup_generated, scratch.children);
            baseline_complete_residual_eq_delta(v, min_add_deg, deltaP, scratch.cup_generated, scratch.children);
            baseline_complete_residual_eq_delta_plus1(v, min_add_deg, deltaP, scratch.cup_generated, scratch.children);

            if (g_cup_mode == CupMode::Compare) {
                collect_children_baseline_full(v, min_add_deg, deltaP, scratch.baseline_children);
                scratch.compare_children = scratch.children;
                if (scratch.compare_children.size() > 1) {
                    std::sort(scratch.compare_children.begin(), scratch.compare_children.end());
                    scratch.compare_children.erase(
                        std::unique(scratch.compare_children.begin(), scratch.compare_children.end()),
                        scratch.compare_children.end()
                    );
                }
                if (scratch.baseline_children.size() > 1) {
                    std::sort(scratch.baseline_children.begin(), scratch.baseline_children.end());
                    scratch.baseline_children.erase(
                        std::unique(scratch.baseline_children.begin(), scratch.baseline_children.end()),
                        scratch.baseline_children.end()
                    );
                }
                cup_stats.compare_parity_checks++;
                assert(scratch.compare_children == scratch.baseline_children);
                finalize_cup_stats(scratch.baseline_children, scratch.cup_generated, cup_ctx);
            } else {
                finalize_cup_stats(scratch.children, scratch.cup_generated, cup_ctx);
            }
        }
    }

    // Deterministic expansion order
    if (scratch.children.size() > 1) {
        std::sort(scratch.children.begin(), scratch.children.end());
        scratch.children.erase(std::unique(scratch.children.begin(), scratch.children.end()), scratch.children.end());
    }
    publish_current_node_to_parent_frame(v, scratch.children, min_add_deg, deltaP, cup_ctx.changed_footprint);
    for (int u : scratch.children) {
        add_to_inside_P(u, &local_cup_frame);
        remove_from_inside_P(u);
    }
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <filename> [--theta <theta>] [--minimum <min>] [--maximum <max>]"
                  << " [--cup-stats|--cup-compare|--cup-assist|--cup-strict-assist|--no-cup]" << std::endl;
        return 1;
    }

    // 1. Re‑assemble the (possibly space‑containing) filename.
    std::vector<std::string> filename_parts;
    int arg_idx = 1;
    for (; arg_idx < argc; ++arg_idx) {
        if (std::strncmp(argv[arg_idx], "--", 2) == 0) break; // first option detected
        filename_parts.emplace_back(argv[arg_idx]);
    }

    if (filename_parts.empty()) {
        std::cerr << "Error: filename not provided." << std::endl;
        return 1;
    }

    std::string filename;
    for (size_t k = 0; k < filename_parts.size(); ++k) {
        if (k) filename += ' ';
        filename += filename_parts[k];
    }

    // 2. Parse optional arguments that follow the filename.
    double theta = 1.0;
    int minimum = 1;
    int maximum = std::numeric_limits<int>::max();
    

    for (; arg_idx < argc; ++arg_idx) {
        std::string arg = argv[arg_idx];
        if (arg == "--theta" && arg_idx + 1 < argc) {
            theta = std::stod(argv[++arg_idx]);
        } else if (arg == "--minimum" && arg_idx + 1 < argc) {
            minimum = std::stoi(argv[++arg_idx]);
        } else if (arg == "--mode" && arg_idx + 1 < argc) {
            int mode = std::stoi(argv[++arg_idx]);
            if (mode == 1) { g_enable_order_bound = false; g_enable_edge_bound = false; g_enable_turan = false; }
            else if (mode == 2) { g_enable_order_bound = true; g_enable_edge_bound = false; g_enable_turan = false; }
            else if (mode == 3) { g_enable_order_bound = true; g_enable_edge_bound = true; g_enable_turan = false; }
            else { g_enable_order_bound = true; g_enable_edge_bound = true; g_enable_turan = true; }
        } else if (arg == "--no-order") {
            g_enable_order_bound = false;
        } else if (arg == "--no-edge") {
            g_enable_edge_bound = false;
        } else if (arg == "--no-turan") {
            g_enable_turan = false;
        } else if (arg == "--cup-stats") {
            g_cup_mode = CupMode::Stats;
        } else if (arg == "--cup-compare") {
            g_cup_mode = CupMode::Compare;
        } else if (arg == "--cup-assist") {
            g_cup_mode = CupMode::Assist;
        } else if (arg == "--cup-strict-assist") {
            g_cup_mode = CupMode::StrictAssist;
        } else if (arg == "--no-cup") {
            g_cup_mode = CupMode::Off;
        } else if (arg == "--no-cup-auto-disable") {
            g_enable_cup_auto_disable = false;
        } else if (arg == "--order") {
            g_enable_order_bound = true;
        } else if (arg == "--edge") {
            g_enable_edge_bound = true;
        } else if (arg == "--turan") {
            g_enable_turan = true;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }
    if (minimum <= 0) {
        std::cerr << "Error: --minimum must be positive." << std::endl;
        return 1;
    }
    if (theta <= 0.0 || theta > 1.0) {
        std::cerr << "Error: --theta must be in the open interval (0,1)." << std::endl;
        return 1;
    }

    int R = std::ceil(1.0 / (1.0 - theta * (minimum - 1) / static_cast<double>(minimum)));
    std::cout << "Computed R: " << R << std::endl;
    std::cout << "Gates => order:" << (g_enable_order_bound?"on":"off")
              << " edge:" << (g_enable_edge_bound?"on":"off")
              << " turan:" << (g_enable_turan?"on":"off")
              << " cup:" << cup_mode_name(g_cup_mode) << std::endl;

    // 4. Derive directory path in a cross‑platform way.
    std::filesystem::path p(filename);
    std::string dir_path = p.parent_path().empty() ? "." : p.parent_path().string();
    std::cout << "Directory path: " << dir_path << std::endl;

    Graph graph;
    // Try to load BBkC binaries from the same directory as the provided .grh path
    {
        std::filesystem::path bdeg = std::filesystem::path(dir_path) / "b_degree.bin";
        std::filesystem::path badj = std::filesystem::path(dir_path) / "b_adj.bin";
        if (std::filesystem::exists(bdeg) && std::filesystem::exists(badj)) {
            std::cout << "Detected BBkC binaries in: " << dir_path << ", loading graph from binaries..." << std::endl;
            graph.read_graph_from_bbkcbinaries(dir_path, filename);
        } else {
            graph.read_graph_from_file(filename);
        }
    }
    
    // --- Core shrink mirroring FPCE roots: compute cores on our CSR,
    //     then build EBBkC CSR as the induced (R-1)-core subgraph (only if Turán enabled).
    graph.compute_core_numbers_from_csr();  // fills graph.core_numbers and graph.degeneracy
    std::cout << "degeneracy: " << graph.degeneracy << std::endl;
    if (g_enable_turan && R >= 3) {
        const int n = graph.total_nodes;
        std::vector<int> keep(n, 0), newid(n, -1), rev;
        int kept = 0;
        for (int u = 0; u < n; ++u) { // iterating over all nodes
            if (u < (int)graph.core_numbers.size() && graph.core_numbers[u] >= (R - 1)) {
                keep[u] = 1; // keep node u if it is in the (R-1)-core
                newid[u] = kept++; // assign a new index to node u
            }
        }
        // Build induced CSR into EBBkC arrays
        std::vector<uint32_t> r_off((size_t)kept + 1u, 0u); // offset array for the induced CSR
        std::vector<int>      r_dst; // destination array for the induced CSR
        r_dst.reserve(graph.csr_neighbors.size());
        rev.assign(kept, -1);
        for (int u = 0; u < n; ++u) { // iterating over all nodes
            if (!keep[u]) continue; // u is not in R-1 core, so skip
            rev[newid[u]] = u; // original vertex u gets newid[u] as its new index
            uint32_t begin = graph.csr_offsets[(size_t)u]; // begin of the neighbor list for node u
            uint32_t end   = graph.csr_offsets[(size_t)u + 1]; // end of the neighbor list for node u
            for (uint32_t e = begin; e < end; ++e) { // iterating over all neighbors of node u
                int v = (int)graph.csr_neighbors[e]; // neighbor of u
                if (v >= 0 && v < n && keep[v]) r_dst.push_back(newid[v]); // add the new index of the neighbor to the destination array
            }
            r_off[(size_t)(newid[u] + 1)] = (uint32_t)r_dst.size(); // set the offset for the next node
        }
        graph.csr_node_off.swap(r_off); // swap the offset array with the induced CSR
        graph.csr_edge_dst.swap(r_dst); // swap the destination array with the induced CSR
        graph.csr_n = (uint32_t)kept; // set the number of nodes in the induced CSR
        graph.csr_m = (uint32_t)graph.csr_edge_dst.size(); // set the number of edges in the induced CSR
        graph.csr_rev_map.swap(rev); // swap the reverse map with the induced CSR
        std::cout << "EBBkC CSR restricted to (R-1)-core: n=" << graph.csr_n
                  << " m=" << graph.csr_m << std::endl;
    }
    
    // If maximum is not specified, use total vertices (existing code below also does this)
    if (maximum == std::numeric_limits<int>::max()) {
        maximum = graph.total_nodes;
    }

    // 1) Compute ξ_G from CSR, then μ(θ, ξ_G)
    if (g_enable_order_bound || g_enable_edge_bound) {
        graph.compute_core_numbers_from_csr();
    }
    int mu = 0; // maximum possible size of pseudo-clique
    if (g_enable_order_bound) {
        mu = graph.compute_order_bound(theta, graph.degeneracy);
        std::cout << "Order bound μ = " << mu << " (degeneracy ξ_G = " << graph.degeneracy << ")\n";
        if (minimum > mu) {
            std::cout << "Pruning by Order Bound: no (ℓ,θ)-pseudo-clique can exist since ℓ (" 
                    << minimum << ") > μ (" << mu << ")\n";
            return 0;
        }
        if (maximum > mu) maximum = mu;
    } else {
        std::cout << "Order bound disabled" << std::endl;
    }

    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum, mu); 

    // Seeding strategy: Turán via EBBkC when enabled and R >= 3; otherwise node-by-node
    if (g_enable_turan && R >= 3 && graph.csr_n > 0) {
        const int n = static_cast<int>(graph.csr_n); // number of nodes in the induced CSR
        const int m = static_cast<int>(graph.csr_m); // number of edges in the induced CSR
        EBBkC_t::list_k_clique_mem_stream_from_csr(
            n,
            m,
            graph.csr_node_off.data(), // offset array for the induced CSR
            graph.csr_edge_dst.data(), // destination array for the induced CSR
            R, /*threads=*/2,
            [&](const std::vector<int>& r_clique){ // callback function for the EBBkC library
                if (!graph.csr_rev_map.empty()) {
                    std::vector<int> mapped;
                    mapped.reserve(r_clique.size()); // reserve space for the mapped clique
                    for (int id : r_clique) mapped.push_back(graph.csr_rev_map[id]); // map the new index to the original index
                    PC.enumerate_with_turan(mapped); // enumerate the clique with the mapped indices
                } else {
                    PC.enumerate_with_turan(r_clique); // enumerate the clique with the original indices
                }
            }
        );
        graph.csr_node_off.clear(); graph.csr_node_off.shrink_to_fit(); // clear the offset array
        graph.csr_edge_dst.clear(); graph.csr_edge_dst.shrink_to_fit(); // clear the destination array
    } else {
        for (int node = 0; node < graph.total_nodes; ++node) { // iterate over all nodes
            PC.add_to_inside_P(node);
            PC.remove_from_inside_P(node);
        }
    }

    std::vector<int> pseudo_clique_counts = PC.get_pseudo_cliques_count();

    // Print the clique sizes
    std::cout << "Maximal pseudo-clique counts:" << std::endl;
    bool found_any = false;
    long long total_found = 0;
    for (size_t sz = 0; sz < pseudo_clique_counts.size(); ++sz) {
        if (pseudo_clique_counts[sz] > 0) {
            std::cout << "Size " << sz << ": " << pseudo_clique_counts[sz] << "\n";
            found_any = true;
            total_found += pseudo_clique_counts[sz];
        }
    }
     
    if (!found_any) {
        std::cout << "(No pseudo-cliques found meeting the criteria)" << std::endl;
    }
 
    std::cout << "Total Maximal Pseudo-Cliques: " << total_found << std::endl;
    std::cout << "\nTotal Iterations: " << PC.get_iter_count() << std::endl;
    std::cout << "Edge-bound prunes saved: " << PC.get_numcalls_saved_by_edge_bound() << std::endl;
    const auto& cup_stats = PC.get_cup_stats();
    std::cout << "CUP mode: " << cup_mode_name(g_cup_mode)
              << " (runtime " << (PC.is_cup_runtime_disabled() ? "disabled" : "active") << ")" << std::endl;
    std::cout << "CUP nodes with previous uncles: " << cup_stats.nodes_with_previous_uncles << std::endl;
    std::cout << "CUP nodes without previous uncles: " << cup_stats.nodes_without_previous_uncles << std::endl;
    std::cout << "CUP assist attempts: " << cup_stats.assist_attempts << std::endl;
    std::cout << "CUP assist successes: " << cup_stats.assist_successes << std::endl;
    std::cout << "CUP fallback validations: " << cup_stats.fallback_validations << std::endl;
    std::cout << "CUP direct accepts: " << cup_stats.direct_accepts << std::endl;
    std::cout << "CUP recovered baseline children: " << cup_stats.baseline_children_recovered << std::endl;
    std::cout << "CUP covered low buckets: " << cup_stats.covered_low_buckets << std::endl;
    std::cout << "CUP covered eq-delta buckets: " << cup_stats.covered_eq_delta_buckets << std::endl;
    std::cout << "CUP covered eq-delta+1 buckets: " << cup_stats.covered_eq_delta_plus1_buckets << std::endl;
    std::cout << "CUP compare parity checks: " << cup_stats.compare_parity_checks << std::endl;
    std::cout << "CUP avg shortlist size: "
              << (cup_stats.assist_attempts
                    ? static_cast<double>(cup_stats.total_shortlist_size) /
                          static_cast<double>(cup_stats.assist_attempts)
                    : 0.0)
              << std::endl;
    std::cout << "CUP avg exact window size: "
              << (PC.get_iter_count()
                    ? static_cast<double>(cup_stats.total_child_window_size) /
                          static_cast<double>(PC.get_iter_count())
                    : 0.0)
              << std::endl;
    std::cout << "CUP avg relevant-uncle count: "
              << (cup_stats.assist_attempts
                    ? static_cast<double>(cup_stats.total_relevant_uncle_count) /
                          static_cast<double>(cup_stats.assist_attempts)
                    : 0.0)
              << std::endl;
    std::cout << "CUP avg changed-footprint: "
              << (cup_stats.assist_attempts
                    ? static_cast<double>(cup_stats.total_changed_footprint) /
                          static_cast<double>(cup_stats.assist_attempts)
                    : 0.0)
              << std::endl;
    std::cout << "CUP auto-disables: " << cup_stats.auto_disables << std::endl;

    return 0;
}

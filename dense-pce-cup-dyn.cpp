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
static bool g_enable_uncle_pruning = true; // Conditional Uncle Pruning (parent-feasible pool)

// Dynamic CUP gating parameters.
//
// Gate 1 — NP-spread check (O(1) before pool build):
//   Only build parent pool when maxNPdeg <= deltaP + g_cup_np_slack.
//   Rationale: parent_pool covers [min_add_deg .. maxNPdeg].
//              baseline scan covers [min_add_deg .. deltaP+1].
//   When maxNPdeg > deltaP+slack the pool is a strict superset of what baseline
//   scans → CUP will always cost more than baseline → skip pool build entirely.
//   slack=0: mathematically exact (pool ⊆ baseline iff maxNPdeg ≤ deltaP+1).
//   slack>0: allows slight overlap for cases where condition-C prunes the pool.
static int g_cup_np_slack = 0;   // set via --cup-slack=N

// Gate 0 — degeneracy ceiling (checked once at startup, before any iteration):
//   On graphs with high degeneracy ξ_G, the NP degree spread is always wide
//   (maxNPdeg >> deltaP+1 at almost every node).  Gates 1 and 2 would fire on
//   nearly every iteration, paying gate-check overhead while CUP fires rarely.
//   If ξ_G > g_cup_max_degeneracy, disable uncle pruning globally for this graph.
//   Data: tech ξ=23 (regression), bio-grid-human ξ=12 (improvement).
//   Default: 15.  Set to INT_MAX to disable Gate 0.
static int g_cup_max_degeneracy = 15;  // set via --cup-max-deg=N

// (g_cup_pool_ratio kept for CLI compatibility but no longer used in hot path)
static double g_cup_pool_ratio = 1.5;

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

    // ---- Conditional Uncle Pruning (CUP) support ----
    // We pass a pointer to the PARENT's feasible-extension pool down the recursion
    // stack. A child can then generate its candidates from this pool under the
    // CUP condition, instead of scanning NP buckets again.
    std::vector<const std::vector<int>*> parent_pool_stack;

    // O(1) neighbor test for the current v* (built only when CUP is active).
    std::vector<uint32_t> nbr_mark;
    uint32_t nbr_token = 1;

    // O(1) membership test for "u is already in parent pool" to shrink the
    // neighbor-supplement scan and avoid duplicates.
    std::vector<uint32_t> pool_mark;
    uint32_t pool_token = 1;

    inline const std::vector<int>* parent_feasible_pool() const {
        return parent_pool_stack.empty() ? nullptr : parent_pool_stack.back();
    }

    inline void bump_token(std::vector<uint32_t>& mark, uint32_t& token) {
        // extremely rare wrap; reset safely
        if (++token == 0) {
            std::fill(mark.begin(), mark.end(), 0);
            token = 1;
        }
    }

    inline void mark_neighbors_of(int x) {
        bump_token(nbr_mark, nbr_token);
        const uint32_t b = graph.csr_offsets[static_cast<size_t>(x)];
        const uint32_t e = graph.csr_offsets[static_cast<size_t>(x) + 1];
        for (uint32_t i = b; i < e; ++i) {
            const uint32_t w = graph.csr_neighbors[i];
            if (w < nbr_mark.size()) nbr_mark[w] = nbr_token;
        }
    }

    inline bool is_neighbor_of_vstar(int u) const {
        return (u >= 0 && static_cast<size_t>(u) < nbr_mark.size() && nbr_mark[static_cast<size_t>(u)] == nbr_token);
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
    unsigned long long numcalls_saved_by_uncle_pruning = 0ULL;

    // Dynamic CUP diagnostics
    // Gate 1: pool build skipped because proxy condition C was false at parent
    unsigned long long numcalls_gate1_skipped = 0ULL;
    // Gate 2: pool was built but too large vs. baseline estimate; fell back to baseline
    unsigned long long numcalls_gate2_fallback = 0ULL;

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

        // CUP: marker arrays for O(1) neighbor + pool membership checks (reused via token)
        nbr_mark.assign(static_cast<size_t>(graph.total_nodes), 0);
        pool_mark.assign(static_cast<size_t>(graph.total_nodes), 0);
        nbr_token = 1;
        pool_token = 1;
        parent_pool_stack.clear();
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

    inline void build_parent_feasible_pool(std::vector<int>& out, int min_add_deg) const {
        out.clear();
        if (min_add_deg < 0) min_add_deg = 0;
        if (maxNPdeg < min_add_deg) return;

        // Micro-opt #1: reserve size using countNP histogram.
        int reserve_n = 0;
        for (int d = min_add_deg; d <= maxNPdeg && d < (int)countNP.size(); ++d) {
            reserve_n += countNP[d];
        }
        if (reserve_n > 0) out.reserve(static_cast<size_t>(reserve_n));

        for (int d = min_add_deg; d <= maxNPdeg && d < (int)headNP.size(); ++d) {
            for (int u = headNP[d]; u != -1; u = nextNP[u]) {
                out.push_back(u);
            }
        }
    }


    void add_to_inside_P(int v, const std::vector<int>* parent_pool = nullptr);
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
    unsigned long long get_numcalls_saved_by_uncle_pruning() const {
        return numcalls_saved_by_uncle_pruning;
    }
    unsigned long long get_numcalls_gate1_skipped() const {
        return numcalls_gate1_skipped;
    }
    unsigned long long get_numcalls_gate2_fallback() const {
        return numcalls_gate2_fallback;
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

void PseudoCliqueEnumerator::add_to_inside_P(int v, const std::vector<int>* parent_pool) {
    parent_pool_stack.push_back(parent_pool);
    add_vertex_internal(v);   // state updates only
    iter(v);
    parent_pool_stack.pop_back();
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

    // CUP: build the parent-feasible pool for THIS node (passed to its children).
    // Pool(P) := { u ∉ P : deg_P(u) ≥ ceil(theta_P) }, i.e. all vertices that are pseudo-clique feasible w.r.t. P.
    // This pool is computed once per node and is re-used by all children for conditional pruning.
    //
    // GATE 1 — NP-spread check (O(1), fires before any pool allocation):
    //
    // parent_pool covers NP vertices with deg in [min_add_deg .. maxNPdeg].
    // Baseline scan covers NP vertices with deg in [min_add_deg .. deltaP+1].
    //
    // Therefore parent_pool ⊆ baseline  iff  maxNPdeg <= deltaP + 1.
    // When maxNPdeg > deltaP + g_cup_np_slack the pool is a STRICT SUPERSET of what
    // baseline would scan, so CUP cannot be cheaper — skip pool build entirely.
    //
    // This is the correct gating condition.  The previous proxy (theta*(k+1) >= degP)
    // was vacuously true for all valid PCs and never fired (0 skips on all benchmarks).
    std::vector<int> feasible_pool;
    bool pool_built = false;
    if (g_enable_uncle_pruning) {
        if (maxNPdeg <= deltaP + 1 + g_cup_np_slack) {
            build_parent_feasible_pool(feasible_pool, min_add_deg);
            pool_built = true;
        } else {
            numcalls_gate1_skipped++;  // maxNPdeg spread too wide; pool would exceed baseline
        }
    }

    std::vector<int> children_vec;
    children_vec.reserve(64);

    auto try_add_child = [&](int u, bool u_adj_v) {
        if (tracks[u].inP) return;
        int deg = static_cast<int>(tracks[u].degNP);
        if (deg < min_add_deg) return;

        // Type 1: deg < δ(P)
        if (deg < deltaP) {
            children_vec.push_back(u);
            return;
        }

        // Type 2: deg = δ(P)
        if (deg == deltaP) {
            if (u < v) {
                children_vec.push_back(u);
                return;
            }
            if (!u_adj_v) return;
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (x < u && !graph.has_edge(u, x)) return;
            }
            children_vec.push_back(u);
            return;
        }

        // Type 3: deg = δ(P)+1
        if (deg == deltaP + 1) {
            if (!(u < v)) return;
            if (!u_adj_v) return;

            // Must connect to all vertices in P[δ(P)]
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (!graph.has_edge(u, x)) return;
            }

            // Must connect to all vertices in P[δ(P)+1] with id < u
            if (deltaP + 1 < (int)headP.size()) {
                for (int x = headP[deltaP + 1]; x != -1; x = nextP[x]) {
                    if (x < u && !graph.has_edge(u, x)) return;
                }
            }
            children_vec.push_back(u);
            return;
        }
    };

    bool used_cup = false;
    const std::vector<int>* parent_pool = (g_enable_uncle_pruning ? parent_feasible_pool() : nullptr);

    if (g_enable_uncle_pruning && parent_pool != nullptr) {
        // GATE 2 — child-level NP-spread check (O(1), replaces the old O(3)-loop):
        //
        // Gate 1 verified the PARENT's NP was flat (maxNPdeg_parent <= deltaP_parent+1).
        // But after inserting vertex u into P to form the child:
        //   - neighbors of u in NP get degNP += 1  →  maxNPdeg can increase by 1
        //   - δ(P') can DROP if u has fewer P-edges than the old minimum
        // Both effects together can widen the child's NP spread even when the parent was flat.
        //
        // Re-apply the same O(1) NP-spread check at the CHILD's current state.
        // If it fails, skip CUP and pay O(1) instead of the old O(3)-loop + fallback.
        // When it passes, parent_pool ≈ child's feasible NP range → no ratio loop needed.
        if (maxNPdeg > deltaP + 1 + g_cup_np_slack) {
            numcalls_gate2_fallback++;
            // child NP spread too wide; pool exceeds child's 3-bucket baseline → skip CUP
        } else {
            // Both Gate 1 (parent) and Gate 2 (child) confirmed flat NP.
            // Now check condition C: θk ≥ deg_{K\{v*(K)}}(v*(K)) + I
            const long double theta_k = (long double)theta * (long double)total_nodes_in_P;
            const int deg_parent_s = static_cast<int>(tracks[v].degP);  // deg_{K\{v}}(v)
            const bool cond0 = (theta_k + 1e-12L >= (long double)deg_parent_s);
            const bool cond1 = (theta_k + 1e-12L >= (long double)(deg_parent_s + 1));

            if (cond0 || cond1) {
                used_cup = true;
                numcalls_saved_by_uncle_pruning++;

                // Mark N(v) once for O(1) adjacency checks to v.
                mark_neighbors_of(v);
                bump_token(pool_mark, pool_token);

                // 1) Scan parent-feasible pool.
                //    - cond1: pool covers ALL candidates.
                //    - only cond0: pool covers all I=0 candidates; I=1 gap handled below.
                for (int u : *parent_pool) {
                    if (u < 0 || u >= graph.total_nodes) continue;
                    pool_mark[static_cast<size_t>(u)] = pool_token;
                    if (tracks[u].inP) continue;
                    const bool u_adj_v = is_neighbor_of_vstar(u);
                    try_add_child(u, u_adj_v);
                }

                // 2) If cond1 is false, add the I=1 gap by scanning neighbors of v.
                if (!cond1) {
                    const uint32_t beg = graph.csr_offsets[static_cast<size_t>(v)];
                    const uint32_t endN = graph.csr_offsets[static_cast<size_t>(v) + 1];
                    for (uint32_t idx = beg; idx < endN; ++idx) {
                        const int u = static_cast<int>(graph.csr_neighbors[idx]);
                        if (u < 0 || u >= graph.total_nodes) continue;
                        if (pool_mark[static_cast<size_t>(u)] == pool_token) continue;
                        if (tracks[u].inP) continue;
                        try_add_child(u, true);
                    }
                }
            }
        } // end Gate 2 child NP-spread check
    } // end if (g_enable_uncle_pruning && parent_pool != nullptr)

    if (!used_cup) {
        // Original bucket-based child generation (baseline).

        // Type 1: deg in [min_add_deg, δ(P)-1]
        for (int deg = std::max(0, min_add_deg); deg <= deltaP - 1; ++deg) {
            if (deg < 0 || deg >= (int)headNP.size()) continue;
            for (int u = headNP[deg]; u != -1; u = nextNP[u]) {
                children_vec.push_back(u);
            }
        }

        // Type 2: deg = δ(P)
        if (deltaP >= 0 && deltaP < (int)headNP.size()) {
            for (int u = headNP[deltaP]; u != -1; u = nextNP[u]) {
                if (tracks[u].degNP < min_add_deg) continue;
                if (u < v) {
                    children_vec.push_back(u);
                } else {
                    if (graph.has_edge(v, u)) {
                        bool ok = true;
                        for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                            if (x < u && !graph.has_edge(u, x)) { ok = false; break; }
                        }
                        if (ok) children_vec.push_back(u);
                    }
                }
            }
        }

        // Type 3: deg = δ(P)+1
        if (deltaP + 1 >= 0 && deltaP + 1 < (int)headNP.size()) {
            for (int u = headNP[deltaP + 1]; u != -1; u = nextNP[u]) {
                if (tracks[u].degNP < min_add_deg) continue;
                if (graph.has_edge(v, u) && v > u) {
                    bool ok = true;
                    for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                        if (!graph.has_edge(u, x)) { ok = false; break; }
                    }
                    if (ok) {
                        for (int x = headP[deltaP + 1]; x != -1; x = nextP[x]) {
                            if (x < u && !graph.has_edge(u, x)) { ok = false; break; }
                        }
                    }
                    if (ok) children_vec.push_back(u);
                }
            }
        }
    }

    // CUP may introduce duplicates; always normalize expansion list deterministically.
    if (children_vec.size() > 1) {
        std::sort(children_vec.begin(), children_vec.end());
        children_vec.erase(std::unique(children_vec.begin(), children_vec.end()), children_vec.end());
    }

    // Local order-bound stop at μ: do not expand past μ
    if (g_enable_order_bound && order_ub > 0 && total_nodes_in_P == order_ub) return;

    // Deterministic expansion order
    for (int u : children_vec) {
        add_to_inside_P(u, (g_enable_uncle_pruning && pool_built) ? &feasible_pool : nullptr);
        remove_from_inside_P(u);
    }
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename> [--theta <theta>] [--minimum <min>] [--maximum <max>]" << std::endl;
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
        } else if (arg == "--no-uncle") {
            g_enable_uncle_pruning = false;
        } else if (arg == "--uncle") {
            g_enable_uncle_pruning = true;
        } else if (arg.rfind("--cup-ratio=", 0) == 0) {
            // e.g. --cup-ratio=1.5  (Gate 2 threshold; large value disables Gate 2)
            g_cup_pool_ratio = std::stod(arg.substr(12));
        } else if (arg.rfind("--cup-slack=", 0) == 0) {
            // e.g. --cup-slack=0  (Gate 1&2 slack; 0=exact, 1=allow 1 extra NP bucket)
            g_cup_np_slack = std::stoi(arg.substr(12));
        } else if (arg.rfind("--cup-max-deg=", 0) == 0) {
            // e.g. --cup-max-deg=15  (Gate 0: disable CUP if ξ_G exceeds this)
            // Use --cup-max-deg=999 to disable Gate 0 entirely.
            g_cup_max_degeneracy = std::stoi(arg.substr(14));
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
              << " uncle:" << (g_enable_uncle_pruning?"on":"off")
              << " cup-ratio:" << g_cup_pool_ratio
              << " cup-slack:" << g_cup_np_slack
              << " cup-max-deg:" << g_cup_max_degeneracy << std::endl;

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

    // GATE 0 — degeneracy ceiling: disable CUP on dense graphs where NP spread
    // is structurally guaranteed to be wide (ξ_G > threshold).
    // On such graphs Gates 1+2 fire on nearly every node, paying overhead for
    // pool builds and child NP-spread checks while CUP activates < 1% of the time.
    if (g_enable_uncle_pruning && graph.degeneracy > g_cup_max_degeneracy) {
        std::cout << "Gate0: degeneracy ξ_G=" << graph.degeneracy
                  << " > threshold " << g_cup_max_degeneracy
                  << " → uncle pruning disabled for this graph." << std::endl;
        g_enable_uncle_pruning = false;
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
    std::cout << "Edge-bound prunes saved:    " << PC.get_numcalls_saved_by_edge_bound() << std::endl;
    std::cout << "Uncle-pruning CUP hits:     " << PC.get_numcalls_saved_by_uncle_pruning() << std::endl;
    std::cout << "  Gate1 skip (parent NP wide): " << PC.get_numcalls_gate1_skipped()
              << "  (parent maxNPdeg > deltaP+1+" << g_cup_np_slack << "; pool build skipped)" << std::endl;
    std::cout << "  Gate2 skip (child NP wide) : " << PC.get_numcalls_gate2_fallback()
              << "  (child maxNPdeg > deltaP+1+" << g_cup_np_slack << " after vertex insertion; CUP skipped)" << std::endl;

    return 0;
}

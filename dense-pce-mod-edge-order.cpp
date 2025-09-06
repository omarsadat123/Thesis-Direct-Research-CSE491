#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>
#include <stack>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <climits>
#include <cstring>
#include <cerrno>
#include <filesystem>
#include <set> 
#include "EBBkC/src/edge_oriented.h"

// Globals required by EBBkC in library mode
int K = 0;
int L = 2;
unsigned long long N = 0ULL;

// Define PATH_MAX if not already defined
// #ifndef PATH_MAX
// #define PATH_MAX 4096
// #endif

std::vector<std::vector<int>> read_clique_file(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening clique file: " << path << std::endl;
        return {};
    }
    std::cout << "Reading clique file: " << path << std::endl;
    std::vector<std::vector<int>> cliques;
    std::string line;
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<int> clique;
        int node;
        while (ss >> node) {
            clique.push_back(node);
        }
        cliques.push_back(clique);
    }
    return cliques;
}

// counting sort with tracklist
class Graph {
public:
    Graph() = default;
    // Compact adjacency: sorted vector per vertex
    std::vector<std::vector<int>> adj_map;

    int degeneracy = 0; // Store ξ_G (maximum coreness)
    int total_nodes = 0 ;
    // CSR view for EBBkC zero-copy entry
    std::vector<uint32_t> csr_node_off;
    std::vector<int> csr_edge_dst;
    uint32_t csr_n = 0;
    uint32_t csr_m = 0;

    void add_edge(int u, int v) {
        if (u >= (int)adj_map.size() || v >= (int)adj_map.size()) {
            adj_map.resize(std::max(u, v) + 1);
        }
        adj_map[u].push_back(v);
        adj_map[v].push_back(u);
    }

    void add_empty_edge(int u) { //placeholder for empty edge
        adj_map[u][-1] = false;        
    }

    void print_graph() const { //prints out each vertex and its list of neighbors
        for (size_t i = 0; i < adj_map.size(); ++i) {
            const auto& neighbors = adj_map[i];
            std::cout << "Node " << i << " : ";
            for (int neighbor : neighbors) {
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }
    }

    // void print_graph() const {
    //     for (const auto& node : adj_map) {
    //         std::vector<int> neighbors;
    //         for (const auto& neighbor : node.second) {
    //             neighbors.push_back(neighbor.first);
    //         }
    //         std::sort(neighbors.begin(), neighbors.end());

    //         std::cout << "Node " << node.first << " : ";
    //         for (int neighbor : neighbors) {
    //             std::cout << neighbor << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    std::string graph_file_path;
    void finalize_adjacency() {
        for (auto& nbrs : adj_map) {
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        }
    }

    bool has_edge(int u, int v) const {
        if (u < 0 || v < 0 || u >= (int)adj_map.size()) return false;
        const auto& nbrs = adj_map[u];
        return std::binary_search(nbrs.begin(), nbrs.end(), v);
    }

    bool read_graph_from_file(const std::string& filename) {
        graph_file_path = filename;   // <- add this line
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        int current_vertex = 0;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line == "[EOF]") break;
            total_nodes++;
            if (line.empty()) { //just an isolated vertex
                // add_empty_edge(current_vertex);
                current_vertex++;
                continue;
            }

            std::stringstream ss(line); //helps split the line of text into individual numbers
            int neighbor;
            while (ss >> neighbor) {
                add_edge(current_vertex, neighbor);
            }
            current_vertex++;
        }
        finalize_adjacency();
        return true;
    }

    // Load graph from BBkC binaries produced by Cohesive_subgraph_book (b_degree.bin, b_adj.bin)
    // dir should be the directory containing the binaries; original_graph_path is kept for downstream path logic
    bool read_graph_from_bbkcbinaries(const std::string& dir, const std::string& original_graph_path) {
        graph_file_path = original_graph_path;

        std::filesystem::path deg_path = std::filesystem::path(dir) / "b_degree.bin";
        std::filesystem::path adj_path = std::filesystem::path(dir) / "b_adj.bin";

        std::ifstream deg_file(deg_path, std::ios::binary);
        std::ifstream adj_file(adj_path, std::ios::binary);
        if (!deg_file.is_open() || !adj_file.is_open()) {
            std::cerr << "Error: Could not open BBkC binaries in directory " << dir << std::endl;
            return false;
        }

        // b_degree.bin layout: [ui tt][ui n][ui m][ui degree[0..n-1]]
        uint32_t tt = 0;
        uint32_t n = 0;
        uint32_t m = 0;
        deg_file.read(reinterpret_cast<char*>(&tt), sizeof(uint32_t));
        deg_file.read(reinterpret_cast<char*>(&n), sizeof(uint32_t));
        deg_file.read(reinterpret_cast<char*>(&m), sizeof(uint32_t));
        if (!deg_file.good()) {
            std::cerr << "Error: Failed reading degree header from " << deg_path.string() << std::endl;
            return false;
        }
        if (tt != sizeof(uint32_t)) {
            std::cerr << "Warning: Unexpected ui size header in b_degree.bin (" << tt << ")" << std::endl;
        }

        std::vector<uint32_t> degrees(n);
        if (n > 0) {
            deg_file.read(reinterpret_cast<char*>(degrees.data()), static_cast<std::streamsize>(n * sizeof(uint32_t)));
            if (!deg_file.good()) {
                std::cerr << "Error: Failed reading degree array from " << deg_path.string() << std::endl;
                return false;
            }
        }

        // Prepare adjacency container
        adj_map.clear();
        adj_map.resize(n);
        total_nodes = static_cast<int>(n);
        // Prepare CSR arrays for EBBkC
        csr_n = n;
        csr_m = m;
        csr_node_off.assign(n + 1, 0);
        for (uint32_t i = 0; i < n; ++i) csr_node_off[i + 1] = csr_node_off[i] + degrees[i];
        if (csr_node_off[n] != m) {
            std::cerr << "Warning: degree prefix sum != m (" << csr_node_off[n] << " vs " << m << ")" << std::endl;
            csr_m = csr_node_off[n];
        }
        csr_edge_dst.resize(csr_m);

        // Read neighbors sequentially from b_adj.bin; length should be m entries
        uint64_t consumed = 0;
        uint64_t fill_idx = 0;
        for (uint32_t u = 0; u < n; ++u) {
            uint32_t du = degrees[u];
            for (uint32_t k = 0; k < du; ++k) {
                uint32_t v;
                adj_file.read(reinterpret_cast<char*>(&v), sizeof(uint32_t));
                if (!adj_file.good()) {
                    std::cerr << "Error: Unexpected EOF in b_adj.bin while reading neighbors (u=" << u << ")" << std::endl;
                    return false;
                }
                // Add as undirected; b_adj already contains both directions, but map insert is idempotent
                add_edge(static_cast<int>(u), static_cast<int>(v));
                consumed++;
                if (fill_idx < csr_edge_dst.size()) csr_edge_dst[fill_idx++] = static_cast<int>(v);
            }
        }

        // Optional sanity check
        if (consumed != m) {
            std::cerr << "Warning: Consumed neighbor count (" << consumed << ") != header m (" << m << ")" << std::endl;
        }
        finalize_adjacency();
        return true;
    }

    std::vector<int> core_numbers;
    void compute_core_numbers() {
        unsigned int n = adj_map.size();
        core_numbers.resize(n);  // Resize the member vector
        std::vector<unsigned int> vert(n), pos(n), deg(n);
        unsigned int md = 0;

        // Calculate degrees and find maximum degree
        for (unsigned int i = 0; i < n; i++) {
            deg[i] = adj_map[i].size();
            if (deg[i] > md) {
                md = deg[i];
            }
        }

        // Initialize bin array
        std::vector<unsigned int> bin(md + 1, 0);
        
        // Count vertices of each degree
        for (unsigned int v = 0; v < n; v++) {
            bin[deg[v]]++;
        }

        // Compute starting positions in bins
        unsigned int start = 0;
        for (unsigned int d = 0; d <= md; d++) {
            unsigned int num = bin[d];
            bin[d] = start;
            start += num;
        }

        // Place vertices in vert array by degree
        for (unsigned int v = 0; v < n; v++) {
            pos[v] = bin[deg[v]];
            vert[pos[v]] = v;
            bin[deg[v]]++;
        }

        // Restore bin starting positions
        for (int d = md; d >= 1; d--) {
            bin[d] = bin[d-1];
        }
        bin[0] = 0;

        // Main core computation
        for (unsigned int i = 0; i < n; i++) {
            unsigned int v = vert[i];
            core_numbers[v] = deg[v];
            
            // Process neighbors
            for (int u : adj_map[v]) {
                if (deg[u] > deg[v]) {
                    unsigned int du = deg[u];
                    unsigned int pu = pos[u];
                    unsigned int pw = bin[du];
                    unsigned int w = vert[pw];
                    
                    // Swap vertices if needed
                    if (u != w) {
                        pos[u] = pw;
                        vert[pu] = w;
                        pos[w] = pu;
                        vert[pw] = u;
                    }
                    
                    bin[du]++;
                    deg[u]--;
                }
            }
        }

        // NEW: set ξ_G = max coreness
        int maxc = 0;
        for (int c : core_numbers) if (c > maxc) maxc = c;
        degeneracy = maxc;   // <-- important
        std::cout << "degeneracy: " << degeneracy << std::endl;
    }

    int compute_order_bound(double theta, int degeneracy) {
        // Use extended precision to avoid downward rounding
        const long double th = static_cast<long double>(theta);
        const long double xi = static_cast<long double>(degeneracy);
        const long double eps = 1e-12L;
    
        // Corollary 1: |S| ≤ ⌊ 2 ξ_G / θ ⌋
        long double b1_ld = std::floor((2.0L * xi) / th + eps);
        int bound1 = static_cast<int>(b1_ld);
    
        // Corollary 2: if θ > ξ_G / (ξ_G + 1),
        // |S| ≤ ⌊ 1 / (1 − ξ_G / ((ξ_G + 1) θ)) ⌋
        int bound2 = std::numeric_limits<int>::max();
        long double thresh = xi / (xi + 1.0L);
        if (th > thresh + eps) {
            long double denom = 1.0L - xi / ((xi + 1.0L) * th);
            // denom should be positive here; still guard against tiny negatives
            if (denom <= 0.0L) denom = eps;
            long double b2_ld = std::floor(1.0L / denom + eps);
            bound2 = static_cast<int>(b2_ld);
        }
    
        return std::min(bound1, bound2);
    }
    
};

class PseudoCliqueEnumerator {
private:
    Graph& graph; //Avoids copying potentially large graph data
    float theta; //Threshold for pseudo-clique detection
    int min_size; //Minimum size of pseudo-clique to consider
    int max_size; //Maximum size of pseudo-clique to consider
    int total_nodes; //Total number of vertices in the graph
    float theta_P; //Current threshold for pseudo-clique detection
    int total_nodes_in_P; //Total number of vertices in the current pseudo-clique P
    int total_edges_in_P; //Total number of edges in the current pseudo-clique P
    int iter_count = 0; //Counts the number of iterations of the algorithm

    // --- Edge-bound state ---
    int lhs_required_edges = 0;             // ceil(theta * l*(l-1)/2)
    std::vector<int> core_hist;             // histogram of core numbers in current P
    int cur_min_core = INT_MAX;             // c(P) = min core in P
    int max_core_seen = 0;                  // size for core_hist

    // --- accounting like fpce.c
    unsigned long long numcalls_saved_by_edge_bound = 0ULL;

    int min_degree_in_P_fast() const {
        if (total_nodes_in_P == 0) return 0;
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            if (!inside_P[d].empty()) return d;       // first non-empty degree bucket
        }
        return 0;
    }
    int current_min_core_in_P_fast() const {
        return (cur_min_core == INT_MAX) ? 0 : cur_min_core;
    }

    std::stack<int> children; //Stack of vertices to be processed

    std::vector<std::vector<int> > inside_P; //Groups vertices currently in pseudo-clique P by their degree within P
    
    std::vector<std::vector<int> > neighbors_and_P; //Groups ALL vertices by their degree with respect to current pseudo-clique P
    std::vector<std::vector<int> > pseudo_cliques; //Stores actual pseudo-cliques found, grouped by size
    //[0,1,2,3,4,[8,7,10,14,4]...]
    std::vector<int> pseudo_cliques_count; //Counts how many pseudo-cliques exist for each size

    // std::unordered_map<int, TrackInfo> tracks;
    std::vector<std::vector<int> > tracks; //Tracks the state of each vertex in the algorithm
    // tracks[0] : is v currently in P?
    // tracks[1] : what's v's degree inside P?
    // tracks[2] : How many neighbors does v have in P?
    // tracks[3] : what is v's position inside its inside_P bucket? [[],[8,7,12,17],...]
    // tracks[4] : what is v's position inside its neighbors_and_P bucket?


    int R; // Size of seed cliques (r-cliques)
    //bool maximal_only = true; // Always enumerate maximal pseudo-cliques
    // for duplicate checking (maximal mode only)
    std::set<std::vector<int>> found_maximal_cliques; // CHANGED: unordered_set to set
    void add_vertex_internal(int v) {
        if (total_nodes_in_P > max_size) {
            return;
        }
        
        //Updating v's info
        tracks[v][0] = true;                             // Mark v as being part of P
        tracks[v][1] = tracks[v][2];                     // Its new "degree_in_P" is its old degree_with_P"
        inside_P[tracks[v][1]].push_back(v);             // Add v to the correct degree bucket in the inside_P list
        tracks[v][3] = inside_P[tracks[v][1]].size() - 1; // Record its position in that bucket  
        
        // track c(P)
        int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0;
        if ((int)core_hist.size() <= cv) core_hist.resize(cv + 1, 0);
        core_hist[cv] += 1;
        if (cv < cur_min_core) cur_min_core = cv;

    
        //Working for v's neighbors
        for (int adj_v : graph.adj_map[v]) {
            // std::cout << "    ite: " << adj_v  << "\n";
    
            //Process neighbors of vertex v that are already in P
            if (tracks[adj_v][0]) { // If neighbor is already in P
                total_edges_in_P += 1; // New edge within P
                //code to move adj_v to a higher degree bucket in inside_P
                int degree = tracks[adj_v][1]; // Get degree of neighbor inside P
                tracks[inside_P[degree].back()][3] = tracks[adj_v][3]; // Update position of neighbor in P
    
                // Swap and pop
                std::swap(inside_P[degree][tracks[adj_v][3]], inside_P[degree].back()); // Swap neighbor with last vertex in degree bucket
                inside_P[degree].pop_back(); // Remove last vertex from degree bucket
    
                inside_P[degree + 1].push_back(adj_v); // Add neighbor to next degree bucket
                tracks[adj_v][3] = inside_P[degree + 1].size() - 1;
                tracks[adj_v][1] = degree + 1;
            }
    
            // Get the current "degree with respect to P" and position for the neighbor adj_v
            int degree = tracks[adj_v][2]; // Get degree of neighbor wrt P
            int pos = tracks[adj_v][4]; // Get position of neighbor in neighbors_and_P
            int last_node = neighbors_and_P[degree].back(); // Get last vertex in degree bucket
    
            // Swap and pop
            std::swap(neighbors_and_P[degree][pos], neighbors_and_P[degree].back());
    
            // std::cout << "       " << last_node << "\n";
            
            tracks[last_node][4] = pos;
            // std::cout << "-----" << "\n";
    
            neighbors_and_P[degree].pop_back();
    
            neighbors_and_P[degree + 1].push_back(adj_v);
            tracks[adj_v][4] = neighbors_and_P[degree + 1].size() - 1;
            tracks[adj_v][2] = degree + 1;
        }
    
        // Update global statistics
        total_nodes_in_P += 1;
        set_theta_P(); //Update theta_P
    
        // Check if current set qualifies as a pseudo-clique
        // if (total_nodes_in_P >= min_size) {
        //     // std::vector<int> pseudo_clique;
        //     // for (const auto& nodes : inside_P) {
        //     //     pseudo_clique.insert(pseudo_clique.end(), nodes.begin(), nodes.end());
        //     // }
        //     // pseudo_cliques.push_back(pseudo_clique);
        //     pseudo_cliques_count[total_nodes_in_P] += 1;
        // }
    }

    bool prune_by_edge_bound() const {
        const int S = min_size;                 // ℓ
        const int P = total_nodes_in_P;         // |P|
        const int t = S - P;                    // vertices still needed

        // Only prune when |P|≥R and |P|<ℓ. (Call-site also checks this.)
        if (P < 1 || t <= 0) return false;

        const int delta       = min_degree_in_P_fast();         // δ(P)
        const int cP          = current_min_core_in_P_fast();   // c(P)
        const double sumE     = static_cast<double>(total_edges_in_P);
        const int LHS         = lhs_required_edges;

        // Lemma 9 shape coded in fpce.c
        const int g = cP - delta - 1;

        if (g >= 0) {
            // max_value = min( c(P), δ(P) + t )
            const int max_value = std::min(cP, delta + t);
            // RHS1 = |E[P]| + max_value*(t + δ + 0.5) - 0.5*(δ+1)*(2δ+1)
            const double rhs1 = std::floor(
                sumE + max_value * (t + delta + 0.5)
                    - 0.5 * (delta + 1) * (2 * delta + 1)
            );
            if (LHS > rhs1) return true;   // EDGE_BOUND1
        } else {
            // c(P) == δ(P): RHS2 = |E[P]| + t*δ
            const double rhs2 = std::floor(sumE + t * static_cast<double>(delta));
            if (LHS > rhs2) return true;   // EDGE_BOUND2
        }

        // Lemma 8 fallback: RHS3 = |E[P]| + t*(δ + (t+1)/2)
        const double rhs3 = std::floor(sumE + t * (delta + (t + 1) * 0.5));
        if (LHS > rhs3) return true;       // EDGE_BOUND3

        return false;
    }



public:
    PseudoCliqueEnumerator(Graph& graph, float theta = 1.0, int min_size = 1, int max_size = -1, int R = 0)
        : graph(graph), theta(theta), min_size(min_size), R(R) {
        
        total_nodes = graph.adj_map.size(); 

        std::cout << "total nodes: " << total_nodes << " "   << "\n";

        this->max_size = max_size != -1 ? max_size : total_nodes; //max_size to total nodes if not specified

        pseudo_cliques.resize(max_size + 1); //Stores actual pseudo-cliques found, grouped by size e.g. pseudo_cliques[2] = {{0,1}, {1,2}, {2,3}} pseudo_cliques[3] = {{0,1,2}, {1,2,3}}

        pseudo_cliques_count.resize(max_size + 1, 0); //Counts how many pseudo-cliques exist for each size

        inside_P.resize(max_size + 1); //Groups vertices currently in pseudo-clique P by their degree within P e.g. inside_P[2] = {{0,1}, {1,2}, {2,3}}

        neighbors_and_P.resize(max_size + 1); //Groups ALL vertices by their degree with respect to current pseudo-clique P e.g. neighbors_and_P[2] = {7, 8} Vertices 7,8 have 2 connections to P
        //initializing the degree-based tracking system
        std::vector<int> keys(graph.total_nodes); //keys = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        std::iota(keys.begin(), keys.end(), 0); //Fills the vector with sequential integers starting from 0


        neighbors_and_P[0] = keys; //All vertices have degree 0 with respect to empty pseudo-clique P

        // neighbors_and_P[0] = get_keys(graph.adj_map);      
        
        int idx = 0;
        for (const auto& v : neighbors_and_P[0]) {
            tracks.push_back({0, -1, 0, -1, idx}); //[is_in_P, degree_in_P, degree_in_NP, pos_in_P, pos_in_NP]
            ++idx;
        }

        theta_P = 0;
        total_nodes_in_P = 0;
        total_edges_in_P = 0;

        // LHS once (fpce.c keeps it global)
        lhs_required_edges = static_cast<int>(
            std::ceil(theta * (min_size * (min_size - 1) / 2.0))
        );

        // Prepare core histogram from graph.core_numbers
        max_core_seen = 0;
        for (int c : graph.core_numbers) max_core_seen = std::max(max_core_seen, c);
        core_hist.assign(max_core_seen + 1, 0);
        cur_min_core = INT_MAX;

    }

    void set_theta_P() { //θ(K) = θ * clq(|K| + 1) - |E[K]|
        theta_P = (theta * total_nodes_in_P * (total_nodes_in_P + 1) / 2.0) - total_edges_in_P;
    }

    void print_all() const { //debugging function, prints the current state of the algorithm
        std::cout << "P: "; //Shows vertices currently in the pseudo-clique P, grouped by their degree within P 
        for (const auto& p : inside_P) {
            for (int v : p) std::cout << v << " ";
            std::cout << "| "; //P: | 0 1 | 2 3 | |
        }
        std::cout << "\nNP: "; //Shows ALL vertices grouped by their total degree with respect to P (includes both vertices in P and outside P)
        for (const auto& np : neighbors_and_P) {
            for (int v : np) std::cout << v << " ";
            std::cout << "| "; //NP: 5 6 | | 0 1 | 2 3 4 | |
        }
        // std::cout << "\ntracks: ";
        // for (const auto& track : tracks) {
        //     std::cout << "\n\n" ;
        //     int key = track.first;
        //     auto& value = track.second;

        //     std::cout << key << ": {is_inside_P: " << value[0]
        //               << ", degree_in_P: " << value[1]
        //               << ", degree_in_NP: " << value[2]
        //               << ", position_inside_P: " << value[3]
        //               << ", position_inside_NP: " << value[4] << "} ";
        // }
        std::cout << "\ntheta_P: " << theta_P << " (" << total_nodes_in_P << ", " << total_edges_in_P << ")\n"; //theta_P: 2.5 (4, 5)
    }

    void add_to_inside_P(int v);
    void remove_from_inside_P(int v);
    void iter(int v);

    std::vector<int> get_pseudo_cliques_count(){
        return pseudo_cliques_count;
    };

    int get_iter_count(){
        return iter_count;
    };

    unsigned long long get_edge_prunes() const { return numcalls_saved_by_edge_bound; }

    // std::vector<int> get_keys(const std::vector<int, std::unordered_map<int, bool> >& map) {
    //     std::vector<int> keys;
    //     for (const auto& pair : map) {
    //         keys.push_back(pair.first);
    //     }
    //     std::sort(keys.begin(), keys.end());
    //     return keys;
    // }


    void enumerate_with_turan(const std::vector<int>& turan_seed);  // Enumerate pseudo-cliques using Turan's theorem
    // Enumerate pseudo-cliques using provided R-clique seeds (in-memory, no file I/O)
    void enumerate_with_seeds(const std::vector<std::vector<int>>& r_cliques);
    // ADD THIS NEW PUBLIC METHOD for reporting unique pseudo-cliques
    void report_if_new(const std::vector<int>& clique) {
        if (clique.size() < this->min_size) return;

        // Create a sorted copy to act as a canonical key
        std::vector<int> key = clique;
        std::sort(key.begin(), key.end());

        // The .second of the insert result is true if the insertion was new
        if (found_maximal_cliques.insert(key).second) {
            pseudo_cliques_count[clique.size()]++;
        }
    }

    // Non-maximal: count only if current set is canonical under seed partition
    // Non-maximal enumeration is disabled
    bool is_current_seed_canonical() const {
        // if (R <= 0 || static_cast<int>(current_seed_sorted_ids.size()) != R) return true;
        // std::vector<int> cur = collect_current_pseudo_clique();
        // if (static_cast<int>(cur.size()) < R) return true;
        // std::sort(cur.begin(), cur.end());
        // for (int i = 0; i < R; ++i) {
        //     if (cur[i] != current_seed_sorted_ids[i]) return false;
        // }
        return true;
    }
    
    // Helper to collect and sort current pseudo-clique P
    std::vector<int> collect_current_pseudo_clique() const {
        std::vector<int> current_p;
        current_p.reserve(total_nodes_in_P);
        for (const auto& bucket : inside_P) {
            if (!bucket.empty()) {
                current_p.insert(current_p.end(), bucket.begin(), bucket.end());
            }
        }
        std::sort(current_p.begin(), current_p.end());
        return current_p;
    }

    // Robust maximality check similar to pce.c
    bool is_current_maximal() const {
        int min_add_deg = static_cast<int>(std::ceil(theta_P));
        if (min_add_deg < 0) min_add_deg = 0;
        for (int d = min_add_deg; d <= total_nodes_in_P; ++d) {
            size_t occ_num = (d >= 0 && d < static_cast<int>(neighbors_and_P.size())) ? neighbors_and_P[d].size() : 0;
            size_t in_p_num = (d >= 0 && d < static_cast<int>(inside_P.size())) ? inside_P[d].size() : 0;
            if (occ_num > in_p_num) return false; // some outside vertex can extend P
        }
        return true;
    }
}; 

void PseudoCliqueEnumerator::enumerate_with_turan(const std::vector<int>& turan_seed) {
    if (turan_seed.empty()) {
        std::cout << "Turán path is empty; nothing to extend." << std::endl;
        return;
    }

    // A) Deterministic order: (deg, id)
    auto sorted_clique = turan_seed;
    std::sort(sorted_clique.begin(), sorted_clique.end(), [&](int a, int b) {
        return std::make_pair(graph.adj_map[a].size(), a) < std::make_pair(graph.adj_map[b].size(), b);
    });

    // B) Insert seed into P
    for (int v : sorted_clique) add_vertex_internal(v);

    // C) v* from δ(P) bucket (min id)
    const int deltaP = min_degree_in_P_fast();
    int v_star = -1;
    if (deltaP >= 0 && deltaP < static_cast<int>(inside_P.size()) && !inside_P[deltaP].empty()) {
        v_star = *std::min_element(inside_P[deltaP].begin(), inside_P[deltaP].end());
    } else {
        v_star = *std::min_element(sorted_clique.begin(), sorted_clique.end()); // safe fallback
    }

    // No extra seed canonicalization here; match enumerate_with_seeds semantics

    // E) Recurse from true v*
    iter(v_star);

    // F) Clean up
    for (auto it = sorted_clique.rbegin(); it != sorted_clique.rend(); ++it) remove_from_inside_P(*it);
}


// Enumerate using in-memory R-clique seeds (no file I/O)
void PseudoCliqueEnumerator::enumerate_with_seeds(const std::vector<std::vector<int>>& r_cliques) {
    if (r_cliques.empty()) {
        std::cout << "No " << R << "-cliques found to extend." << std::endl;
        return;
    }

    for (const auto& r_clique : r_cliques) {
        // A) Order seed vertices deterministically and insert them into P
        auto sorted_clique = r_clique;
        std::sort(sorted_clique.begin(), sorted_clique.end(), [&](int a, int b) {
            return std::make_pair(graph.adj_map[a].size(), a) < std::make_pair(graph.adj_map[b].size(), b);
        });

        for (int vertex : sorted_clique) {
            add_vertex_internal(vertex);
        }

        // B) Compute v* = min-id among vertices with degree δ(P) inside P
        int v_star = -1;
        const int deltaP = min_degree_in_P_fast(); // δ(P)

        if (deltaP >= 0 && deltaP < static_cast<int>(inside_P.size()) && !inside_P[deltaP].empty()) {
            // Choose the smallest ID among the δ(P) bucket
            v_star = *std::min_element(inside_P[deltaP].begin(), inside_P[deltaP].end());
        } else {
            // Fallback (shouldn't happen if we just inserted a seed): smallest ID in the seed
            v_star = *std::min_element(sorted_clique.begin(), sorted_clique.end());
        }

        // C) Launch recursion from the true v*
        iter(v_star);

        // D) Clean up P for the next seed
        for (auto it = sorted_clique.rbegin(); it != sorted_clique.rend(); ++it) {
            remove_from_inside_P(*it);
        }
    }
}


void PseudoCliqueEnumerator::add_to_inside_P(int v) { //Add vertex v to the current pseudo-clique P and update all tracking structures in O(degree(v)) time.
    // std::cout << "adding: " << v  << "\n";

    add_vertex_internal(v); // This updates the state

    // Count all pseudo-cliques (non-maximal mode) if canonical for this seed
    // if (!maximal_only && total_nodes_in_P >= min_size) {
    //     if (is_current_seed_canonical() && static_cast<size_t>(total_nodes_in_P) < pseudo_cliques_count.size()) {
    //         pseudo_cliques_count[total_nodes_in_P] += 1;
    //     }
    // }
    // Non-maximal enumeration removed

    iter(v); //Find and process the children of this NEW pseudo-clique
}

void PseudoCliqueEnumerator::remove_from_inside_P(int v) {
        tracks[v][0] = false; //Mark v as not in P
        int degree = tracks[v][1]; //Get degree of v in P
        int pos = tracks[v][3]; //Get position of v in P

        if (!inside_P[degree].empty()) { //If v is in P, just a check
            int last_node = inside_P[degree].back(); //Get last vertex in degree bucket
            std::swap(inside_P[degree][pos], inside_P[degree].back()); //Swap v with last vertex in degree bucket
            inside_P[degree].pop_back(); //Remove last vertex from degree bucket
            tracks[last_node][3] = pos; //Update position of last vertex in degree bucket
        }

        // untrack c(P)
        int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0;
        if (cv >= 0 && cv < (int)core_hist.size()) {
            core_hist[cv] -= 1;
            if (cv == cur_min_core && core_hist[cv] == 0) {
                while (cur_min_core < (int)core_hist.size() && core_hist[cur_min_core] == 0)
                    cur_min_core++;
                if (cur_min_core >= (int)core_hist.size()) cur_min_core = INT_MAX;
            }
        }


        //update all of v's neighbors' info
        for (int adj_u : graph.adj_map[v]) {
            // Process neighbors of vertex v that are also in P
            if (tracks[adj_u][0]) {
                total_edges_in_P--;
                int degree = tracks[adj_u][1];
                int pos = tracks[adj_u][3];

                if (!inside_P[degree].empty()) {
                    int last_node = inside_P[degree].back();
                    std::swap(inside_P[degree][pos], inside_P[degree].back());
                    inside_P[degree].pop_back();
                    tracks[last_node][3] = pos;
                }

                inside_P[degree - 1].push_back(adj_u);
                tracks[adj_u][3] = inside_P[degree - 1].size() - 1;
                tracks[adj_u][1]--;
            }

            // Process neighbors of vertex v that are not in P
            int degree_np = tracks[adj_u][2]; // Get degree of neighbor wrt P
            int pos_np = tracks[adj_u][4]; // Get position in neighbors_and_P

            if (!neighbors_and_P[degree_np].empty()) {
                int last_node = neighbors_and_P[degree_np].back();
                std::swap(neighbors_and_P[degree_np][pos_np], neighbors_and_P[degree_np].back());
                neighbors_and_P[degree_np].pop_back();
                tracks[last_node][4] = pos_np;
            }

            neighbors_and_P[degree_np - 1].push_back(adj_u); // Add neighbor to lower degree bucket
            tracks[adj_u][4] = neighbors_and_P[degree_np - 1].size() - 1; // Update position of neighbor in neighbors_and_P
            tracks[adj_u][2]--; // Update degree of neighbor in neighbors_and_P
        }

        total_nodes_in_P--;
        set_theta_P();
}

void PseudoCliqueEnumerator::iter(int v) { //the "sophisticated process"
    // std::cout << "iter: " << v << "\n";
    iter_count++;

    // EDGE bound: only when R ≤ |P| < ℓ (fpce.c behavior)
    if (total_nodes_in_P >= R && total_nodes_in_P < min_size) {
        if (prune_by_edge_bound()) {
            numcalls_saved_by_edge_bound++;
            return; // prune this subtree
        }
    }

    // A) Use a LOCAL children container (do NOT use a member/shared one)
    std::vector<int> children;
    children.reserve(256); // or a heuristic you like

    // B) Integer threshold for "add-degree in P" using ceil(theta_P)
    const int min_add_deg = static_cast<int>(std::ceil(theta_P));


    int c = 0; // Counter for number of children found
    int v_star_degree = tracks[v][2]; // Get degree of v in P, v is v*
    
    // Child type 1: Lower Degree Vertices
    for (int deg = std::max(0, static_cast<int>(theta_P)); deg < v_star_degree; ++deg) {
        for (int u : neighbors_and_P[deg]) {
            if (tracks[u][0] || tracks[u][2] < min_add_deg) {
                continue;
            }
            c++;
            children.push_back(u);
        }
    }

    // Child type 2: Same Degree Vertices
    for (int u : neighbors_and_P[v_star_degree]) {
        if (tracks[u][0] || tracks[u][2] < min_add_deg) { // Skip if u is already in P or has degree less than theta_P
            continue;
        }
        if (u < v) { // Tie-break using vertex index
            children.push_back(u);
            c++;
        } 
        else if (graph.has_edge(v, u)) { //if u>v, practical implementation of the main condition from Lemma 4: u <_K l(u, K).
            bool valid = true;
            for (int x : inside_P[v_star_degree]) { //iterating over same degree vertices
                if (!graph.has_edge(u, x) && x < u) { //if u and x have edge and x is lexicographically smaller than u
                    valid = false;
                    break;
                }
            }
            if (valid) {
                children.push_back(u);
                c++;

            }
        }
    }

    // ---- Child type 3: vertices with degree == δ(P)+1 in P (correct gating) ----
    // v_star_degree is δ(P). We target vertices whose degree in P is δ(P)+1.
    {
        const int deltaP    = v_star_degree;
        const int targetDeg = deltaP + 1;

        // Bounds/availability check
        if (targetDeg >= 0 && targetDeg < static_cast<int>(inside_P.size())) {

            // Iterate candidates that currently have degree == δ(P)+1 inside P
            // (Assumes neighbors_and_P[deg] enumerates vertices of that degree in P’s neighborhood)
            for (int u : neighbors_and_P[targetDeg]) {

                // Skip if 'u' already in P (or otherwise marked) OR doesn't meet the add-degree threshold
                if (tracks[u][0] || tracks[u][2] < min_add_deg) continue;

                // Tie-breaks relative to v*:
                //  - u must be adjacent to v*,
                //  - and enforce canonical lex order u < v (equivalently: v > u)
                if (!graph.has_edge(v, u) || !(v > u)) continue;

                bool valid = true;

                // (1) u must be adjacent to all vertices in the δ(P) bucket (the min-degree vertices in P)
                //     This ensures we only consider "next-layer" candidates that fully connect to the
                //     current δ(P) core as required by the expansion rules.
                for (int x : inside_P[deltaP]) {
                    if (!graph.has_edge(u, x)) { valid = false; break; }
                }

                // (2) Lemma 4-style lex constraint within the (δ(P)+1) bucket:
                //     For any x in the same bucket with x < u, u must be adjacent to x.
                //     This preserves the original ordering rule and prevents duplicates.
                if (valid) {
                    for (int x : inside_P[targetDeg]) {
                        if (x < u && !graph.has_edge(u, x)) { valid = false; break; }
                    }
                }

                if (!valid) continue;

                // Passed all checks: add as a Type-3 child
                children.push_back(u);
                ++c;
            }
        }
    }


    // Report maximal pseudo-cliques when requested (pce.c style check)
    if (/*maximal_only &&*/ total_nodes_in_P >= min_size && is_current_maximal()) {
        auto current_p = collect_current_pseudo_clique();
        report_if_new(current_p);
    }

    // Iterate over children and add/remove them from inside_P
    while (c > 0) {

        add_to_inside_P(children.back());
        remove_from_inside_P(children.back());
        children.pop_back();
        c--;
    }
}

// In dense-pce.cpp

// void PseudoCliqueEnumerator::iter(int v_star) {
//     iter_count++;
    
//     // 1. Calculate necessary parameters for the current state P
//     int p_size = total_nodes_in_P;
//     double min_edges_for_child = theta * p_size * (p_size + 1) / 2.0 - total_edges_in_P;
//     int min_degree_in_P = (p_size > 0) ? inside_P.begin()->size() > 0 ? 0 : std::find_if(inside_P.begin(), inside_P.end(), [](const auto& bucket){ return !bucket.empty(); }) - inside_P.begin() : 0;

//     // 2. Identify all potential candidates
//     std::vector<int> candidates;
//     // A candidate `u` must have enough neighbors in P to satisfy the density requirement.
//     for (int deg = std::ceil(min_edges_for_child); deg <= p_size; ++deg) {
//         for (int u : neighbors_and_P[deg]) {
//             if (!tracks[u][0]) { // Ensure u is not already in P
//                 candidates.push_back(u);
//             }
//         }
//     }

//     // 3. Filter candidates based on the MD Ordering rule (Lemma 5)
//     int c = 0;
//     for (int u : candidates) {
//         // This check implements the core of Property 1 / Lemma 5.
//         // A child is valid only if its connection degree into P is less than the minimum degree
//         // of P, OR if it's equal and it comes lexicographically before v_star.
//         // This is the essential rule that builds the unique reverse search tree.
//         if (tracks[u][2] < min_degree_in_P + 1 || (tracks[u][2] == min_degree_in_P + 1 && u < v_star)) {
//             children.push(u);
//             c++;
//         }
//     }

//     // 4. Report if maximal, otherwise recurse (same as your current code)
//     if (c == 0) { 
//         std::vector<int> current_p;
//         for (const auto& bucket : inside_P) {
//             if (!bucket.empty()) {
//                 current_p.insert(current_p.end(), bucket.begin(), bucket.end());
//             }
//         }
//         report_if_new(current_p);
//     }

//     while (c > 0) {
//         add_to_inside_P(children.top());
//         remove_from_inside_P(children.top());
//         children.pop();
//         c--;
//     }
// }


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]<< " <filename> [--theta <theta>] [--minimum <min>] [--maximum <max>]" << std::endl;
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
    int minimum   = 1;
    int maximum   = std::numeric_limits<int>::max();

    for (; arg_idx < argc; ++arg_idx) {
        std::string arg = argv[arg_idx];
        if (arg == "--theta" && arg_idx + 1 < argc) {
            theta = std::stod(argv[++arg_idx]);
        } else if (arg == "--minimum" && arg_idx + 1 < argc) {
            minimum = std::stoi(argv[++arg_idx]);
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

    // 3. Compute R and show diagnostics.
    int R = std::ceil(1.0 / (1.0 - theta * (minimum - 1) / static_cast<double>(minimum)));
    std::cout << "Computed R: " << R << std::endl;

    // 4. Derive directory path in a cross‑platform way.
    std::filesystem::path p(filename);
    std::string dir_path = p.parent_path().empty() ? "." : p.parent_path().string();
    std::cout << "Directory path: " << dir_path << std::endl;


    // 6. Read the graph and enumerate pseudo‑cliques.
    Graph graph;
    // Try to load BBkC binaries from the same directory as the provided .grh path
    bool graph_loaded = false;
    {
        std::filesystem::path bdeg = std::filesystem::path(dir_path) / "b_degree.bin";
        std::filesystem::path badj = std::filesystem::path(dir_path) / "b_adj.bin";
        if (std::filesystem::exists(bdeg) && std::filesystem::exists(badj)) {
            std::cout << "Detected BBkC binaries in: " << dir_path << ", loading graph from binaries..." << std::endl;
            graph_loaded = graph.read_graph_from_bbkcbinaries(dir_path, filename);
        } else {
            graph_loaded = graph.read_graph_from_file(filename);
        }
    }
    if (!graph_loaded) {
        std::cerr << "Aborting: graph could not be loaded." << std::endl;
        return 1;
    }
    graph.compute_core_numbers();

    // Order Bound here with the *real* degeneracy
    int mu = graph.compute_order_bound(theta, graph.degeneracy);
    if (minimum > mu) {
        std::cout << "Pruning by Order Bound: No (l,θ)-pseudo-cliques exist as l ("
                << minimum << ") > μ (" << mu << ")\n";
        return 0; 
    }
    else{
        std::cout<<"NOT Pruning by Order Bound: (l,θ)-pseudo-cliques exist as l (" << minimum << ") <= μ (" << mu << ")" << std::endl;
    }

    if (maximum == std::numeric_limits<int>::max()) {
        maximum = static_cast<int>(graph.adj_map.size());
    }
    // Clamp to order bound μ to cap |P| and reduce memory (bucket arrays sized to max_size+1)
    maximum = std::min(maximum, mu);
 
    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum, R);

    // 5. Stream R-clique seeds via EBBkC and extend on the fly (O(1) memory)
    if (!graph.csr_node_off.empty() && !graph.csr_edge_dst.empty()) {
        EBBkC_t::list_k_clique_mem_stream_from_csr(
            graph.csr_n,
            graph.csr_m,
            graph.csr_node_off.data(),
            graph.csr_edge_dst.data(),
            R, 2,
            [&](const std::vector<int>& r_clique){ PC.enumerate_with_turan(r_clique); }
        );
    } else {
        EBBkC_t::list_k_clique_mem_stream(
            dir_path.c_str(), R, 2,
            [&](const std::vector<int>& r_clique){ PC.enumerate_with_turan(r_clique); }
        );
    }

    // 7. Display results.
    std::vector<int> pseudo_clique_counts = PC.get_pseudo_cliques_count();
    std::cout << "Maximal pseudo-clique counts:" << std::endl;
    bool found_any = false;
    unsigned long long total_pseudo_cliques = 0ULL;
    for (size_t sz = 0; sz < pseudo_clique_counts.size(); ++sz) {
        if (pseudo_clique_counts[sz] > 0) {
            std::cout << "Size " << sz << ": " << pseudo_clique_counts[sz] << "\n";
            found_any = true;
        }
        total_pseudo_cliques += static_cast<unsigned long long>(pseudo_clique_counts[sz]);
    }
     
    if (!found_any) {
        std::cout << "(No pseudo-cliques found meeting the criteria)" << std::endl;
    }
    std::cout << "Total pseudo-cliques: " << total_pseudo_cliques << std::endl;
 
    std::cout << "\nTotal Iterations: " << PC.get_iter_count() << std::endl;
    std::cout << "#calls of saved by EDGE bound = "<< PC.get_edge_prunes() << "\n";

    return 0;
}
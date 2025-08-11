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

// Define PATH_MAX if not already defined
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

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
    // std::unordered_map<int, std::unordered_map<int, bool> > adj_map;
    std::vector<std::unordered_map<int, bool>> adj_map; //unordered_set

    int total_nodes = 0 ;

    void add_edge(int u, int v) {
        if (u >= adj_map.size() || v >= adj_map.size()) {
           adj_map.resize(std::max(u, v) + 1);
            }
        adj_map[u][v] = true; 
        adj_map[v][u] = true;
    }

    void add_empty_edge(int u) { //placeholder for empty edge
        adj_map[u][-1] = false;        
    }

    void print_graph() const { //prints out each vertex and its list of neighbors
        for (size_t i = 0; i < adj_map.size(); ++i) {
            const auto& neighbors_map = adj_map[i]; //unordered_map
            std::vector<int> neighbors;

            for (const auto& neighbor : neighbors_map) {
                neighbors.push_back(neighbor.first); //{2, true}
            }

            std::sort(neighbors.begin(), neighbors.end());

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
    void read_graph_from_file(const std::string& filename) {
        graph_file_path = filename;   // <- add this line
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
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
            // std::string value;
            // while (std::getline(ss, value, ' ')) {
            //     int neighbor = std::stoi(value); //string into integer
            int neighbor;
            while (ss >> neighbor) {
                add_edge(current_vertex, neighbor);
            }
            current_vertex++;
        }
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
            for (const auto& neighbor : adj_map[v]) {
                unsigned int u = neighbor.first;
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
        for (const auto& adj_v_ite : graph.adj_map[v]) {
            
            int adj_v = adj_v_ite.first;
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

        //Solves: Vector Reallocation
        int reserve_size = total_nodes / 4 ;
        for (int i = 0; i < max_size; ++i) {
            inside_P[i].reserve(reserve_size); //add up to reserve_size elements later
            neighbors_and_P[i].reserve(reserve_size); 

        } 
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


    void enumerate_with_turan();  // Enumerate pseudo-cliques using Turan's theorem
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

void PseudoCliqueEnumerator::enumerate_with_turan() {
    // 1. Build path to the clique file
    std::string graph_path = graph.graph_file_path;
    size_t last_slash = graph_path.find_last_of("/\\");
    std::string dir = (last_slash == std::string::npos) ? "." : graph_path.substr(0, last_slash);
    std::cout << "dir: " << dir << std::endl;
    std::string clique_file = dir + "/cliques_K" + std::to_string(R);
    
    // 2. Load the r-clique seeds
    auto r_cliques = read_clique_file(clique_file);
    if (r_cliques.empty()) {
        std::cout << "No " << R << "-cliques found to extend." << std::endl;
        return;
    }

    // 3. Process each r-clique seed
    for (const auto& r_clique : r_cliques) {
        auto sorted_clique = r_clique;
        std::sort(sorted_clique.begin(), sorted_clique.end(), [&](int a, int b) {
            return std::make_pair(graph.adj_map[a].size(), a) < std::make_pair(graph.adj_map[b].size(), b);
        });

        // ======================= START: NEW DEBUG LINE =======================
        //std::cout << "[DEBUG] Executing search branch starting from seed r-clique: ";
        // for(int v : sorted_clique) { std::cout << v << " "; }
        // std::cout << std::endl;
        // ======================== END: NEW DEBUG LINE ========================
        
        // --- A. Seed the enumerator's state ---
        for (int vertex : sorted_clique) {
            add_vertex_internal(vertex); // Use the new helper that doesn't recurse
        }
        // Non-maximal enumeration removed
        // // Establish canonical seed ids: ascending order of the seed vertices by id
        // current_seed_sorted_ids = r_clique;
        // std::sort(current_seed_sorted_ids.begin(), current_seed_sorted_ids.end());

        // // Count the seeded R-clique itself in non-maximal mode if canonical
        // if (!maximal_only && total_nodes_in_P >= min_size) {
        //     if (is_current_seed_canonical() && static_cast<size_t>(total_nodes_in_P) < pseudo_cliques_count.size()) {
        //         pseudo_cliques_count[total_nodes_in_P] += 1;
        //     }
        // }

        // --- B. Kick-off the recursive search for extensions ---
        // int last_vertex = *std::max_element(r_clique.begin(), r_clique.end());
        int last_vertex = sorted_clique.back(); // Use last added vertex, not max_element
        iter(last_vertex);

        // --- C. Backtrack to reset the state for the next seed ---
        // Must remove vertices in the reverse order of addition.
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
        for (const auto& adj_v : graph.adj_map[v]) {
            // Process neighbors of vertex v that are also in P
            if (tracks[adj_v.first][0]) {
                total_edges_in_P--;
                int degree = tracks[adj_v.first][1];
                int pos = tracks[adj_v.first][3];

                if (!inside_P[degree].empty()) {
                    int last_node = inside_P[degree].back();
                    std::swap(inside_P[degree][pos], inside_P[degree].back());
                    inside_P[degree].pop_back();
                    tracks[last_node][3] = pos;
                }

                inside_P[degree - 1].push_back(adj_v.first);
                tracks[adj_v.first][3] = inside_P[degree - 1].size() - 1;
                tracks[adj_v.first][1]--;
            }

            // Process neighbors of vertex v that are not in P
            int degree_np = tracks[adj_v.first][2]; // Get degree of neighbor adj_v wrt P
            int pos_np = tracks[adj_v.first][4]; // Get position of neighbor adj_v in neighbors_and_P

            if (!neighbors_and_P[degree_np].empty()) {
                int last_node = neighbors_and_P[degree_np].back();
                std::swap(neighbors_and_P[degree_np][pos_np], neighbors_and_P[degree_np].back());
                neighbors_and_P[degree_np].pop_back();
                tracks[last_node][4] = pos_np;
            }

            neighbors_and_P[degree_np - 1].push_back(adj_v.first); // Add neighbor to lower degree bucket
            tracks[adj_v.first][4] = neighbors_and_P[degree_np - 1].size() - 1; // Update position of neighbor in neighbors_and_P
            tracks[adj_v.first][2]--; // Update degree of neighbor in neighbors_and_P
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


    int c = 0; // Counter for number of children found
    int v_star_degree = tracks[v][2]; // Get degree of v in P, v is v*
    
    // Child type 1: Lower Degree Vertices
    for (int deg = std::max(0, static_cast<int>(theta_P)); deg < v_star_degree; ++deg) {
        for (int u : neighbors_and_P[deg]) {
            if (tracks[u][0] || tracks[u][2] < theta_P) {
                continue;
            }
            c++;
            children.push(u);
        }
    }

    // Child type 2: Same Degree Vertices
    for (int u : neighbors_and_P[v_star_degree]) {
        if (tracks[u][0] || tracks[u][2] < theta_P) { // Skip if u is already in P or has degree less than theta_P
            continue;
        }
        if (u < v) { // Tie-break using vertex index
            children.push(u);
            c++;
        } 
        else if (graph.adj_map[v].count(u)) { //if u>v, practical implementation of the main condition from Lemma 4: u <_K l(u, K).
            bool valid = true;
            for (int x : inside_P[v_star_degree]) { //iterating over same degree vertices
                if (!graph.adj_map[u].count(x) && x < u) { //if u and x have edge and x is lexicographically smaller than u
                    valid = false;
                    break;
                }
            }
            if (valid) {
                children.push(u);
                c++;

            }
        }
    }

    // Child type 3: Higher Degree Vertices (restore original ordering rule)
    if (v_star_degree + 1 >= inside_P[v_star_degree].size()) {
        for (int u : neighbors_and_P[v_star_degree + 1]) {
            if (tracks[u][0] || tracks[u][2] < theta_P) {
                continue;
            }
            if (graph.adj_map[v].count(u) && v > u) {
                bool valid = true;
                for (int x : inside_P[v_star_degree]) {
                    if (!graph.adj_map[u].count(x)) {
                        valid = false;
                        break;
                    }
                }
                for (int x : inside_P[v_star_degree + 1]) { // lemma 4
                    if (!graph.adj_map[u].count(x) && x < u) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    children.push(u);
                    c++;
                }
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

        add_to_inside_P(children.top());
        remove_from_inside_P(children.top());
        children.pop();
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

    /* ------------------------------------------------------------------
     * 1. Re‑assemble the (possibly space‑containing) filename.
     *    Every argument up to—but NOT including—the first that starts
     *    with "--" is considered part of the path.
     * ----------------------------------------------------------------*/
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

    /* ------------------------------------------------------------------
     * 2. Parse optional arguments that follow the filename.
     * ----------------------------------------------------------------*/
    double theta = 1.0;
    int minimum   = 1;
    int maximum   = std::numeric_limits<int>::max();
    //bool maximal_only = true;

    for (; arg_idx < argc; ++arg_idx) {
        std::string arg = argv[arg_idx];
        if (arg == "--theta" && arg_idx + 1 < argc) {
            theta = std::stod(argv[++arg_idx]);
        } else if (arg == "--minimum" && arg_idx + 1 < argc) {
            minimum = std::stoi(argv[++arg_idx]);
        } 
        // else if (arg == "--maximum" && arg_idx + 1 < argc) {
        //     maximum = std::stoi(argv[++arg_idx]);
        // } 
        //else if (arg == "--maximal") {
            // no-op, always maximal mode}
        else {
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

    /* ------------------------------------------------------------------
     * 3. Compute R and show diagnostics.
     * ----------------------------------------------------------------*/
    int R = std::ceil(1.0 / (1.0 - theta * (minimum - 1) / static_cast<double>(minimum)));
    std::cout << "Computed R: " << R << std::endl;

    /* ------------------------------------------------------------------
     * 4. Derive directory path in a cross‑platform way.
     * ----------------------------------------------------------------*/
    std::filesystem::path p(filename);
    std::string dir_path = p.parent_path().empty() ? "." : p.parent_path().string();
    std::cout << "Directory path: " << dir_path << std::endl;

    /* ------------------------------------------------------------------
     * 5. Launch BBkC.
     * ----------------------------------------------------------------*/
    std::string cmd = "./BBkC e \"" + dir_path + "\" " + std::to_string(R) + " " + std::to_string(2);
    std::cout << "Executing: " << cmd << std::endl;
    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Error: BBkC failed to execute (return code: " << ret << ")" << std::endl;
        return 1;
    }

    /* ------------------------------------------------------------------
     * 6. Read the graph and enumerate pseudo‑cliques.
     *    (Assumes Graph and PseudoCliqueEnumerator provide the same API
     *     as in your original code.)
     * ----------------------------------------------------------------*/
    Graph graph;
    graph.read_graph_from_file(filename);
    graph.compute_core_numbers();

    if (maximum == std::numeric_limits<int>::max()) {
        maximum = static_cast<int>(graph.adj_map.size());
    }

    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum, R);
    // for (size_t node = 0; node < graph.adj_map.size(); ++node) {
    //     PC.add_to_inside_P(node);
    //     PC.remove_from_inside_P(node);
    // }
    PC.enumerate_with_turan();

    /* ------------------------------------------------------------------
     * 7. Display results.
     * ----------------------------------------------------------------*/
     std::vector<int> pseudo_clique_counts = PC.get_pseudo_cliques_count();
     std::cout << (/*maximal_only ?*/ "Maximal pseudo-clique counts:") << std::endl;
     bool found_any = false;
     for (size_t sz = 0; sz < pseudo_clique_counts.size(); ++sz) {
         if (pseudo_clique_counts[sz] > 0) {
             std::cout << "Size " << sz << ": " << pseudo_clique_counts[sz] << "\n";
             found_any = true;
         }
     }
     
     if (!found_any) {
         std::cout << "(No pseudo-cliques found meeting the criteria)" << std::endl;
     }
 
     std::cout << "\nTotal Iterations: " << PC.get_iter_count() << std::endl;
     std::cout << "#calls of PCE_iter saved by EDGE bound = "<< PC.get_edge_prunes() << "\n";

    return 0;
}
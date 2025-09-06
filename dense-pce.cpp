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
#include <set>
#include <cmath>



// counting sort with tracklist
class Graph {
public:
    Graph() = default;
    // std::unordered_map<int, std::unordered_map<int, bool> > adj_map;
    std::vector<std::unordered_map<int, bool>> adj_map;

    int total_nodes = 0 ;

    void add_edge(int u, int v) {
        if (u >= adj_map.size() || v >= adj_map.size()) {
           adj_map.resize(std::max(u, v) + 1);
            }
        adj_map[u][v] = true;
        adj_map[v][u] = true;
    }

    void add_empty_edge(int u) {
        adj_map[u][-1] = false;
        
    }

    void print_graph() const {
        for (size_t i = 0; i < adj_map.size(); ++i) {
            const auto& neighbors_map = adj_map[i];
            std::vector<int> neighbors;

            for (const auto& neighbor : neighbors_map) {
                neighbors.push_back(neighbor.first);
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

    void read_graph_from_file(const std::string& filename) {
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
            if (line.empty()) {
                current_vertex++;
                continue;
            }

            std::stringstream ss(line);
            int neighbor;
            while (ss >> neighbor) {
                add_edge(current_vertex, neighbor);
            }
            current_vertex++;
        }
    }
};

class PseudoCliqueEnumerator {
public:
    PseudoCliqueEnumerator(Graph& graph, float theta = 1.0, int min_size = 1, int max_size = -1, bool maximal_only = false)
        : graph(graph), theta(theta), min_size(min_size), maximal_only(maximal_only) {
        
        total_nodes = graph.adj_map.size(); 

        std::cout << "total nodes: " << total_nodes << " "   << "\n";

        this->max_size = max_size != -1 ? max_size : total_nodes;

        pseudo_cliques.resize(max_size + 1);
        pseudo_cliques_count.resize(max_size + 1, 0);

        inside_P.resize(max_size + 1);
        neighbors_and_P.resize(max_size + 1);

        int reserve_size = total_nodes / 4 ;
        for (int i = 0; i < max_size; ++i) {
            inside_P[i].reserve(reserve_size); 
            neighbors_and_P[i].reserve(reserve_size); 

        }
        std::vector<int> keys(graph.total_nodes);
        std::iota(keys.begin(), keys.end(), 0);


        neighbors_and_P[0] = keys;

        // neighbors_and_P[0] = get_keys(graph.adj_map);
        
        
        

        int idx = 0;
        for (const auto& v : neighbors_and_P[0]) {
            tracks.push_back({0, -1, 0, -1, idx});
            ++idx;
        }

        theta_P = 0;
        total_nodes_in_P = 0;
        total_edges_in_P = 0;
    }

    void set_theta_P() {
        theta_P = (theta * total_nodes_in_P * (total_nodes_in_P + 1) / 2.0) - total_edges_in_P;
    }

    void print_all() const {
        std::cout << "P: ";
        for (const auto& p : inside_P) {
            for (int v : p) std::cout << v << " ";
            std::cout << "| ";
        }
        std::cout << "\nNP: ";
        for (const auto& np : neighbors_and_P) {
            for (int v : np) std::cout << v << " ";
            std::cout << "| ";
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
        std::cout << "\ntheta_P: " << theta_P << " (" << total_nodes_in_P << ", " << total_edges_in_P << ")\n";
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

private:
    // Collect current pseudo-clique P as a sorted vector
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

    // Maximality check following pce.c logic:
    // No vertex outside P has degree into P >= ceil(theta_P)
    bool is_current_maximal() const {
        int min_add_deg = static_cast<int>(std::ceil(theta_P));
        if (min_add_deg < 0) min_add_deg = 0;
        // Compare histogram counts: all vertices vs vertices in P
        for (int d = min_add_deg; d <= total_nodes_in_P; ++d) {
            size_t occ_num = (d >= 0 && d < static_cast<int>(neighbors_and_P.size())) ? neighbors_and_P[d].size() : 0;
            size_t in_p_num = (d >= 0 && d < static_cast<int>(inside_P.size())) ? inside_P[d].size() : 0;
            if (occ_num > in_p_num) return false; // there exists an outside vertex extendable
        }
        return true;
    }


    Graph& graph;
    float theta;
    int min_size;
    int max_size;
    int total_nodes;
    float theta_P;
    int total_nodes_in_P;
    int total_edges_in_P;
    int iter_count = 0;
    bool maximal_only = false;

    std::stack<int> children;

    std::vector<std::vector<int> > inside_P;
    std::vector<std::vector<int> > neighbors_and_P;
    std::vector<std::vector<int> > pseudo_cliques;
    std::vector<int> pseudo_cliques_count;

    // std::unordered_map<int, TrackInfo> tracks;
    std::vector<std::vector<int> > tracks;

    // To avoid duplicate reporting of maximal pseudo-cliques
    std::set<std::vector<int>> reported_maximal_pseudo_cliques;

    


    // std::vector<int> get_keys(const std::vector<int, std::unordered_map<int, bool> >& map) {
    //     std::vector<int> keys;
    //     for (const auto& pair : map) {
    //         keys.push_back(pair.first);
    //     }
    //     std::sort(keys.begin(), keys.end());
    //     return keys;
    // }
}; 

void PseudoCliqueEnumerator::add_to_inside_P(int v) {
    // std::cout << "adding: " << v  << "\n";

    if (total_nodes_in_P > max_size) {
        return;
    }
    
    tracks[v][0] = true;
    tracks[v][1] = tracks[v][2];
    inside_P[tracks[v][1]].push_back(v);
    tracks[v][3] = inside_P[tracks[v][1]].size() - 1;
    

    for (const auto& adj_v_ite : graph.adj_map[v]) {
        
        int adj_v = adj_v_ite.first;
        // std::cout << "    ite: " << adj_v  << "\n";


        if (tracks[adj_v][0]) {
            total_edges_in_P += 1;
            int degree = tracks[adj_v][1];
            tracks[inside_P[degree].back()][3] = tracks[adj_v][3];

            // Swap and pop
            std::swap(inside_P[degree][tracks[adj_v][3]], inside_P[degree].back());
            inside_P[degree].pop_back();

            inside_P[degree + 1].push_back(adj_v);
            tracks[adj_v][3] = inside_P[degree + 1].size() - 1;
            tracks[adj_v][1] = degree + 1;
        }

        // Update degree_in_NP
        int degree = tracks[adj_v][2];
        int pos = tracks[adj_v][4];
        int last_node = neighbors_and_P[degree].back();


        
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

    total_nodes_in_P += 1;
    set_theta_P();

    // Check if current set qualifies as a pseudo-clique
    // In maximal-only mode, we defer counting until we know it's maximal (leaf state)
    if (!maximal_only && total_nodes_in_P >= min_size) {
        pseudo_cliques_count[total_nodes_in_P] += 1;
    }

    iter(v);
}


void PseudoCliqueEnumerator::remove_from_inside_P(int v) {
        tracks[v][0] = false;
        int degree = tracks[v][1];
        int pos = tracks[v][3];

        if (!inside_P[degree].empty()) {
            int last_node = inside_P[degree].back();
            std::swap(inside_P[degree][pos], inside_P[degree].back());
            inside_P[degree].pop_back();
            tracks[last_node][3] = pos;
        }

        for (const auto& adj_v : graph.adj_map[v]) {
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

            int degree_np = tracks[adj_v.first][2];
            int pos_np = tracks[adj_v.first][4];

            if (!neighbors_and_P[degree_np].empty()) {
                int last_node = neighbors_and_P[degree_np].back();
                std::swap(neighbors_and_P[degree_np][pos_np], neighbors_and_P[degree_np].back());
                neighbors_and_P[degree_np].pop_back();
                tracks[last_node][4] = pos_np;
            }

            neighbors_and_P[degree_np - 1].push_back(adj_v.first);
            tracks[adj_v.first][4] = neighbors_and_P[degree_np - 1].size() - 1;
            tracks[adj_v.first][2]--;
        }

        total_nodes_in_P--;
        set_theta_P();
}

void PseudoCliqueEnumerator::iter(int v) {
    // std::cout << "iter: " << v << "\n";
    iter_count++;
    int c = 0;
    int v_star_degree = tracks[v][2];
    
    // Child type 1
    for (int deg = std::max(0, static_cast<int>(theta_P)); deg < v_star_degree; ++deg) {
        for (int u : neighbors_and_P[deg]) {
            if (tracks[u][0] || tracks[u][2] < theta_P) {
                continue;
            }
            c++;
            children.push(u);
        }
    }

    for (int u : neighbors_and_P[v_star_degree]) {
        if (tracks[u][0] || tracks[u][2] < theta_P) {
            continue;
        }
        if (u < v) {
            children.push(u);
            c++;

        } 
        else if (graph.adj_map[v].count(u)) {
            bool valid = true;
            for (int x : inside_P[v_star_degree]) {
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
                for (int x : inside_P[v_star_degree + 1]) {
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

    // Maximality check following pce.c: report current P if maximal
    if (maximal_only && total_nodes_in_P >= min_size && is_current_maximal()) {
        std::vector<int> current_p = collect_current_pseudo_clique();
        if (!current_p.empty()) {
            if (reported_maximal_pseudo_cliques.insert(current_p).second) {
                if (static_cast<size_t>(total_nodes_in_P) < pseudo_cliques_count.size()) {
                    pseudo_cliques_count[total_nodes_in_P] += 1;
                }
            }
        }
    }

    // Iterate over children and add/remove them from inside_P
    while (c > 0) {

        add_to_inside_P(children.top());
        remove_from_inside_P(children.top());
        children.pop();
        c--;
    }
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename> [--theta <theta>] [--minimum <min>] [--maximum <max>]" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    double theta = 1.0;
    int minimum = 1;
    int maximum = std::numeric_limits<int>::max();
    bool maximal_only = false;

    // Parse optional arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--theta" && i + 1 < argc) {
            theta = std::stod(argv[++i]);
        } else if (arg == "--minimum" && i + 1 < argc) {
            minimum = std::stoi(argv[++i]);
        } else if (arg == "--maximum" && i + 1 < argc) {
            maximum = std::stoi(argv[++i]);
        } else if (arg == "--maximal") {
            maximal_only = true;
        }
    }

    Graph graph;

    graph.read_graph_from_file(argv[1]);

    // graph.print_graph();

    std::cout << theta << "\n";
    
    // If maximum is not specified, use the total number of vertices
    if (maximum == std::numeric_limits<int>::max()) {
        maximum = graph.adj_map.size();
    }

    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum, maximal_only);
    
    // std::vector<int> nodes;
    // for (const auto& pair : graph.adj_map) {
    //     nodes.push_back(pair.first);
    // }
    // std::sort(nodes.begin(), nodes.end());

    for (size_t node = 0; node < graph.adj_map.size(); ++node) {
        // std::cout << "Starting node: " << node << std::endl;
        PC.add_to_inside_P(node);
        // std::cout << "removing node: " << node << std::endl;
        PC.remove_from_inside_P(node);
    }

    std::vector<int> pseudo_clique_counts = PC.get_pseudo_cliques_count();

    // Print the clique sizes
    std::cout << (maximal_only ? "Maximal pseudo-clique counts:" : "Pseudo-clique counts:") << std::endl;
    for (size_t sz = 0; sz < pseudo_clique_counts.size(); ++sz) {
        if (pseudo_clique_counts[sz] > 0) {
            std::cout << "Size " << sz << ": " << pseudo_clique_counts[sz] << "\n";
        }
    }
    std::cout << std::endl << "Total Iterations: " << PC.get_iter_count() << "\n";

    return 0;
}

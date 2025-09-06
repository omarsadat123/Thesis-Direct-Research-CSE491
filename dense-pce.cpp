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
            if (line.empty()) { //just an isolated vertex
                // add_empty_edge(current_vertex);
                current_vertex++;
                continue;
            }

            std::stringstream ss(line); //helps split the line of text into individual numbers
            std::string value;
            while (std::getline(ss, value, ' ')) {
                int neighbor = std::stoi(value); //string into integer
                add_edge(current_vertex, neighbor);
            }
            current_vertex++;
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

public:
    PseudoCliqueEnumerator(Graph& graph, float theta = 1.0, int min_size = 1, int max_size = -1)
        : graph(graph), theta(theta), min_size(min_size) {
        
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

    // std::vector<int> get_keys(const std::vector<int, std::unordered_map<int, bool> >& map) {
    //     std::vector<int> keys;
    //     for (const auto& pair : map) {
    //         keys.push_back(pair.first);
    //     }
    //     std::sort(keys.begin(), keys.end());
    //     return keys;
    // }
}; 

void PseudoCliqueEnumerator::add_to_inside_P(int v) { //Add vertex v to the current pseudo-clique P and update all tracking structures in O(degree(v)) time.
    // std::cout << "adding: " << v  << "\n";

    if (total_nodes_in_P > max_size) {
        return;
    }
    
    //Updating v's info
    tracks[v][0] = true;                             // Mark v as being part of P
    tracks[v][1] = tracks[v][2];                     // Its new "degree_in_P" is its old degree_with_P"
    inside_P[tracks[v][1]].push_back(v);             // Add v to the correct degree bucket in the inside_P list
    tracks[v][3] = inside_P[tracks[v][1]].size() - 1; // Record its position in that bucket    

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
    if (total_nodes_in_P >= min_size) {
        // std::vector<int> pseudo_clique;
        // for (const auto& nodes : inside_P) {
        //     pseudo_clique.insert(pseudo_clique.end(), nodes.begin(), nodes.end());
        // }
        // pseudo_cliques.push_back(pseudo_clique);
        pseudo_cliques_count[total_nodes_in_P] += 1;
    }

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

    // Child type 3: Higher Degree Vertices
    if (v_star_degree + 1 >= inside_P[v_star_degree].size()) {
        for (int u : neighbors_and_P[v_star_degree + 1]) {
            if (tracks[u][0] || tracks[u][2] < theta_P) {//if u already in P or has degree less than theta_P
                continue;
            }
            if (graph.adj_map[v].count(u) && v > u) { //For Type 3, v is connected to u and we need v > u (lexicographic ordering)
                bool valid = true;
                for (int x : inside_P[v_star_degree]) {
                    if (!graph.adj_map[u].count(x)) {
                        valid = false;
                        break;
                    }
                }
                for (int x : inside_P[v_star_degree + 1]) { //lemma 4
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

    // Parse optional arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--theta" && i + 1 < argc) {
            theta = std::stod(argv[++i]);
        } else if (arg == "--minimum" && i + 1 < argc) {
            minimum = std::stoi(argv[++i]);
        } else if (arg == "--maximum" && i + 1 < argc) {
            maximum = std::stoi(argv[++i]);
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

    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum);
    
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
    std::cout << "Pseudo-clique counts:" << std::endl;
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
    std::cout << std::endl << "Total Iterations: " << PC.get_iter_count() << "\n";

    return 0;
}
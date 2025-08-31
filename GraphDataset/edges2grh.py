import sys
import os
import collections

def convert_edges_to_grh(edges_file_path, grh_file_path):
    """
    Converts a graph from an edge list file to a .grh adjacency list format.

    Args:
        edges_file_path (str): The path to the input .edges file.
        grh_file_path (str): The path for the output .grh file.
    """
    adj_list = collections.defaultdict(list)
    max_vertex_id = -1

    print(f"Reading edges from '{edges_file_path}'...")
    try:
        with open(edges_file_path, 'r') as f:
            for line in f:
                # Skip empty or comment lines
                if not line.strip() or line.startswith('#'):
                    continue
                
                parts = line.strip().split()
                try:
                    # Assumes undirected graph, connect both ways
                    u, v = int(parts[0]), int(parts[1])
                    adj_list[u].append(v)
                    adj_list[v].append(u)
                    
                    # Track the highest vertex ID
                    current_max = max(u, v)
                    if current_max > max_vertex_id:
                        max_vertex_id = current_max
                except (ValueError, IndexError):
                    print(f"Warning: Skipping malformed line: '{line.strip()}'")

    except FileNotFoundError:
        print(f"Error: The file '{edges_file_path}' was not found.")
        return

    print(f"Graph processed. Max vertex ID found: {max_vertex_id}.")
    print(f"Writing .grh file to '{grh_file_path}'...")

    with open(grh_file_path, 'w') as f:
        # Ensure file has lines for all vertices from 0 to max_vertex_id
        for i in range(max_vertex_id + 1):
            # Get neighbors if vertex exists, otherwise empty list
            neighbors = adj_list.get(i, [])
            
            # Sort neighbors for consistent output and convert to string
            neighbors_str = ' '.join(map(str, sorted(neighbors)))
            
            f.write(neighbors_str + '\n')

    print(f"Conversion complete! ✨ File saved as '{grh_file_path}'.")


# --- Main execution block ---
if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python edges2grh.py <path_to_your_edges_file>")
        sys.exit(1) # Exit the script indicating an error

    # Get the input file path from the command line
    input_file = sys.argv[1]

    # Automatically determine the output filename
    # This takes the input filename, removes its extension, and adds '.grh'
    base_name = os.path.splitext(input_file)[0]
    output_file = base_name + '.grh'
    
    # Run the conversion function
    convert_edges_to_grh(input_file, output_file)
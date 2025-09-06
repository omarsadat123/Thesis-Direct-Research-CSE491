#!/bin/bash


# Output file
timestamp=$(date +"%Y%m%d_%H%M%S")
output_file="output__log_${timestamp}.txt"

# Clear the output file if it already exists
> "$output_file"


for graph_file in $(ls  density_graphs/*.grh); do
    # Process each .grh file in the density_graphs directory
    echo "Processing file: $graph_file"
    
    # Check if the file exists
    if [[ -f "$graph_file" ]]; then
        echo "Running for m=$m with file $graph_file" | tee -a "$output_file"
        
        # Run the command and measure the time
        { time ./dense-pce "$graph_file" --maximum 300 --theta 0.95; } 2>&1 | tee -a "$output_file"

        # Run the command and measure the time
        { time ./pce C -l 1 -u 100 "$graph_file" 0.95; } 2>&1 | tee -a "$output_file"
    else
        echo "File not found: $graph_file" | tee -a "$output_file"
    fi
done
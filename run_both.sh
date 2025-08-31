#!/bin/bash

# Set the range of values for m
start_m=10     # Start value for m
end_m=100      # End value for m
step_m=10      # Step size for m
n=1000

# Output file
output_file="output_both_log_${n}.txt"

# Clear the output file if it already exists
> "$output_file"

# Loop over the range of values
for m in $(seq $start_m $step_m $end_m); do
    # Update the filename with the current value of m
    graph_file="synth-graphs-${n}/scale_free_graph_m_${m}.grh"
    
    # Check if the file exists
    if [[ -f "$graph_file" ]]; then
        echo "Running for m=$m with file $graph_file" | tee -a "$output_file"
        
        # Run the command and measure the time
        { time ./cs3 "$graph_file" --maximum 300 --theta 0.95; } 2>&1 | tee -a "$output_file"

        # Run the command and measure the time
        { time ./pce C -l 1 -u 100 "$graph_file" 0.95; } 2>&1 | tee -a "$output_file"
    else
        echo "File not found: $graph_file" | tee -a "$output_file"
    fi
done
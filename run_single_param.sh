#!/bin/bash

# Output file
timestamp=$(date +"%Y%m%d_%H%M%S")
output_file="output__log_${timestamp}.txt"

# Clear the output file if it already exists
> "$output_file"

# Ensure script runs from its own directory
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

# Optional parallelism: use -j N to run up to N graphs concurrently (default 1 = sequential)
PARALLEL_JOBS=1
# Algorithm parameters with defaults
MIN_SIZE=10
THETA=0.9
# Graph directory with default
GRAPH_DIR="testGraphs"

# Show usage if help requested
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -d DIR    Directory containing graph subdirectories (default: testGraphs)"
    echo "  -j N      Number of parallel jobs (default: 1)"
    echo "  -l SIZE   Minimum pseudo-clique size (default: 10)"
    echo "  -t THETA  Density threshold (default: 0.9)"
    echo "  -h        Show this help message"
}

while getopts "j:l:t:d:h" opt; do
	case $opt in
		j) PARALLEL_JOBS="$OPTARG" ;;
		l) MIN_SIZE="$OPTARG" ;;
		t) THETA="$OPTARG" ;;
		d) GRAPH_DIR="$OPTARG" ;;
		h) show_usage; exit 0 ;;
		*) echo "Invalid option: -$OPTARG" >&2; show_usage; exit 1 ;;
	esac
done
shift $((OPTIND - 1))

# Validate graph directory exists
if [[ ! -d "$GRAPH_DIR" ]]; then
    echo "Error: Graph directory '$GRAPH_DIR' not found" | tee -a "$output_file"
    echo "Please specify a valid directory with -d option" | tee -a "$output_file"
    exit 1
fi

# Avoid oversubscription inside libraries when running multiple processes
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export TBB_NUM_THREADS="${TBB_NUM_THREADS:-1}"

# Per-graph logs directory
logs_dir="logs/run_${timestamp}"
mkdir -p "$logs_dir"

# Function to extract key metrics from time output
extract_metrics() {
    local time_output="$1"
    local program_name="$2"
    local log_file="$3"
    local user_time=$(echo "$time_output" | sed -n 's/.*User time (seconds): //p')
    local sys_time=$(echo "$time_output" | sed -n 's/.*System time (seconds): //p')
    local elapsed_time=$(echo "$time_output" | sed -n 's/.*Elapsed (wall clock) time (h:mm:ss or m:ss): //p')
    local max_memory_kb=$(echo "$time_output" | sed -n 's/.*Maximum resident set size (kbytes): //p')
    local total_cpu="N/A"
    if [[ -n "$user_time" && -n "$sys_time" ]]; then
        total_cpu=$(echo "$user_time + $sys_time" | bc -l 2>/dev/null || echo "N/A")
    fi
    local peak_mem_display="$max_memory_kb"
    if [[ "$max_memory_kb" =~ ^[0-9]+$ ]]; then
        peak_mem_display=$(echo "scale=2; $max_memory_kb / 1024" | bc -l 2>/dev/null)
        peak_mem_display="${peak_mem_display} MB"
    fi
    echo "=== $program_name Performance Summary ===" | tee -a "$log_file"
    echo "CPU Time (user+sys): ${total_cpu}s" | tee -a "$log_file"
    echo "Real Time: ${elapsed_time}" | tee -a "$log_file"
    echo "Peak Memory: ${peak_mem_display}" | tee -a "$log_file"
    echo "" | tee -a "$log_file"
}

# Run a command streaming program output, while capturing time -v to a temp file
run_and_report() {
    local name="$1"; shift
    local log_file="$1"; shift
    local tmp_file
    tmp_file=$(mktemp)
    echo "Running $name" | tee -a "$log_file"
    /usr/bin/time -v -o "$tmp_file" "$@" 2> >(tee -a "$log_file" >&2) | tee -a "$log_file"
    local time_output
    time_output="$(cat "$tmp_file")"
    rm -f "$tmp_file"
    extract_metrics "$time_output" "$name" "$log_file"
}

# The base directory containing subfolders with .grh files
graph_directory="$GRAPH_DIR"

# Executable
integrated_exec="./build_integrated/dense-pce-ab-integrated"

# Verify executable exists
if [[ ! -x "$integrated_exec" ]]; then
	echo "Error: $integrated_exec not found or not executable. Build it first." | tee -a "$output_file"
	exit 1
fi

handle_dir() {
	local dir="$1"
	local found_grh=$(echo "$dir"/*.grh)

	if [[ -n "$found_grh" ]]; then
		local subdir_name
		subdir_name="$(basename "$dir")"
		local graph_base
		graph_base="$(basename "$found_grh" .grh)"
		local log_file
		log_file="${logs_dir}/${subdir_name}__${graph_base}.log"

		echo "Processing directory: $dir" | tee -a "$log_file"
		echo "Using graph: $found_grh" | tee -a "$log_file"
		echo "Parameters: l=${MIN_SIZE}, theta=${THETA}" | tee -a "$log_file"
		echo "" | tee -a "$log_file"

		run_and_report "dense-pce-integrated" "$log_file" "$integrated_exec" "$found_grh" --minimum "$MIN_SIZE" --theta "$THETA"

		# The run for 'fpce' has been removed from this script.

		echo "----------------------------------------" | tee -a "$log_file"
	else
		echo "No .grh file found in: $dir" | tee -a "$output_file"
	fi
}

# Log run parameters
echo "Run parameters:" | tee -a "$output_file"
echo "  Graph directory: $GRAPH_DIR" | tee -a "$output_file"
echo "  Parallel jobs: $PARALLEL_JOBS" | tee -a "$output_file"
echo "  Minimum size: $MIN_SIZE" | tee -a "$output_file"
echo "  Theta: $THETA" | tee -a "$output_file"
echo "  Logs directory: $logs_dir" | tee -a "$output_file"
echo "" | tee -a "$output_file"

# Iterate each subdirectory in graph directory and process, optionally in parallel
for dir in "$graph_directory"/*; do
	if [[ -d "$dir" ]]; then
		if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
			handle_dir "$dir" &
			while [[ $(jobs -r -p | wc -l) -ge $PARALLEL_JOBS ]]; do
				wait -n
			done
		else
			handle_dir "$dir"
		fi
	fi
done

# Wait for all background jobs to finish (if any)
wait

echo "Script finished. Output logged to $output_file and per-graph logs under $logs_dir"
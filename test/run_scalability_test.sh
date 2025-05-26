#!/bin/bash

EXECUTABLE=main
OUTPUT_FILE="performance_results.txt"

echo "MPI Scalability Test Results" > "$OUTPUT_FILE"
echo "============================" >> "$OUTPUT_FILE"

# Option flags
VERBOSE=false
POSTPROCESS=false

# Parse options
while getopts ":vp" opt; do
  case $opt in
    v)
      VERBOSE=true
      ;;
    p)
      POSTPROCESS=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

shift $((OPTIND - 1))

# Loop over core counts
for cores in 1 2 4
do
    echo "Running with $cores cores..."
    echo -e -n "\n### $cores cores ### \n" >> "$OUTPUT_FILE"

    # Time the execution
    START=$(date +%s.%N)

    # Run the program with or without -v
    if [ "$VERBOSE" = true ]; then
        mpiexec --bind-to none -n $cores "$EXECUTABLE" -v >> "$OUTPUT_FILE" 2>&1
    else
        mpiexec --bind-to none -n $cores "$EXECUTABLE" >> "$OUTPUT_FILE" 2>&1
    fi



    END=$(date +%s.%N)

    # Calculate execution time
    ELAPSED=$(echo "$END - $START" | bc)
    echo "Execution time: 0$ELAPSED seconds" >> "$OUTPUT_FILE"

done


echo "Done. Results written to $OUTPUT_FILE."

if [ "$POSTPROCESS" = true ]; then
    echo "Running Mayavi post-processing..."
    mayavi2 -d output.vtk -m Surface -m Outline
fi
#!/bin/bash

# Script: run_all_tests.sh
# Purpose: Enter each folder in current directory, run CTurboBFM with input.ini,
#          report success/failure.

# Store current directory
BASE_DIR=$(pwd)

# Initialize counters
SUCCESS=0
FAILURE=0

# Loop over all directories in current folder
for dir in */ ; do
    # Remove trailing slash
    dir="${dir%/}"
    
    echo "=== Running test in $dir ==="
    
    # Check if input.ini exists
    if [[ ! -f "$dir/input.ini" ]]; then
        echo "  input.ini not found, skipping."
        continue
    fi
    
    # Enter directory
    cd "$dir" || { echo "  Failed to enter $dir"; continue; }
    
    # Run solver
    "$BASE_DIR"/../bin/turbobfm input.ini > log_test.txt
    EXIT_CODE=$?
    
    # Check exit code
    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "  SUCCESS"
        ((SUCCESS++))
    else
        echo "  FAILURE (exit code $EXIT_CODE)"
        ((FAILURE++))
    fi
    
    # Return to base directory
    cd "$BASE_DIR"
done

# Summary
echo "==============================="
echo "Total directories: $((SUCCESS + FAILURE))"
echo "Succeeded: $SUCCESS"
echo "Failed: $FAILURE"
echo "==============================="

echo "NOTE: a successful test only means that the solver ran without errors and the residuals decreased by 4 orders of magnitude."
echo "It does NOT guarantee that the results are physically correct."

if [[ $FAILURE -eq 0 ]]; then
    exit 0
else
    exit 1
fi

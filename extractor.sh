#!/bin/bash

# base_name: sequence pair name (e.g., human-denisovan)
# matrix: BLOSUM matrix name (e.g., BLOSUM62)
# This script extracts alignment scores, percentage identity, and estimated E-values from metrics files generated by sw.py program.

# function to check metrics file with original or reversed sequence pair names
get_metrics_file() {
    local base_name="$1"
    local matrix="$2"
    local original_file="OUTPUT/${base_name}_${matrix}_output/${base_name}_${matrix}_metrics.txt"
    
    # check if the original file exists
    if [ -f "$original_file" ]; then
        echo "$original_file"
        return
    fi
    
    # try reversed version if original not found (human-denisovan <-> denisovan-human)
    if [[ "$base_name" == *"-"* ]]; then
        local reversed_name=$(echo "$base_name" | awk -F'-' '{print $2"-"$1}')
        local reversed_file="OUTPUT/${reversed_name}_${matrix}_output/${reversed_name}_${matrix}_metrics.txt"
        
        if [ -f "$reversed_file" ]; then
            echo "$reversed_file"
            return
        fi
    fi
    
    # return empty if neither file exists
    echo ""
}
# one arg: base_name (sequence pair)
if [ "$#" -eq 1 ]; then
    base_name="$1"
    
    # extract metrics for all BLOSUM matrices of a given sequence pair
    echo "Alignment Score"
    for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
        metrics_file=$(get_metrics_file "$base_name" "$matrix")
        if [ -n "$metrics_file" ]; then
            align_score=$(grep -i "alignment score" "$metrics_file" | awk -F": " '{print $2}')
            echo "[${matrix}] ${align_score}"
        fi
    done

    echo ""
    echo "E-value"
    for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
        metrics_file=$(get_metrics_file "$base_name" "$matrix")
        if [ -n "$metrics_file" ]; then
            e_value=$(grep -i "estimated e-value" "$metrics_file" | awk -F": " '{print $2}')
            echo "[${matrix}] ${e_value}"
        fi
    done

    echo ""
    echo "Percentage Identity"
    for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
        metrics_file=$(get_metrics_file "$base_name" "$matrix")
        if [ -n "$metrics_file" ]; then
            perc_identity=$(grep -i "percentage identity" "$metrics_file" | awk -F": " '{print $2}')
            echo "[${matrix}] ${perc_identity}"
        fi
    done

# two arg: base_name (sequence pair) and BLOSUM matrix
elif [ "$#" -eq 2 ]; then
    base_name="$1"
    # convert matrix name to uppercase for case insensitivity
    matrix_name=$(echo "$2" | tr '[:lower:]' '[:upper:]')
    
    # check if it's a valid BLOSUM matrix format
    if [[ ! "$matrix_name" =~ ^BLOSUM(45|50|62|80)$ ]]; then
        echo "Error: Invalid BLOSUM matrix. Please use BLOSUM45, BLOSUM50, BLOSUM62, or BLOSUM80."
        exit 1
    fi
    
    # extract the metrics file path (original or reversed)
    metrics_file=$(get_metrics_file "$base_name" "$matrix_name")
    if [ -z "$metrics_file" ]; then
        echo "Error: Metrics file not found for '$base_name' with $matrix_name."
        exit 1
    fi

    align_score=$(grep -i "alignment score" "$metrics_file" | awk -F": " '{print $2}')
    e_value=$(grep -i "estimated e-value" "$metrics_file" | awk -F": " '{print $2}')
    perc_identity=$(grep -i "percentage identity" "$metrics_file" | awk -F": " '{print $2}')

    echo "Alignment Score: ${align_score}"
    echo "E-value: ${e_value}"
    echo "Percentage Identity: ${perc_identity}"

# zero arg:
# show all results grouped by sequence pairs
else

    # find all output directories and extract unique base names
    output_dirs=$(find OUTPUT -maxdepth 1 -type d -name "*_output" | sort)
    
    if [ -z "$output_dirs" ]; then
        # If no output directories found, show help message
        echo "No alignment output directories found."
        echo "Usage:"
        echo "  ./extractor.sh                          # Show results for all sequence pairs"
        echo "  ./extractor.sh <sequence-pair>          # Show results for all BLOSUM matrices"
        echo "  ./extractor.sh <sequence-pair> <matrix> # Show results for specific BLOSUM matrix"
        echo ""
        echo "Examples:"
        echo "  ./extractor.sh denisovan-human          # All matrices for denisovan-human"
        echo "  ./extractor.sh human-denisovan          # Works the same as above (reversible)"
        echo "  ./extractor.sh denisovan-human BLOSUM80 # Only BLOSUM80 for denisovan-human"
        echo "  ./extractor.sh denisovan-human blosum45 # Case insensitive matrix name"
        exit 0
    fi
    
    # extract unique base pairs (before first underscore) without using associative arrays
    base_pairs=$(for dir in $output_dirs; do
        dir_name=$(basename "$dir" | sed 's/_output$//')
        echo "$dir_name" | cut -d'_' -f1
    done | sort -u)
    
    # list alignment scores for all base pairs across all matrices
    echo "Alignment Score"
    for base_pair in $base_pairs; do
        echo "($base_pair)"
        for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
            metrics_file="OUTPUT/${base_pair}_${matrix}_output/${base_pair}_${matrix}_metrics.txt"
            if [ -f "$metrics_file" ]; then
                align_score=$(grep -i "alignment score" "$metrics_file" | awk -F": " '{print $2}')
                echo "[${matrix}] ${align_score}"
            fi
        done
        echo ""
    done
    
    echo "E-value"
    for base_pair in $base_pairs; do
        echo "($base_pair)"
        for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
            metrics_file="OUTPUT/${base_pair}_${matrix}_output/${base_pair}_${matrix}_metrics.txt"
            if [ -f "$metrics_file" ]; then
                e_value=$(grep -i "estimated e-value" "$metrics_file" | awk -F": " '{print $2}')
                echo "[${matrix}] ${e_value}"
            fi
        done
        echo ""
    done

    echo "Percentage Identity"
    for base_pair in $base_pairs; do
        echo "($base_pair)"
        for matrix in "BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80"; do
            metrics_file="OUTPUT/${base_pair}_${matrix}_output/${base_pair}_${matrix}_metrics.txt"
            if [ -f "$metrics_file" ]; then
                perc_identity=$(grep -i "percentage identity" "$metrics_file" | awk -F": " '{print $2}')
                echo "[${matrix}] ${perc_identity}"
            fi
        done
        echo ""
    done
fi
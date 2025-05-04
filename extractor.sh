#!/bin/bash

if [ "$#" -eq 1 ]; then
    # Use the provided folder base name to extract metrics from that folder.
    folder_base="$1"
    metrics_file="${folder_base}_output/${folder_base}_metrics.txt"
    if [ ! -f "$metrics_file" ]; then
        echo "Error: Metrics file '$metrics_file' does not exist."
        exit 1
    fi

    align_score=$(grep -i "alignment score" "$metrics_file" | awk -F": " '{print $2}')
    e_value=$(grep -i "estimated e-value" "$metrics_file" | awk -F": " '{print $2}')

    echo "Alignment Score : ${align_score}"
    echo "E-value: ${e_value}"
else
    # No parameter provided: process all *_metrics.txt in any *_output folder.
    echo "Alignment Score:"
    for metrics_file in ./*_output/*_metrics.txt; do
        # Extract folder base from the filename (assuming file is like folder_base_metrics.txt)
        folder_base=$(basename "$metrics_file" _metrics.txt)
        align_score=$(grep -i "alignment score" "$metrics_file" | awk -F": " '{print $2}')
        echo "[${folder_base}] : ${align_score}"
    done

    echo ""
    echo "E-value:"
    for metrics_file in ./*_output/*_metrics.txt; do
        folder_base=$(basename "$metrics_file" _metrics.txt)
        e_value=$(grep -i "estimated e-value" "$metrics_file" | awk -F": " '{print $2}')
        echo "[${folder_base}] : ${e_value}"
    done
fi
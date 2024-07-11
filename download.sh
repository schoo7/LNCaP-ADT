#!/bin/bash

# Set the initial parameters
input_file="/home/simon/Documents/ids.csv"
params_template="nf-params.json"
work_dir="/media/simon/work/work"
download_method="sratools"
chunk_size=100

# Get the total number of lines in the input file
total_lines=$(wc -l < "$input_file")

# Calculate the number of chunks
num_chunks=$(( (total_lines + chunk_size - 1) / chunk_size ))

for ((i=0; i<num_chunks; i++)); do
    start_line=$(( i * chunk_size + 1 ))
    end_line=$(( start_line + chunk_size - 1 ))
    
    # Create the new CSV file for the current chunk
    chunk_file="/home/simon/Documents/ids_$((i+1)).csv"
    sed -n "${start_line},${end_line}p" "$input_file" > "$chunk_file"
    
    # Update the nf-params.json file
    chunk_outdir="/media/simon/data/ctpc_$((i+1))"
    params_file="nf-params_$((i+1)).json"
    jq --arg input "$chunk_file" --arg outdir "$chunk_outdir" --arg method "$download_method" \
    '.input = $input | .outdir = $outdir | .download_method = $method' "$params_template" > "$params_file"
    
    # Run the nextflow command
    nextflow run nf-core/fetchngs -r 1.12.0 -profile docker -work-dir "$work_dir" -resume -params-file "$params_file"
done


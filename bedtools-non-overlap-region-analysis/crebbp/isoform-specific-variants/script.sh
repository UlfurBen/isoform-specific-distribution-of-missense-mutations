#!/bin/bash

# Define the folder containing the BED files
BED_FOLDER="crebbp"

# Change to the directory containing the BED files
cd "$BED_FOLDER"

# Get a list of all .bed files in the directory
BED_FILES=(*.bed)

# Loop over each file as the primary file (file1)
for file1 in "${BED_FILES[@]}"; do
    # Create a temporary file to store intermediate results
    temp_file=$(mktemp)
    
    # Initialize the intermediate result with the content of file1
    cp "$file1" "$temp_file"
    
    # Loop over each file as the secondary file (file2)
    for file2 in "${BED_FILES[@]}"; do
        # Skip if file1 and file2 are the same
        if [[ "$file1" != "$file2" ]]; then
            # Perform bedtools intersect and update the intermediate result
            bedtools intersect -v -a "$temp_file" -b "$file2" > "${temp_file}_new"
            mv "${temp_file}_new" "$temp_file"
        fi
    done
    
    # Save the final unique result to a file
    mv "$temp_file" "${file1%.bed}_uniq.bed"
done

# Cleanup temporary files if any are left
rm -f "${temp_file}_new"

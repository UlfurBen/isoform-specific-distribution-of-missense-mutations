#!/bin/bash

# Define the folder containing the BED files
BED_FOLDER="crebbp"
EXPECTED_FIELDS=26  # Adjust this number to the expected number of fields

# Change to the directory containing the BED files
cd "$BED_FOLDER" || { echo "Directory not found: $BED_FOLDER"; exit 1; }

# Get a list of all .bed files in the directory
BED_FILES=(*.bed)

# Step 1: Clean the BED files to ensure proper tab separation
echo "Cleaning BED files..."
for file in "${BED_FILES[@]}"; do
    # Replace any sequence of spaces and tabs with a single tab
    awk '{$1=$1}1' OFS='\t' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done

# Step 2: Verify and correct the number of fields in each line
echo "Verifying and correcting field counts..."
for file in "${BED_FILES[@]}"; do
    echo "Processing file: $file"
    awk -v expected_fields=$EXPECTED_FIELDS -F'\t' '
    {
        if (NF != expected_fields) {
            print "Warning: Line " NR " has " NF " fields instead of " expected_fields
        } else {
            print $0
        }
    }' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done

# Step 3: Find unique intervals using bedtools intersect
echo "Finding unique intervals..."
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

# Cleanup temporary files
rm -f "${temp_file}"

echo "Processing complete."

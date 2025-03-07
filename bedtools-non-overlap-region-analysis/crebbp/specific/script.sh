cd crebbp
BED_FILES=(*.bed)

# Step 3: Find unique intervals using bedtools intersect
echo "Finding unique intervals..."

# Create a directory to store the unique output files
output_dir="unique_intervals"
mkdir -p "$output_dir"

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
    
    # Save the final unique result to the output directory
    mv "$temp_file" "${output_dir}/${file1%.bed}_uniq.bed"
done

echo "Processing complete. Unique files are saved in the '$output_dir' directory."

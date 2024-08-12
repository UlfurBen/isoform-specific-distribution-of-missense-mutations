#!/bin/bash

# Define the path to the gene names with ENST values CSV file
CSV_FILE=~/gene_names_with_enst.csv

# Define the path to the "genes" folder
GENES_FOLDER="benign-genes"

# Read the CSV file line by line
while IFS=, read -r GENE_NAME ENST_ID; do
    # Navigate to the gene's folder
    GENE_FOLDER="${GENES_FOLDER}/${GENE_NAME}"
    
    if [ -d "$GENE_FOLDER" ]; then
        cd "$GENE_FOLDER" || exit
        
        # Get a list of .bed files
        BED_FILES=(*.bed)
        
        # Check and remove trailing tabs in all .bed files
        for BED_FILE in "${BED_FILES[@]}"; do
            if [ -f "$BED_FILE" ]; then
                # Remove trailing tabs from each line in the BED file
                sed -i 's/\t$//' "$BED_FILE"
            fi
        done
        
        # Check if there are more than one .bed files in the folder
        if [ "${#BED_FILES[@]}" -gt 1 ]; then
            # Delete the .bed file that matches the ENST ID
            BED_FILE="${ENST_ID}.bed"
            if [ -f "$BED_FILE" ]; then
                echo "Deleting $BED_FILE in $GENE_FOLDER"
                rm "$BED_FILE"
            fi
            
            # Update the list of remaining .bed files
            BED_FILES=(*.bed)
        fi

        # Check if there are two or more .bed files remaining
        if [ "${#BED_FILES[@]}" -ge 2 ]; then
            # Create a directory to store the unique output files
            OUTPUT_DIR="unique_intervals"
            mkdir -p "$OUTPUT_DIR"
            
            # Loop over each bed file to find unique intervals
            for FILE1 in "${BED_FILES[@]}"; do
                # Create a temporary file to store intermediate results
                TEMP_FILE=$(mktemp)
                
                # Initialize the intermediate result with the content of file1
                cp "$FILE1" "$TEMP_FILE"
                
                # Loop over each file as the secondary file (FILE2)
                for FILE2 in "${BED_FILES[@]}"; do
                    # Skip if FILE1 and FILE2 are the same
                    if [[ "$FILE1" != "$FILE2" ]]; then
                        # Perform bedtools intersect and update the intermediate result
                        bedtools intersect -v -a "$TEMP_FILE" -b "$FILE2" > "${TEMP_FILE}_new"
                        mv "${TEMP_FILE}_new" "$TEMP_FILE"
                    fi
                done
                
                # Save the final unique result to the output directory
                mv "$TEMP_FILE" "${OUTPUT_DIR}/${FILE1%.bed}_uniq.bed"
            done
            
            echo "Unique intervals for ${GENE_NAME} have been saved in ${OUTPUT_DIR}/"
        else
            echo "Not enough bed files in ${GENE_NAME} to perform intersection."
        fi
        
        # Return to the root directory
        cd - > /dev/null
    else
        echo "Gene folder $GENE_FOLDER does not exist."
    fi
done < "$CSV_FILE"

echo "All processing complete."

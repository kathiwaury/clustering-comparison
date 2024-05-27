#!/bin/bash

# Check if the user provided a file containing PDB IDs
if [ $# -ne 1 ]; then
    echo "Usage: $0 <pdb_id_file>"
    exit 1
fi

# Create folder to store PDB files
mkdir -p pdb_files

# Loop through each PDB ID in the file
while IFS= read -r pdb_id; do
    # Retrieve PDB file
    wget -q -O "pdb_files/${pdb_id}.pdb" "https://files.rcsb.org/download/${pdb_id}.pdb"
    
    # Check if retrieval was successful
    if [ $? -eq 0 ]; then
        echo "Successfully retrieved ${pdb_id}.pdb"
    else
        echo "Failed to retrieve ${pdb_id}.pdb"
    fi
done < "$1"

echo "All PDB files retrieved and saved to pdb_files folder."


#!/bin/bash
# Script to run inferCNV analysis pipeline
#PBS -N tianjin_funcs            # Job name for PBS scheduler
#PBS -l walltime=123:45:06       # Maximum wall time for job execution (123 hours, 45 minutes, 6 seconds)
#PBS -l nodes=1:ppn=1            # Request 1 node and 1 processor per node
#PBS -o /home/bowang/tj.log      # Output log file for standard output (stdout)
#PBS -j oe                       # Merge standard error (stderr) into the standard output (stdout) log file
#PBS -V                          # Export all current environment variables to the job

# Set project identifier
PROJ=prj301                     # Project variable (identifier for organizing related files or tasks)

# Increase stack size limit to prevent stack overflow during execution
ulimit -s 32768                 # Set the stack size limit to 32 MB

# Define a function to run inferCNV analysis
indercnv() {
    # Define function arguments
    local raw_counts_matrix=$1  # Argument 1: Path to the raw counts matrix file
    local annotations_file=$2  # Argument 2: Path to the annotations file
    local gene_order_file=$3   # Argument 3: Path to the gene order file
    local ref_group_names=$4   # Argument 4: Reference group names (R-compatible format)
    local out_dir=$5           # Argument 5: Directory to store outputs
    local atlas=$6             # Argument 6: (Optional) Path to atlas file, not used in current logic

    # Create output directory if it does not exist
    mkdir -p $out_dir

    # Preprocess annotations file: replace hyphens (-) with dots (.)
    if [ ! -f label.2.txt ]; then
        awk -v OFS='\t' '{gsub("-",".",$1);print}' label.txt > label.2.txt
    fi

    # Run inferCNV R script with specified arguments
    R CMD BATCH --no-save --no-restore \
        "--args $raw_counts_matrix $annotations_file $gene_order_file $ref_group_names $out_dir" \
        $path_R/infercnv.R $out_dir/infercnv.log

}

# Main script execution starts here
# Change directory to the input folder for the experiment
cd $wd/exp1/input

# Prepare reference groups for R script
# Format content of ref.txt to an R-compatible vector (e.g., c('group1','group2'))
ref=$(cat ref.txt | sed "s/,/','/g" | sed "s/^/'/" | sed "s/$/'/")

# Print formatted reference groups (for debugging or logging)
echo $ref

# Run the inferCNV pipeline using the `indercnv` function
# Arguments:
#   1. counts.mat: Raw counts matrix
#   2. label.2.txt: Processed annotations file
#   3. seurat.symbol.anno.txt: Gene order file
#   4. "c($ref)": Reference group names in R-compatible format
#   5. $od/exp1/input: Output directory
indercnv counts.mat label.2.txt seurat.symbol.anno.txt "c($ref)" $od/exp1/input
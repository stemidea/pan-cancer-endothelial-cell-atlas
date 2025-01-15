#!/bin/bash
#PBS -N tianjin_funcs  # Job name is set as 'tianjin_funcs'
#PBS -l walltime=123:45:06  # Max allowed execution time (123 hours, 45 minutes, 6 seconds)
#PBS -l nodes=1:ppn=1  # Specifies resources: 1 node, 1 processor per node
#PBS -o /home/bowang/tj.log  # Path to the output log file
#PBS -j oe  # Merges standard output and error streams
#PBS -V  # Exports environment variables to the job

### FUNCTION: count2loom
# Converts single-cell RNA-seq count data to a `loom` file format using velocyto.
count2loom(){
    local sampleID=$1  # The sample ID (directory containing raw data)
    local gtf=$2  # Optional: Path to the GTF file
    local transcriptome=$3  # Optional: Transcriptome annotation file

    # Use default paths if not provided
    if [ -z $gtf ]; then
        gtf=$gtf_human  # Default GTF for human
    fi
    if [ -z $transcriptome ]; then
        transcriptome=$transcriptome_hs/genes/genes.gtf  # Default transcriptome
    fi

    # Run velocyto to generate the loom file (uses 16 threads for speed)
    # Logs output to $sampleID.count2loom.log
    velocyto run10x -m $gtf $wd/count/$sampleID $transcriptome -@ 16 > $sampleID.count2loom.log
}

### Execution: Call count2loom for multiple samples
# Converts raw count data for each sample to loom format.
# Add more sample IDs as needed...
# List of sample IDs
samples=("GY_20230719_T" "HFQ_20230810_T" "JY_20230726_T" "KSF_20230607_N" "KSF_20230607_T" "MAL_20230905_T" 
         "QL_20230817_T" "SM_20230830_T" "SXM_20230714_T" "WJD_20230724_T" "XZ_20230712_T" 
         "ZGQ_20230904_T" "ZYD_20230810_T")

# Convert count data for all samples
for sample in "${samples[@]}"; do
    count2loom $sample
done

### FUNCTION: looms2h5ad
# Merges multiple loom files into a single h5ad file for downstream analysis.
looms2h5ad(){
    local looms=$1  # Comma-separated list of loom file paths
    local sample_obs=$2  # Metadata CSV containing sample observations
    local out=$3  # Output name for the merged h5ad file

    # Run the Python script to perform the merging
    # Logs the output to $wd/velocyto/looms2h5ad.log
    python3 $path_Python/szzy_merge_looms2h5ad.py $looms $sample_obs $out > $wd/velocyto/looms2h5ad.log 2>&1
}

### Execution: Merge loom files into a single h5ad
# Collects all loom files in the `count` directory and merges them.
looms=$(ls $wd/count/*/velocyto/*.loom -1 | tr '\n' ',' | sed 's/,$//')  # Generates a comma-separated list of loom files
looms2h5ad $looms cellID_obs.csv $PROJ  # Merges them into $PROJ.h5ad

### FUNCTION: velo_py
# Runs RNA velocity analysis on the merged h5ad file and generates velocity plots.
velo_py(){
    local h5ad=$1  # Input h5ad file (merged loom data)
    local anno=$2  # Annotation file (e.g., cell subtypes)
    local outpre=$3  # Prefix for output files

    # Run the Python script to perform RNA velocity analysis
    python3 $path_Python/velo.py $h5ad cellID_obs.csv cell_embeddings.csv clusters.csv $anno $outpre > velocyto.log 2>&1
}

### Execution: RNA velocity analysis
# Example call for RNA velocity analysis on the merged h5ad file.
velo_py $PROJ.h5ad subtype ATLAS

# Alternatively, this line can be used to execute the velocity analysis.
# python3 $path_Python/velo.py $PROJ.h5ad cellID_obs.csv cell_embeddings.csv clusters.csv subtype ATLAS > velocyto.log 2>&1


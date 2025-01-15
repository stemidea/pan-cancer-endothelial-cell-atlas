#!/bin/bash
#PBS -N tianjin_funcs                 # Job name for the PBS scheduler
#PBS -l walltime=123:45:06            # Maximum runtime for the job (123 hours, 45 minutes, 6 seconds)
#PBS -l nodes=1:ppn=1                 # Requests 1 node with 1 processor per node
#PBS -o /home/bowang/tj.log           # Redirect standard output to the specified log file
#PBS -j oe                            # Merge standard error into the standard output log
#PBS -V                               # Export current environment variables to the job environment

# Activate the PySCENIC conda environment
conda activate pyscenic               # Activate the Conda environment containing PySCENIC
PROJ=prj301                           # Project identifier (used for organizing outputs)

mkdir -p $od                          # Ensure the output directory ($od) exists
cd $od                                # Change to the output directory

# Define paths to custom Python scripts and SCENIC resources
export path_python=$home/Prog/Python  # Path to the directory containing custom Python programs
Ref=$root/home/software/SCENIC        # Base directory for SCENIC-related reference files

# Define SCENIC resources for human (hg38) and mouse (mm9)
ref_hg38=$Ref/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather    # Human genome reference file for SCENIC
motif_hg38=$Ref/motifs-v9-nr.hgnc-m0.001-o0.0.tbl                    # Human motif annotations
tfs_hg38=$Ref/hs_hgnc_tfs.txt                                        # Human transcription factor list

ref_mm9=$Ref/mm9-tss-centered-10kb-7species.mc9nr.feather            # Mouse genome reference file for SCENIC
motif_mm9=$Ref/motifs-v9-nr.mgi-m0.001-o0.0.tbl                      # Mouse motif annotations
tfs_mm9=$Ref/mm_mgi_tfs.txt                                          # Mouse transcription factor list
pyScenic3in1() {
    local csv=$1      # Argument 1: Input CSV file containing filtered gene expression counts
    local sp=$2       # Argument 2: Species identifier ('hs' for human, 'mm' for mouse)

    # Default to human ('hs') if the species is not specified
    if [ -z $sp ]; then
        sp=hs
    fi

    # Assign species-specific resources based on the species identifier
    if [ "$sp" == "hs" ]; then
        ref=$ref_hg38
        motif=$motif_hg38
        tfs=$tfs_hg38
    else
        ref=$ref_mm9
        motif=$motif_mm9
        tfs=$tfs_mm9
    fi

    # Step 1: Convert CSV to Loom
    python $root/home/bowang/Prog/Python/csv2loom.py $csv $(basename $csv).sample.loom

    # Step 2: Run PySCENIC pipeline with species-specific resources
    pyScenic $(basename $csv) $ref $motif $tfs

    # Step 3: Export results using R
    R CMD BATCH --no-save --no-restore \
        "--args csv='$csv'" \
        $root/home/bowang/Prog/R/pySCENIC.R pySCENIC.R.$(basename $csv).log
}
# Change to the directory containing filtered gene expression data
cd $wd/exp3/filtered_counts

# Loop through all filtered gene expression CSV files in the directory
for csv_file in downsampled_atlas_*_filtered_counts.csv; do
    echo "Processing $csv_file..."
    pyScenic3in1 $csv_file
done


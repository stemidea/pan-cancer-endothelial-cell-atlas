# Load the Seurat library
library(Seurat)

# Create a directory to store Seurat objects
dir.create("./seurat/")  
# This command creates a new folder named "seurat" in the current working directory to save the output files.

# Get a list of subdirectories inside the "./data" directory
my.dir <- list.dirs(path = "./data", recursive = F)  
# `list.dirs()` lists all directories (non-recursively). Each directory is expected to contain raw gene expression data.

# Extract the names of the subdirectories without the "./data/" prefix
my.atlas <- gsub("./data/", "", my.dir)  
# `gsub()` removes "./data/" from the paths, leaving only the names of the subdirectories (e.g., "sample1", "sample2").

# Loop through each directory to process its data
for (i in 1:length(my.dir)){  
  # Iterates over each directory in `my.dir`

  # Read the 10X Genomics data from the directory
  atlas.data <- Read10X(data.dir=my.dir[i])  
  # `Read10X()` reads the gene expression matrix (e.g., barcodes, features, and matrix files) from the specified directory.

  # Create a Seurat object from the expression matrix
  atlas <- CreateSeuratObject(counts = atlas.data, project = my.atlas[i])  
  # `CreateSeuratObject()` initializes a Seurat object with the gene expression matrix (`counts`).
  # The `project` parameter assigns a name to the project, based on the subdirectory name (`my.atlas[i]`).

  # (Optional) Calculate the percentage of mitochondrial gene content
  # atlas[["percent.mt"]] <- PercentageFeatureSet(atlas, pattern = "^MT-")  
  # Uncomment this line if you want to calculate mitochondrial gene percentage. 
  # This is typically used for quality control in single-cell RNA-seq analysis.
  # "^MT-" identifies mitochondrial genes in human datasets (adjust pattern for other species).

  # Save the Seurat object to the "seurat" directory
  saveRDS(atlas, file = paste("./seurat/", my.atlas[i], ".rds", sep = ""))  
  # `saveRDS()` saves the Seurat object to a file named after the subdirectory (e.g., "sample1.rds").
  # The files are stored in the "./seurat/" directory for later use.
}


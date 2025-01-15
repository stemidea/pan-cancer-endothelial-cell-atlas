# Load required libraries
library(CytoTRACE)  # For computing differentiation potential from single-cell data
library(Seurat)     # For single-cell data analysis and visualization
library(tidyverse)  # For data manipulation and visualization

# Set output directory for results
outs <- "./outs.cytoTRACE/"
out.dir <- outs

# Load Seurat object containing single-cell data
atlas <- readRDS("atlas.rds")

# Inspect metadata to understand available annotations and clustering
colnames(atlas@meta.data)       # Check the column names in the metadata
unique(atlas@meta.data$seurat_clusters)  # Check unique cluster identifiers

# Extract the raw count matrix from the Seurat object
mat <- as.matrix(GetAssayData(atlas, slot = "counts"))  
# The matrix will be used as input for CytoTRACE, representing gene expression counts.

# Set a random seed for reproducibility
set.seed(123)

# Run CytoTRACE to calculate differentiation potential for each cell
results <- CytoTRACE(mat)  
# `results` contains differentiation potential scores and associated gene information.

# Extract cluster metadata for associating differentiation potential with clusters
meta <- atlas@meta.data %>%
  select(seurat_clusters)  # Keep only the cluster information

# Prepare phenotype vector from cluster metadata
pheno <- as.vector(meta$seurat_clusters)  # Convert cluster metadata to a vector
names(pheno) <- rownames(meta)            # Assign cell barcodes as names for the phenotype vector

# Extract UMAP embeddings for visualization
umap <- as.data.frame(Embeddings(object = atlas[["umap"]]))  
# UMAP embeddings are used for plotting differentiation potential in 2D space.

# Visualize CytoTRACE results with specific genes and UMAP embeddings
# VWF: An endothelial marker; ACTA2: A marker for smooth muscle cells.
plotCytoTRACE(results, phenotype = pheno, emb = umap, gene = "VWF")  
plotCytoTRACE(results, phenotype = pheno, emb = umap, gene = "ACTA2")

# Visualize top 10 genes contributing to differentiation potential
plotCytoGenes(results, numOfGenes = 10)  

# Save CytoTRACE results for further analysis or sharing
saveRDS(results, "cytotrace.rds")  


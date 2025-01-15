# Load required libraries
library(Seurat)  # For single-cell analysis
library(dplyr)   # For data manipulation
library(ggplot2) # For data visualization
library(RColorBrewer) # For color palettes in plots

# Parse command-line arguments
args = commandArgs()
## Command-line usage example: Rscript cluster.R 1 20 res1mt20
args[1]  # First argument (typically the script name)
args[2]  # Unused in the script
args[6]  # Resolution for clustering
args[7]  # Unused in the script
args[8]  # Output directory name (e.g., "res1mt20")

# Set the working directory
setwd("./")  # Use the current directory
wd <- "./"  # Define working directory variable for later use

# Create output directory
dir.create(paste("./", args[8], sep = ""))
out_dir <- paste("./", args[8], "/", sep = "")
out_dir  # Print the output directory path

# Load the Seurat object
atlas <- readRDS(file = "atlas.combined.rds")  
# Load the integrated single-cell dataset

# Inspect metadata
meta <- atlas@meta.data  
# Metadata for the single cells, stored in the Seurat object

# Set the cell identities to the original sample identifiers
Idents(atlas) <- atlas[["orig.ident"]]

# Generate QC metrics plot
pdf(file = paste(out_dir, "QCmetrics.pdf", sep = ""), width = 15, height = 5, useDingbats = FALSE)
VlnPlot(atlas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)  
# Violin plots for QC metrics: number of features, counts, and mitochondrial gene percentages
dev.off()

# Perform scaling and dimensional reduction
atlas <- ScaleData(atlas, verbose = FALSE)  # Scale data to prepare for PCA
atlas <- RunPCA(atlas, features = VariableFeatures(object = atlas))  # Run PCA on variable features
atlas <- RunUMAP(atlas, dims = 1:30)  # Compute UMAP for visualization (30 dimensions)
atlas <- RunTSNE(atlas, dims = 1:30)  # Compute t-SNE for visualization (30 dimensions)

# Identify neighbors and clusters
atlas <- FindNeighbors(atlas, reduction = "pca", dims = 1:20)  # Find nearest neighbors using PCA
atlas <- FindClusters(atlas, resolution = as.numeric(args[6]))  # Cluster cells using specified resolution
table(atlas@meta.data$orig.ident, atlas@meta.data$seurat_clusters)  
# Show the distribution of clusters across original samples

# Visualization of clusters using UMAP and t-SNE
pdf(paste(out_dir, "umap.pdf", sep =""), width = 12, height = 5, useDingbats = FALSE)
p1 <- DimPlot(atlas, reduction = "umap", group.by = "orig.ident")  # UMAP colored by sample
p2 <- DimPlot(atlas, reduction = "umap", label = TRUE, repel = TRUE)  # UMAP with cluster labels
p1 + p2  # Combine the two plots
dev.off()

pdf(paste(out_dir, "tsne.pdf", sep =""), width = 12, height = 5, useDingbats = FALSE)
p1 <- DimPlot(atlas, reduction = "tsne", group.by = "orig.ident")  # t-SNE colored by sample
p2 <- DimPlot(atlas, reduction = "tsne", label = TRUE, repel = TRUE)  # t-SNE with cluster labels
p1 + p2  # Combine the two plots
dev.off()

# Split UMAP plots by sample and cluster
pdf(paste(out_dir, "split.umap.pdf", sep =""), width = 10, height = 5, useDingbats = FALSE)
DimPlot(atlas, reduction = "umap", split.by = "orig.ident")  # Split UMAP by sample
DimPlot(atlas, reduction = "tsne", split.by = "orig.ident")  # Split t-SNE by sample
dev.off()

pdf(paste(out_dir, "split.sc.umap.pdf", sep =""), width = 30, height = 3, useDingbats = FALSE)
DimPlot(atlas, reduction = "umap", split.by = "seurat_clusters")  # Split UMAP by clusters
dev.off()

# Calculate cluster averages and visualize as heatmap
DefaultAssay(atlas) <- "RNA"
Idents(atlas) <- atlas[["seurat_clusters"]]  # Set cell identities to clusters
cluster.averages <- AverageExpression(atlas, return.seurat = T)  
# Compute average expression per cluster

features.markers <- c("CDT1", "GMNN", "TOP2A", "PCNA", "PTPRC", "CD2", 
                      "CD3E", "CD3D", "CD3G", "CD4", "CD8A", "CD79A", 
                      "CD79B", "MS4A1", "NCR1", "KIR3DL2", "GZMB", 
                      "KLRD1", "CD74", "HLA-DRA", "CD14", "FCGR3A", 
                      "FCN1", "CD163", "C1QA", "C1QB", "C1QC", "ITGAX", 
                      "ITGAM", "S100A8", "S100A9", "MAFB", "VSIG4", 
                      "EPCAM", "KRT5", "UPK3B", "KRT8", "KLK3", "PECAM1", 
                      "VWF", "CD34", "ACTA2", "CD36", "TAGLN", "TPM2", 
                      "MMP2", "SFRP2", "COL11A1", "CD44", "ICAM1", 
                      "GFAP", "MBP", "ANG", "SPARCL1", "COL6A1", 
                      "IFITM1", "GNG11", "COL18A1", "IGFBP4", "MCAM", 
                      "SPARC", "COL4A2", "COL6A2", "MGP", "CD248", 
                      "RASD2", "PLXDC1", "ARHGEF17", "ADGRA2", "TNS3", 
                      "ANTXR1", "MRC2")  # Selected markers for visualization

pdf(file = paste(out_dir, "cluster_averages.pdf", sep = ""), width = 10, height = 10, useDingbats = FALSE)
DoHeatmap(cluster.averages, features = features.markers, size = 3, draw.lines = F)  
# Heatmap of average expression for selected markers
dev.off()

# Save the updated Seurat object
saveRDS(atlas, paste(out_dir, "atlas.rds", sep = ""))  


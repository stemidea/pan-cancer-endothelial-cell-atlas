# Load necessary libraries
library(SingleR)   # For automated cell type annotation
library(celldex)   # For pre-built reference datasets for SingleR
library(Seurat)    # For single-cell data processing and visualization
library(tidyverse) # For data manipulation and visualization

# Create output directory for storing results
out.dir <- "./outs.singleR/"
dir.create(out.dir)  # Create directory if it doesn't already exist

# Load the Human Primary Cell Atlas reference data for annotation
hpca <- readRDS("./HumanPrimaryCellAtlasData.rds")  
# hpca contains cell type labels and expression profiles for reference

# Load the single-cell dataset
OBJ <- readRDS("atlas.rds")  # Load the Seurat object with single-cell data

# Ensure that Seurat cluster IDs are factors
OBJ@meta.data$seurat_clusters <- OBJ@meta.data$seurat_clusters[drop = TRUE]  

# Set the default assay to "RNA" for annotation
DefaultAssay(OBJ) <- "RNA"

# Extract expression data for SingleR
sce_for_SingleR <- GetAssayData(OBJ, slot = "data")  
# Extract normalized expression data for all cells

# Retrieve cluster annotations from metadata
clusters <- OBJ@meta.data$seurat_clusters  

# Perform cell type annotation using the main labels in the reference
pred.hpca.main <- SingleR(
  test = sce_for_SingleR,             # Test dataset (cells to annotate)
  ref = hpca,                         # Reference dataset
  labels = hpca$label.main,           # Main cell type labels from the reference
  clusters = clusters,                # Clusters to annotate
  assay.type.test = "logcounts",      # Assay type for test dataset
  assay.type.ref = "logcounts"        # Assay type for reference dataset
)

# Perform cell type annotation using the fine-grained labels in the reference
pred.hpca.fine <- SingleR(
  test = sce_for_SingleR,
  ref = hpca,
  labels = hpca$label.fine,           # Fine-grained cell type labels
  clusters = clusters,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

# Create a data frame to store annotations for each cluster
cellType <- data.frame(
  ClusterID = levels(OBJ@meta.data$seurat_clusters), 
  HPCA.main = pred.hpca.main$labels, # Main labels predicted by SingleR
  HPCA.fine = pred.hpca.fine$labels  # Fine-grained labels predicted by SingleR
)

# Map predicted annotations to the Seurat metadata
OBJ@meta.data$HPCA.main <- cellType[match(clusters, cellType$ClusterID), 'HPCA.main']
OBJ@meta.data$HPCA.fine <- cellType[match(clusters, cellType$ClusterID), 'HPCA.fine']

# Visualize the main annotations on UMAP
p <- DimPlot(OBJ, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "HPCA.main")
pdf(paste(out.dir, "main.pdf", sep = ""), width = 8, height = 6, useDingbats = FALSE)
print(p)  # Save the UMAP plot with main annotations
dev.off()

# Visualize the fine-grained annotations on UMAP
p <- DimPlot(OBJ, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "HPCA.fine")
pdf(paste(out.dir, "fine.pdf", sep = ""), width = 9, height = 6, useDingbats = FALSE)
print(p)  # Save the UMAP plot with fine annotations
dev.off()

# Retrieve metadata from the Seurat object
meta <- OBJ@meta.data  

# Count the number of cells per Seurat cluster and HPCA main annotation
df <- meta %>%
  group_by(seurat_clusters, HPCA.main) %>%
  count()  # Count the number of cells in each cluster-label combination

# Save the annotation results as a CSV file
write.csv(df, paste(out.dir, "singleR.anno.sc.csv", sep = ""), row.names = FALSE)

# Save the updated Seurat object with SingleR annotations
saveRDS(OBJ, "atlas.rds")  


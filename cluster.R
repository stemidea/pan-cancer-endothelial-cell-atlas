# Libraries ---------------------------------------------------------------
# Load required libraries for analysis and visualization
library(Seurat)
library(reshape2)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(formattable)

# Output Directory --------------------------------------------------------
# Define and create the output directory for results
out_dir <- "./outs.cluster/"
dir.create(out_dir)

# Load Data ---------------------------------------------------------------
# Load the integrated Seurat object
merged_seurat <- readRDS("atlas.rds")

# Display the dimensions of the dataset
dim(merged_seurat)

# Add Metadata ------------------------------------------------------------
# Calculate mitochondrial and immunoglobulin content percentages
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
merged_seurat[["percent.ig"]] <- PercentageFeatureSet(merged_seurat, pattern = "^IG")

# Plot Percentages --------------------------------------------------------
# Plot distribution of immunoglobulin content
hist(merged_seurat@meta.data$percent.ig, breaks = 1000)

# Set dataset as identity class
Idents(merged_seurat) <- merged_seurat[["dataset"]]

# Save violin plots for mitochondrial and immunoglobulin percentages
pdf(paste(out_dir, "percent.mt.pdf", sep = ""), width = 18, height = 5, useDingbats = F)
VlnPlot(merged_seurat, features = "percent.mt", pt.size = 0) + NoLegend()
dev.off()

pdf(paste(out_dir, "percent.ig.pdf", sep = ""), width = 18, height = 5, useDingbats = F)
VlnPlot(merged_seurat, features = "percent.ig", pt.size = 0) + NoLegend()
dev.off()

# Save violin plots for RNA metrics
pdf(paste(out_dir, "gene.pdf", sep = ""), width = 18, height = 5, useDingbats = F)
VlnPlot(merged_seurat, features = "nFeature_RNA", pt.size = 0) + NoLegend()
dev.off()

pdf(paste(out_dir, "umi.pdf", sep = ""), width = 18, height = 5, useDingbats = F)
VlnPlot(merged_seurat, features = "nCount_RNA", pt.size = 0) + NoLegend()
dev.off()

# Filter Data -------------------------------------------------------------
# Filter cells based on quality metrics
dim(merged_seurat)
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nCount_RNA < 50000 & percent.mt < 25 & percent.ig < 2)
dim(merged_seurat)

# Remove specific datasets
merged_seurat <- subset(merged_seurat, subset = dataset != "GSE131928")
merged_seurat <- subset(merged_seurat, subset = dataset != "GSE227718")
dim(merged_seurat)

# Exclude Cells Expressing Marker Genes -----------------------------------
# Count marker gene expression in cells
markers <- c("PTPRC", "EPCAM", "HBB", "HBA2", "JCHAIN", "CD79A", "CD79B")
for (gene in markers) {
  merged_seurat[[paste0("count_", gene)]] <- merged_seurat[['RNA']]@counts[gene, ]
}

# Subset Seurat object to exclude cells expressing these markers
merged_seurat <- subset(merged_seurat, subset = count_PTPRC == 0 & 
                        count_EPCAM == 0 & 
                        count_HBB == 0 & 
                        count_HBA2 == 0 & 
                        count_JCHAIN == 0 & 
                        count_CD79A == 0 & 
                        count_CD79B == 0)
dim(merged_seurat)

# Variable Features -------------------------------------------------------
# Identify and save variable features
mv.features <- VariableFeatures(object = merged_seurat)
write.csv(mv.features, paste(out_dir, "mv.features.csv", sep = ""), row.names = F)

# Exclude IG genes
mv.features <- mv.features[!grepl("^IG", mv.features)]

# Data Scaling and PCA ----------------------------------------------------
# Scale data and run PCA using variable features
merged_seurat <- ScaleData(merged_seurat, assay = "RNA", features = mv.features)
merged_seurat <- RunPCA(merged_seurat, features = mv.features)

# Clustering and UMAP -----------------------------------------------------
# Find neighbors and clusters, and run UMAP
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.3)
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)

# Save UMAP plot
pdf(file = paste(out_dir, "ump.pdf", sep = ""), width = 5, height = 4, useDingbats = FALSE)
DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters", label = T)
dev.off()

# Marker Analysis ---------------------------------------------------------
# Find cluster markers
atlas.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers and create heatmaps for top genes
write.csv(atlas.markers, file = paste(out_dir, "atlas.marker.csv", sep = ""))
top10 <- atlas.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file = paste(out_dir, "heatmap_top10_markers.pdf", sep = ""), width = 12, height = 18, useDingbats = FALSE)
DoHeatmap(object = AverageExpression(merged_seurat, return.seurat = T), features = top10$gene, size = 3, disp.max = 1.4, draw.lines = F) + 
  scale_fill_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0)
dev.off()

# Dataset Representation --------------------------------------------------
# Create normalized cluster representation heatmap
table_data <- table(merged_seurat@meta.data$dataset, merged_seurat@meta.data$seurat_clusters)
long_data <- as.data.frame.matrix(table_data) %>%
  melt(id.vars = "Dataset") %>%
  rename(Dataset = Var1, Cluster = Var2, Count = value) %>%
  group_by(Dataset) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

p <- ggplot(long_data, aes(x = Dataset, y = Cluster, fill = Percentage)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Normalized Heatmap of Clusters per Dataset", x = "Dataset", y = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = paste(out_dir, "representation_datasets.pdf", sep = ""), width = 15, height = 10, useDingbats = FALSE)
print(p)
dev.off()

# Save Processed Data -----------------------------------------------------
# Save the final Seurat object
saveRDS(merged_seurat, paste(out_dir, "atlas.rds", sep = ""))


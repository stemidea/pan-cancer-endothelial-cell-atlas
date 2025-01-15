# Load necessary libraries
library(Seurat)         # For single-cell data analysis
library(harmony)        # For batch correction and integration
library(tidyverse)      # For data manipulation
library(RColorBrewer)   # For color palettes

# Define output directory
out.dir <- "./outs.predict/"
dir.create(out.dir, showWarnings = FALSE)  # Create output directory if it doesn't exist

# Load integrated Seurat object and query dataset
atlas <- readRDS("atlas.rds")  # Load integrated Seurat object (reference dataset)
new_data <- readRDS("query.rds")  # Load new dataset to be integrated

# Load marker variable features, excluding immunoglobulin genes
mv.features <- read.csv("mv.features.csv")
mv.features <- mv.features$x
mv.features <- mv.features[!grepl("^IG", mv.features)]  # Remove features starting with "IG"

# Step 1: Normalization and feature scaling for the query dataset
new_data <- NormalizeData(new_data)  # Normalize gene expression
new_data <- FindVariableFeatures(new_data)  # Identify highly variable features
new_data <- ScaleData(new_data, features = mv.features)  # Scale the specified features

# Step 2: PCA on the query dataset using the specified features
new_data <- RunPCA(new_data, features = mv.features)

# Confirm PCA dimensions
cat("Atlas PCA Dimensions: ", dim(Embeddings(atlas, reduction = "pca")), "\n")
cat("New Data PCA Dimensions: ", dim(Embeddings(new_data, reduction = "pca")), "\n")

# Step 3: Find transfer anchors between the reference and query datasets
anchors <- FindTransferAnchors(reference = atlas, query = new_data, dims = 1:30)

# Step 4: Transfer cell type annotations from the reference to the query dataset
predictions <- TransferData(
  anchorset = anchors,
  reference = atlas,
  query = new_data,
  refdata = atlas@meta.data$celltype.2,  # Cell type annotations from the reference
  dims = 1:30,
  weight.reduction = "pcaproject"
)

# Step 5: Check prediction results
print(head(predictions))  # Preview predictions
table(predictions@meta.data$predicted.id)  # Frequency of predicted cell types

# Step 6: Set predicted identities and plot violin plots for selected features
predictions@meta.data$predicted.id <- factor(
  predictions@meta.data$predicted.id, 
  levels = c("EC0:VEC.1_PLCG2", "EC1:CEC.1_ESM1", "EC2:LEC", "EC3:EEC_EGFR", "EC4:AEC", 
             "EC5:VEC.2_CLU", "EC6:MEC.1_RGS5", "EC7:UEC.1_IG", "EC8:MEC.2_DCN", 
             "EC9:UEC.2_AAtransporter", "EC10:CEC.2_FCN3", "EC11:CEC.3_CXCR4", 
             "EC12:VEC.3_FTO", "EC13:CEC.4_CRABP2", "EC14:PEC.1", "EC15:CEC.5_KLK3", 
             "EC16:CEC.6_FOSB", "EC17:PEC.2")
)

# Generate violin plots for selected genes
pdf(file = paste(out.dir, "vlnplot.pdf", sep = ""), width = 10, height = 6)
Idents(predictions) <- predictions[["predicted.id"]]
selected_genes <- c("PECAM1", "EPCAM", "EGFR", "VWF", "CLDN5", "DCN", "RGS5", 
                    "PLXDC1", "TGFB1", "ACTB", "PFKFB3", "MKI67", "GMNN", "PCNA", 
                    "PROX1", "FYVE1", "NPY1R", "NPY2R", "NPY4R", "NPY5R", 
                    "nFeature_RNA", "nCount_RNA")
lapply(selected_genes, function(gene) VlnPlot(predictions, features = gene))
dev.off()

# Step 7: Identify marker genes for each cluster
atlas.markers <- FindAllMarkers(
  predictions, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)
write.csv(atlas.markers, file = paste(out.dir, "atlas.marker.csv", sep = ""))

# Select top markers
top10 <- atlas.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- atlas.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top5 <- atlas.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Step 8: Generate heatmaps for top markers
cluster.averages <- AverageExpression(predictions, return.seurat = TRUE)

pdf(file = paste(out.dir, "heatmap_top5_markers.pdf", sep = ""), width = 6.5, height = 10, useDingbats = FALSE)
DoHeatmap(
  object = cluster.averages,
  features = top5$gene,
  size = 3,
  disp.max = 1.4,
  disp.min = -1.5
) + scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 0
)
dev.off()

pdf(file = paste(out.dir, "heatmap_top10_markers.pdf", sep = ""), width = 6.5, height = 18, useDingbats = FALSE)
DoHeatmap(
  object = cluster.averages,
  features = top10$gene,
  size = 3,
  disp.max = 1.4,
  disp.min = -1.5
) + scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 0
)
dev.off()

# Step 9: Save predictions and metadata
saveRDS(predictions, paste(out.dir, "predictions.rds", sep = ""))


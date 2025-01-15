# Suppress startup messages from libraries
suppressPackageStartupMessages({
  library(slingshot)           # For lineage inference
  library(SingleCellExperiment) # Single-cell data manipulation
  library(RColorBrewer)        # For color palettes
  library(scales)              # Scaling functions
  library(viridis)             # Viridis color maps
  library(UpSetR)              # Set visualization
  library(pheatmap)            # Heatmaps
  library(msigdbr)             # Gene sets database
  library(fgsea)               # Gene set enrichment analysis
  library(knitr)               # Reporting tools
  library(ggplot2)             # Advanced plotting
  library(gridExtra)           # Grid-based plots
  library(tradeSeq)            # Tools for analyzing Slingshot trajectories
  library(loomR)               # For handling loom files
  library(scater)              # Quality control and visualization
  library(Seurat)              # Single-cell data management
  library(tidyverse)           # Data manipulation and visualization
})

# Set working directory and create output folder
out.dir <- "./outs.slingshot/"
dir.create(out.dir, showWarnings = FALSE)

# Load the Seurat object
atlas <- readRDS("./atlas.rds")

# Extract UMAP embeddings and add them to metadata
umapCoord <- as.data.frame(Embeddings(object = atlas[["umap"]]))
atlas <- AddMetaData(atlas, metadata = umapCoord)

# Set cell type as identity
Idents(atlas) <- atlas[["celltype.2"]]

# Convert Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(atlas, assay = "RNA")

# Perform lineage inference with Slingshot
ss <- slingshot(sce, 
                reducedDim = 'UMAP', 
                clusterLabels = colData(sce)$ident, 
                start.clus = 'EC6:EndoMT.1_RGS5', 
                approx_points = 150)

# Summarize pseudotime for a specific lineage
summary(ss$slingPseudotime_1)

# Handling pseudotime values for multiple lineages
lineage_pseudotime <- lapply(1:10, function(i) {
  pseudotime_var <- paste0("slingPseudotime_", i)
  temp <- ss[[pseudotime_var]]
  temp[is.na(temp)] <- 0
  return(temp)
})

# Combine pseudotime values across lineages
ss$slingPseudotime <- apply(do.call(cbind, lineage_pseudotime), 1, max)

# Summarize combined pseudotime
summary(ss$slingPseudotime)

# Create pseudotime trajectory plots
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(ss$slingPseudotime_1, breaks = 100)]

# Save lineage trajectory plot
pdf(file.path(out.dir, "lineages.pdf"), width = 5, height = 5, useDingbats = FALSE)
plot(reducedDims(ss)$UMAP, 
     col = plotcol, 
     pch = 16, 
     asp = 1, 
     xlim = c(-15, 10), 
     ylim = c(-10, 10))
lines(SlingshotDataSet(ss), lwd = 2, col = 'black')
dev.off()

# Generate UMAP plot with pseudotime
pdf(file.path(out.dir, "trajectory.pdf"), width = 8, height = 5, useDingbats = FALSE)
df <- as.data.frame(reducedDims(ss)$UMAP)
df <- data.frame(df, cluster = colData(sce)$ident)
ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = cluster), shape = 21, size = 3) +
  theme_classic() +
  labs(title = "UMAP Trajectory") +
  xlim(-10, 10) + ylim(-8, 8)
dev.off()

# Plot combined pseudotime using ggplot2
pdf(file.path(out.dir, "time_combined.pdf"), width = 8, height = 5, useDingbats = FALSE)
df <- as.data.frame(reducedDims(ss)$UMAP)
df <- data.frame(df, st = ss$slingPseudotime)
ggplot(data = df, aes(x = UMAP_1, y = UMAP_2, color = st)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Combined Pseudotime") +
  xlim(-10, 10) + ylim(-10, 10)
dev.off()

# Generate pseudotime plot for each lineage
for (i in 1:10) {
  pseudotime_var <- paste0("slingPseudotime_", i)
  if (pseudotime_var %in% colnames(colData(ss))) {
    df <- as.data.frame(reducedDims(ss)$UMAP)
    df$st <- colData(ss)[[pseudotime_var]]
    
    pdf(file.path(out.dir, paste0("time_lineage_", i, ".pdf")), width = 5.5, height = 5, useDingbats = FALSE)
    ggplot(data = df, aes(x = UMAP_1, y = UMAP_2, color = st)) +
      geom_point(size = 3) +
      scale_color_viridis(option = "D") +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(title = paste("Pseudotime for Lineage", i)) +
      xlim(-10, 10) + ylim(-10, 10)
    dev.off()
  }
}

# Save combined pseudotime plot with lineages highlighted
pdf(file.path(out.dir, "lineages_highlighted.pdf"), width = 5, height = 5, useDingbats = FALSE)
plot(reducedDims(ss)$UMAP, 
     asp = 1, 
     pch = 16, 
     xlab = "UMAP-1", 
     ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(ss$slingPseudotime, breaks = 100)], 
     xlim = c(-8, 8), ylim = c(-8, 8), 
     xaxt = "n", yaxt = "n")
axis(1, at = c(-10, -5, 0, 5, 10))
axis(2, at = c(-10, -5, 0, 5, 10))
lines(SlingshotDataSet(ss), type = 'lineages', show.constraints = TRUE)
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), file.path(out.dir, "sessionInfo.txt"))


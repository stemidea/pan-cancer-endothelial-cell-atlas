# Load necessary libraries
library(CellChat)  # For cell-cell communication analysis
library(patchwork) # For combining ggplot-based plots
library(Seurat)    # For single-cell data handling
library(NMF)       # For matrix factorization
library(ggalluvial) # For alluvial diagrams
options(stringsAsFactors = FALSE) # Prevent automatic conversion of strings to factors

# Set working directory and create output folder
setwd("./")
wd <- "./"
out_dir <- "./outs.cellchat/"  # Output folder for results
dir.create(out_dir, showWarnings = FALSE)  # Create output directory if it doesn't exist

# Load input data (Seurat object)
atlas <- readRDS(file = "atlas.rds")  # Load pre-processed single-cell data
meta <- atlas@meta.data  # Extract metadata for further processing

# Display unique identifiers in the metadata
unique(atlas@meta.data$orig.ident)

## OPTIONAL: Create a subset of the data
# atlas <- subset(atlas, subset = group.1 != "Normal")

# Define a variable of interest to use as cell identity
unique(atlas@meta.data$celltype.2)
atlas[["celltype.2"]] <- factor(
  atlas@meta.data$celltype.2,
  levels = c("Epithelial_cells", "Monocyte", "Macrophage", "Neutrophils", "DC", 
             "B_cell", "T_cells",  "NK_cell", "Smooth_muscle_cells", "EC0:VEC.1_PLCG2", 
             "EC1:CEC.1_ESM1", "EC2:LEC", "EC3:Mixture", "EC4:AEC", "EC5:VEC.2_CLU", 
             "EC6:EndoMT.1_RGS5", "EC7:UEC.1_IG", "EC8:EndoMT.2_DCN", 
             "EC9:UEC.2_AAtransporter", "EC10:CEC.2_FCN3", "EC11:CEC.3_CXCR4", 
             "EC12:VEC.3_FTO", "EC13:CEC.4_CRABP2", "EC14:PEC.1", 
             "EC15:CEC.5_KLK3", "EC17:PEC.2")
)

# Set cell identities for downstream analysis
Idents(atlas) <- atlas@meta.data$celltype.2

# Set default assay
DefaultAssay(atlas) <- "RNA"

# Initialize CellChat object
cellchat <- createCellChat(object = atlas, group.by = "ident")  # Group by cell identity
groupSize <- as.numeric(table(cellchat@idents))  # Store group sizes for visualization

# Load human CellChat database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB  # Assign database to CellChat object

# Preprocess data for CellChat
cellchat <- subsetData(cellchat)  # Subset data to remove uninformative genes
cellchat <- identifyOverExpressedGenes(cellchat)  # Identify overexpressed genes
cellchat <- identifyOverExpressedInteractions(cellchat)  # Identify ligand-receptor pairs
cellchat <- projectData(cellchat, PPI.human)  # Map data to protein-protein interactions
cellchat <- computeCommunProb(cellchat)  # Compute communication probabilities
cellchat <- computeCommunProbPathway(cellchat)  # Compute pathway-level probabilities
cellchat <- aggregateNet(cellchat)  # Aggregate the interaction network
cellchat <- netAnalysis_computeCentrality(cellchat)  # Compute network centrality

# Generate interaction plots
pdf(file = paste0(out_dir, "circleplot_ni.pdf"), width = 6, height = 6, useDingbats = FALSE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

# Visualize individual pathways in a circular layout
pdf(file = paste0(out_dir, "individual_circle.pdf"), width = 15, height = 15, useDingbats = FALSE)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Aggregate plots for specific pathways (e.g., VEGF)
pathways.show <- c("VEGF")
vertex.receiver <- seq(1, 4)

# Hierarchical layout
pdf(file = paste0(out_dir, "aggregate_hierachy.pdf"), width = 18, height = 10, useDingbats = FALSE)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
dev.off()

# Circle plot
pdf(file = paste0(out_dir, "aggregate_circle.pdf"), width = 8, height = 8, useDingbats = FALSE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
pdf(file = paste0(out_dir, "aggregate_chord.pdf"), width = 8, height = 8, useDingbats = FALSE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

# Heatmap of signaling
pdf(file = paste0(out_dir, "aggregate_heatmap.pdf"), width = 8, height = 8, useDingbats = FALSE)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

# Contribution of individual ligand-receptor pairs
pdf(file = paste0(out_dir, "net_contribution.pdf"), width = 8, height = 8, useDingbats = FALSE)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

# Save CellChat object and session info
saveRDS(cellchat, file = paste0(out_dir, "cellchat.rds"))
writeLines(capture.output(sessionInfo()), paste0(out_dir, "sessionInfo.txt"))


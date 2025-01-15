# Load necessary libraries
library(Seurat)  # For single-cell data analysis and processing
library(monocle) # For pseudotime trajectory inference

# Create output directory for Monocle2 results
out.dir <- "./outs.monocle2/"
dir.create(out.dir, showWarnings = FALSE)

# Load Seurat object
OBJ <- readRDS("atlas.rds")

# Check metadata column names to explore available annotations
colnames(OBJ@meta.data)

# Set 'subtype' metadata to 'tissue' column
OBJ[["subtype"]] <- OBJ[["tissue"]]

# Convert Seurat object counts to a sparse matrix for Monocle
data <- as(as.matrix(OBJ@assays$RNA@counts), 'sparseMatrix')

# Prepare phenotype (cell metadata) and feature data (gene metadata)
pd <- new('AnnotatedDataFrame', data = OBJ@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# Create Monocle2 CellDataSet
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Detect genes expressed in at least 10 cells
HSMM <- detectGenes(HSMM, min_expr = 0.5)
detect_genes <- fData(HSMM)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

# Perform differential gene testing to find ordering genes
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~subtype + nCount_RNA + percent.mt", 
                                      reducedModelFormulaStr = "~nCount_RNA + percent.mt", 
                                      cores = 1)

# Select ordering genes based on significant p-values and cell expression
ordering_genes <- row.names(subset(diff_test_res, 
                                   pval < 1e-10 & qval < 1e-10 & 
                                   num_cells_expressed > round(ncol(data) * 0.1)))
HSMM <- setOrderingFilter(HSMM, ordering_genes)

# Reduce dimensionality and order cells in pseudotime
HSMM <- reduceDimension(HSMM, max_components = 2)
HSMM <- orderCells(HSMM)

# Save trajectory plots
pdf(file.path(out.dir, "ordering1.pdf"), 6, 6, useDingbats = FALSE)
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "subtype")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "Pseudotime") + scale_color_viridis_c()
dev.off()

# Save Monocle object
saveRDS(HSMM, file = file.path(out.dir, "monocle.rds"))

# Define root state based on a specific subtype
GM_state <- function(cds) {
  if (length(unique(pData(cds)$State)) > 1) {
    T0_counts <- table(pData(cds)$State, pData(cds)$subtype)[, "Normal"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return(1)
  }
}
HSMM <- orderCells(HSMM, root_state = GM_state(HSMM))

# Save rooted trajectory plots
pdf(file.path(out.dir, "rooted_ordering2.pdf"), 6, 6, useDingbats = FALSE)
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "subtype")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "Pseudotime") + scale_color_viridis_c()
dev.off()

# Perform branch analysis
branch_points <- c(1, 2, 3)  # Define branch points to analyze
for (bp in branch_points) {
  BEAM_res <- BEAM(HSMM, branch_point = bp, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$pval, BEAM_res$qval), ]
  BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval", "num_cells_expressed")]
  
  # Remove mitochondrial, ribosomal protein genes
  BEAM_res <- BEAM_res[grep("(^MT-)|(^RPL)|(^RPS)", rownames(BEAM_res), invert = TRUE), ]
  
  # Select genes for heatmap
  heatgenes <- row.names(subset(BEAM_res, pval < 1e-10 & qval < 1e-10 & num_cells_expressed > round(ncol(data) * 0.1)))
  pall <- plot_genes_branched_heatmap(HSMM[heatgenes, ],
                                      branch_point = bp,
                                      num_clusters = 3,
                                      cores = 1,
                                      use_gene_short_name = TRUE,
                                      show_rownames = TRUE,
                                      return_heatmap = TRUE)
  
  # Save branch-specific heatmap and gene information
  pdf(file.path(out.dir, paste0("branch", bp, "-v2.pdf")), 6, 6, useDingbats = FALSE)
  print(pall$ph_res)
  dev.off()
  
  geneclass <- data.frame(Gene = rownames(pall$BranchA_exprs),
                          annotation_row = pall$annotation_row,
                          pval = BEAM_res[rownames(pall$BranchA_exprs), "pval"],
                          qval = BEAM_res[rownames(pall$BranchA_exprs), "qval"],
                          num_cells_expressed = BEAM_res[rownames(pall$BranchA_exprs), "num_cells_expressed"])
  write.table(geneclass, file.path(out.dir, paste0("branched", bp, ".geneclass.txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Save session information
writeLines(capture.output(sessionInfo()), file.path(out.dir, "sessionInfo.txt"))


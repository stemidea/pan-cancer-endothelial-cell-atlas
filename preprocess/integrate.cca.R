# Load necessary libraries
library(Seurat)  # For single-cell data analysis and integration
library(dplyr)   # For data manipulation

# Set the working directory
setwd("./")  # Sets the working directory to the current directory

# Step 1: Read in individual Seurat objects from the "./seurat" folder
my.seurat <- list.files(path = "./seurat")  # List all files in the "seurat" directory
atlas.list <- list()  # Initialize an empty list to store Seurat objects

for (i in 1:length(my.seurat)) {  # Loop through each file in the directory
  atlas.list[[i]] <- readRDS(file = paste("./seurat/", my.seurat[i], sep = ""))  
  # Load each Seurat object and store it in the list
}

# Step 2: Pre-process each Seurat object in the list
atlas.list <- lapply(X = atlas.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)  # Normalize gene expression values
  x <- FindVariableFeatures(x, verbose = FALSE)  # Identify variable features
  return(x)  # Return the processed object
})

# Step 3: Identify features for integration
features <- SelectIntegrationFeatures(object.list = atlas.list)  
# Select features shared across datasets to use for integration

# Step 4: Scale data and run PCA for each dataset
atlas.list <- lapply(X = atlas.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)  # Scale selected features
  x <- RunPCA(x, features = features, verbose = FALSE)  # Perform PCA for dimensionality reduction
  return(x)  # Return the processed object
})

# Step 5: Identify anchors for integration
atlas.anchors <- FindIntegrationAnchors(
  object.list = atlas.list,        # List of pre-processed Seurat objects
  anchor.features = features,     # Features selected for integration
  dims = 1:50                     # Number of dimensions to use for integration
)

# Optional alternative: Use RPCA-based integration for large datasets
# anchors <- FindIntegrationAnchors(
#   object.list = atlas.list,
#   reference = c(1, 2),          # Specify reference datasets for integration
#   reduction = "rpca",           # Use reciprocal PCA for anchor identification
#   dims = 1:50
# )

# Step 6: Integrate the datasets
atlas <- IntegrateData(
  anchorset = atlas.anchors,  # Anchor set identified in the previous step
  dims = 1:50                 # Number of dimensions to use for integration
)

# Step 7: Save the integrated Seurat object
saveRDS(atlas, file = "./atlas.combined.rds")  
# Save the integrated Seurat object to an RDS file for future use


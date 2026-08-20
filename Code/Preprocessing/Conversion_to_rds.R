# Purpose: Converting the H5AD files to RDS format.

# Loading required libraries
library(Seurat)
library(reticulate)
library(anndata)

# Load the h5ad file from scvi, scArches, or Scanorama
dataset <- "" # Specify the datatset
data <- read_h5ad(paste0("./",dataset,"_X_corrected.h5ad")) 

# Extract the corrected gene expression matrix and create a Seurat object
expr <- as(as.matrix(data$X), "sparseMatrix")
colnames(expr) <- data$var_names
rownames(expr) <- data$obs_names
data <- CreateSeuratObject(counts = t(expr), meta.data = data$obs)

# Saving the data
saveRDS(data, paste0("./",dataset,"_X_corrected.RDS"))

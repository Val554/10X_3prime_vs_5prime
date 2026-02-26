# Purpose: Apply limma batch correction within each subjects/individuals and save result.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(limma)

# Read the raw preprocessed data for each dataset (Seurat object)
data <- readRDS(paste0(getwd(),"/prep_FINAL.rds"))

# Split the data by patient
seurat_list <- SplitObject(data, split.by = "new_id")

# Perform normalization for all subjects/individuals
for (i in seq_along(seurat_list)) {
  seurat_obj <- seurat_list[[i]]
  
  # Normalize and Extract the normalized count matrix
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  counts <- GetAssayData(seurat_obj, layer = "data") # Log-normalized counts
  
  # Get batch information (e.g., 5' vs 3') from metadata
  batch <- seurat_obj@meta.data$assay
  
  start.time <- Sys.time()  
  # Run limma 
  corrected_counts <- removeBatchEffect(
    x = as.matrix(counts),
    batch = batch
  )
  end.time <- Sys.time()
  print(end.time - start.time)
  
  # Converting the dense matrix to sparse
  threshold <- 1e-5  # Adjust based on inspection
  corrected_counts[corrected_counts < threshold & corrected_counts > -threshold] <- 0
  corrected_counts <- Matrix(corrected_counts, sparse = TRUE)
  
  # Update the corrected counts back to the Seurat object
  seurat_list[[i]] <- SetAssayData(seurat_obj, layer = "data", new.data = corrected_counts)
}

# Merge the Seurat objects into one
data <- merge(seurat_list[[1]], y = seurat_list[-1])

# Saving the final object for downstream analysis
saveRDS(data, paste0(getwd(),"/limma_corrected.RDS"))

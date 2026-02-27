# Purpose: Apply mnnCorrect batch correction within each patient/individual in parallel and save result.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(batchelor)
library(Matrix)
library(foreach)
library(doParallel)

# Read the raw preprocessed data for each dataset (Seurat object)
data <- readRDS(paste0(getwd(),"/prep_FINAL.rds"))

# Split the data by patient
seurat_list <- SplitObject(data, split.by = "new_id")

# Register parallel backend with a suitable number of cores
cl <- makeCluster(2)  # Adjust the number of cores as needed
registerDoParallel(cl)

# Parallelized processing
seurat_list <- foreach(i = seq_along(seurat_list), .combine = c, .packages = c("Seurat", "batchelor", "Matrix")) %dopar% {
  # Extract the Seurat object for the current patient
  seurat_obj <- seurat_list[[i]]
  
  # Normalize and extract the normalized count matrix
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  
  counts <- as.matrix(GetAssayData(seurat_obj, layer = "data"))  # Log-normalized counts
  
  # Get batch information (e.g., 5' vs 3') from metadata
  batch <- seurat_obj@meta.data$assay
  # Get highly variable features
  hvf <- VariableFeatures(seurat_obj)
  
  start <- Sys.time()  
  # Run mnnCorrect
  mnn_result <- mnnCorrect(
    counts,
    batch = batch,
    subset.row = hvf,
    correct.all = TRUE
  )
  end <- Sys.time()
  print(end - start)
  
  # Extract corrected counts
  corrected_counts <- assay(mnn_result, "corrected")
  
  # Convert dense matrix to sparse matrix
  threshold <- 1e-5  # Adjust based on inspection
  corrected_counts[corrected_counts < threshold & corrected_counts > -threshold] <- 0
  corrected_counts <- Matrix(corrected_counts, sparse = TRUE)
  
  # Update the corrected counts back to the Seurat object
  SetAssayData(seurat_obj, layer = "data", new.data = corrected_counts)
}

# Stop the cluster
stopCluster(cl)

# Merge the Seurat objects into one
data <- merge(seurat_list[[1]], y = seurat_list[-1])

# Saving the final object for downstream analysis
saveRDS(data, paste0(getwd(),"/mnnCorrect_corrected.RDS"))

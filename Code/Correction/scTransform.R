# Purpose: Apply scTransform within each subjects/individuals and save result.

# Loading required libraries
library(Seurat)
library(patchwork)
library(sctransform)

# Read the raw preprocessed data for each dataset (Seurat object)
data <- readRDS(paste0(getwd(),"/prep_FINAL.rds"))

# Split the data by patient
seurat_list <- SplitObject(data, split.by = "new_id")

# Load the highly variable features (HVF) list, which is the same for all patients and was generated in the previous step. 
# OR load the HVF list with Bias Genes
# OR any list of genes you choose to evaluate
hvf <- readRDS(paste0(getwd(),"/hvf.rds"))

# Perform normalization for all subjects/individuals
for (i in seq_along(seurat_list)) {
  start <- Sys.time()
  seurat_list[[i]] <- SCTransform(seurat_list[[i]], vars.to.regress = "assay", verbose = TRUE, residual.features = hvf)
  end <- Sys.time()
  print(end - start)
}

# Merge the Seurat objects into one
data <- merge(seurat_list[[1]], y = seurat_list[-1])

# Saving the final object for downstream analysis
saveRDS(data, paste0(getwd(),"/scTransformRegression_corrected.RDS"))

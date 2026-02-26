# Purpose: Apply M3Drop batch correction within each patient/individual and save result.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(M3Drop)

# Read the raw preprocessed data for each dataset (Seurat object)
data <- readRDS(paste0(getwd(),"/prep_FINAL.rds"))

# Subset the 3' and 5' data and split it by patient
m_tr_3 <- subset(data, subset = assay == "10x 3' v2")
patient_list_3 <- SplitObject(m_tr_3, split.by = "new_id")
m_tr_5 <- subset(data, subset = assay == "10x 5' v1")
patient_list_5 <- SplitObject(m_tr_5, split.by = "new_id")

# Create the M3Drop normalization function
norm <- function(object) {
  for (object_name in names(object)) {
    print(paste("Loop for ", object_name))
    print(object[[object_name]])
    object_temp <- object[[object_name]]
    
    # Extract the normalized count matrix
    expression_matrix <- as.matrix(GetAssayData(object_temp, layer = "counts"))
    start <- Sys.time()    
    
    # Converting the expression matrix data into M3Drop format, as per developers' instructions
    count_mat <- NBumiConvertData(expression_matrix, is.counts=TRUE)
    
    # Get the residuals using the NBumi model, which is the default for UMI data. This will return a matrix of Pearson residuals.
    pearson <- NBumiPearsonResiduals(count_mat, fits = NULL)
    end <- Sys.time()
    print(end - start)
    
    # Set the normalized expression matrix back to the Seurat object
    object[[object_name]] <- CreateSeuratObject(pearson, meta.data = object_temp@meta.data)
  }
  return(object)
}

# Perform M3Drop normalization for each subject's 3' and 5' data and merge the results
m_tr_3 <- norm(patient_list_3)
m_tr_3 <- merge(m_tr_3[[1]], y = m_tr_3[-1])

m_tr_5 <- norm(patient_list_5)
m_tr_5 <- merge(m_tr_5[[1]], y = m_tr_5[-1])

# Merge the Seurat objects into one
data <- merge(m_tr_3, m_tr_5)

# Saving the final object for downstream analysis
saveRDS(data, paste0(getwd(),"/M3Drop_corrected.RDS"))

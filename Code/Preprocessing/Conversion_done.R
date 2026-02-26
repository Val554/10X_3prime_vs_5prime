# Purpose: Converting the RDS files to h5ad format and merged datasets for training and testing of the scArches models.

# Loading required libraries
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(zellkonverter)

# Read the raw preprocessed data for each dataset (Seurat object)
data1 <- readRDS("DS1_prep_FINAL.rds")
data2 <- readRDS("DS2_prep_FINAL.rds") 
data3 <- readRDS("DS3_prep_FINAL.rds") 
data4 <- readRDS("DS4_prep_FINAL.rds") 
data5 <- readRDS("DS5_prep_FINAL.rds") 
data6 <- readRDS("DS6_prep_FINAL.rds") 

# Clean up the meta data:
keep <- c("nCount_RNA", "nFeature_RNA", "new_id", "cell_type_ontology_term_id", "assay", "disease", "tissue")
data1@meta.data <- data1@meta.data[, keep]
data1@meta.data$assay <- ifelse(data1@meta.data$assay == "10x 3' v2", "3X", "5X")

data2@meta.data <- data2@meta.data[, keep]
data2@meta.data$assay <- ifelse(data2@meta.data$assay == "10x 3' v2", "3X", "5X")

data3@meta.data <- data3@meta.data[, keep]
data3@meta.data$assay <- ifelse(data3@meta.data$assay == "10x 3' v2", "3X", "5X")

data4@meta.data <- data4@meta.data[, keep]
data4@meta.data$assay <- ifelse(data4@meta.data$assay == "10x 3' v3", "3X", "5X")

data5@meta.data <- data5@meta.data[, keep]
data5@meta.data$assay <- ifelse(data5@meta.data$assay == "10x 3' v2", "3X", "5X")

data6@meta.data <- data6@meta.data[, keep]
data6@meta.data$assay <- ifelse(data6@meta.data$assay == "10x 3' v2", "3X", "5X")


# Merge different combinations of the datasets. 
# The merged object will be used for training, while the individual datasets will be used for testing and validation.
merged_1 <- merge(data2, y = c(data3, data4, data5, data6), add.cell.ids = c("DS2", "DS3", "DS4", "DS5", "DS6"))
merged_1 <- JoinLayers(merged_1)
merged_2 <- merge(data1, y = c(data3, data4, data5, data6), add.cell.ids = c("DS1", "DS3", "DS4", "DS5", "DS6"))
merged_2 <- JoinLayers(merged_2)
merged_3 <- merge(data1, y = c(data2, data4, data5, data6), add.cell.ids = c("DS1", "DS2", "DS4", "DS5", "DS6"))
merged_3 <- JoinLayers(merged_3)
merged_5 <- merge(data1, y = c(data2, data3, data4, data6), add.cell.ids = c("DS1", "DS2", "DS3", "DS4", "DS6"))
merged_5 <- JoinLayers(merged_5)

# Save the training datasets in h5ad format 
# An example for one dataset, repeat for the rest.
sce <- as.SingleCellExperiment(merged_1, assay = "RNA")
writeH5AD(sce, "merged_1_prep_FINAL.h5ad")

# Save all the datasets in h5ad format independently for testing
# An example for one dataset, repeat for the rest.
data1_sub <- subset(data1, subset = assay == "3X") # Used for testing
data1_sub_5 <- subset(data1, subset = assay == "5X") # Actual 5' used for comparison

sce <- as.SingleCellExperiment(data1_sub, assay = "RNA")
writeH5AD(sce, "data1_3_prep_FINAL.h5ad")
saveRDS(data1_sub_5, "./data1_5_prep_FINAL.rds")



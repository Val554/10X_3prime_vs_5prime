# Purpose: Scale (where needed) the corrected (or log transformed) data, run UMAP, calculate mixing metrics and Jaccard Index/Overlap Score for each subject/individual and save results.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(presto)

# Load the highly variable features (HVF) list, which is the same for all patients and was generated in the previous step. 
hvf <- readRDS(paste0(getwd(),"/hvf.rds"))

# Create the Jaccard Index/ Overlap Score function
jaccard_index <- function(set1, set2) {
  if (length(set1) == 0 & length(set2) == 0) return(NA)  # Handle empty sets
  length(intersect(set1, set2)) / min(length(set1), length(set2)) # Overlap score, min of the lengths of two sets
}

# Loading the corrected object after a chosen normalization/batch correction - X.
obj <- readRDS(paste0(getwd(),"/X_corrected.RDS"))

VariableFeatures(obj) <- hvf

#--- Additional preprocessing step for the uncorrected (log normalized) data.
obj <- split(obj, f = obj@meta.data$new_id)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, nfeatures = 3000) # does it on the raw counts
obj <- JoinLayers(obj)
#---

#--- Alternative to ScaleData for M3Drop instead due to their internal scaling
obj <- SetAssayData(obj, layer = "scale.data", new.data = as.matrix(GetAssayData(obj, layer = "counts")))
obj <- SetAssayData(obj, layer = "data", new.data = as.matrix(GetAssayData(obj, layer = "counts")))
#---

# Split the object by subject/individual
obj <- SplitObject(obj, split.by = "new_id") 

# Setup the variables to store results
jaccard_scores <- list()
find_var_f <- list()
mixing_metrics_assay <- numeric(length(names(obj)))
names(mixing_metrics_assay) <- names(obj)
mixing_metrics_cell <- numeric(length(names(obj)))
names(mixing_metrics_cell) <- names(obj)
raw_markers <- list()

# Run the set of function for each subject/individual
for (i in names(obj)) {
  print(i)
  
  # a] Standard UMAP 
  obj[[i]] <- ScaleData(obj[[i]]) # ComBat, fastMNN, limma, mnnCorrect, Z-transformation; NOT for M3Drop or scTransform (SCT assay)
  obj[[i]] <- RunPCA(obj[[i]], npcs = 40)
  obj[[i]] <- RunUMAP(obj[[i]], dims = 1:20, reduction = "pca")
  
  # b] Mixing Metrics
  # By assay (3' vs 5')
  mix_assay <- MixingMetric(
    obj[[i]],
    grouping.var = "assay",
    reduction = "pca",
    dims = 1:20,
    k = 5,
    max.k = 150,
    eps = 0,
    verbose = TRUE
  )
  
  # By cell type
  mix_cell <- MixingMetric(
    obj[[i]],
    grouping.var = "cell_type_ontology_term_id",
    reduction = "pca",
    dims = 1:20,
    k = 5,
    max.k = 150,
    eps = 0,
    verbose = TRUE
  )
  
  # Storing the mixing metric results, we take 150 - metric to have higher values indicate better mixing
  mixing_metrics_assay[[i]] <- median(150 - mix_assay)
  mixing_metrics_cell[[i]] <- median(150 - mix_cell)
  
  # c] Jaccard Index / Overlap Score calculations
  
  Idents(obj[[i]]) <- obj[[i]]$cell_type_ontology_term_id
  
  # Find markers for 3' and 5' assays separately. The `fc.slot` = "scale.data" ensures the effect size is calculated as the average difference.
  mark_3 <- FindAllMarkers(subset(obj[[i]], subset = assay == "10x 3' v2"), # Set the proper version of the 3' assay
                           group.by = "cell_type_ontology_term_id",
                           logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5, slot = "data", fc.slot = "scale.data")
  
  mark_5 <- FindAllMarkers(subset(obj[[i]], subset = assay == "10x 5' v1"), # Set the proper version of the 5' assay
                           group.by = "cell_type_ontology_term_id", 
                           logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5, slot = "data", fc.slot = "scale.data")
  
  # Subset the significant genes
  sig_genes_3 <- subset(mark_3, p_val_adj < 0.05)
  sig_genes_5 <- subset(mark_5, p_val_adj < 0.05)
  
  # Summary of the number of significant genes per cluster for 3' and 5' assays, stored in the `find_var_f` list
  find_var_f[[paste0(i, "_3")]] <- sig_genes_3 %>%
    group_by(cluster) %>%
    summarise(count = n()) %>%
    arrange(count)
  find_var_f[[paste0(i, "_5")]] <- sig_genes_5 %>%
    group_by(cluster) %>%
    summarise(count = n()) %>%
    arrange(count)
  
  # Calculate the Jaccard Index / Overlap Score for each cluster between the 3' and 5' assays. 
  # We will store these scores in a named vector where names are the cell clusters.
  score <- numeric(length(intersect(sig_genes_3$cluster, sig_genes_5$cluster)))
  names(score) <- intersect(sig_genes_3$cluster, sig_genes_5$cluster)
  for (cell in intersect(sig_genes_3$cluster, sig_genes_5$cluster)) {
    score[[cell]] <- jaccard_index(sig_genes_3$gene[sig_genes_3$cluster == cell],
                                   sig_genes_5$gene[sig_genes_5$cluster == cell])
  }
  
  # Store the raw markers for potential downstream analysis
  raw_markers[[paste0(i, "_3")]] <- mark_3
  raw_markers[[paste0(i, "_5")]] <- mark_5
  
  # Convert to a named vector
  jaccard_scores[[i]] <- score
}

saveRDS(find_var_f, paste0(getwd(),"/X_FAF.rds"))
saveRDS(raw_markers, paste0(getwd(),"/X_raw_markers.rds"))
saveRDS(mixing_metrics_assay, paste0(getwd(),"/X_MM_assay.rds"))
saveRDS(mixing_metrics_cell, paste0(getwd(),"/X_MM_cell.rds"))
saveRDS(obj, paste0(getwd(),"/Scaled_X_corrected.rds"))

# Summarize the Jaccard/Overlap scores by calculating the average score across clusters 
# for each subject/individual and save the result in a data frame.
df <- map_dfr(jaccard_scores, enframe, .id = "Sample")
df <- df %>%
  group_by(name) %>%
  summarise(Average_Jaccard_Score = mean(value, na.rm = TRUE))
df$norm <- "X"
saveRDS(df, paste0(getwd(),"/X_JS_sum.rds"))


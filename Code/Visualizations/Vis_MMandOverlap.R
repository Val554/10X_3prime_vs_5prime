# Purpose: Summarizing the Mixing Metrics and Jaccard/Overlap Scores for the comparison of the 3' and 5' assays for each batch correction method. 
# This is done for further visualizations and comparisons across the different techniques. 

# Loading required libraries
library(reshape2) 
library(ggplot2)
library(dplyr)

# Repeat each of these steps for each dataset

# ----- 1] Mixing Metrics ----

## ---- a] by assay ----

tis <- "" # Specify the tissue for dataset 2
LN <- readRDS(paste0("./LogNorm_MM_assay",tis,".rds"))
MM_df <- data.frame(Individual = names(LN), Log_Normalize = LN)
CB <- readRDS(paste0("./ComBat_MM_assay",tis,".rds"))
MM_df$ComBat <- CB
mnn <- readRDS(paste0("./mnn_MM_assay",tis,".rds"))
MM_df$mnnCorrect <- mnn
fast <- readRDS(paste0("./fast_MM_assay",tis,".rds"))
MM_df$fastMNN <- fast
scanorama <- readRDS(paste0("./scanorama_MM_assay",tis,".rds"))
MM_df$Scanorama <- scanorama
limma <- readRDS(paste0("./limma_MM_assay",tis,".rds"))
MM_df$limma <- limma
Z_t <- readRDS(paste0("./Z_MM_assay",tis,".rds"))
MM_df$Z_transform <- Z_t
scvi <- readRDS(paste0("./scvi_MM_assay",tis,".rds"))
MM_df$SCVI <- scvi
sca <- readRDS(paste0("./scArches_MM_assay",tis,".rds"))
MM_df$scArches <- sca
SC <- readRDS(paste0("./SC_MM_assay",tis,".rds"))
MM_df$scTransform <- SC
SC_s <- readRDS(paste0("./SC_split_MM_assay",tis,".rds"))
MM_df$scTransform_split <- SC_s
M3 <- readRDS(paste0("./m3drop_MM_assay",tis,".rds"))
MM_df$M3Drop <- M3

# Reshape data for plotting
melted_data <- reshape2::melt(MM_df, id.vars = "Individual")
saveRDS(melted_data, paste0("./FINAL_MM_assay",tis,".rds"))

## ---- b] by cell type ----

tis <- ""
LN <- readRDS(paste0("./LogNorm_MM_cell",tis,".rds"))
MM_df <- data.frame(Individual = names(LN), Log_Normalize = LN)
CB <- readRDS(paste0("./ComBat_MM_cell",tis,".rds"))
MM_df$ComBat <- CB
mnn <- readRDS(paste0("./mnn_MM_cell",tis,".rds"))
MM_df$mnnCorrect <- mnn
fast <- readRDS(paste0("./fast_MM_cell",tis,".rds"))
MM_df$fastMNN <- fast
scanorama <- readRDS(paste0("./scanorama_MM_cell",tis,".rds"))
MM_df$Scanorama <- scanorama
limma <- readRDS(paste0("./limma_MM_cell",tis,".rds"))
MM_df$limma <- limma
Z_t <- readRDS(paste0("./Z_MM_cell",tis,".rds"))
MM_df$Z_transform <- Z_t
scvi <- readRDS(paste0("./scvi_MM_cell",tis,".rds"))
MM_df$SCVI <- scvi
sca <- readRDS(paste0("./scArches_MM_cell",tis,".rds"))
MM_df$scArches <- sca
SC <- readRDS(paste0("./SC_MM_cell",tis,".rds"))
MM_df$scTransform <- SC
SC_s <- readRDS(paste0("./SC_split_MM_cell",tis,".rds"))
MM_df$scTransform_split <- SC_s
M3 <- readRDS(paste0("./m3drop_MM_cell",tis,".rds"))
MM_df$M3Drop <- M3

# Reshape data for plotting
melted_data <- reshape2::melt(MM_df, id.vars = "Individual")
saveRDS(melted_data, paste0("./FINAL_MM_cell",tis,".rds"))

# ----- 2] Overlap Scores ----

tis <- "" # Specify the tissue for dataset 2
LN <- readRDS(paste0("./LogNorm_JS_sum",tis,".rds"))
JS_df <- LN
CB <- readRDS(paste0("./ComBat_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, CB)
mnn <- readRDS(paste0("./mnn_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, mnn)
Fast <- readRDS(paste0("./fast_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, Fast)
scanorama <- readRDS(paste0("./scanorama_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, scanorama)
limma <- readRDS(paste0("./limma_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, limma)
Z_t <- readRDS(paste0("./Z_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, Z_t)
scvi <- readRDS(paste0("./scvi_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, scvi)
sca <- readRDS(paste0("./scArches_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, sca)
SC <- readRDS(paste0("./SC_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, SC)
SC_s <- readRDS(paste0("./SC_split_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, SC_s)
M3 <- readRDS(paste0("./m3drop_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, M3)

saveRDS(JS_df, paste0("./FINAL_JS_sum",tis,".rds"))


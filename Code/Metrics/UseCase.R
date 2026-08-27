# Purpose: Evaluating the performance of different batch correction techniques in terms of their ability to recover the true DEGs between two cell types with partial overlap across the 3' and 5' assays.
# Additionally, evaluating the general differences between the 3' and 5' assays 

# Loading required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(presto)
library(MAST)

# ------ Partial cell type overlap ------

## ------- Loading the data --------
### --- DS1 ---
data <- readRDS("./prep_FINAL_v5.rds")
data <- subset(data, new_id %in% c("F29_45P", "F30_45P", "F38_45P",  "P1_CD3P")) # "F41_45P", "F45_45P" - no Tregs
data@meta.data <- droplevels(data@meta.data)

cell_type1 <- "regulatory T cell"
cell_type2 <- "CD8-positive, alpha-beta T cell"
###

### --- DS6  ---
data <- readRDS("./prep_FINAL.rds")

cell_type1 <- "dendritic cell"
cell_type2 <- "B cell"
###

### --- DS2lv ---
data <- readRDS("./prep_FINAL_v5.rds")
data <- subset(data, new_id %in% c("F32_CD45P_liver", "F34_CD45P_liver", "F38_CD45P_liver",  "F41_CD45P_liver", "F45_CD45P_liver")) # CD45N not included
data@meta.data <- droplevels(data@meta.data)

cell_type1 <- "dendritic cell"
cell_type2 <- "natural killer cell"
###

### Next Step
assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

data_3 <- subset(data, subset = assay == assay_3)
data_5 <- subset(data, subset = assay == assay_5)

## ------ Identify the Gold Truth DEGs --------
# DS1 - ~1,500 to 3,000 cells -> ~ 700 cells (25%), ~ 1,500 cells (50%), ~ 2,200 cells (75%)
# DS6 - ~3,000 to 3,750 cells -> ~750 to 900 cells (25%)
# DS2lv - ~ 3,000 to 3,750 cells -> ~750 to 900 cells (25%)

Idents(data_3) <- "cell_type_ontology_term_id"
Idents(data_5) <- "cell_type_ontology_term_id"

### By individual 
data_3 <- SplitObject(data_3, split.by = "new_id")
data_5 <- SplitObject(data_5, split.by = "new_id")

raw_markers <- list()
truth <- list()

for(i in names(data_3)) {
  print(i)
  obj_3 <- data_3[[i]]
  obj_5 <- data_5[[i]]
  
  obj_3 <- NormalizeData(obj_3, normalization.method = "LogNormalize", scale.factor = 10000)
  obj_5 <- NormalizeData(obj_5, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Set the assays into appreopriate layer to calculate the difference (not FC)
  obj_3 <- SetAssayData(obj_3, layer = "scale.data", new.data = GetAssayData(obj_3, layer = "data"))
  obj_5 <- SetAssayData(obj_5, layer = "scale.data", new.data = GetAssayData(obj_5, layer = "data"))
  
  # Get the DEGs
  gold3 <- FindMarkers(obj_3, ident.1 = cell_type1, ident.2 = cell_type2, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5, fc.slot = "scale.data")
  gold5 <- FindMarkers(obj_5, ident.1 = cell_type1, ident.2 = cell_type2, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5, fc.slot = "scale.data")

  # Match the order of the genes in two gold standards:
  order <- match(rownames(gold3), rownames(gold5))
  gold5 <- gold5[order,]
  
  if(!(all(rownames(gold5) == rownames(gold3)))) {
    stop("Gene names do not match between 5 and 3 prime for ", i)
  } # validate
  
  raw_markers[[paste0(i, "_3")]] <- gold3
  raw_markers[[paste0(i, "_5")]] <- gold5
  
  # Set the desired thresholds 
  truth[[i]] <- data.frame(
    gene = rownames(gold3),
    dir_sig = ifelse((gold3$p_val_adj < 0.05 & gold5$p_val_adj < 0.05) & (sign(gold3$avg_diff) == sign(gold5$avg_diff)),
                     ifelse(gold3$avg_diff >= 0,"up","down"), 
                     "ns"),
    dir_gen = ifelse(sign(gold3$avg_diff) == sign(gold5$avg_diff),
                     ifelse(gold3$avg_diff >= 0,"up","down"), 
                     "opp_dir")
  )
}

saveRDS(raw_markers, "./DS_indiv_truth_raw_markers_intersect_noFC.rds") # Specify the datasets
saveRDS(truth, "./DS_indiv_truth_labels_intersect_noFC.rds") # Specify the datasets

# Then INTERSECTION 
truth <- data.frame(
  gene = rownames(gold3),
  dir_sig = ifelse((gold3$p_val_adj < 0.05 & gold5$p_val_adj < 0.05) & (abs(gold3$avg_log2FC) > 0.25 & abs(gold5$avg_log2FC) > 0.25) & sign(gold3$avg_log2FC) == sign(gold5$avg_log2FC),
                   ifelse(gold3$avg_log2FC > 0,"up","down"), 
                   "ns"),
  dir_gen = ifelse(sign(gold3$avg_log2FC) == sign(gold5$avg_log2FC),
                   ifelse(gold3$avg_log2FC > 0,"up","down"), 
                   "opp_dir")
)

saveRDS(raw_markers, "./DS6_mer_truth_raw_markers_intersect.rds")
saveRDS(truth, "./DS6_mer_truth_labels_intersect.rds")


## ------- Introduce artificial sparsification --------
# DS1:
cell_3 <- c("CD4-positive, alpha-beta T cell", "alpha-beta T cell", "CD8-positive, alpha-beta T cell", "gamma-delta T cell", "T-helper 17 cell", "double negative thymocyte", "thymocyte") # list of cell types for 3' 25% DS1
cell_5 <- c("CD4-positive, alpha-beta T cell", "alpha-beta T cell", "regulatory T cell", "T cell", "double-positive, alpha-beta thymocyte") # list of cell types for 5' 25% DS1

# DS6:
cell_3 <- c("granulocyte monocyte progenitor cell", "B cell progenitors", "dendritic cell", "early lymphoid progenitor", "macrophage", "promonocyte", "RBCs") 
cell_5 <- c("granulocyte monocyte progenitor cell", "B cell progenitors", "B cell", "mast cell", "megakaryocyte", "osteoclast") 

# DS2 liver:
cell_3 <- c("B cell progenitors", "B cell", "hematopoietic progenitor cell", "innate lymphoid cell", "T cell", "dendritic cell", "endothelial cell", "granulocyte") 
cell_5 <- c("B cell progenitors", "B cell", "hematopoietic progenitor cell", "innate lymphoid cell", "T cell", "natural killer cell") 

# Subset the desired cell types for 3' and 5' assays for each dataset.
data_3 <- subset(data, subset = cell_type_ontology_term_id %in% cell_3 & assay == assay_3) 
data_5 <- subset(data, subset = cell_type_ontology_term_id %in% cell_5 & assay == assay_5)

data_3@meta.data <- droplevels(data_3@meta.data)
data_5@meta.data <- droplevels(data_5@meta.data)

## ------- Merge the data --------
merged_data <- merge(data_3, y = data_5)
merged_data <- JoinLayers(merged_data)
saveRDS(merged_data, "./DS_merged_data_scarce_cell_types_25.rds") # Specify the dataset

## ------- Proceed with the correction/log normalization using the appropriate .R files

## ------ Load the integrated data --------

# -----> For 100%
# Scaled_LogNorm, Scaled_fastMNN, Scaled_scanorama, Scaled_limma, Scaled_m3drop, Scaled_scvi
norm_data <- readRDS("./Scaled_X.rds") # Specify the correction technique - X
norm_data <- norm_data[c("F29_45P", "F30_45P", "F38_45P",  "P1_CD3P")] # for DS1
norm_data <- norm_data[c("F32_CD45P_liver", "F34_CD45P_liver", "F38_CD45P_liver",  "F41_CD45P_liver")] # for DS2lv

# ----> For the 25%
norm_data <- readRDS("./DS1_merged_list_scTransf_regr_assay_split_25.rds")

Idents(norm_data) <- "cell_type_ontology_term_id" # for DS2 and 6
# For M3Drop, scanorama, and scvi
norm_data <- SetAssayData(norm_data, layer = "scale.data", new.data = GetAssayData(norm_data, layer = "counts"))
# For fastMNN, limma, and logN
norm_data <- SetAssayData(norm_data, layer = "scale.data", new.data = GetAssayData(norm_data, layer = "data"))

# Now the data is ready for the DE analysis, but we need to split it by individual first to be able to compare the results with the gold standard.
norm_data_list <- SplitObject(norm_data, split.by = "new_id")

raw_markers_integ <- list()
test <- list()

# By individual 
for(i in names(norm_data_list)) {
  print(i)
  obj <- norm_data_list[[i]]
  
  #obj <- SetAssayData(obj, layer = "scale.data", new.data = GetAssayData(obj, layer = "counts")) # For 100% -> M3, scanorama, scvi
  #obj <- SetAssayData(obj, layer = "scale.data", new.data = GetAssayData(obj, layer = "data")) # For 100% -> logN, fastMNN, limma
  
  Idents(obj) <- "cell_type_ontology_term_id" 

  if(!(cell_type1 %in% levels(Idents(obj)))) {
    print("Error")
  }
  if(!(cell_type2 %in% levels(Idents(obj)))) {
    print("Error")
  }
  
  integ_markers <- FindMarkers(obj, ident.1 = cell_type1, ident.2 = cell_type2, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5, slot = "scale.data", fc.slot = "scale.data")

  raw_markers_integ[[i]] <- integ_markers
  
  # Categorize
  test[[i]] <- data.frame(
    gene = rownames(integ_markers),
    dir_sig = ifelse((integ_markers$p_val_adj < 0.05),
                     ifelse(integ_markers$avg_diff >= 0,"up","down"), 
                     "ns"),
    dir_gen = ifelse(integ_markers$avg_diff >= 0,"up","down")
  )
}

# For 25%
saveRDS(raw_markers_integ, "./DS_X_raw_markers_integ_25_noFC.rds") # Specify the dataset and the correction technique - X
saveRDS(test, "./DS_X_test_labels_25_noFC.rds")
# For 100%
saveRDS(raw_markers_integ, "./DS_X_raw_markers_integ_100_noFC.rds")
saveRDS(test, "./DS_X_test_labels_100_noFC.rds")


## ------ Evaluation of the results --------

#### ---------- Create common truth from per indiv truths ---------- ####
ds1 <- readRDS("./DS_indiv_truth_labels_intersect_noFC.rds") # Specify the dataset

three_up   <- Reduce(intersect, lapply(ds1, \(x) x$gene[x$dir_sig == "up"])) 
three_down <- Reduce(intersect, lapply(ds1, \(x) x$gene[x$dir_sig == "down"]))

rownames(ds1[[1]]) <- ds1[[1]]$gene
rownames(ds1[[2]]) <- ds1[[2]]$gene
rownames(ds1[[3]]) <- ds1[[3]]$gene
rownames(ds1[[4]]) <- ds1[[4]]$gene

stopifnot(setequal(rownames(ds1[[1]]), rownames(ds1[[2]])))
stopifnot(setequal(rownames(ds1[[1]]), rownames(ds1[[3]])))
stopifnot(setequal(rownames(ds1[[1]]), rownames(ds1[[4]])))

order <- match(rownames(ds1[[1]]), rownames(ds1[[2]]))
ds1[[2]] <- ds1[[2]][order,]

order <- match(rownames(ds1[[1]]), rownames(ds1[[3]]))
ds1[[3]] <- ds1[[3]][order,]

order <- match(rownames(ds1[[1]]), rownames(ds1[[4]]))
ds1[[4]] <- ds1[[4]][order,]

## For DS1 and 2lv
new <- data.frame(
  gene = ds1[[1]]$gene,
  dir_sig = ifelse(
    ds1[[1]]$gene %in% three_up, "up",
    ifelse(ds1[[1]]$gene %in% three_down, "down", "ns")),
  dir_gen = ifelse(
    ds1[[1]]$dir_gen == "up" & ds1[[2]]$dir_gen == "up" & ds1[[3]]$dir_gen == "up" & ds1[[4]]$dir_gen == "up", "up",
    ifelse(
      ds1[[1]]$dir_gen == "down" & ds1[[2]]$dir_gen == "down" & ds1[[3]]$dir_gen == "down" & ds1[[4]]$dir_gen == "down", "down",
      "opp_dir"
    )
  )
)

# For DS6
new <- data.frame(
  gene = ds1[[1]]$gene,
  dir_sig = ifelse(
    ds1[[1]]$gene %in% three_up, "up",
    ifelse(ds1[[1]]$gene %in% three_down, "down", "ns")),
  dir_gen = ifelse(
    ds1[[1]]$dir_gen == "up" & ds1[[2]]$dir_gen == "up" & ds1[[3]]$dir_gen == "up", "up",
    ifelse(
      ds1[[1]]$dir_gen == "down" & ds1[[2]]$dir_gen == "down" & ds1[[3]]$dir_gen == "down", "down",
      "opp_dir"
    )
  )
)

new$dir_sig[new$dir_gen == "opp_dir" & new$dir_sig != "ns"] <- "ns"

saveRDS(new, "./clean_sparse/DS6_intersect_truth_labels_intersect_noFC.rds")


# ---- MCC Scores ----

techno <- c("logN", "fast", "scanorama", "limma", "M3Drop", "scvi")

DS <- "DS1"
n <- "100" # 25 or 100
indiv <- "intersect" 

for(tec in techno) {
  print(tec)
  truth <- readRDS(paste0("./",DS,"_",indiv,"_truth_labels_intersect_noFC.rds"))
  test <- readRDS(paste0("./",DS,"_",tec,"_test_labels_",n,"_noFC.rds"))

  # Here, we need two lists: truth and test
  mcc <- numeric(length(names(test)))
  names(mcc) <- names(test)

  confusion_matrix_data <- list()

  for(i in names(test)) {
    print(i)

    truth_i <- truth # one truth
    test_i <- test[[i]]
  
    # Match the order of the genes in two lists:
    common <- intersect(truth_i$gene, test_i$gene)
    test_i <- test_i[test_i$gene %in% common,]
    truth_i <- truth_i[truth_i$gene %in% common,]
    test_i <- test_i[match(truth_i$gene, test_i$gene),]
    if(!all(test_i$gene == truth_i$gene)) {
      stop("Gene names do not match between truth and test for ", i)
    }
  
   # Calculate MCC score:
   ## Step 1: Exclude genes that are significant in the test but not in truth if it’s in the same direction
   exclude <- ifelse((test_i$dir_sig != "ns") & (truth_i$dir_sig == "ns") & (test_i$dir_gen == truth_i$dir_gen), TRUE, FALSE)
   sum(exclude)
   test_i <- test_i[!exclude,]
   truth_i <- truth_i[!exclude,]
   if(!all(test_i$gene == truth_i$gene)) {
     stop("Gene names do not match between truth and test for ", i)
   }
  
   # Calculate the confusion matrix
  
   # Significant in both and in the same direction
   TP <- sum(truth_i$dir_sig != "ns" & truth_i$dir_sig == test_i$dir_sig)
   # Opposite direction disregarding the level of significance 
   FP <- sum((truth_i$dir_gen == "up" & test_i$dir_gen == "down") | (truth_i$dir_gen == "down" & test_i$dir_gen == "up")) 
   # Significant in the gold standard but not in integrated
   FN <- sum((truth_i$dir_gen == test_i$dir_gen) & (truth_i$dir_sig != "ns" & test_i$dir_sig == "ns"))
   # Not significant in both
   TN <- sum((truth_i$dir_gen == test_i$dir_gen) & (truth_i$dir_sig == "ns" & test_i$dir_sig == "ns"))
  
   # Check:
   # Significant in both and in the same direction
   TPg <- truth_i$gene[truth_i$dir_sig != "ns" & truth_i$dir_sig == test_i$dir_sig]
   # Opposite direction disregarding the level of significance 
   FPg <- truth_i[(truth_i$dir_gen == "up" & test_i$dir_gen == "down") | (truth_i$dir_gen == "down" & test_i$dir_gen == "up"),]
   # Significant in the gold standard but not in integrated
   FNg <- truth_i$gene[(truth_i$dir_gen == test_i$dir_gen) & (truth_i$dir_sig != "ns" & test_i$dir_sig == "ns")]
   # Not significant in both
   TNg <- truth_i$gene[(truth_i$dir_gen == test_i$dir_gen) & (truth_i$dir_sig == "ns" & test_i$dir_sig == "ns")]
   intersect(TNg, FNg)
  
   # Total counts
   tot <- TP+TN+FP+FN
   opp <- sum(truth_i$dir_gen == "opp_dir")
   tot + opp
   nrow(truth_i) 
  
   if(tot + opp != nrow(truth_i)) {
     stop("Total counts do not match for ", i)
   }
  
   TP <- as.numeric(TP)
   FP <- as.numeric(FP)
   TN <- as.numeric(TN)
   FN <- as.numeric(FN)
   MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
   MCC * 100
  
   mcc[[i]] <- round(MCC * 100, 3)
  
   confusion_matrix_data[[i]] <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
                                  dimnames = list(Truth = c("Positive", "Negative"),
                                                  Integrated = c("Positive", "Negative")))
  }

  # Save the MCC scores
  saveRDS(mcc, paste0("./mcc_scores/",DS,"_",tec,"_MCC_",n,"_",indiv,"_noFC.rds"))
  saveRDS(confusion_matrix_data, paste0("./mcc_scores/",DS,"_",tec,"_confusion_",n,"_",indiv,"_noFC.rds"))
}


# ---- Look at 3' vs 5' differences -----

###### --- Find markers by subject/individual -----
data <- readRDS("./Scaled_logNorm.rds") # Do it for all datasets

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.0001,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

gene_list <- table(unlist(genes))

saveRDS(dge_list, "./DS_3_vs_5_sign_by_indiv.rds") # Specify the dataset and the tissue, if needed
saveRDS(genes, "./DS_3_vs_5_sign_by_indiv_genes.rds")

###### Combine all ####
ds1 <- readRDS("./DS1_3_vs_5_sign_by_indiv_genes.rds")
ds2bm <- readRDS("./DS2_3_vs_5_sign_by_indiv_genes_bm.rds")
ds2sp <- readRDS("./DS2_3_vs_5_sign_by_indiv_genes_sp.rds")
ds2kd <- readRDS("./DS2_3_vs_5_sign_by_indiv_genes_kd.rds")
ds2lv <- readRDS("./DS2_3_vs_5_sign_by_indiv_genes_lv.rds")
ds2th <- readRDS("./DS2_3_vs_5_sign_by_indiv_genes_th.rds")
ds3 <- readRDS("./DS3_3_vs_5_sign_by_indiv_genes.rds")
ds4 <- readRDS("./DS4_3_vs_5_sign_by_indiv_genes.rds")
ds5 <- readRDS("./DS5_3_vs_5_sign_by_indiv_genes.rds")
ds6 <- readRDS("./DS6_3_vs_5_sign_by_indiv_genes.rds")

all <- c(ds1, ds2bm, ds2sp, ds2kd, ds2lv, ds2th, ds3, ds4, ds5, ds6)
genes <- list()
for (i in 1:length(all)) {
  all[[i]] <- all[[i]][all[[i]]$p_val_adj < 1e-100,] # attempted different p-value thresholds to get a reasonable number of genes for the table
  genes[[i]] <- rownames(all[[i]])
}
gene_table <- table(unlist(genes))

saveRDS(gene_table, "./DS_all_3_vs_5_sign_by_indiv_1e-100.rds")



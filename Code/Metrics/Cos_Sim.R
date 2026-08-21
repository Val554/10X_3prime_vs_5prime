# Purpose: Create random sets of genes with matching expression profiles to the bias genes. Calculate Cosine Similarity for all sets.

# Loading required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(presto)

# Establishing helper functions #
my_rowMeans <- function(x) {
  if (!is.null(ncol(x))) {
    if (ncol(x) > 1) {
      return(Matrix::rowMeans(x))
    }
  }
  return(x);
}

my_rowSums <- function(x) {
  if (!is.null(ncol(x))) {
    if (ncol(x) > 1) {
      return(Matrix::rowSums(x))
    }
  }
  return(x);
}

my_rowVars <- function(x) {
  if (!is.null(nrow(x))) {
    if (nrow(x) > 1) {
      center <- Matrix::rowMeans(x)
      n <- ncol(x)
      vars <- n/(n-1)*(Matrix::rowMeans(x^2) - center^2)
      return(vars)
    }
  }
  return(rep(0, length(x)));
}

group_rowmeans <- function(MAT, group_labs, type=c("mean","sum", "var")) {
  d <- split(seq(ncol(MAT)), group_labs);
  if (type[1] == "mean") {
    if(nrow(MAT) > 1) {
      mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
    } else {
      mus <- sapply(d, function(group) mean(MAT[,group])) # only one row
    }
  } 
  if (type[1] == "sum") {
    if (nrow(MAT) > 1) {
      mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
    } else {
      mus <- sapply(d, function(group) sum(MAT[,group])) # only one row
    }
  } 
  if (type[1] == "var") {
    if (nrow(MAT) > 1) {
      mus <- sapply(d, function(group) my_rowVars(MAT[,group]))
    } else {
      mus <- sapply(d, function(group) stats::var(MAT[,group])) # only one row
    }
  }
  return(mus);
}

# Setting a function to compute the cosine similarity score
compute_metrics <- function(vec1, vec2) {
  if(length(vec1) != length(vec2)) {
    stop("Vectors must be of the same length")
  }
  if(all(vec1 == 0) || all(vec2 == 0)) {
    return(NA)  # Return NA if either vector is all zeros
  }
  cosine_similarity <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  return(cosine_similarity)
}

# Load the data
data <- readRDS("./Scaled_LogNorm.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
num <- 800
set <- "1"

data_merged <- merge(data[[1]], y = data[-1])
data_merged <- JoinLayers(data_merged)

common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50
#common_ind <- c(names(gene_table[gene_table %in% c(3:35)])) # 600: 100

indiv_genes <- list()

# For the selected list of genes calculate the proportions of lowly/highly expressed genes for each set (800/600)
for(indiv_name in names(data)) {
  print(indiv_name)
  indiv <- data[[indiv_name]]
  indiv <- indiv[(rownames(indiv) %in% common_ind), ]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  assay <- assay_list[[1]]
  assay <- SplitObject(assay, split.by = "cell_type_ontology_term_id")
  
  genes_3 <- list()
  for(j in names(assay)) {
    cell_type <- assay[[j]]
    if (dim(cell_type)[2] < 5) {
      next
    }
    expr <- SeuratObject::GetAssayData(cell_type, layer = "counts")
    genes_3[[j]] <- rownames(expr[rowSums(expr) > ceiling(0.05 * ncol(expr)),]) 
  }
  
  assay <- assay_list[[2]]
  assay <- SplitObject(assay, split.by = "cell_type_ontology_term_id")
  
  genes_5 <- list()
  for(j in names(assay)) {
    cell_type <- assay[[j]]
    if (dim(cell_type)[2] < 5) {
      next
    }
    expr <- GetAssayData(cell_type, layer = "counts")
    genes_5[[j]] <- rownames(expr[rowSums(expr) > ceiling(0.05 * ncol(expr)),]) 
  }
  
  genes_3 <- Reduce(intersect, genes_3)
  genes_5 <- Reduce(intersect, genes_5)
  
  common_genes <- intersect(genes_3, genes_5)
  indiv_genes[[indiv_name]] <- common_genes
}

high_genes <- Reduce(union, indiv_genes)
low_genes <- common_ind[!common_ind %in% high_genes]

# Compute the highly expressed gene list - random
for(indiv_name in names(data)) {
  print(indiv_name)
  indiv <- data[[indiv_name]]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  assay <- assay_list[[1]]
  assay <- SplitObject(assay, split.by = "cell_type_ontology_term_id")
  
  genes_3 <- list()
  for(j in names(assay)) {
    cell_type <- assay[[j]]
    if (dim(cell_type)[2] < 5) {
      next
    }
    expr <- SeuratObject::GetAssayData(cell_type, layer = "counts")
    genes_3[[j]] <- rownames(expr[rowSums(expr) > ceiling(0.05 * ncol(expr)),]) 
  }
  
  assay <- assay_list[[2]]
  assay <- SplitObject(assay, split.by = "cell_type_ontology_term_id")
  
  genes_5 <- list()
  for(j in names(assay)) {
    cell_type <- assay[[j]]
    if (dim(cell_type)[2] < 5) {
      next
    }
    expr <- GetAssayData(cell_type, layer = "counts")
    genes_5[[j]] <- rownames(expr[rowSums(expr) > ceiling(0.05 * ncol(expr)),]) 
  }
  
  genes_3 <- Reduce(intersect, genes_3)
  genes_5 <- Reduce(intersect, genes_5)
  
  common_genes <- intersect(genes_3, genes_5)
  indiv_genes[[indiv_name]] <- common_genes
}

common_genes <- Reduce(union, indiv_genes)
common_genes <- common_genes[!(common_genes %in% common_ind)]
inx <- sample(1:length(common_genes), replace = FALSE, length(common_ind))
high_rand <- common_genes[inx]

low_rand_g <- rownames(data_merged)[!rownames(data_merged) %in% common_genes] # for lowly
low_rand_g <- low_rand_g[!(low_rand_g %in% common_ind)]
inx <- sample(1:length(low_rand_g), replace = FALSE, length(common_ind))
low_rand <- low_rand_g[inx]

gen_rand_g <- rownames(data_merged)[!rownames(data_merged) %in% common_ind] # for general
inx <- sample(1:length(gen_rand_g), replace = FALSE, length(common_ind))
gen_rand <- gen_rand_g[inx]

# create a list of random genes with the same proportions of high and low genes that is in the actual gene list
inx_h <- sample(1:length(common_genes), replace = FALSE, length(high_genes))
inx_l <- sample(1:length(low_rand_g), replace = FALSE, length(low_genes))
prop_rand <- c(common_genes[inx_h], low_rand_g[inx_l])

inx_h <- sample(1:length(common_genes), replace = FALSE, length(high_genes))
inx_l <- sample(1:length(low_rand_g), replace = FALSE, length(low_genes))
prop_rand2 <- c(common_genes[inx_h], low_rand_g[inx_l]) 

inx_h <- sample(1:length(common_genes), replace = FALSE, length(high_genes))
inx_l <- sample(1:length(low_rand_g), replace = FALSE, length(low_genes))
prop_rand3 <- c(common_genes[inx_h], low_rand_g[inx_l])


##### Cosine similarity calculations  ####


# ------ Random gene list proportional ----

### Random List 1 
result_list <- list()
cell_count_list <- list()
cell_names <- list()

for(indiv_name in names(data)) {
  print(indiv_name)
  # Split by assay
  indiv <- data[[indiv_name]]
  indiv <- ScaleData(indiv, features = prop_rand)
  indiv <- indiv[(rownames(indiv) %in% prop_rand), ]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  # Step 3a: Calculate average expression by cell type for 3' assay
  
  # 3'
  assay <- assay_list[[1]] # made sure previously that 3' is the first assay in the factor
  assay@meta.data <- droplevels(assay@meta.data) # drop unused levels if a factor
  
  expr_3 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  # check that the order is the same, crash if false!!!
  if(!identical(colnames(expr_3), rownames(assay@meta.data))) {
    stop("Error: The column names in the 3 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  #print(table(cell_type)) # save these for reference, might not keep; 5 is too few, exclude from the analysis
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_3")]] <- x
  
  # Calculate average gene expression for each cell type in the dataset
  averaged_expression_3 <- group_rowmeans(expr_3, cell_type, "mean") 
  
  # Step 3b: Calculate average expression by cell type for 5' assay
  
  # 5'
  assay <- assay_list[[2]]
  assay@meta.data <- droplevels(assay@meta.data) 
  
  expr_5 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  if(!identical(colnames(expr_5), rownames(assay@meta.data))) {
    stop("Error: The column names in the 5 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_5")]] <- x
  
  # Calculating average gene expression for each gene by cell_type
  averaged_expression_5 <- group_rowmeans(expr_5, cell_type, "mean")
  
  # Check which cell type match between the 3' and 5'
  # Only get those that > 5 per cell type and match between 3' and 5'
  cell_counts <- cell_count_list[[paste0(indiv_name, "_3")]]
  cells_3 <- names(cell_counts[cell_counts > 5])
  cell_counts <- cell_count_list[[paste0(indiv_name, "_5")]]
  cells_5 <- names(cell_counts[cell_counts > 5])
  
  ident_cell_types <- intersect(cells_3, cells_5)
  cell_names[[indiv_name]] <- ident_cell_types
  
  # Only include genes that are matching
  idx <- intersect(rownames(averaged_expression_3), rownames(averaged_expression_5))
  averaged_expression_3 <- averaged_expression_3[idx,] # subsetting only matching genes
  averaged_expression_5 <- averaged_expression_5[idx,]
  
  # Order matching genes
  averaged_expression_3 <- averaged_expression_3[match(rownames(averaged_expression_5), 
                                                       rownames(averaged_expression_3)), ]
  
  # Create an named vector for the output
  metrics <- numeric(length(ident_cell_types))
  names(metrics) <- ident_cell_types
  
  # For matching cell types between 3' and 5' calculate the metrics
  for(i in ident_cell_types) {
    #i <- "CL:0000084"
    print(i)
    metrics[i] <- compute_metrics(averaged_expression_3[,i], averaged_expression_5[,i])
    
  }
  
  # Results for each patient are stored in a list of lists
  result_list[[indiv_name]] <- metrics
}

result_list[["IndivNames"]] <- names(data)
result_list[["CellTypes"]] <- cell_names
result_list[["AllCellTypes"]] <- unique(unlist(result_list$CellTypes))
result_list[["CellCounts"]] <- cell_count_list

result_list[["Genes"]] <- prop_rand

saveRDS(result_list, paste0("./DS",set,"_cos_",num,"_rand_prop_SCALED.rds"))

#### Random List 2
result_list <- list()
cell_count_list <- list()
cell_names <- list()

for(indiv_name in names(data)) {
  print(indiv_name)
  # Split by assay
  indiv <- data[[indiv_name]]
  indiv <- ScaleData(indiv, features = prop_rand2)
  indiv <- indiv[(rownames(indiv) %in% prop_rand2), ]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  # Step 3a: Calculate average expression by cell type for 3' assay
  
  # 3'
  assay <- assay_list[[1]] # made sure previously that 3' is the first assay in the factor
  assay@meta.data <- droplevels(assay@meta.data) # drop unused levels if a factor
  
  expr_3 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  # check that the order is the same, crash if false!!!
  if(!identical(colnames(expr_3), rownames(assay@meta.data))) {
    stop("Error: The column names in the 3 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  #print(table(cell_type)) # save these for reference, might not keep; 5 is too few, exclude from the analysis
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_3")]] <- x
  
  # Calculate average gene expression for each cell type in the dataset
  averaged_expression_3 <- group_rowmeans(expr_3, cell_type, "mean") 
  
  # Step 3b: Calculate average expression by cell type for 5' assay
  
  # 5'
  assay <- assay_list[[2]]
  assay@meta.data <- droplevels(assay@meta.data) 
  
  expr_5 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  if(!identical(colnames(expr_5), rownames(assay@meta.data))) {
    stop("Error: The column names in the 5 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_5")]] <- x
  
  # Calculating average gene expression for each gene by cell_type
  averaged_expression_5 <- group_rowmeans(expr_5, cell_type, "mean")
  
  # Check which cell type match between the 3' and 5'
  # Only get those that > 5 per cell type and match between 3' and 5'
  cell_counts <- cell_count_list[[paste0(indiv_name, "_3")]]
  cells_3 <- names(cell_counts[cell_counts > 5])
  cell_counts <- cell_count_list[[paste0(indiv_name, "_5")]]
  cells_5 <- names(cell_counts[cell_counts > 5])
  
  ident_cell_types <- intersect(cells_3, cells_5)
  cell_names[[indiv_name]] <- ident_cell_types
  
  # Only include genes that are matching
  idx <- intersect(rownames(averaged_expression_3), rownames(averaged_expression_5))
  averaged_expression_3 <- averaged_expression_3[idx,] # subsetting only matching genes
  averaged_expression_5 <- averaged_expression_5[idx,]
  
  # Order matching genes
  averaged_expression_3 <- averaged_expression_3[match(rownames(averaged_expression_5), 
                                                       rownames(averaged_expression_3)), ]
  
  # Create an named vector for the output
  metrics <- numeric(length(ident_cell_types))
  names(metrics) <- ident_cell_types
  
  # For matching cell types between 3' and 5' calculate the metrics
  for(i in ident_cell_types) {
    print(i)
    metrics[i] <- compute_metrics(averaged_expression_3[,i], averaged_expression_5[,i])
    
  }
  
  # Results for each patient are stored in a list of lists
  result_list[[indiv_name]] <- metrics
}

result_list[["IndivNames"]] <- names(data)
result_list[["CellTypes"]] <- cell_names
result_list[["AllCellTypes"]] <- unique(unlist(result_list$CellTypes))
result_list[["CellCounts"]] <- cell_count_list

result_list[["Genes"]] <- prop_rand2

saveRDS(result_list, paste0("./DS",set,"_cos_",num,"_rand_prop_2_SCALED.rds"))

### Random List 3
result_list <- list()
cell_count_list <- list()
cell_names <- list()

for(indiv_name in names(data)) {
  print(indiv_name)
  # Split by assay
  indiv <- data[[indiv_name]]
  indiv <- ScaleData(indiv, features = prop_rand3)
  indiv <- indiv[(rownames(indiv) %in% prop_rand3), ]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  # Step 3a: Calculate average expression by cell type for 3' assay
  
  # 3'
  assay <- assay_list[[1]] # made sure previously that 3' is the first assay in the factor
  assay@meta.data <- droplevels(assay@meta.data) # drop unused levels if a factor
  
  expr_3 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  # check that the order is the same, crash if false!!!
  if(!identical(colnames(expr_3), rownames(assay@meta.data))) {
    stop("Error: The column names in the 3 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  #print(table(cell_type)) # save these for reference, might not keep; 5 is too few, exclude from the analysis
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_3")]] <- x
  
  # Calculate average gene expression for each cell type in the dataset
  averaged_expression_3 <- group_rowmeans(expr_3, cell_type, "mean") 
  
  # Step 3b: Calculate average expression by cell type for 5' assay
  
  # 5'
  assay <- assay_list[[2]]
  assay@meta.data <- droplevels(assay@meta.data) 
  
  expr_5 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  if(!identical(colnames(expr_5), rownames(assay@meta.data))) {
    stop("Error: The column names in the 5 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_5")]] <- x
  
  # Calculating average gene expression for each gene by cell_type
  averaged_expression_5 <- group_rowmeans(expr_5, cell_type, "mean")
  
  # Check which cell type match between the 3' and 5'
  # Only get those that > 5 per cell type and match between 3' and 5'
  cell_counts <- cell_count_list[[paste0(indiv_name, "_3")]]
  cells_3 <- names(cell_counts[cell_counts > 5])
  cell_counts <- cell_count_list[[paste0(indiv_name, "_5")]]
  cells_5 <- names(cell_counts[cell_counts > 5])
  
  ident_cell_types <- intersect(cells_3, cells_5)
  cell_names[[indiv_name]] <- ident_cell_types
  
  # Only include genes that are matching
  idx <- intersect(rownames(averaged_expression_3), rownames(averaged_expression_5))
  averaged_expression_3 <- averaged_expression_3[idx,] # subsetting only matching genes
  averaged_expression_5 <- averaged_expression_5[idx,]
  
  # Order matching genes
  averaged_expression_3 <- averaged_expression_3[match(rownames(averaged_expression_5), 
                                                       rownames(averaged_expression_3)), ]
  
  # Create an named vector for the output
  metrics <- numeric(length(ident_cell_types))
  names(metrics) <- ident_cell_types
  
  # For matching cell types between 3' and 5' calculate the metrics
  for(i in ident_cell_types) {
    print(i)
    metrics[i] <- compute_metrics(averaged_expression_3[,i], averaged_expression_5[,i])
    
  }
  
  # Results for each patient are stored in a list of lists
  result_list[[indiv_name]] <- metrics
}

result_list[["IndivNames"]] <- names(data)
result_list[["CellTypes"]] <- cell_names
result_list[["AllCellTypes"]] <- unique(unlist(result_list$CellTypes))
result_list[["CellCounts"]] <- cell_count_list

result_list[["Genes"]] <- prop_rand3

saveRDS(result_list, paste0("./DS",set,"_cos_",num,"_rand_prop_3_SCALED.rds"))


### General lists
set <- "6"

num <- 600
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-100.rds")
common_ind <- c(names(gene_table[gene_table %in% c(3:35)])) # 600: 100

num <- 800
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50


# gene list general
result_list <- list()
cell_count_list <- list()
cell_names <- list()
for(indiv_name in names(data)) {
  print(indiv_name)
  # Split by assay
  indiv <- data[[indiv_name]]
  indiv <- ScaleData(indiv, features = common_ind)
  indiv <- indiv[(rownames(indiv) %in% common_ind), ]
  
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  # Step 3a: Calculate average expression by cell type for 3' assay
  
  # 3'
  assay <- assay_list[[1]] # made sure previously that 3' is the first assay in the factor
  assay@meta.data <- droplevels(assay@meta.data) # drop unused levels if a factor
  
  expr_3 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  # check that the order is the same, crash if false!!!
  if(!identical(colnames(expr_3), rownames(assay@meta.data))) {
    stop("Error: The column names in the 3 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  #print(table(cell_type)) # save these for reference, might not keep; 5 is too few, exclude from the analysis
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_3")]] <- x
  
  # Calculate average gene expression for each cell type in the dataset
  averaged_expression_3 <- group_rowmeans(expr_3, cell_type, "mean") 
  
  # Step 3b: Calculate average expression by cell type for 5' assay
  
  # 5'
  assay <- assay_list[[2]]
  assay@meta.data <- droplevels(assay@meta.data) 
  
  expr_5 <- GetAssayData(assay, layer = "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  if(!identical(colnames(expr_5), rownames(assay@meta.data))) {
    stop("Error: The column names in the 5 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_5")]] <- x
  
  # Calculating average gene expression for each gene by cell_type
  averaged_expression_5 <- group_rowmeans(expr_5, cell_type, "mean")
  
  # Check which cell type match between the 3' and 5'
  # Only get those that > 5 per cell type and match between 3' and 5'
  cell_counts <- cell_count_list[[paste0(indiv_name, "_3")]]
  cells_3 <- names(cell_counts[cell_counts > 5])
  cell_counts <- cell_count_list[[paste0(indiv_name, "_5")]]
  cells_5 <- names(cell_counts[cell_counts > 5])
  
  ident_cell_types <- intersect(cells_3, cells_5)
  cell_names[[indiv_name]] <- ident_cell_types
  
  # Only include genes that are matching
  idx <- intersect(rownames(averaged_expression_3), rownames(averaged_expression_5))
  averaged_expression_3 <- averaged_expression_3[idx,] # subsetting only matching genes
  averaged_expression_5 <- averaged_expression_5[idx,]
  
  # Order matching genes
  averaged_expression_3 <- averaged_expression_3[match(rownames(averaged_expression_5), 
                                                       rownames(averaged_expression_3)), ]
  
  # Create an named vector for the output
  metrics <- numeric(length(ident_cell_types))
  names(metrics) <- ident_cell_types
  
  # For matching cell types between 3' and 5' calculate the metrics
  for(i in ident_cell_types) {
    print(i)
    metrics[i] <- compute_metrics(averaged_expression_3[,i], averaged_expression_5[,i])
    
  }
  
  # Results for each patient are stored in a list of lists
  result_list[[indiv_name]] <- metrics
}

result_list[["IndivNames"]] <- names(data)
result_list[["CellTypes"]] <- cell_names
result_list[["AllCellTypes"]] <- unique(unlist(result_list$CellTypes))
result_list[["CellCounts"]] <- cell_count_list

result_list[["Genes"]] <- common_ind

saveRDS(result_list, paste0("./DS",set,"_cos_",num,"_SCALED.rds"))




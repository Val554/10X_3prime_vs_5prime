# Purpose: Calculating the metrics for the comparison of the 3' and 5' assays for each batch correction method. 
# The metrics include correlation, cosine similarity, MSE, JSD, Euclidean distance, Manhattan distance, mean and variance of gene expression for each cell type. 
# The results are stored in a list of lists for each individual and cell type and saved as RDS files for further analysis and visualization.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(parallel) 
library(philentropy)

# Read the corrected data for each dataset and batch correction technique (Seurat object)
object <- readRDS("./X_corrected.RDS") 

#--- Additional step for the uncorrected data (log normalized)
object <- split(object, f = object@meta.data$new_id)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, nfeatures = 3000) 
saveRDS(VariableFeatures(object), "./hvf.rds") # Save the highly variable features for later use.
object <- JoinLayers(object)
#---


# Function to compute all relevant metrics 
compute_metrics <- function(vec1, vec2) {
  correlation <- cor(vec1, vec2, use = "pairwise.complete.obs")
  cosine_similarity <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  
  scale <- min(c(min(vec1), min(vec2)))
  vec1_adj <- vec1 - scale + 1e-6
  vec1_prob <- vec1_adj / sum(vec1_adj)
  vec2_adj <- vec2 - scale + 1e-6
  vec2_prob <- vec2_adj / sum(vec2_adj)
  JSD <- JSD(rbind(vec1_prob, vec2_prob))
  
  MSE <- mean((vec1 - vec2)^2)
  Euc <- sqrt(sum((vec1 - vec2)^2))  # Euclidean Distance
  Manh <- sum(abs(vec1 - vec2))      # Manhattan Distance
  mean_3 <- mean(vec1)
  var_3 <- var(vec1)
  mean_5 <- mean(vec2)
  var_5 <- var(vec2)
  
  return(c(correlation = correlation, cosine_similarity = cosine_similarity,
           MSE = MSE, JSD = JSD, Euc = Euc, Manh = Manh, mean_3 = mean_3, mean_5 = mean_5, var_3 = var_3,
           var_5 = var_5))
}

# --- Helper Functions

# Calculate the row means for a matrix, if there is only one column, return the original vector
my_rowMeans <- function(x) {
  if (!is.null(ncol(x))) {
    if (ncol(x) > 1) {
      return(Matrix::rowMeans(x))
    }
  }
  return(x);
}

# Calculate the row sums for a matrix, if there is only one column, return the original vector
my_rowSums <- function(x) {
  if (!is.null(ncol(x))) {
    if (ncol(x) > 1) {
      return(Matrix::rowSums(x))
    }
  }
  return(x);
}

# Calculate the row variances for a matrix, if there is only one column, return a vector of zeros with the same length as the number of rows
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

# Universal function to perform either mean, sum, or variance calculations for a matrix grouped by a factor. 
group_rowmeans <- function(MAT, group_labs, type=c("mean", "sum", "var")) {
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

# --- Begin the metrics calculations

# List of genes on which the evaluation is performed
genes_for_analysis <- readRDS("./hvf.rds") # Load the list of highly variable features (hvfs) for the respective dataset or the list of bias genes 

# Subset the Seurat object to include only the hvfs or bias genes for the corrected data, as the evaluation is performed on these features.
object <- subset(object, features = genes_for_analysis) # Either this or the alternative below for the uncorrected data

# Split the object by individual and set up the empty lists to store results
indiv_list <- SplitObject(object, split.by = "new_id")

#--- Alternative for the uncorrected data (log normalized)
for (i in seq_along(indiv_list)) {
  indiv_list[[i]] <- ScaleData(indiv_list[[i]])
  indiv_list[[i]] <- subset(indiv_list[[i]], features = genes_for_analysis)
}
#---

# Set up empty lists to store results
result_list <- list()
cell_count_list <- list()
cell_names <- list()

# Calculate metrics per each subject/individual
for(indiv_name in names(indiv_list)) {
  print(indiv_name)
  # Split by assay
  indiv <- indiv_list[[indiv_name]]
  assay_list <- SplitObject(indiv, split.by = "assay")
  
  # Calculate average expression by cell type for 3' assay
  assay <- assay_list[[1]] # Made sure previously that 3' is the first assay in the factor
  assay@meta.data <- droplevels(assay@meta.data) # Drop unused levels if a factor
  
  expr_3 <- GetAssayData(assay, layer = "data") # Choose the necessary layer with corrected/uncorrected data ("counts", "data", "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id

  # Check that the order is the same, produce the error if false
  if(!identical(colnames(expr_3), rownames(assay@meta.data))) {
    stop("Error: The column names in the 3 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  # Count the number of cells per cell type in the dataset and store it in a list for later use. Only include those with more than 5 cells per cell type for the metrics calculations.
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_3")]] <- x
  
  # Calculate average gene expression for each gene by cell_type
  averaged_expression_3 <- group_rowmeans(expr_3, cell_type, "mean") 
  
  # Calculate average expression by cell type for 5' assay
  assay <- assay_list[[2]]
  assay@meta.data <- droplevels(assay@meta.data) 
  
  expr_5 <- GetAssayData(assay, layer = "data") # Choose the necessary layer with corrected/uncorrected data ("counts", "data", "scale.data")
  cell_type <- assay@meta.data$cell_type_ontology_term_id
  
  # Check that the order is the same, produce the error if false
  if(!identical(colnames(expr_5), rownames(assay@meta.data))) {
    stop("Error: The column names in the 5 prime matrix do not match the order of rownames in the Seurat object!")
  }
  
  # Count the number of cells per cell type in the dataset and store it in a list for later use. Only include those with more than 5 cells per cell type for the metrics calculations.
  x <- as.vector(table(cell_type))
  names(x) <- names(table(cell_type))
  cell_count_list[[paste0(indiv_name, "_5")]] <- x
  
  # Calculating average gene expression for each gene by cell_type
  averaged_expression_5 <- group_rowmeans(expr_5, cell_type, "mean")
  
  # Check which cell types match between the 3' and 5'
  # Only get those that have > 5 cells per cell type and match between 3' and 5'
  cell_counts <- cell_count_list[[paste0(indiv_name, "_3")]]
  cells_3 <- names(cell_counts[cell_counts > 5])
  cell_counts <- cell_count_list[[paste0(indiv_name, "_5")]]
  cells_5 <- names(cell_counts[cell_counts > 5])
  
  ident_cell_types <- intersect(cells_3, cells_5)
  cell_names[[indiv_name]] <- ident_cell_types
  
  # Only include genes that are matching
  idx <- intersect(rownames(averaged_expression_3), rownames(averaged_expression_5))
  averaged_expression_3 <- averaged_expression_3[idx,] 
  averaged_expression_5 <- averaged_expression_5[idx,]
  
  # Order matching genes
  averaged_expression_3 <- averaged_expression_3[match(rownames(averaged_expression_5), 
                                                       rownames(averaged_expression_3)), ]
  
  # Create an empty matrix for the output
  metrics <- matrix(0, length(ident_cell_types), 10, 
                    dimnames = list(colnames(averaged_expression_3[,ident_cell_types]),
                                    c("corr", "cos", "MSE", "JSD", "Euc", "Manh",
                                      "mean_3", "mean_5", "var_3", "var_5")))
  
  # For matching cell types between 3' and 5' calculate the metrics
  for(i in ident_cell_types) {
    print(i)
    metrics[i,] <- compute_metrics(averaged_expression_3[,i], averaged_expression_5[,i])
  }
  
  # Results for each patient are stored in a list of lists
  result_list[[paste0(indiv_name, "_metrics")]] <- metrics
}

# Store the components of the results in the unified list for further analysis
result_list[["IndivNames"]] <- names(indiv_list)
result_list[["CellTypes"]] <- cell_names
result_list[["AllCellTypes"]] <- levels(as.factor(object@meta.data$cell_type_ontology_term_id))
result_list[["CellCounts"]] <- cell_count_list


# Save the results
saveRDS(result_list, "./result_X_hvf.rds") # Or if the results were calculated on the bias genes, save as result_X_bias.rds


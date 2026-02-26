# Purpose: Apply Z-transformetion within each subjects/individuals and save result.

# Loading required libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)

# Read the raw preprocessed data for each dataset (Seurat object)
data <- readRDS(paste0(getwd(),"/prep_FINAL.rds"))

# Create a list of top 96 housekeeping genes as identified by Lin et. al. 2019
genes <- c("ENSG00000143761", "ENSG00000136819", "ENSG00000122565", "ENSG00000106153", "ENSG00000175216", "ENSG00000173207", "ENSG00000135940",
           "ENSG00000111775", "ENSG00000115944", "ENSG00000009307", "ENSG00000204435", "ENSG00000088205", "ENSG00000135829", "ENSG00000117395",
           "ENSG00000114942", "ENSG00000175390", "ENSG00000084623", "ENSG00000100129", "ENSG00000110321", "ENSG00000120705", "ENSG00000057608",
           "ENSG00000146066", "ENSG00000122566", "ENSG00000165119", "ENSG00000099783", "ENSG00000153187", "ENSG00000187522", "ENSG00000109971",
           "ENSG00000132305", "ENSG00000108829", "ENSG00000101367", "ENSG00000112118", "ENSG00000242485", "ENSG00000053372", "ENSG00000121579",
           "ENSG00000196531", "ENSG00000253506", "ENSG00000132780", "ENSG00000115053", "ENSG00000186010", "ENSG00000181163", "ENSG00000108518",
           "ENSG00000115762", "ENSG00000099817", "ENSG00000100902", "ENSG00000131467", "ENSG00000187514", "ENSG00000108774", "ENSG00000075785",
           "ENSG00000099901", "ENSG00000134453", "ENSG00000147274", "ENSG00000067560", "ENSG00000106399", "ENSG00000198755", "ENSG00000174748",
           "ENSG00000108107", "ENSG00000156482", "ENSG00000130255", "ENSG00000148303", "ENSG00000161016", "ENSG00000110700", "ENSG00000105372",
           "ENSG00000008988", "ENSG00000138326", "ENSG00000118181", "ENSG00000083845", "ENSG00000137154", "ENSG00000171490", "ENSG00000160633",
           "ENSG00000132432", "ENSG00000183431", "ENSG00000115524", "ENSG00000087365", "ENSG00000169976", "ENSG00000136824", "ENSG00000100028",
           "ENSG00000139343", "ENSG00000143977", "ENSG00000133226", "ENSG00000136450", "ENSG00000116754", "ENSG00000112081", "ENSG00000115875",
           "ENSG00000165283", "ENSG00000023734", "ENSG00000055070", "ENSG00000120948", "ENSG00000070814", "ENSG00000054118", "ENSG00000136527",
           "ENSG00000163811", "ENSG00000136758", "ENSG00000166913", "ENSG00000108953", "ENSG00000065548")

# Subset the 3' and 5' data and split it by patient
Z_tr_3 <- subset(data, subset = assay == "10x 3' v2")
patient_list_3 <- SplitObject(Z_tr_3, split.by = "new_id")
Z_tr_5 <- subset(data, subset = assay == "10x 5' v1")
patient_list_5 <- SplitObject(Z_tr_5, split.by = "new_id")

# Create the Z-transform normalization function
norm <- function(object) {
  for (object_name in names(object)) {
    start <- Sys.time()
    print(paste("Loop for ", object_name))
    print(object[[object_name]])
    object_temp <- object[[object_name]]
    object_temp <- NormalizeData(object_temp, normalization.method = "LogNormalize", scale.factor = 10000)
    
    # Extract the expression data for the specified genes
    expression_data <- FetchData(object_temp, vars = genes)
    # Calculate the mean expression
    mean_expr <- apply(expression_data, 1, FUN = mean)
    mean_ref <- mean(mean_expr)
    
    # Calculate the variance of the expression
    var_ref <- var(mean_expr)
    
    # Extract the expression matrix and perform Z-transformation
    expression_matrix <- as.matrix(GetAssayData(object_temp, layer = "data"))
    normalized_matrix <- sweep(expression_matrix, 2, mean_ref, FUN = "-")  # Subtract the mean
    normalized_matrix <- sweep(normalized_matrix, 2, sqrt(var_ref), FUN = "/")  # Divide by the standard deviation
    
    # Converting the dense matrix to sparse
    threshold <- 1e-5  # Adjust based on inspection
    normalized_matrix[normalized_matrix < threshold & normalized_matrix > -threshold] <- 0
    normalized_matrix <- Matrix(normalized_matrix, sparse = TRUE)
    
    end <- Sys.time()
    print(end - start)

    # Set the normalized expression matrix back to the Seurat object
    object[[object_name]] <- SetAssayData(object_temp, layer = "data", new.data = normalized_matrix)
  }
  return(object)
}

# Perform Z-transform normalization for each subject's 3' and 5' data and merge the results
Z_tr_3 <- norm(patient_list_3)
Z_tr_3 <- merge(Z_tr_3[[1]], y = Z_tr_3[-1])

Z_tr_5 <- norm(patient_list_5)
Z_tr_5 <- merge(Z_tr_5[[1]], y = Z_tr_5[-1])

# Merge the Seurat objects into one
data <- merge(Z_tr_3, Z_tr_5)

# Saving the final object for downstream analysis
saveRDS(data, paste0(getwd(),"/Z_corrected.RDS"))

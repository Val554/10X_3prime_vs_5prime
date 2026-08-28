# Purpose: Preprocessing the datasets for the comparison of the 3' and 5' assays.

# Loading required libraries
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)


### ---- DATASET 1 ------ 
# Sample id - donor_id
# 3' or 5' - assay

ds <- readRDS("./prep.rds")

#-------- 1] Number of Cells/sample, stacked by seq method ------ 
### Can see which samples have both methods and isolate them

# Number of cells per sample
cell_count_per_sample <- ds@meta.data %>%
  group_by(new_id, assay) %>%
  summarise(nCells = n(), .groups = "drop") %>%
  arrange(assay)
ggplot(cell_count_per_sample, aes(x = new_id, y = nCells, fill = assay)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Cells per Sample", x = "Sample ID", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#-------- Isolate the fitting subjects ------ 
ds@meta.data <- unite(ds@meta.data, "new_id", donor_id, sort, sep = "_")

samples <- c("F29_45P", "F30_45P", "F38_45P", "F41_45N", "F41_45P", "F45_45N", "F45_45P", "P1_CD3N", "P1_CD3P")
ds <- subset(ds, new_id %in% samples)
ds@meta.data <- droplevels(ds@meta.data) # drop unused levels if needed

#-------- 2] Mean RNA and Features Count/sample, stacked by seq method ------ 

RNA_per_donor_assay <- ds@meta.data %>%
  group_by(new_id, assay) %>%
  summarize(mean_RNA = mean(nCounts_RNA, na.rm = TRUE)) %>%
  ungroup()
ggplot(RNA_per_donor_assay, aes(x = new_id, y = mean_RNA, fill = assay)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean RNA Count by Assay",
       x = "Assay", y = "Mean RNA Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

features_per_donor_assay <- ds@meta.data %>%
  group_by(assay, new_id) %>%
  summarize(mean_nFeatures = mean(nFeatures_RNA, na.rm = TRUE)) %>%
  ungroup()
ggplot(features_per_donor_assay, aes(x = new_id, y = mean_nFeatures, fill = assay)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Feature Count by Assay for Selected Donors",
       x = "Assay", y = "Mean Feature Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

#-------- 3] Number of Cells/sample, stacked by seq method by cell type ------ 
cells_by_cell_type_and_technology <- ds@meta.data %>%
  group_by(cell_type, assay, new_id) %>%
  summarize(nCells = n()) %>%
  ungroup()
ggplot(cells_by_cell_type_and_technology, aes(x = cell_type, y = nCells, fill = assay)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Cells by Cell Type and Sequencing Technology for Selected Donors",
       x = "Cell Type", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid( ~ new_id, scales = "free_y")

#-------- Preprocess ------ 

split_ds <- SplitObject(ds, split.by = "new_id")

clean_seurat <- function(seurat_obj) {
  
  # Ensure that the metadata is accessible
  meta_data <- seurat_obj@meta.data
  
  # Group the metadata by donor_id and cell_type, then count the number of cells per group
  cell_counts <- meta_data %>%
    group_by(new_id, cell_type, assay) %>%
    summarise(cell_count = n(), .groups = 'drop')
  
  # Create a vector of cell types to keep (those with >= 50 cells per donor)
  not_valid_cell_types <- cell_counts %>%
    filter(cell_count < 50) %>%
    group_by(new_id, cell_type) %>%
    # Keep cell types only if they meet the condition in all assays
    filter(n_distinct(assay) == n_distinct(assay[cell_count < 50])) %>%
    pull(cell_type) %>%
    unique()
  
  # Filter the Seurat object to only include cells from valid cell types
  seurat_obj_clean <- subset(seurat_obj, subset = cell_type %!in% not_valid_cell_types)
  return(seurat_obj_clean)
}

for(i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
  split_ds[[i]]@meta.data <- droplevels(split_ds[[i]]@meta.data)
  print(levels(as.factor(split_ds[[i]]@meta.data$cell_type)))
}

ds_clean <- merge(split_ds[[1]], y = split_ds[-1])
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

ds_clean@meta.data[c("cell_type", "cell_type_ontology_term_id")]

saveRDS(ds_clean, "./prep_FINAL.rds")


### ---- DATASET 2 ------ 
# Sample id - donor_id
# 3' or 5' - assay

ds <- readRDS("./dataset2.rds")

# Visualizations are done similarly to dataset 1.

#-------- Isolate the fitting subjects ------ 
samples <- c("F32", "F34", "F38", "F41", "F45")
ds <- subset(ds, (donor_id %in% samples) & (Sort_id %in% c("CD45P", "CD45N"))) 
ds@meta.data <- droplevels(ds@meta.data) # drop unused levels if needed

######### Isolate tissues with both technologies
# Create multiple object for different tissues and sets of donors

spleen_ds <- subset(ds, donor_id %in% c("F45") & tissue == "spleen")
spleen_ds@meta.data <- droplevels(spleen_ds@meta.data)

liver_ds <- subset(ds, donor_id %in% c("F32", "F34", "F38", "F41", "F45") & tissue == "liver")
liver_ds@meta.data <- droplevels(liver_ds@meta.data)

kidney_ds <- subset(ds, donor_id %in% c("F41") & tissue == "kidney")
kidney_ds@meta.data <- droplevels(kidney_ds@meta.data)

thymus_ds <- subset(ds, donor_id %in% c("F38", "F41", "F45") & tissue == "thymus")
thymus_ds@meta.data <- droplevels(thymus_ds@meta.data)

bm_ds <- subset(ds, donor_id %in% c("F41", "F45") & tissue == "bone marrow")
bm_ds@meta.data <- droplevels(bm_ds@meta.data)

#-------- Preprocess (by tissue) ------ 

# Reorder assay levels such that 3 is first and 5 is second:
desired_order <- c("10x 3' v2", "10x 5' v1")
ds@meta.data$assay <- factor(ds@meta.data$assay, levels = desired_order)

# Drop all unused levels in the meta data after the subset
spleen_ds@meta.data <- droplevels(spleen_ds@meta.data)
spleen_ds@meta.data$new_id <- spleen_ds@meta.data$donor_id
liver_ds@meta.data <- droplevels(liver_ds@meta.data)
liver_ds@meta.data$new_id <- liver_ds@meta.data$donor_id
kidney_ds@meta.data <- droplevels(kidney_ds@meta.data)
kidney_ds@meta.data$new_id <- kidney_ds@meta.data$donor_id
thymus_ds@meta.data <- droplevels(thymus_ds@meta.data)
thymus_ds@meta.data$new_id <- thymus_ds@meta.data$donor_id
bm_ds@meta.data <- droplevels(bm_ds@meta.data)
bm_ds@meta.data$new_id <- bm_ds@meta.data$donor_id

saveRDS(spleen_ds, "./prep_spleen.rds")
saveRDS(liver_ds, "./prep_liver.rds")
saveRDS(kidney_ds, "./prep_kidney.rds")
saveRDS(thymus_ds, "./prep_thymus.rds")
saveRDS(bm_ds, "./prep_bm.rds")

# Further cleanup of the cell types is needed for the metrics calculations. 
# We remove the cell types with < 50 cells per subject/ assay, and then merge some of the cell types together to have more cells per cell type.

ds <- readRDS("./prep_tissue.rds") # Choose tissue as needed
ds@meta.data$new_id <- paste0(ds@meta.data$donor_id, "_",ds@meta.data$Sort_id)

filtered_cells <- c("unknown", "erythrocyte", "megakaryocyte","myofibroblast cell", "epithelial cell",
                    "glial cell", "neuron", "osteoclast", "skeletal muscle satellite cell")

ds <- subset(ds, cell_type %!in% filtered_cells)

ds@meta.data <- ds@meta.data %>%
  mutate(
    merged_cell_type = case_when(
      cell_type %in% c("B-1 B cell", "B-2 B cell", "immature B cell", "mature B cell") ~ "B cell",
      cell_type %in% c("CD8-alpha-alpha-positive, alpha-beta intraepithelial T cell", "double negative thymocyte",
                       "double-positive, alpha-beta thymocyte",
                       "naive thymus-derived CD4-positive, alpha-beta T cell",
                       "naive thymus-derived CD8-positive, alpha-beta T cell",
                       "regulatory T cell", "T cell") ~ "T cell",
      cell_type %in% c("common dendritic progenitor", "dendritic cell", "plasmacytoid dendritic cell",
                       "pre-conventional dendritic cell") ~ "dendritic cell",
      cell_type %in% c("common myeloid progenitor", "early lymphoid progenitor",
                       "granulocyte monocyte progenitor cell",
                       "hematopoietic multipotent progenitor cell",
                       "hematopoietic stem cell", "megakaryocyte-erythroid progenitor cell") ~ "hematopoietic progenitor cell",
      cell_type %in% c("fraction A pre-pro B cell", "large pre-B-II cell",
                       "late pro-B cell", "pro-B cell",
                       "small pre-B-II cell") ~ "B cell progenitors",
      cell_type %in% c("granulocyte", "myelocyte",
                       "neutrophil", "promyelocyte") ~ "granulocyte",
      cell_type %in% c("group 3 innate lymphoid cell", "group 2 innate lymphoid cell", "innate lymphoid cell") ~ "innate lymphoid cell",
      cell_type %in% c("Kupffer cell", "macrophage",
                       "monocyte", "promonocyte") ~ "monocyte/macrophage",
      cell_type %in% c("smooth muscle cell", "vascular associated smooth muscle cell") ~ "smooth muscle cell",
      TRUE ~ cell_type  # Default case for any unmatched cell types
    )
  )

split_ds <- SplitObject(ds, split.by = "new_id")

for (i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
}

ds_clean <- merge(split_ds[[1]], y = split_ds[-1])
ds_clean@meta.data <- droplevels(ds_clean@meta.data)
ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$merged_cell_type

saveRDS(ds_clean, "./prep_spleen_FINAL.rds")
saveRDS(ds_clean, "./prep_liver_FINAL.rds")
saveRDS(ds_clean, "./prep_kidney_FINAL.rds")
saveRDS(ds_clean, "./prep_thymus_FINAL.rds")
saveRDS(ds_clean, "./prep_bm_FINAL.rds")

### ---- DATASET 3 ------ 

# Sample id - donor_id
# 3' or 5' - assay

# Load the data for the different cell types separately, as they are stored in different files. 
# Then merge them together for the visualizations and preprocessing.
ds <- readRDS("./hepatocyte1.rds")
ds$cell <- "Hepatocyte-1"
ds2 <- readRDS("./hepatocyte2.rds")
ds2$cell <- "Hepatocyte-2"
B <- readRDS("./Bcells.rds")
B$cell <- "B cell"
chola <- readRDS("./chola.rds")
chola$cell <- "Cholangiocyte"
endo <- readRDS("./endothelial.rds")
endo$cell <- "Endothelial"
macro <- readRDS("./macrophage.rds")
macro$cell <- "Macrophage"
nk <- readRDS("./NK.rds")
nk$cell <- "NK cells"
stell <- readRDS("./stellate.rds")
stell$cell <- "Stellate cell"

ds <- merge(ds, y = c(ds2, B, chola, endo, macro, nk, stell))

# Visualizations are done similarly to dataset 1.

#-------- Isolate the fitting subjects ------ 
samples <- c("C58", "C59", "C61", "C64", "C70")
ds <- subset(ds, donor_id %in% samples) 
ds@meta.data <- droplevels(ds@meta.data) 

#-------- Preprocess ------ 
levels(as.factor(ds@meta.data$donor_id))
ds <- subset(ds, donor_id %in% c("C58","C59","C64"))
ds@meta.data <- droplevels(ds@meta.data)

# Cleaning the cell types:
filtered_cells <- rownames(ds@meta.data[!(ds@meta.data$cell_type %in% 
                                            c("unknown", "vascular associated smooth muscle cell", 
                                              "erythroblast", "erythroid lineage cell", 
                                              "B cell", "mature B cell", 
                                              "naive B cell", "plasmablast", "lymphocyte",
                                              "intrahepatic cholangiocyte", "liver dendritic cell",
                                              "hepatic stellate cell")), ])
ds <- subset(ds, cells = filtered_cells)

ds@meta.data <- ds@meta.data %>%
  mutate(
    merged_cell_type = case_when(
      cell_type %in% c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "natural killer cell") ~ "T/NK cell",
      cell_type %in% c("centrilobular region hepatocyte", "midzonal region hepatocyte",
                       "periportal region hepatocyte", "hepatocyte") ~ "hepatocyte",
      cell_type %in% c("classical monocyte", "Kupffer cell", "monocyte",
                       "macrophage") ~ "monocyte/macrophage",
      cell_type %in% c("endothelial cell of artery", "endothelial cell of hepatic sinusoid",
                       "endothelial cell of pericentral hepatic sinusoid",
                       "endothelial cell of periportal hepatic sinusoid",
                       "vein endothelial cell") ~ "endothelial",
      TRUE ~ cell_type  # Default case for any unmatched cell types
    )
  )


split_ds <- SplitObject(ds, split.by = "donor_id")

for(i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
  levels(as.factor(split_ds[[i]]@meta.data$merged_cell_type))
}

ds_clean<- merge(split_ds[[1]], y = split_ds[-1])
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

ds_clean@meta.data$new_id <- ds_clean@meta.data$donor_id
ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$merged_cell_type

saveRDS(ds_clean, "./prep_FINAL_v5.rds")

### ---- DATASET 4 ------ 
# Sample id - donor_id
# 3' or 5' - assay

ds3_1 <- readRDS("./10X_3-rep1.rds")
ds3_2 <- readRDS("./10X_3-rep2.rds")
ds5_1 <- readRDS("./10X_5-rep1.rds")
ds5_2 <- readRDS("./10X_5-rep2.rds")

ds <- merge(ds3_1, y = c(ds3_2, ds5_1, ds5_2))

# Visualizations are done similarly to dataset 1.

#-------- Preprocess ------ 
ds <- readRDS("./merged.rds")
ds <- subset(ds, cell_type_ontology_term_id != "unknown")
ds@meta.data <- droplevels(ds@meta.data)

# Create a new_id for futire processing 
ds@meta.data$new_id <- "RG1237"
saveRDS(ds, "./prep_FINAL_v5.rds")

### ---- DATASET 5 ------ 
# Sample id - donor_id
# 3' or 5' - assay

ds <- readRDS("./data.rds")
ds <- subset(ds, donor_id %in% c("Leader_Meard_2021_522"))
ds@meta.data <- droplevels(ds@meta.data) 

# Visualizations are done similarly to dataset 1.

#-------- Preprocess ------ 

# Merging cell types:
ds@meta.data <- ds@meta.data %>%
  mutate(
    merged_cell_type = case_when(
      cell_type %in% c("dendritic cell", "conventional dendritic cell", "CD1c-positive myeloid dendritic cell") ~ "dendritic cell",
      cell_type %in% c("epithelial cell of lung", "club cell",
                       "multiciliated epithelial cell", "pulmonary alveolar type 1 cell",
                       "pulmonary alveolar type 2 cell") ~ "lung epithelial cell",
      cell_type %in% c("bronchus fibroblast of lung",
                       "fibroblast of lung") ~ "fibroblasts",
      cell_type %in% c("smooth muscle cell", "pericyte") ~ "smooth muscle/pericyte cell",
      cell_type %in% c("endothelial cell of lymphatic vessel", "capillary endothelial cell",
                       "pulmonary artery endothelial cell", "vein endothelial cell") ~ "endothelial cell",
      TRUE ~ cell_type  # Default case for any unmatched cell types
    )
  )

ds@meta.data$new_id <- paste0(ds@meta.data$donor_id, "_", ds@meta.data$origin)

split_ds <- SplitObject(ds, split.by = "new_id")

for(i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
  levels(as.factor(split_ds[[i]]@meta.data$merged_cell_type))
}

ds_clean<- merge(split_ds[[1]], y = split_ds[-1])
ds_clean <- JoinLayers(ds_clean)
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$merged_cell_type
levels(as.factor(ds_clean@meta.data$cell_type_ontology_term_id))
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

saveRDS(ds_clean, "./prep_FINAL.rds")

### ---- DATASET 6 ------ 
# Sample id - donor_id
# 3' or 5' - assay

ds <- readRDS("./data.RDS")
ds <- subset(ds, donor_id %in% c("F30", "F41", "F45"))

# Visualizations are done similarly to dataset 1.

#-------- Preprocess ------ 

# Merging cell types:
ds@meta.data <- ds@meta.data %>%
  mutate(
    merged_cell_type = case_when(
      cell_type %in% c("precursor B cell", "immature B cell", "naive B cell") ~ "B cell",
      cell_type %in% c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
                       "regulatory T cell", "mature NK T cell",
                       "CD16-negative, CD56-bright natural killer cell, human",
                       "natural killer cell",
                       "immature natural killer cell") ~ "T/NK cell",
      cell_type %in% c("common dendritic progenitor", "dendritic cell") ~ "dendritic cell",
      cell_type %in% c("fraction A pre-pro B cell", "pro-B cell") ~ "B cell progenitors",
      cell_type %in% c("myelocyte", "neutrophil",
                       "promyelocyte", "eosinophil", "basophil") ~ "granulocyte",
      cell_type %in% c("endothelial cell of sinusoid", "endothelial tip cell", "endothelial cell") ~ "endothelial cell",
      cell_type %in% c("muscle cell", "muscle precursor cell", "precursor cell",
                       "chondrocyte", "fibroblast", "myofibroblast cell") ~ "stromal/mesenchymal lineage",
      cell_type %in% c("primitive red blood cell", "erythrocyte") ~ "RBCs",
      TRUE ~ cell_type  # Default case for any unmatched cell types
    )
  )

split_ds <- SplitObject(ds, split.by = "donor_id")

for(i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
  levels(as.factor(split_ds[[i]]@meta.data$merged_cell_type))
}

ds_clean<- merge(split_ds[[1]], y = split_ds[-1])
ds_clean@meta.data <- droplevels(ds_clean@meta.data)
ds_clean <- JoinLayers(ds_clean)

ds_clean@meta.data$new_id <- ds_clean@meta.data$donor_id
ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$merged_cell_type

# Clean up unbalanced cell types
ds_clean <- subset(ds_clean, new_id != "F30" | 
                     cell_type_ontology_term_id %!in% c("CD14-positive monocyte","granulocyte",
                                                        "megakaryocyte", "progenitor cell"))
ds_clean <- subset(ds_clean, new_id != "F41" | 
                     cell_type_ontology_term_id %!in% c("CD14-positive monocyte","granulocyte", "RBCs"))
ds_clean <- subset(ds_clean, new_id != "F45" | 
                     cell_type_ontology_term_id %!in% c("CD14-positive monocyte","granulocyte",
                                                        "promonocyte"))
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

saveRDS(ds_clean, "./prep_FINAL.rds")

## ----- Validation Datasets ----- ##

### ---- DATASET 1 ------ 
# Sample id - donor_id
# 3' or 5' - assay
# Cell Type - broad_cell_type

ds <- readRDS("./data.RDS")
ds <- subset(ds, assay %in% c("10x 3' v3", "10x 5' v1"))

# Visualizations are done similarly to dataset 1.

#-------- Preprocess ------ 

split_ds <- SplitObject(ds, split.by = "donor_id")

for(i in 1:length(split_ds)) {
  split_ds[[i]] <- clean_seurat(split_ds[[i]])
  levels(as.factor(split_ds[[i]]@meta.data$merged_cell_type))
}

ds_clean<- merge(split_ds[[1]], y = split_ds[-1])
ds_clean@meta.data <- droplevels(ds_clean@meta.data)
ds_clean <- JoinLayers(ds_clean)

ds_clean@meta.data$new_id <- ds_clean@meta.data$donor_id
ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$broad_cell_type
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

saveRDS(ds_clean, "./prep_FINAL.rds")

### ---- DATASET 2 ------ 
# Sample id - donor_id
# 3' or 5' - assay
# Cell Type - broad_cell_type

ds <- readRDS("./data.RDS")

ds$sort_tis <- paste0(ds$tissue, "_", ds$Sort_id)
ds$don_sort_tis <- paste0(ds$donor_id, "_", ds$sort_tis)
samples <- c("F45_kidney_CD45N", "F45_skin of body_CD45P", "F38_skin of body_CD45P", "F45_liver_CD45P")
ds <- subset(ds, (don_sort_tis %in% samples)) 
ds@meta.data <- droplevels(ds@meta.data) # drop unused levels if needed

# Visualizations are done similarly to dataset 1.

#-------- Preprocess ------ 
# Preprocessing and cell type merging is done similarly to dataset 2
#...

ds_clean@meta.data$new_id <- ds_clean@meta.data$donor_id
ds_clean@meta.data$orig_cell_type_ontology_term_id <- ds_clean@meta.data$cell_type_ontology_term_id
ds_clean@meta.data$cell_type_ontology_term_id <- ds_clean@meta.data$merged_cell_type
ds_clean@meta.data <- droplevels(ds_clean@meta.data)

saveRDS(subset(ds_clean, tissue == "skin of body"), "./prep_FINAL_skin.rds")
saveRDS(subset(ds_clean, tissue == "kidney"), "./prep_FINAL_kidney.rds")
saveRDS(subset(ds_clean, tissue == "liver"), "./prep_FINAL_liver.rds")


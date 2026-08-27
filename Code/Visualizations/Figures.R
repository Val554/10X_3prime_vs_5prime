# Purpose: Plotting all the figures 

# Loading required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(presto)
library(ggbreak)
library(Seurat)
library(biomaRt)
library(forecast)

### ----- Figure 1: ----- 

#### A) ---- CDFs for the 3' vs 5' ----
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-100.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_000001.rds")
gene_table <- table(gene_table)
gene_table <- gene_table[length(gene_table):1]

cum_counts <- cumsum(gene_table)

# Plot in reverse order: from high to low
df <- data.frame(
  x = as.numeric(names(cum_counts)),
  y = cum_counts
)

# Lineplot 
ggplot(df, aes(x = x, y = y)) +
  geom_step(color = "#2C7BE5", size = 1.2) +
  scale_x_reverse(breaks = df$x) +
  labs(
    title = "Cumulative Distribution (CDF) of Counts",
    x = "Value",
    y = "Cumulative Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Barplot with a break 

ggplot(df, aes(x = x, y = y)) +
  geom_col(fill = "#9BC3E6", width = 0.85,color = "black") +
  scale_x_reverse(breaks = df$x) +
  labs(
    title = "Cumulative Distribution (CDF) of Counts",
    x = "Value",
    y = "Cumulative Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.title       = element_blank(),
    #axis.text.x       = element_text(angle = 90,face = "bold", size = 18),
    panel.grid.major.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border    = element_blank(),
    axis.line       = element_line(color = "black", linewidth = 1)  # <-- axes restored
  ) #+
# scale_y_break() +
#scale_y_break(c(1750,3500, 3750,9000)) +# for 50
#scale_y_break(c(1500,5500)) + # for 100
#scale_y_continuous(
#  breaks = seq(0, 10000, by = 500)
#)


#### B) ---- UMAPs ----

# PREPARE DATA
data <- readRDS("./Scaled_LogNorm.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-100.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_000001.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50
common_ind <- c(names(gene_table[gene_table %in% c(3:35)])) # 600 :100
common_ind <- c(names(gene_table[gene_table %in% c(15:35)])) # 3000: 0.000001

data_merged <- merge(data[[1]], y = data[-1])
data_merged <- JoinLayers(data_merged)

data_merged <- RunPCA(data_merged)
ElbowPlot(data_merged, ndims = 50)
data_merged <- RunUMAP(data_merged, dims = 1:30)

# Reduce the size of the object for easier plotting and save it for later use.
data_merged1 <- SetAssayData(data_merged, layer = "scale.data", new.data = NULL)
data_merged1 <- SetAssayData(data_merged1, layer = "data", new.data = NULL)
data_merged1 <- data_merged1[(rownames(data_merged1) %in% common_ind[50:100]), ]
saveRDS(data_merged1, "./Scaled_LogNorm_PLOT.rds")

data_merged1 <- data_merged[!(rownames(data_merged) %in% common_ind), ]

data_merged1 <- split(data_merged1, f = data_merged1@meta.data$new_id)
data_merged1 <- NormalizeData(data_merged1)
data_merged1 <- FindVariableFeatures(data_merged1, selection.method = "vst", nfeatures = 3000)
data_merged1 <- JoinLayers(data_merged1)

data_merged_list <- SplitObject(data_merged1, split.by = "new_id")
data_merged_list <- lapply(data_merged_list, function(x){
  x <- ScaleData(x)
})

data_merged1 <- merge(data_merged_list[[1]], y = data_merged_list[-1])
data_merged1 <- JoinLayers(data_merged1)

data_merged1 <- RunPCA(data_merged1)
data_merged1 <- RunUMAP(data_merged1, dims = 1:30)

data_merged1 <- SetAssayData(data_merged1, layer = "scale.data", new.data = NULL)
data_merged1 <- SetAssayData(data_merged1, layer = "data", new.data = NULL)
data_merged1 <- data_merged1[(rownames(data_merged1) %in% c("ENSG00000144290")), ]
saveRDS(data_merged1, "./Scaled_LogNorm_PLOT_no600.rds")
saveRDS(data_merged1, "./Scaled_LogNorm_PLOT_no3000.rds")


#### ------ PLOT ----######
original <- readRDS("./Scaled_LogNorm_PLOT.rds")
original <- readRDS("./Scaled_LogNorm_PLOT_no600.rds") # Adjust the number as needed

DimPlot(original, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "right"                # remove legend
  )
DimPlot(original, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

DimPlot(original, reduction = "umap", group.by = "cell_type_ontology_term_id")  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

### Mixing Metrics for original
original <- merge(original[[1]], y = original[-1])
original <- JoinLayers(original)
original <- RunPCA(original)
original <- RunUMAP(original, dims = 1:30)

mix_assay <- MixingMetric(
  original,
  grouping.var = "assay",
  reduction = "pca",
  dims = 1:30,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_assay)

mix_assay <- MixingMetric(
  original,
  grouping.var = "assay",
  reduction = "umap",
  dims = 1:2,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_assay)

mix_cell <- MixingMetric(
  original,
  grouping.var = "cell_type_ontology_term_id",
  reduction = "umap",
  dims = 1:2,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_cell)

reduced <- readRDS("./Scaled_LogNorm10000_PLOT_no3000.rds") # Adjust the number as needed

DimPlot(reduced, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4")) +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

DimPlot(reduced, reduction = "umap", group.by = "cell_type_ontology_term_id")  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

### Mixing Metrics for reduced

mix_assay <- MixingMetric(
  reduced,
  grouping.var = "assay",
  reduction = "umap",
  dims = 1:2,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_assay)

mix_assay <- MixingMetric(
  reduced,
  grouping.var = "assay",
  reduction = "pca",
  dims = 1:30,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_assay)

mix_cell <- MixingMetric(
  reduced,
  grouping.var = "cell_type_ontology_term_id",
  reduction = "umap",
  dims = 1:2,
  k = 5,
  max.k = 150,
  eps = 0,
  verbose = TRUE
)
median(150 - mix_cell)

#### C) ---- Cos Sim ----
set <- "" # Specify the dataset number for the file names
result_list <- readRDS(paste0("./DS",set,"_cos_",800,"_SCALED.rds"))
result_list_rand <- readRDS(paste0("./DS",set,"_cos_",800,"_rand_prop_SCALED.rds"))
result_list_rand2 <- readRDS(paste0("./DS",set,"_cos_",800,"_rand_prop_2_SCALED.rds"))
result_list_rand3 <- readRDS(paste0("./DS",set,"_cos_",800,"_rand_prop_3_SCALED.rds"))

len <- length(result_list[["AllCellTypes"]])
# For low vs high:
cell_avg <- list(
  Gene_List = setNames(rep(0, len), result_list[["AllCellTypes"]]),
  Random = setNames(rep(0, len), result_list[["AllCellTypes"]]),
  Random_2 = setNames(rep(0, len), result_list[["AllCellTypes"]]),
  Random_3 = setNames(rep(0, len), result_list[["AllCellTypes"]])
)
for(indiv_name in result_list[["IndivNames"]]) {
  for(cell in result_list[["CellTypes"]][[indiv_name]]) {
    cell_avg$Gene_List[[cell]] <- cell_avg$Gene_List[[cell]] + result_list[[indiv_name]][cell]
    cell_avg$Random[[cell]] <- cell_avg$Random[[cell]] + result_list_rand[[indiv_name]][cell]
    cell_avg$Random_2[[cell]] <- cell_avg$Random_2[[cell]] + result_list_rand2[[indiv_name]][cell]
    cell_avg$Random_3[[cell]] <- cell_avg$Random_3[[cell]] + result_list_rand3[[indiv_name]][cell]
  }
}

all_cell_types <- unlist(lapply(result_list[["CellTypes"]], unique))
cell_type_counts <- table(all_cell_types)

# Divide each metric value by its corresponding frequency
adjusted_metrics <- lapply(cell_avg, function(metric) {
  metric / cell_type_counts[names(metric)]  # Divide by corresponding frequencies
})

# Convert the list to a data frame for ggplot
data_df <- adjusted_metrics %>%
  purrr::imap_dfr(~ data.frame(name = names(.x), value = .x, vector = .y))

data_df$vector <- factor(data_df$vector, levels = c("Gene_List", "Random", "Random_2", "Random_3"))

col <- scale_fill_manual(values = c(
  "Gene_List" = "#1f78b4",
  "Random" = "#b2df8a",
  "Random_2" = "#cab2d6",
  "Random_3" = "#fdbf6f"
))

ggplot(data_df, aes(x = vector, y = value.Freq, fill = vector)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(x = "Genes", y = "Scaled Score", shape = "Cell Type") +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9,11,12,13, 14,15)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  scale_y_continuous(
    breaks = seq(-0.75, 1, by = 0.25)) +
  guides(fill = "none")


### ----- Figure 2: ----- 

split_df <- readRDS("./updated_FINAL_%_all_DSandTechniques_to_baseline.rds")
split_df <- readRDS("./800_FINAL_%_all_DSandTechniques_to_baseline.rds")

## Plotting
split_df$corr$norm <- factor(split_df$corr$norm, levels = c("ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                            "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$cos$norm <- factor(split_df$cos$norm, levels = c("ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                          "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$Euc$norm <- factor(split_df$corr$norm, levels = c("ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$MSE$norm <- factor(split_df$corr$norm, levels = c("ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$JSD$norm <- factor(split_df$corr$norm, levels = c("ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))

pref <- "800_"

col <- scale_fill_manual(values = c(
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))
# a) Corr
ggplot(split_df$corr, aes(x = norm, y = average_impr, fill = norm)) +
  geom_boxplot(outliers = FALSE, notch = TRUE, notchwidth = 0.8) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Correlation Coefficient", x = "Normalization", y = "Relative Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

# b) Cos
ggplot(split_df$cos, aes(x = norm, y = average_impr, fill = norm)) +
  geom_boxplot(outliers = FALSE, notch = TRUE, notchwidth = 0.8) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Cosine Similarity", x = "Normalization", y = "Relative Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  guides(fill = "none")

# c) MSE
ggplot(split_df$MSE, aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Mean Squared Error", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  guides(fill = "none")  +
  scale_y_break(c(16000,48000))

ggplot(subset(split_df$MSE, subset = norm !="Z-transform" & (norm !="SCVI") & norm !="scArches"), aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Mean Squared Error", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

# d) JSD
ggplot(split_df$JSD, aes(x = norm, y = average_impr0, fill = norm)) +
  geom_boxplot(outliers = FALSE,  notch = TRUE, notchwidth = 0.8) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Jensen-Shannon Divergence", x = "Normalization", y = "Relative Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

# e) Euc
ggplot(split_df$Euc, aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE, notch = TRUE, notchwidth = 0.8) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Euclidean Distance", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(subset(split_df$Euc, subset = norm !="Z-transform" ), aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Euclidean Distance", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

# f) Manh
ggplot(split_df$Manh, aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Manhattan Distance", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(subset(split_df$Manh, subset = norm !="Z-transform"), aes(x = norm, y = average_per, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "Manhattan Distance", x = "Normalization", y = "Percent Change (%)", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")


### ----- Figure 3: ----- 

#### A) ---- UMAPs ----

# DS1: F29_45P
# DS2sp: F45_CD45P_spleen
# DS3: C64
# DS4: RG1237
# DS5: Leader_Merad_2021_522_tumor_primary
# DS6: F30

bold_cols <- c(
  "#E41A1C", "#ff7f00", "#377EB8", "#984EA3",
  "#4DAF4A", "#2F5B73", "#fb9a99", "#F781BF", "#fdbf6f"
)

col <- scale_fill_manual(values = c(
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

obj <- readRDS("./Scaled_X_indiv_PLOT.rds") # This is done for all of the normalization techniques - X

indiv <- "" # Specify

DimPlot(obj[[indiv]], reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    #axis.line = element_line(color = "black", linewidth = 1),
    axis.line = element_blank(), # no axis!!!
    legend.position = "right"                # remove legend
  )
DimPlot(obj, reduction = "umap", group.by = "cell_type_ontology_term_id")  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "right"                # remove legend
  )+
  scale_color_manual(values = bold_cols)

DimPlot(obj[[indiv]], reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

DimPlot(obj[[indiv]], reduction = "umap", group.by = "cell_type_ontology_term_id")  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )+
  scale_color_manual(values = bold_cols)


#### B) ---- MM ----
tis <- "" # Specify

# By Assay
melted_data <- readRDS(paste0("./FINAL_MM_assay",tis,".rds"))
DS <- "" # Specify

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")


# By Cell
melted_data <- readRDS(paste0("./FINAL_MM_cell",tis,".rds"))

ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

#### Supp) ---- Jaccard ----
tis <- "" # Specify

JS_df <- readRDS(paste0("./FINAL_JS_sum",tis,".rds"))
DS <- "" # Specify

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

ggplot(JS_df, aes(x = norm, y = Average_Jaccard_Score, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12,13)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(JS_df, aes(x = norm, y = Average_Jaccard_Score, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12,13)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

### ----- Figure 4: ----- 

#### A) ---- UseCase ----

indiv <- "_intersect_noFC" # "" for merged, "_indiv" for indiv, "_union" or "_intersect" => we used intersect OR NEW: "_intersect_noFC"
DS <- "" # Specify

log <- readRDS(paste0("./",DS,"_logN_MCC_25",indiv,".rds"))
fast <- readRDS(paste0("./",DS,"_fast_MCC_25",indiv,".rds"))
scanorama <- readRDS(paste0("./",DS,"_scanorama_MCC_25",indiv,".rds"))
limma <- readRDS(paste0("./",DS,"_limma_MCC_25",indiv,".rds"))
m3drop <- readRDS(paste0("./",DS,"_M3Drop_MCC_25",indiv,".rds"))
scvi <- readRDS(paste0("./",DS,"_scvi_MCC_25",indiv,".rds"))

# Combine into a list
lst <- list(
  Log_Normalize = log,
  limma = limma,
  fastMNN = fast,
  Scanorama = scanorama,
  SCVI = scvi,
  M3Drop = m3drop
)

# Convert to long format
df <- map_dfr(names(lst), \(n) {
  tibble(
    group = n,
    name = names(lst[[n]]),
    value = lst[[n]]
  )
})

if (DS =="DS2lv"){
  df$name[df$name == "F32_CD45P_liver"] <- "F32_CD45P"
  df$name[df$name == "F34_CD45P_liver"] <- "F34_CD45P"
  df$name[df$name == "F38_CD45P_liver"] <- "F38_CD45P"
  df$name[df$name == "F41_CD45P_liver"] <- "F41_CD45P"
}

df$group <- factor(df$group, levels = c("Log_Normalize", "limma", "fastMNN", "SCVI", "M3Drop"))

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "limma" = "#1f78b4",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "SCVI" = "#fdbf6f",
  "M3Drop" = "#e31a1c"
))

ggplot(df,  aes(x = group, y = value, fill = group)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "MCC", x = "Normalization", y = "MCC Score (%)", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  scale_y_continuous(
    limits = c(-2, 100),
    breaks = seq(0, 100, 25)
  ) +
  guides(fill = "none")


log <- readRDS(paste0("./",DS,"_logN_MCC_100",indiv,".rds"))
fast <- readRDS(paste0("./",DS,"_fast_MCC_100",indiv,".rds"))
scanorama <- readRDS(paste0("./",DS,"_scanorama_MCC_100",indiv,".rds"))
limma <- readRDS(paste0("./",DS,"_limma_MCC_100",indiv,".rds"))
m3drop <- readRDS(paste0("./",DS,"_M3Drop_MCC_100",indiv,".rds"))
scvi <- readRDS(paste0("./",DS,"_scvi_MCC_100",indiv,".rds"))

# Combine into a list
lst <- list(
  Log_Normalize = log,
  limma = limma,
  fastMNN = fast,
  Scanorama = scanorama,
  SCVI = scvi,
  M3Drop = m3drop
)

# Convert to long format
df <- map_dfr(names(lst), \(n) {
  tibble(
    group = n,
    name = names(lst[[n]]),
    value = lst[[n]]
  )
})

if (DS =="DS2lv"){
  df$name[df$name == "F32_CD45P_liver"] <- "F32_CD45P"
  df$name[df$name == "F34_CD45P_liver"] <- "F34_CD45P"
  df$name[df$name == "F38_CD45P_liver"] <- "F38_CD45P"
  df$name[df$name == "F41_CD45P_liver"] <- "F41_CD45P"
}

df$group <- factor(df$group, levels = c("Log_Normalize", "limma", "fastMNN", "SCVI", "M3Drop"))

ggplot(df,  aes(x = group, y = value, fill = group)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid") + 
  labs(title = "MCC", x = "Normalization", y = "MCC Score (%)", shape = "Individual") +
  #scale_fill_brewer(type = "qual", palette = 3) +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25)
  ) +
  guides(fill = "none")


#### Expression ##### 
DS <- "" # Specify

truth <- readRDS(paste0("./",DS,"_intersect_truth_labels_intersect_noFC.rds"))

fast <- readRDS(paste0("./",DS,"_merged_fast_25.rds"))
limma <- readRDS(paste0("./",DS,"_merged_limma_25.rds"))
logN <- readRDS(paste0("./",DS,"_merged_logN_25.rds"))
M3 <- readRDS(paste0("./",DS,"_merged_M3Drop_25.rds"))
scvi <- readRDS(paste0("./",DS,"_merged_scvi_25.rds"))
scanorama <- readRDS(paste0("./",DS,"_merged_scanorama_25.rds"))

fast_MM <- readRDS(paste0("./",DS,"_fast_test_labels_25_noFC.rds"))
limma_MM <- readRDS(paste0("./",DS,"_limma_test_labels_25_noFC.rds"))
logN_MM <- readRDS(paste0("./",DS,"_logN_test_labels_25_noFC.rds"))
M3_MM <- readRDS(paste0("./",DS,"_M3Drop_test_labels_25_noFC.rds"))
scvi_MM <- readRDS(paste0("./",DS,"_scvi_test_labels_25_noFC.rds"))
scanorama_MM <- readRDS(paste0("./",DS,"_scanorama_test_labels_25_noFC.rds"))

# Summary of DE direction
df <- data.frame(method = "Truth", direction = c("up", "down"), value = c(sum(truth$dir_sig == "up"), sum(truth$dir_sig == "down")), indiv = "all")

#fast
table(fast_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(fast_MM, function(x) {
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("fastMNN",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

#limma
table(limma_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(limma_MM, function(x) {
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("limma",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

#logN
table(logN_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(logN_MM, function(x) {
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("Log Normalize",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

#M3
table(M3_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(M3_MM, function(x) {
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("M3Drop",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

#scvi
table(scvi_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(scvi_MM, function(x) {
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("SCVI",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

#scanorama
table(Scanorama_MM[[1]]$dir_sig)[c(1,3)]
# -> Prop of up/down between the truth and individual in each method
sum <- lapply(Scanorama_MM, function(x) {
  print(table(x$dir_sig))
  res <- table(x$dir_sig)[c(1,3)]
  return(res)
}
)
sum <- unlist(sum)
new_rows <- data.frame(method = rep("Scanorama",length(sum)), direction = c(rep(c("down","up"),length(sum)/2)), value = sum, indiv = names(sum))
df <- rbind(df, new_rows)

df$individual <- gsub("\\..*", "", df$indiv) # optional 
df$method <- factor(df$method, levels = c("Truth", "Log Normalize", "limma", "fastMNN", "Scanorama", "SCVI", "M3Drop"))
df$direction <- factor(df$direction, levels = c("up", "down"))

saveRDS(df, paste0("./",DS,"_PLOT_UP_DOWN_noFC.rds"))

df <- readRDS(paste0("./",DS,"_PLOT_UP_DOWN_noFC.rds"))

ggplot(df, aes(x = direction, y = value, fill = method)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_point(aes(shape = individual, group = method),
             position = position_dodge(width = 0.8), 
             size = 2.2, color = "black", alpha = 1.2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4)) +
  scale_fill_manual(
    values = c(
      "Truth" = "#6a3d9a",
      "fastMNN" = "#33a02c",
      "limma" = "#1f78b4",
      "Log Normalize" = "lightgrey",
      "M3Drop" = "#e31a1c",
      "SCVI" = "#fdbf6f"
    )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", size = 16),
    legend.title = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  labs(
    x = "Direction",
    y = "Gene Count",
    fill = "Method",
    shape = "Individual"
  ) +
  scale_y_continuous(
    limits = c(0, 26000),
    breaks = seq(0, 25000, 5000)
  )


### ----- Figure 5: ----- 

#### Venn Diagram #####
truth <- readRDS("./DS1_intersect_truth_labels_intersect_noFC.rds") 
markers <- truth$gene[truth$dir_sig == "up" | truth$dir_sig == "down"]

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

length(intersect(markers, common_ind))

con <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
ens_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = intersect(markers, common_ind),
                   mart = con)

truth <- readRDS("./DS2lv_intersect_truth_labels_intersect_noFC.rds")
markers2 <- truth$gene[truth$dir_sig == "up" | truth$dir_sig == "down"]

length(intersect(markers2, common_ind)) 

con <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
ens_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = intersect(markers2, common_ind),
                   mart = con)

intersect(markers, markers2) 

x <- intersect(intersect(markers, common_ind), intersect(markers6, common_ind)) 

truth <- readRDS("./DS6_intersect_truth_labels_intersect_noFC.rds")
markers6 <- truth$gene[truth$dir_sig == "up" | truth$dir_sig == "down"]

length(intersect(markers6, common_ind))

#### Ridge Plots ####

##### ------ DS1 - Treg vd CD8 -----

DS <- "DS1"

fast <- readRDS(paste0("./",DS,"_merged_fast_25.rds"))
limma <- readRDS(paste0("./",DS,"_merged_limma_25.rds"))
Scanorama <- readRDS(paste0("./",DS,"_merged_scanorama_25.rds"))
logN <- readRDS(paste0("./",DS,"_merged_logN_25.rds"))
M3 <- readRDS(paste0("./",DS,"_merged_M3Drop_25.rds"))
scvi <- readRDS(paste0("./",DS,"_merged_scvi_25.rds"))

# Original:
original <- readRDS("./Scaled_LogNorm.rds")
original <- merge(original[[1]], original[-1])
original <- JoinLayers(original)

###### FOXP3 UP ####
list <- c("logN_MM", "fast_MM", "limma_MM", "M3_MM", "scvi_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000049768", ])
  }
}


# 3'
levels(as.factor(original@meta.data$assay))
Idents(original) <- "cell_type"

original <- RenameIdents(
  original,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#A0B3E3",
    "Regulatory T cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
# check if this gene is in the list:
for(i in 1:length(fast_MM)){
  print(fast_MM[[i]][fast_MM[[i]]$gene == "ENSG00000049768", ])
}
#
Idents(fast) <- "cell_type"
fast <- RenameIdents(
  fast,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(fast) <- factor(Idents(fast),
                           levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(fast, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
# check if this gene is in the list:

Idents(limma) <- "cell_type"
limma <- RenameIdents(
  limma,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(limma) <- factor(Idents(limma),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(limma, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type"
logN <- RenameIdents(
  logN,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(logN) <- factor(Idents(logN),
                        levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(logN, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type"
M3 <- RenameIdents(
  M3,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(M3) <- factor(Idents(M3),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(M3, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type"
scvi <- RenameIdents(
  scvi,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(scvi) <- factor(Idents(scvi),
                     levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(scvi, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type"
Scanorama <- RenameIdents(
  Scanorama,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(Scanorama, features = "ENSG00000049768", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))


###### CD8A DOWN; GZMB DOWN #####

list <- c("logN_MM", "fast_MM", "limma_MM", "M3_MM", "scvi_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000153563", ])
  }
}

# 3'
levels(as.factor(original@meta.data$assay))
Idents(original) <- "cell_type"

original <- RenameIdents(
  original,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#A0B3E3",
    "Regulatory T cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN

Idents(fast) <- "cell_type"
fast <- RenameIdents(
  fast,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(fast, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
# check if this gene is in the list:

Idents(limma) <- "cell_type"
limma <- RenameIdents(
  limma,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(limma, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type"
logN <- RenameIdents(
  logN,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(logN, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type"
M3 <- RenameIdents(
  M3,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(M3, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type"
scvi <- RenameIdents(
  scvi,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(scvi, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type"
Scanorama <- RenameIdents(
  Scanorama,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("CD8, αβ T cell", "Regulatory T cell"))

RidgePlot(Scanorama, features = "ENSG00000153563", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

###### CD74 OVERLAPPING ####
list <- c("logN_MM", "limma_MM", "fast_MM", "scvi_MM", "M3_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000019582", ])
  }
}

# 3'

levels(as.factor(original@meta.data$assay))
Idents(original) <- "cell_type"

original <- RenameIdents(
  original,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#A0B3E3",
    "Regulatory T cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
Idents(fast) <- "cell_type"
fast <- RenameIdents(
  fast,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(fast, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
# check if this gene is in the list:

Idents(limma) <- "cell_type"
limma <- RenameIdents(
  limma,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(limma, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type"
logN <- RenameIdents(
  logN,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(logN, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell")) + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type"
M3 <- RenameIdents(
  M3,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(M3, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type"
scvi <- RenameIdents(
  scvi,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(scvi, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type"
Scanorama <- RenameIdents(
  Scanorama,
  "CD8-positive, alpha-beta T cell" = "CD8, αβ T cell",
  "regulatory T cell" = "Regulatory T cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Regulatory T cell", "CD8, αβ T cell"))

RidgePlot(Scanorama, features = "ENSG00000019582", idents = c("CD8, αβ T cell", "Regulatory T cell"), layer = "counts") + #FOXP3
  scale_fill_manual(values = c(
    "CD8, αβ T cell" = "#FF9B99",
    "Regulatory T cell" = "#375DBE"
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))


##### ------ DS6 - B cell vs dendr -----
DS <- "DS6"

fast <- readRDS(paste0("./",DS,"_merged_fast_25.rds"))
limma <- readRDS(paste0("./",DS,"_merged_limma_25.rds"))
Scanorama <- readRDS(paste0("./",DS,"_merged_scanorama_25.rds"))
logN <- readRDS(paste0("./",DS,"_merged_logN_25.rds"))
M3 <- readRDS(paste0("./",DS,"_merged_M3Drop_25.rds"))
scvi <- readRDS(paste0("./",DS,"_merged_scvi_25.rds"))

#original:
original <- readRDS("./Scaled_LogNorm.rds")
original <- merge(original[[1]], original[-1])
original <- JoinLayers(original)

###### CD79A DOWN #####
# check the significance
list <- c("logN_MM", "fast_MM", "limma_MM", "M3_MM", "scvi_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000105369", ])
  }
}
###

Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("B cell", "Dendritic cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000105369", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000105369", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN

Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell"
)

Idents(fast) <- factor(Idents(fast),
                           levels = c("B cell", "Dendritic cell"))

RidgePlot(fast, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell"
)

Idents(limma) <- factor(Idents(limma),
                       levels = c("B cell", "Dendritic cell"))

RidgePlot(limma, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell"
)

Idents(logN) <- factor(Idents(logN),
                        levels = c("B cell", "Dendritic cell"))

RidgePlot(logN, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell"
)

Idents(M3) <- factor(Idents(M3),
                       levels = c("B cell", "Dendritic cell"))

RidgePlot(M3, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell"
)

Idents(scvi) <- factor(Idents(scvi),
                     levels = c("B cell", "Dendritic cell"))

RidgePlot(scvi, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("B cell", "Dendritic cell"))

RidgePlot(Scanorama, features = "ENSG00000105369", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

###### CST3  UP #####
# check the significance
list <- c("logN_MM", "fast_MM", "limma_MM", "M3_MM", "scvi_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000101439", ])
  }
}
###

Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Dendritic cell", "B cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000101439", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000101439", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("Dendritic cell", "B cell"))

RidgePlot(fast, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("Dendritic cell" ,"B cell"))

RidgePlot(limma, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("Dendritic cell", "B cell"))

RidgePlot(logN, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("Dendritic cell", "B cell"))

RidgePlot(M3, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("Dendritic cell","B cell"))

RidgePlot(scvi, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Dendritic cell","B cell"))

RidgePlot(Scanorama, features = "ENSG00000101439", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

###### CD99 OVERLAPPING #####
# check the significance
list <- c("logN_MM", "fast_MM", "limma_MM", "M3_MM", "scvi_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000002586", ])
  }
}

###

Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Dendritic cell", "B cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000002586", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000002586", idents = c("B cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN

Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("Dendritic cell", "B cell"))

RidgePlot(fast, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("Dendritic cell" ,"B cell"))

RidgePlot(limma, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("Dendritic cell", "B cell"))

RidgePlot(logN, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("Dendritic cell", "B cell"))

RidgePlot(M3, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("Dendritic cell","B cell"))

RidgePlot(scvi, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Dendritic cell","B cell"))

RidgePlot(Scanorama, features = "ENSG00000002586", idents = c("B cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "B cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

##### ------ DS2lv - dendr vs NK -----
DS <- "DS2lv"

# check the significance
list <- c("logN_MM", "limma_MM", "fast_MM", "scvi_MM", "M3_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000134539", ])
  }
}
###

fast <- readRDS(paste0("./",DS,"_merged_fast_25.rds"))
limma <- readRDS(paste0("./",DS,"_merged_limma_25.rds"))
Scanorama <- readRDS(paste0("./",DS,"_merged_scanorama_25.rds"))
logN <- readRDS(paste0("./",DS,"_merged_logN_25.rds"))
M3 <- readRDS(paste0("./",DS,"_merged_M3Drop_25.rds"))
scvi <- readRDS(paste0("./",DS,"_merged_scvi_25.rds"))

#original:
original <- readRDS("./Scaled_LogNorm_liver.rds")
original <- merge(original[[1]], original[-1])
original <- JoinLayers(original)

###### KLRD1 DOWN - ENSG00000134539 #####
Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(fast) <- factor(Idents(fast),
                           levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(fast, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(limma) <- factor(Idents(limma),
                       levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(limma, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(logN) <- factor(Idents(logN),
                        levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(logN, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(M3) <- factor(Idents(M3),
                       levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(M3, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(scvi) <- factor(Idents(scvi),
                     levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(scvi, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Natural Killer cell", "Dendritic cell"))

RidgePlot(Scanorama, features = "ENSG00000134539", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

###### HLA-DRA UP - ENSG00000204287 #####
# check the significance
list <- c("logN_MM", "limma_MM", "fast_MM", "scvi_MM", "M3_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000204287", ])
  }
}
###

Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(fast, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(limma, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(logN, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(M3, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(scvi, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(Scanorama, features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

###### HLA-DRB1 OVERLAPPING #####
# check the significance
list <- c("logN_MM", "limma_MM", "fast_MM", "scvi_MM", "M3_MM")
# check if this gene is in the list:
for(j in list) {
  obj <- get(j)
  print(j)
  for(i in 1:length(obj)){
    print(obj[[i]][obj[[i]]$gene == "ENSG00000196126", ])
  }
}
###

Idents(original) <- "cell_type_ontology_term_id"
original <- RenameIdents(
  original,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(original) <- factor(Idents(original),
                           levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(subset(original, assay == "10x 3' v2"), features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#FF221F" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

RidgePlot(subset(original, assay == "10x 5' v1"), features = "ENSG00000204287", idents = c("Natural Killer cell", "Dendritic cell"), layer = "data") + #FOXP3
  scale_fill_manual(values = c(
    "Dendritic cell" = "#A0B3E3",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# FastMNN
Idents(fast) <- "cell_type_ontology_term_id"
fast <- RenameIdents(
  fast,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(fast) <- factor(Idents(fast),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(fast, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Limma
Idents(limma) <- "cell_type_ontology_term_id"
limma <- RenameIdents(
  limma,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(limma) <- factor(Idents(limma),
                        levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(limma, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Log Normalize
Idents(logN) <- "cell_type_ontology_term_id"
logN <- RenameIdents(
  logN,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(logN) <- factor(Idents(logN),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(logN, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell")) + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# M3Drop
Idents(M3) <- "cell_type_ontology_term_id"
M3 <- RenameIdents(
  M3,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(M3) <- factor(Idents(M3),
                     levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(M3, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# SCVI
Idents(scvi) <- "cell_type_ontology_term_id"
scvi <- RenameIdents(
  scvi,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(scvi) <- factor(Idents(scvi),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(scvi, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

# Scanorama
Idents(Scanorama) <- "cell_type_ontology_term_id"
Scanorama <- RenameIdents(
  Scanorama,
  "dendritic cell" = "Dendritic cell",
  "natural killer cell" = "Natural Killer cell"
)

Idents(Scanorama) <- factor(Idents(Scanorama),
                       levels = c("Dendritic cell", "Natural Killer cell"))

RidgePlot(Scanorama, features = "ENSG00000196126", idents = c("Natural Killer cell", "Dendritic cell"), layer = "counts") + 
  scale_fill_manual(values = c(
    "Dendritic cell" = "#FF9B99",
    "Natural Killer cell" = "#375DBE" 
  )) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black", size = 12.5),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_line(color = "black"))

### ----- Sup Fig 1: ----

# This is done for all of the datasets
ds <- readRDS("./prep_FINAL.rds")

cells_by_cell_type_and_technology <- ds@meta.data %>%
  group_by(cell_type, assay, new_id) %>%
  summarize(nCells = n()) %>%
  ungroup()

ggplot(cells_by_cell_type_and_technology, aes(x = cell_type, y = nCells, fill = assay)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cell Type", y = "Number of Cells") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1, size = 12),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  ) +
  facet_wrap( ~ new_id, scales = "free_y",  nrow = 3)+ 
  scale_fill_manual(
    values = c("#e31a1c", "#1f78b4")) +
  scale_y_continuous(
    limits = c(0, 1500),
    breaks = seq(0, 1500, 250))


### ----- Sup Fig 2: ----

# Prepare the data
data <- readRDS("./Scaled_LogNorm.rds")
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

data_merged <- merge(data[[1]], y = data[-1])
data_merged <- JoinLayers(data_merged)

data_merged <- RunPCA(data_merged)
ElbowPlot(data_merged, ndims = 50)
data_merged <- RunUMAP(data_merged, dims = 1:30)

data_merged1 <- SetAssayData(data_merged, layer = "scale.data", new.data = NULL)
data_merged1 <- SetAssayData(data_merged1, layer = "data", new.data = NULL)
data_merged1 <- data_merged1[(rownames(data_merged1) %in% common_ind[50:100]), ]
saveRDS(data_merged1, "./Scaled_LogNorm_PLOT.rds")

data_merged <- data_merged[!(rownames(data_merged) %in% common_ind), ]

data_merged <- split(data_merged, f = data_merged@meta.data$new_id)
data_merged <- NormalizeData(data_merged)
data_merged <- FindVariableFeatures(data_merged, selection.method = "vst", nfeatures = 3000)
data_merged <- JoinLayers(data_merged)

data_merged_list <- SplitObject(data_merged, split.by = "new_id")
data_merged_list <- lapply(data_merged_list, function(x){
  x <- ScaleData(x)
})

data_merged <- merge(data_merged_list[[1]], y = data_merged_list[-1])
data_merged <- JoinLayers(data_merged)

data_merged <- RunPCA(data_merged)
ElbowPlot(data_merged, ndims = 50)
data_merged <- RunUMAP(data_merged, dims = 1:30)

data_merged1 <- SetAssayData(data_merged, layer = "scale.data", new.data = NULL)
data_merged1 <- SetAssayData(data_merged1, layer = "data", new.data = NULL)
data_merged1 <- data_merged1[(rownames(data_merged1) %in% c("ENSG00000144290")), ]
saveRDS(data_merged1, "./Scaled_LogNorm_PLOT_no800.rds")

#### Plot 
original <- readRDS("./Scaled_LogNorm_PLOT.rds")

DimPlot(original, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "right"                # remove legend
  )

DimPlot(original, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4"))  +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

reduced <- readRDS("./Scaled_LogNorm_PLOT_no800.rds")

DimPlot(reduced, reduction = "umap", group.by = "assay", cols = c("#e31a1c", "#1f78b4")) +
  ggtitle(NULL) +                    # remove the title
  theme(
    axis.title = element_blank(),    # remove axis labels
    axis.text = element_blank(),     # remove axis tick labels
    axis.ticks = element_blank(),    # remove axis ticks
    plot.title = element_blank(),     # remove plot title (extra safe)
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "none"                # remove legend
  )

### ----- Sup Fig 5: ----

# DS1
data <- readRDS("./Scaled_LogNorm10000.rds") # Do the same for all datasets

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

ds1 <- data.frame(
  individual = paste0(names(dge_list),"_", "Park"),
  num_significant_genes = sapply(dge_list, nrow)
)

# DS2
# - bone marrow
data <- readRDS("./Scaled_LogNorm10000_bone marrow.rds") # Do it for all datasets
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Suo"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- liver
data <- readRDS("./Scaled_LogNorm10000_liver.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Suo"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- spleen
data <- readRDS("./Scaled_LogNorm10000_spleen.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Suo"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- thymus
data <- readRDS("./Scaled_LogNorm10000_thymus.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Suo"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- kidney
data <- readRDS("./Scaled_LogNorm10000_kidney.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Suo"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

# DS3
data <- readRDS("./Scaled_LogNorm10000.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Andrews"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

# DS4
data <- readRDS("./Scaled_LogNorm10000.rds")
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v3"
assay_5 <- "10x 5' v2"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Simone"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

# DS5
data <- readRDS("./Scaled_LogNorm10000.rds") 
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list)),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

# DS6
data <- readRDS("./Scaled_LogNorm10000.rds") # Do it for all datasets
levels(as.factor(data[[1]]@meta.data$assay))

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

dge_list <- list()
genes <- list()

for (i in names(data)) {
  obj <- data[[i]]
  obj <- subset(obj, features = common_ind)
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0, return.thresh = 1.5)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Jardine"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

## Validation 1 - kidney
data <- readRDS("./Validation_prep_FINAL.rds")
levels(as.factor(data@meta.data$assay))

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

data <- subset(data, features = common_ind)

assay_3 <- "10x 3' v3"
assay_5 <- "10x 5' v1"

dge_list <- list()
genes <- list()

data <- split(data, f = data@meta.data$new_id)
data <- NormalizeData(data)
data <- JoinLayers(data)

data <- SplitObject(data, split.by = "new_id")

for (i in names(data)) {
  #i <- 1
  obj <- data[[i]]
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0.00, min.cells.feature = 0, min.cells.group = 0)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Validation_1"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

## Validation 2 - kidney, liver, skin
data <- readRDS("./Validation_prep_FINAL_kidney.rds")
levels(as.factor(data@meta.data$assay))

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

data <- subset(data, features = common_ind)

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

dge_list <- list()
genes <- list()

data <- split(data, f = data@meta.data$new_id)
data <- NormalizeData(data)
data <- JoinLayers(data)

data <- SplitObject(data, split.by = "new_id")

for (i in names(data)) {
  #i <- 1
  obj <- data[[i]]
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0.00, min.cells.feature = 0, min.cells.group = 0)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Validation_2 (kidney)"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- liver
data <- readRDS("./Validation_prep_FINAL_liver.rds")
levels(as.factor(data@meta.data$assay))

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

data <- subset(data, features = common_ind)

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

dge_list <- list()
genes <- list()

data <- split(data, f = data@meta.data$new_id)
data <- NormalizeData(data)
data <- JoinLayers(data)

data <- SplitObject(data, split.by = "new_id")

for (i in names(data)) {
  #i <- 1
  obj <- data[[i]]
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0.00, min.cells.feature = 0, min.cells.group = 0)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Validation_2 (liver)"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

#-- skin
data <- readRDS("./Validation_prep_FINAL_skin.rds")
levels(as.factor(data@meta.data$assay))

gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- c(names(gene_table[gene_table %in% c(5:35)])) # 800: 50

data <- subset(data, features = common_ind)

assay_3 <- "10x 3' v2"
assay_5 <- "10x 5' v1"

dge_list <- list()
genes <- list()

data <- split(data, f = data@meta.data$new_id)
data <- NormalizeData(data)
data <- JoinLayers(data)

data <- SplitObject(data, split.by = "new_id")

for (i in names(data)) {
  #i <- 1
  obj <- data[[i]]
  Idents(obj) <- "assay"
  DGE <- FindMarkers(obj, ident.1 = assay_3, ident.2 = assay_5, logfc.threshold = 0, min.pct = 0.00, min.cells.feature = 0, min.cells.group = 0)
  sign <- DGE[DGE$p_val_adj < 0.05,] # Change the threshold to get a reasonable number of genes for the table
  
  dge_list[[i]] <- sign
  genes[[i]] <- rownames(sign)
}

temp <- data.frame(
  individual = paste0(names(dge_list),"_", "Validation_2 (skin)"),
  num_significant_genes = sapply(dge_list, nrow)
)

ds1 <- rbind(ds1, temp)

saveRDS(ds1, "./All_DS_and_Valid_num_of_biased_sign_an_05.rds")

## Plot:
cell_nums <- c(5832,
               5075,
               4381,
               2864,
               2603,
               2653,
               3147,
               1684,
               3203,
               1386,
               3640,
               6594,
               7173,
               2457,
               2287,
               1656,
               3757,
               7857,
               7969,
               7231,
               9956,
               4668,
               15232,
               4906,
               3973,
               3894,
               6006,
               2723,
               14887,
               4744,
               29944,
               5441,
               6314,
               4978,
               3027,
               3564,
               2007,
               3377,
               4280,
               3390,
               11280,
               3574,
               6099)

ds1$cell_nums <- cell_nums
saveRDS(ds1, "./All_DS_and_Valid_num_of_biased_sign_an_05.rds")


## Plot
ggplot(
  ds1,
  aes(
    x = num_significant_genes,
    y = cell_nums,
    color = color,
    shape = color
  )
) +
  geom_point(size = 2.5, alpha = 0.85) +
  
  scale_color_manual(
    values = c(
      "Original" = "darkgrey",
      "Validation" = "#E62828"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Original" = 16,
      "Validation" = 17
    )
  ) +
  
  labs(
    x = "Number of significant genes",
    y = "Number of cells",
    color = NULL,
    shape = NULL
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    
    legend.position = "top",
    legend.text = element_text(size = 11),
    
    plot.margin = margin(8, 12, 8, 8)
  )


### -- Sup Fig 8: ----
split_df <- readRDS("./HVF_FINAL_all_DSandTechniques_NOT_to_baseline.rds")
split_df <- readRDS("./800_FINAL_all_DSandTechniques_NOT_to_baseline.rds")

## Plotting
split_df$corr$norm <- factor(split_df$corr$norm.1, levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                            "Z-transform", "SCVI", "scArches", "Log_Normalize_Scale", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$cos$norm <- factor(split_df$cos$norm.1, levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                          "Z-transform", "SCVI", "scArches", "Log_Normalize_Scale", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$Euc$norm <- factor(split_df$corr$norm.1, levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "Log_Normalize_Scale", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$MSE$norm <- factor(split_df$corr$norm.1, levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "Log_Normalize_Scale", "M3Drop", "scTransform_v5", "scTransform_v5_split"))
split_df$JSD$norm <- factor(split_df$corr$norm.1, levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                                           "Z-transform", "SCVI", "scArches", "Log_Normalize_Scale", "M3Drop", "scTransform_v5", "scTransform_v5_split"))

pref <- "HVF_" # or "800_"
col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "Log_Normalize_Scale" = "lightgrey",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))
# a) Corr
ggplot(split_df$corr, aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = "Correlation Coefficient", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none") +
  scale_y_continuous(
    #limits = c(0, 1),
    breaks = seq(0, 1, 0.1))

# b) Cos
ggplot(split_df$cos, aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Cosine Similarity", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  guides(fill = "none")+
  scale_y_continuous(
    #limits = c(0, 1),
    breaks = seq(0, 1, 0.1))

# c) MSE
ggplot(split_df$MSE, aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Mean Squared Error", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  ) +
  guides(fill = "none") #+
  #scale_y_break(c(15000,53000))

ggplot(subset(split_df$MSE, subset = norm !="Z-transform" & (norm !="SCVI") & norm !="scArches"), aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Mean Squared Error", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")


# d) JSD
ggplot(split_df$JSD, aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Jensen-Shannon Divergence", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")


# e) Euc
ggplot(split_df$Euc, aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Euclidean Distance", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(subset(split_df$Euc, subset = norm !="Z-transform" ), aes(x = norm, y = average, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = DS), width = 0.2, size = 1.2, alpha = 1) + # Add jittered dots, colored by name
  labs(title = "Euclidean Distance", x = "Normalization", y = "Value", shape = "Dataset") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

### ----- Sup Fig 10: ----

### MM for cell type

tis <- "" # Specify the tissue if dataset 2
ds <- "" # Specify the dataset
DS <- "" # Specify the dataset

melted_data <- readRDS(paste0("./FINAL_MM_cell",tis,".rds"))

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = NULL, height = 0, size = 1, alpha = 1) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(melted_data, aes(x = variable, y = value, height = 0, fill = variable)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = Individual), width = 0.3, height = 0, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Mixing Value", shape = "Individual") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

### ----- Sup Fig 11: ----
### Overlap Scores:

tis <- "" # Specify the tissue if dataset 2
ds <- "" # Specify
DS <- "" # Specify

LN <- readRDS(paste0("./LogNorm_JS_sum",tis,".rds"))
JS_df <- LN
CB <- readRDS(paste0("./ComBat_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, CB)
mnn <- readRDS(paste0("./mnn_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, mnn)
Fast <- readRDS(paste0("./fast_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, Fast)
Scanorama <- readRDS(paste0("./scanorama_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, Scanorama)
limma <- readRDS(paste0("./limma_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, limma)
Z_t <- readRDS(paste0("./Z_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, Z_t)
scvi <- readRDS(paste0("./",ds,"scvi_JS_sum.rds"))
JS_df <- rbind(JS_df, scvi)
scA <- readRDS(paste0("./",ds,"scArches_JS_sum.rds"))
JS_df <- rbind(JS_df, scA)
M3 <- readRDS(paste0("./m3drop_JS_sum",tis,".rds"))
JS_df <- rbind(JS_df, M3)
SC <- readRDS(paste0("./SC_JS_sum",tis,"_2.rds"))
JS_df <- rbind(JS_df, SC)
SC_s <- readRDS(paste0("./SC_split_JS_sum",tis,"_2.rds"))
JS_df <- rbind(JS_df, SC_s)

JS_df$norm  <- factor(JS_df$norm , levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                              "Z-transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z-transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

ggplot(JS_df, aes(x = norm, y = Average_Jaccard_Score, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12, 13)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

ggplot(JS_df, aes(x = norm, y = Average_Jaccard_Score, fill = norm)) +
  geom_boxplot(outliers = FALSE) + # Boxplot without outliers
  geom_jitter(aes(shape = name), width = 0.3, size = 1, alpha = 1.2) + # Add jittered dots, colored by name
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12,13)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none") +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0, 1, 0.1))

## Number of significant features per cell type:

ds <- ""
tis <- ""
DS <- ""

Log_Normalize <- readRDS(paste0("./LogNorm_raw_markers",tis,".rds"))
ComBat <- readRDS(paste0("./ComBat_raw_markers",tis,".rds"))
limma <- readRDS(paste0("./limma_raw_markers",tis,".rds"))
mnnCorrect <- readRDS(paste0("./mnn_raw_markers",tis,".rds"))
fastMNN <- readRDS(paste0("./fast_raw_markers",tis,".rds"))
Scanorama <- readRDS(paste0("./scanorama_raw_markers",tis,".rds"))
Z_transform <- readRDS(paste0("./Z_raw_markers",tis,".rds"))
SCVI <- readRDS(paste0("./",DS,"_scvi_raw_markers.rds"))
scArches <- readRDS(paste0("./",DS,"_scArches_raw_markers.rds"))
M3Drop <- readRDS(paste0("./m3drop_raw_markers",tis,".rds"))
scTransform_v5 <- readRDS(paste0("./SC_raw_markers",tis,"_2.rds"))
scTransform_v5_split <- readRDS(paste0("./SC_split_raw_markers",tis,"_2.rds"))

list <- unique(sub("_[^_]+$", "", names(ComBat)))
tech <- c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama", "Z_transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split")

jaccard_scores <- list()

for (t in tech) {
  raw_markers <- get(t)

  for (i in list) {
  
    sig_genes_3 <- subset(raw_markers[[paste0(i, "_3")]], p_val_adj < 0.05 )
    sig_genes_5 <- subset(raw_markers[[paste0(i, "_5")]], p_val_adj < 0.05 )

    score <- numeric(length(intersect(sig_genes_3$cluster, sig_genes_5$cluster)))
    names(score) <- intersect(sig_genes_3$cluster, sig_genes_5$cluster)
    for (cell in intersect(sig_genes_3$cluster, sig_genes_5$cluster)) {
      score[[cell]] <- length(intersect(sig_genes_3$gene[sig_genes_3$cluster == cell], 
                                        sig_genes_5$gene[sig_genes_5$cluster == cell]))
    }
    jaccard_scores[[i]] <- score
  }
  df <- map_dfr(jaccard_scores, enframe, .id = "Sample")
  df <- df %>%
    group_by(name) %>%
    summarise(Average_Score = mean(value, na.rm = TRUE))
  df$norm <- t
  
  if(exists("final")){
    final <- rbind(final, df)
  } else {
    final <- df
  }
}

final$norm  <- factor(final$norm , levels = c("Log_Normalize", "ComBat", "limma", "mnnCorrect", "fastMNN", "Scanorama",
                                              "Z_transform", "SCVI", "scArches", "M3Drop", "scTransform_v5", "scTransform_v5_split"))

col <- scale_fill_manual(values = c(
  "Log_Normalize" = "lightgrey",
  "ComBat" = "#a6cee3",
  "limma" = "#1f78b4",
  "mnnCorrect" = "#b2df8a",
  "fastMNN" = "#33a02c",
  "Scanorama" = "#5e8011",
  "Z_transform" = "#fb9a99",
  "SCVI" = "#fdbf6f",
  "scArches" = "#ff7f00",
  "M3Drop" = "#e31a1c",
  "scTransform_v5" = "#cab2d6",
  "scTransform_v5_split" = "#6a3d9a"
))

ggplot(final, aes(x = norm, y = Average_Score, fill = norm)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_point(aes(shape = name, group = norm),
             position = position_dodge(width = 0.8), 
             size = 2.2, color = "black", alpha = 1.2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12, 13)) +
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  theme_minimal() +
    theme(
      axis.text.y = element_text(color = "black", size = 16),
      legend.title = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line.y = element_line(color = "black", linewidth = 1),
      axis.line.x = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black", linewidth = 1)
    )+
  guides(fill = "none")

ggplot(final, aes(x = norm, y = Average_Score, fill = norm)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_point(aes(shape = name, group = norm),
             position = position_dodge(width = 0.8), 
             size = 2.2, color = "black", alpha = 1.2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11,12, 13)) +
  labs(title = DS, x = "Normalization", y = "Overlap Score", shape = "Cell Type") +
  col +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black", size = 16),
    legend.title = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line.y = element_line(color = "black", linewidth = 1),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )+
  guides(fill = "none")

### ----- Sup Fig 12: ----

# DS1
obj <- readRDS("./prep_FINAL.rds")
obj <- subset(obj, new_id %in% c("F29_45P", "F30_45P", "F38_45P",  "P1_CD3P")) # "F41_45P", "F45_45P" - no Tregs
obj@meta.data <- droplevels(obj@meta.data)
cell_type1 <- "regulatory T cell"
cell_type2 <- "CD8-positive, alpha-beta T cell"
truth <- readRDS("./DS1_intersect_truth_labels_intersect_noFC.rds")

obj <- split(obj, f = obj@meta.data$new_id)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- JoinLayers(obj)
obj <- subset(obj, cell_type %in% c(cell_type1, cell_type2)) 
obj@meta.data <- droplevels(obj@meta.data)

con <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = truth$gene,
                   mart = con)

truth <- merge(truth, ens_genes, by.x = "gene", by.y = "ensembl_gene_id")
gold_heat_up <- truth[truth$dir_sig == "up",]

gold_heat_down <- truth[truth$dir_sig == "down",]
gold_heat_down <- gold_heat_down[order(gold_heat_down$external_gene_name), ]

set.seed(2)
DoHeatmap(obj, features = c(gold_heat_up$gene[46:27], gold_heat_down$gene[20:1]), label = FALSE,  group.by = "cell_type", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) +  
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[1:20], gold_heat_up$external_gene_name[27:46])) +
  theme(legend.position = "none", axis.text.y = element_text(color = "black", size = 12, face = "italic"))
DoHeatmap(obj, features = c(gold_heat_up$gene[46:27], gold_heat_down$gene[20:1]), label = FALSE,  group.by = "cell_type", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) +   
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[1:20], gold_heat_up$external_gene_name[27:46])) +
  theme(legend.position = "right", axis.text.y = element_text(color = "black", size = 12, face = "italic"))


# DS2lv
obj <- readRDS("./prep_FINAL.rds")
obj <- subset(obj, new_id %in% c("F32_CD45P_liver", "F34_CD45P_liver", "F38_CD45P_liver",  "F41_CD45P_liver", "F45_CD45P_liver")) # CD45N not included
obj@meta.data <- droplevels(obj@meta.data)
cell_type1 <- "dendritic cell"
cell_type2 <- "natural killer cell"
truth <- readRDS("./DS2lv_intersect_truth_labels_intersect_noFC.rds")

obj <- split(obj, f = obj@meta.data$new_id)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- JoinLayers(obj)
obj <- subset(obj, cell_type_ontology_term_id %in% c(cell_type1, cell_type2)) 
obj@meta.data <- droplevels(obj@meta.data)

con <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = truth$gene,
                   mart = con)

truth <- merge(truth, ens_genes, by.x = "gene", by.y = "ensembl_gene_id")
gold_heat_up <- truth[truth$dir_sig == "up",]

gold_heat_down <- truth[truth$dir_sig == "down",]

DoHeatmap(obj, features = c(gold_heat_up$gene[1:20], gold_heat_down$gene[1:20]), label = FALSE,  group.by = "cell_type_ontology_term_id", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) +
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[20:1], gold_heat_up$external_gene_name[20:1])) +
  theme(legend.position = "none", axis.text.y = element_text(color = "black", size = 12, face = "italic"))
DoHeatmap(obj, features = c(gold_heat_up$gene[1:20], gold_heat_down$gene[1:20]), label = FALSE,  group.by = "cell_type_ontology_term_id", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) +
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[20:1], gold_heat_up$external_gene_name[20:1])) +
  theme(legend.position = "right", axis.text.y = element_text(color = "black", size = 12, face = "italic"))

# DS6
obj <- readRDS("./prep_FINAL.rds")
cell_type1 <- "dendritic cell"
cell_type2 <- "B cell"
truth <- readRDS("./DS6_intersect_truth_labels_intersect_noFC.rds")

obj <- split(obj, f = obj@meta.data$new_id)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- JoinLayers(obj)
obj <- subset(obj, cell_type_ontology_term_id %in% c(cell_type1, cell_type2)) # cell_type_ontology_term_id
obj@meta.data <- droplevels(obj@meta.data)

con <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = truth$gene,
                   mart = con)

truth <- merge(truth, ens_genes, by.x = "gene", by.y = "ensembl_gene_id")
gold_heat_up <- truth[truth$dir_sig == "up",]

gold_heat_down <- truth[truth$dir_sig == "down",]

DoHeatmap(obj, features = c(gold_heat_up$gene[1:20], gold_heat_down$gene[1:20]), label = FALSE,  group.by = "cell_type_ontology_term_id", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) +  
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[20:1], gold_heat_up$external_gene_name[20:1])) +
  theme(legend.position = "none", axis.text.y = element_text(color = "black", size = 12, face = "italic"))
DoHeatmap(obj, features = c(gold_heat_up$gene[1:20], gold_heat_down$gene[1:20]), label = FALSE,  group.by = "cell_type_ontology_term_id", slot = "data", group.colors = c("#a6cee3", "#cab2d6"), draw.lines = TRUE) + 
  scale_fill_gradientn(colors = c("#F5F5F5", "#E57710", "#AB590D")) + 
  scale_y_discrete(labels = c(gold_heat_down$external_gene_name[20:1], gold_heat_up$external_gene_name[20:1])) +
  theme(legend.position = "right", axis.text.y = element_text(color = "black", size = 12, face = "italic"))




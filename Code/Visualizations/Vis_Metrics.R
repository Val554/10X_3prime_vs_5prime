# Purpose: Summarizing the metrics for the comparison of the 3' and 5' assays for each batch correction method. 
# This is done for further visualizations and comparisons across the different techniques. 

# Loading required libraries
library(reshape2) 
library(ggplot2)
library(dplyr)

# Load one object at a time
### for HVFs
result_list <- readRDS("./result_mnnCorrect_hvf.rds")
result_list <- readRDS("./result_limma_hvf.rds")
result_list <- readRDS("./result_ComBat_hvf.rds")
result_list <- readRDS("./result_M3Drop_hvf.rds")
result_list <- readRDS("./result_z_trans_hvf.rds")
result_list <- readRDS("./result_fastMNN_hvf.rds")
result_list <- readRDS("./result_scanorama_hvf.rds")
result_list <- readRDS("./result_logNorm_hvf.rds")
result_list <- readRDS("./result_logNorm_scaled_hvf.rds")
result_list <- readRDS("./result_scTransf_hvf.rds")
result_list <- readRDS("./result_scTransf_split_hvf.rds")
result_list <- readRDS("./result_scvi_hvf.rds")
result_list <- readRDS("./result_scArch_hvf.rds")

### for 800
ds <- ""
result_list <- readRDS(paste0("./DS",ds,"_result_scvi_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_scArch_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_mnnCorrect_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_limma_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_ComBat_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_M3Drop_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_z_transf_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_fastMNN_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_scanorama_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_scTransf_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_scTransf_split_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_logNorm_800.rds"))
result_list <- readRDS(paste0("./DS",ds,"_result_logNorm_scaled_800.rds"))

# --------- Average each metrics over the subjects for each cell type ----

len <- length(result_list[["AllCellTypes"]])
cell_avg <- list(
  corr = setNames(rep(0, len), result_list[["AllCellTypes"]]),
  cos = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  MSE = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  JSD = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  Euc = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  Manh = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  mean_3 = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  mean_5 = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  var_3 = setNames(rep(0, len),result_list[["AllCellTypes"]]),
  var_5 = setNames(rep(0, len),result_list[["AllCellTypes"]])
)
for(indiv_name in result_list[["IndivNames"]]) {
  for(cell in result_list[["CellTypes"]][[indiv_name]]) {
    cell_avg$corr[[cell]] <- cell_avg$corr[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"corr"]
    cell_avg$cos[[cell]] <- cell_avg$cos[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"cos"]
    cell_avg$MSE[[cell]] <- cell_avg$MSE[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"MSE"]
    cell_avg$JSD[[cell]] <- cell_avg$JSD[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"JSD"]
    cell_avg$Euc[[cell]] <- cell_avg$Euc[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"Euc"]
    cell_avg$Manh[[cell]] <- cell_avg$Manh[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"Manh"]
    cell_avg$mean_3[[cell]] <- cell_avg$mean_3[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"mean_3"]
    cell_avg$mean_5[[cell]] <- cell_avg$mean_5[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"mean_5"]
    cell_avg$var_3[[cell]] <- cell_avg$var_3[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"var_3"]
    cell_avg$var_5[[cell]] <- cell_avg$var_5[[cell]] + result_list[[paste0(indiv_name, "_metrics")]][cell,"var_5"]
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

ds <- "" #updates for HVFs and 800 for biased genes

saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_scvi_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_scArch_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_mnnCorrect_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_limma_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_ComBat_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_M3Drop_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_Z_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_fastMNN_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_scanorama_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_scT_v5_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_scT_v5_split_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_logNorm_metrics_meanIndiv.rds"))
saveRDS(data_df, paste0("./Visual Objects/DS",ds,"_logNormScaled_metrics_meanIndiv.rds"))


# --------- Average across all datasets for each cell type by normalization technique ----

ds <- "" #updated_ for HVFs and 800_ for biased genes

mnn <- readRDS(paste0("./Visual Objects/",ds,"mnnCorrect_metrics_meanIndiv.rds"))
limma <- readRDS(paste0("./Visual Objects/",ds,"limma_metrics_meanIndiv.rds"))
CB <- readRDS(paste0("./Visual Objects/",ds,"ComBat_metrics_meanIndiv.rds"))
M3 <- readRDS(paste0("./Visual Objects/",ds,"M3Drop_metrics_meanIndiv.rds"))
Z <- readRDS(paste0("./Visual Objects/",ds,"Z_metrics_meanIndiv.rds"))
fast <- readRDS(paste0("./Visual Objects/",ds,"fastMNN_metrics_meanIndiv.rds"))
scanorama <- readRDS(paste0("./Visual Objects/",ds,"scanorama_metrics_meanIndiv.rds"))
log <- readRDS(paste0("./Visual Objects/",ds,"logNorm_metrics_meanIndiv.rds"))
logS <- readRDS(paste0("./Visual Objects/",ds,"logNormScaled_metrics_meanIndiv.rds"))
scT5 <- readRDS(paste0("./Visual Objects/",ds,"scT_v5_metrics_meanIndiv.rds"))
scT5_s <- readRDS(paste0("./Visual Objects/",ds,"scT_v5_split_metrics_meanIndiv.rds"))
scvi_2 <- readRDS(paste0("./Visual Objects/",ds,"scvi_metrics_meanIndiv.rds"))
sca_2 <- readRDS(paste0("./Visual Objects/",ds,"scArch_metrics_meanIndiv.rds"))

do_avg <- function(tech){
  tech <- tech %>%
    group_by(vector) %>%
    summarise(average = mean(value.Freq, na.rm = TRUE)) %>%
    ungroup()
  return(tech)
}

mnn <- do_avg(mnn)
mnn$norm <- "mnnCorrect"
limma <- do_avg(limma)
limma$norm <- "limma"
CB <- do_avg(CB)
CB$norm <- "ComBat"
M3 <- do_avg(M3)
M3$norm <- "M3Drop"
Z <- do_avg(Z)
Z$norm <- "Z-transform"
fast <- do_avg(fast)
fast$norm <- "fastMNN"
scanorama <- do_avg(scanorama)
scanorama$norm <- "Scanorama"
log <- do_avg(log)
log$norm <- "Log_Normalize"
logS <- do_avg(logS)
logS$norm <- "Log_Normalize_Scale"
scT5 <- do_avg(scT5)
scT5$norm <- "scTransform_v5"
scT5_s <- do_avg(scT5_s)
scT5_s$norm <- "scTransform_v5_split"
scvi_2 <- do_avg(scvi_2)
scvi_2$norm <- "SCVI"
sca_2 <- do_avg(sca_2)
sca_2$norm <- "scArches"

# Combine all techniques into one data frame
all_techniques <- rbind(
  data.frame(norm = "ComBat", CB),
  data.frame(norm = "limma", limma),
  data.frame(norm = "mnnCorrect", mnn),
  data.frame(norm = "fastMNN", fast),
  data.frame(norm = "Scanorama", scanorama),
  data.frame(norm = "Z-transform", Z),
  data.frame(norm = "SCVI", scvi_2),
  data.frame(norm = "scArches", sca_2),
  data.frame(norm = "logNormalize", log),
  data.frame(norm = "logNormalizeScaled", logS),
  data.frame(norm = "scTransform_v5", scT5),
  data.frame(norm = "scTransform_v5_split", scT5_s),
  data.frame(norm = "M3Drop", M3)
)
saveRDS(all_techniques, paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric",tis,".rds"))

## ---- Normalize each technique to its respective reference ----
normalize_technique <- function(technique, reference) {
  
  technique$average_per <- ((technique$average - reference$average) / reference$average) * 100
  
  ## -1 to 1
  # Compute relative improvement for correlation (bounded between [-1,1])
  technique$average_impr <- ifelse(
    technique$average > reference$average,
    ((technique$average - reference$average) / (1 - reference$average)) * 100,  # Improvement Case
    ((technique$average - reference$average) / (reference$average + 1)) * 100    # Worsening Case
  )
  ## 0, 1
  # Compute relative improvement for metrics bounded between [0,1]
  technique$average_impr0 <- ifelse(
    technique$average > reference$average,
    ((technique$average - reference$average) / (1 - reference$average)) * 100,  # Improvement Case
    ((technique$average - reference$average) / reference$average) * 100         # Worsening Case
  )
  return(technique)
}

mnn <- normalize_technique(mnn, log)
limma <- normalize_technique(limma, log)
CB <- normalize_technique(CB, log)
M3 <- normalize_technique(M3, logS)
Z <- normalize_technique(Z, log)
fast <- normalize_technique(fast, log)
scanorama <- normalize_technique(scanorama, log)
scT5 <- normalize_technique(scT5, logS)
scT5_s <- normalize_technique(scT5_s, logS)
scvi_2 <- normalize_technique(scvi_2, log)
sca_2 <- normalize_technique(sca_2, log)


# Combine all normalized techniques into one data frame
together <- rbind(CB, limma, mnn, fast, scanorama, Z, scvi_2, sca_2, scT5, scT5_s, M3)
together$DS <- "" # The dataset number

loc <- "" #specify folder
ds <- paste0(loc, "updated_") # for HVFs or "800_" for biased genes
saveRDS(together, paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline",tis,".rds"))

# --------- Merged all of the datasets for plotting (when summarized to the baseline) ----

# Get all of the files from their respective directories for each dataset

## For HVFs
loc <- "" #specify folder
ds <- paste0(loc, "updated_")
DS1 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_spleen <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline_spleen.rds"))
DS2_liver <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline_liver.rds"))
DS2_kidney <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline_kidney.rds"))
DS2_thymus <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline_thymus.rds"))
DS2_bm <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline_bone marrow.rds"))
DS3 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS4 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS5 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS6 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))

## ALTERNATIVE: For 800 genes
ds <- paste0(loc, "800_") 
DS1 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_spleen <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_liver <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_kidney <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_thymus <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS2_bm <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS3 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS4 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS5 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))
DS6 <- readRDS(paste0("./Visual Objects/",ds,"FINAL_all_techniques_to_baseline.rds"))


split_df <- split(rbind(DS1, DS2_spleen, DS2_liver, DS2_kidney, DS2_thymus, DS2_bm, DS3, DS4, DS5, DS6), 
                  rbind(DS1, DS2_spleen, DS2_liver, DS2_kidney, DS2_thymus, DS2_bm, DS3, DS4, DS5, DS6)$vector)


saveRDS(split_df, paste0("./",ds,"FINAL_%_all_DSandTechniques_to_baseline.rds"))

# --------- Merged all of the datasets for plotting (when NOT summarized to the baseline) ----

## For HVFs
ds <- "" #updated_ for HVFs and 800_ for biased genes
DS1 <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric.rds"))
DS1$DS <- "DS1"
DS2_spleen <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric_spleen.rds"))
DS2_spleen$DS <- "DS2sp"
DS2_liver <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric_liver.rds"))
DS2_liver$DS <- "DS2lv"
DS2_kidney <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric_kidney.rds"))
DS2_kidney$DS <- "DS2kd"
DS2_thymus <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric_thymus.rds"))
DS2_thymus$DS <- "DS2th"
DS2_bm <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric_bone marrow.rds"))
DS2_bm$DS <- "DS2bm"
DS3 <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric.rds"))
DS3$DS <- "DS3"
DS4 <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric.rds"))
DS4$DS <- "DS4"
DS5 <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric.rds"))
DS5$DS <- "DS5"
DS6 <- readRDS(paste0("./Visual Objects/",ds,"all_techniques_averages_per_metric.rds"))
DS6$DS <- "DS6"

split_df <- split(rbind(DS1, DS2_spleen, DS2_liver, DS2_kidney, DS2_thymus, DS2_bm, DS3, DS4, DS5, DS6), 
                  rbind(DS1, DS2_spleen, DS2_liver, DS2_kidney, DS2_thymus, DS2_bm, DS3, DS4, DS5, DS6)$vector)

saveRDS(split_df, "./HVF_FINAL_all_DSandTechniques_NOT_to_baseline.rds")
# OR for 800 genes
saveRDS(split_df, "./800_FINAL_all_DSandTechniques_NOT_to_baseline.rds")



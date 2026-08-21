# Purpose: Characterization of the recurrent genes showing 3' vs 5' protocol-associated bias.
# Includes gene detection, functional enrichment, transcript features, expression patterns, HVF overlap, and recurrence/directionality analyses.

# Loading required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(scales)
library(rtracklayer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ComplexHeatmap)
library(circlize)
library(grid)

set.seed(1234)

### ----- 1] Define the recurrent biased-gene set -----

# Load the recurrent 3' vs 5' gene table generated in UseCase.R
# 5-35 individuals corresponds to the recurrent gene set used in the manuscript.
gene_table <- readRDS("./DS_all_3_vs_5_sign_by_indiv_1e-50.rds")
common_ind <- names(gene_table[gene_table %in% 5:35])
common_ind <- unique(sub("\\..*$", "", common_ind))

# Read the preprocessed datasets (Seurat objects)
# These filenames follow the convention used in Conversion.R.
data1 <- readRDS("./DS1_prep_FINAL.rds")
data2 <- readRDS("./DS2_prep_FINAL.rds")
data3 <- readRDS("./DS3_prep_FINAL.rds")
data4 <- readRDS("./DS4_prep_FINAL.rds")
data5 <- readRDS("./DS5_prep_FINAL.rds")
data6 <- readRDS("./DS6_prep_FINAL.rds")

data_list <- list(
  DS1 = data1,
  DS2 = data2,
  DS3 = data3,
  DS4 = data4,
  DS5 = data5,
  DS6 = data6
)

# Cell-type metadata column used in each dataset
celltype_columns <- c(
  DS1 = "cell_type",
  DS2 = "cell_type_ontology_term_id",
  DS3 = "cell_type_ontology_term_id",
  DS4 = "cell_type",
  DS5 = "cell_type_ontology_term_id",
  DS6 = "cell_type_ontology_term_id"
)

# --- Helper functions
strip_ensembl_version <- function(x) {
  sub("\\..*$", "", x)
}

extract_gene_vector <- function(x) {
  if (is.character(x)) return(x)
  if (is.data.frame(x) && "Genes" %in% colnames(x)) return(x$Genes)
  if (is.list(x) && !is.null(x$Genes)) return(x$Genes)
  stop("The RDS object must contain a character vector or a 'Genes' column.")
}

save_plot <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename = filename, plot = plot, width = width, height = height, units = "in", dpi = 300)
}

prefix_list_names <- function(x, prefix) {
  if (is.null(names(x))) names(x) <- paste0("sample_", seq_along(x))
  names(x) <- paste(prefix, names(x), sep = "::")
  x
}

### ----- 2] Establish the expressed-gene universe across datasets -----

get_detected_genes <- function(seu, min_cell_fraction = 0.01) {
  counts <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts")
  detected <- Matrix::rowSums(counts > 0) > ceiling(min_cell_fraction * ncol(counts))
  strip_ensembl_version(rownames(counts)[detected])
}

detected_gene_sets <- lapply(data_list, function(seu) {
  get_detected_genes(seu, min_cell_fraction = 0.01)
})

detected_genes <- sort(unique(Reduce(union, detected_gene_sets)))

# Reference gene universe from the same GTF/reference family used for alignment.
gtf <- rtracklayer::import("./genes.gtf") # Same GRCh38 GTF used for quantification

gene_universe <- as.data.frame(gtf) %>%
  filter(type == "gene") %>%
  transmute(
    gene_id = strip_ensembl_version(gene_id),
    gene_name = gene_name,
    gene_type = gene_type
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

not_detected_genes <- setdiff(gene_universe$gene_id, detected_genes)

detection_summary <- tibble(
  n_reference_genes = nrow(gene_universe),
  n_detected_genes = sum(gene_universe$gene_id %in% detected_genes),
  percent_detected = 100 * n_detected_genes / n_reference_genes,
  n_not_detected = n_reference_genes - n_detected_genes
)

write.csv(detection_summary, "./gene_detection_summary.csv", row.names = FALSE)

not_detected_df <- gene_universe %>%
  filter(gene_id %in% not_detected_genes) %>%
  count(gene_type, name = "n") %>%
  mutate(
    percent = 100 * n / sum(n),
    label = paste0(scales::comma(n), " (", round(percent, 1), "%)"),
    inside = n > 0.08 * max(n),
    highlight = if_else(gene_type == "protein_coding", "protein_coding", "other"),
    gene_type = reorder(gene_type, n)
  )

p_not_detected <- ggplot(not_detected_df, aes(x = n, y = gene_type)) +
  geom_col(aes(fill = highlight), color = "black", width = 0.85) +
  scale_fill_manual(
    values = c(protein_coding = "#F3A153", other = "#95AAD3"),
    guide = "none"
  ) +
  geom_text(
    data = subset(not_detected_df, inside),
    aes(label = label),
    hjust = 1.05,
    fontface = "bold",
    size = 4
  ) +
  geom_text(
    data = subset(not_detected_df, !inside),
    aes(label = label),
    hjust = -0.1,
    fontface = "bold",
    size = 4
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of genes", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", face = "bold", size = 12),
    axis.title = element_text(color = "black", face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_not_detected, "not_detected_gene_types.png", width = 6.5, height = 6.5)

### ----- 3] Functional enrichment of undetected and recurrent biased genes -----

undetected_protein_coding <- gene_universe %>%
  filter(gene_id %in% not_detected_genes, gene_type == "protein_coding") %>%
  pull(gene_id)

ego_undetected <- enrichGO(
  gene = undetected_protein_coding,
  universe = gene_universe$gene_id[gene_universe$gene_type == "protein_coding"],
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

write.csv(
  as.data.frame(ego_undetected),
  "./not_detected_protein_coding_enrichGO.csv",
  row.names = FALSE
)

# The manuscript groups enriched terms into broad biological categories for a
# compact summary plot. This is a deterministic keyword-based categorization.
ego_categories <- as.data.frame(ego_undetected) %>%
  filter(p.adjust < 0.05, FoldEnrichment > 2) %>%
  mutate(
    category = case_when(
      str_detect(Description, regex("olfactory|smell|chemical stimulus|chemosensory", ignore_case = TRUE)) ~ "Olfactory / Chemosensory",
      str_detect(Description, regex("taste", ignore_case = TRUE)) ~ "Taste",
      str_detect(Description, regex("cilium|cilia|flagellum|flagellar|axoneme", ignore_case = TRUE)) ~ "Cilia / Flagella",
      str_detect(Description, regex("sperm|spermatid|spermatogenesis|fertilization|germ cell|reproduction|gonadotropin|acrosome|meiot|gamete", ignore_case = TRUE)) ~ "Reproduction / Spermatogenesis",
      str_detect(Description, regex("neuron|neuronal|neurotransmitter|synap|axon|dendrite|neuropeptide|postsynaptic|nervous|catecholamine|drinking|dopamine|serotonin|feeding|eating|startle|G protein-coupled|GPCR|glutamate receptor|AMPA|GABA|gamma-aminobutyric|behavior|neuroendocrine", ignore_case = TRUE)) ~ "Neuronal / Synaptic",
      str_detect(Description, regex("development|morphogenesis|pattern|regionalization|diencephalon|embryo|sex|cell fate specification|spinal cord|muscle cell fate|cell commitment", ignore_case = TRUE)) ~ "Development / Morphogenesis",
      str_detect(Description, regex("potassium|sodium|nerve|calcium|ion transport|membrane potential|channel|arachidonate|hormone", ignore_case = TRUE)) ~ "Ion transport / Excitability",
      str_detect(Description, regex("visual|visible|vision|phototransduction|light stimulus|eye|retina|photoreceptor", ignore_case = TRUE)) ~ "Vision / Eye",
      str_detect(Description, regex("keratin|epiderm|skin|intermediate filament", ignore_case = TRUE)) ~ "Skin / Epidermis",
      str_detect(Description, regex("adhesion|cadherin|filament|membrane", ignore_case = TRUE)) ~ "Cell movement / Adhesion",
      TRUE ~ "Other"
    )
  )

category_summary <- ego_categories %>%
  count(category, name = "n_pathways") %>%
  filter(category != "Other") %>%
  arrange(n_pathways) %>%
  mutate(category = factor(category, levels = category))

p_go_categories <- ggplot(category_summary, aes(x = n_pathways, y = category)) +
  geom_col(fill = "#F3A153", color = "black", width = 0.85) +
  geom_text(aes(label = n_pathways), hjust = -0.1, size = 4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Number of enriched GO terms", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", face = "bold", size = 12),
    axis.title = element_text(color = "black", face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_go_categories, "not_detected_gene_GO_categories.png", width = 6.5, height = 6.5)

ego_biased <- enrichGO(
  gene = common_ind,
  universe = detected_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL"
)

p_biased_go <- dotplot(ego_biased) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 12, color = "black"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_biased_go, "common_ind_GO.png", width = 7, height = 6.5)
write.csv(as.data.frame(ego_biased), "./common_ind_enrichGO.csv", row.names = FALSE)

### ----- 4] GC content and longest protein-coding transcript length -----

edb <- EnsDb.Hsapiens.v86
genome <- BSgenome.Hsapiens.UCSC.hg38

# Matched random gene sets used for the cosine-similarity comparison in Figures.R
rand1 <- readRDS("./DS2bm_cos_800_rand_prop_SCALED.rds")
rand2 <- readRDS("./DS2bm_cos_800_rand_prop_2_SCALED.rds")
rand3 <- readRDS("./DS2bm_cos_800_rand_prop_3_SCALED.rds")

random_gene_sets <- list(
  `Random 1` = unique(strip_ensembl_version(extract_gene_vector(rand1))),
  `Random 2` = unique(strip_ensembl_version(extract_gene_vector(rand2))),
  `Random 3` = unique(strip_ensembl_version(extract_gene_vector(rand3)))
)

gene_groups <- c(list(Biased = common_ind), random_gene_sets)

# Length of the transcript-end sequence used for GC-content calculations.
# The final manuscript analysis used 250 bp. Change this value to repeat the
# 500-bp or 1-kb sensitivity analyses.
gc_width <- 250

gc_at_transcript_end <- function(exons_gr, genome, end = c("3prime", "5prime"), width = 250) {
  end <- match.arg(end)

  seqlevelsStyle(exons_gr) <- seqlevelsStyle(genome)[1]
  exons_gr <- sort(exons_gr, ignore.strand = FALSE)

  strand_value <- unique(as.character(strand(exons_gr)))
  if (length(strand_value) != 1 || strand_value == "*") {
    return(NA_real_)
  }

  exon_sequences <- Biostrings::getSeq(genome, exons_gr)
  spliced <- do.call(Biostrings::xscat, as.list(exon_sequences))

  if (strand_value == "-") {
    spliced <- Biostrings::reverseComplement(spliced)
  }

  n <- min(width, length(spliced))
  region <- if (end == "5prime") {
    Biostrings::subseq(spliced, start = 1, width = n)
  } else {
    Biostrings::subseq(spliced, start = length(spliced) - n + 1, width = n)
  }

  Biostrings::letterFrequency(region, letters = "GC", as.prob = TRUE)[1]
}

summarize_gene_features <- function(gene_id, edb, genome, gc_width = 250) {
  gene_id <- strip_ensembl_version(gene_id)

  txs <- ensembldb::transcripts(
    edb,
    filter = AnnotationFilter::GeneIdFilter(gene_id),
    columns = c("tx_id", "tx_biotype")
  )

  if (length(txs) == 0) {
    return(tibble(
      gene_id = gene_id,
      selected_tx_id = NA_character_,
      n_protein_coding_tx = 0L,
      transcript_length_bp = NA_real_,
      n_exons = NA_integer_,
      GC_5prime = NA_real_,
      GC_3prime = NA_real_
    ))
  }

  protein_coding_tx <- txs$tx_id[txs$tx_biotype == "protein_coding"]
  if (length(protein_coding_tx) == 0) {
    return(tibble(
      gene_id = gene_id,
      selected_tx_id = NA_character_,
      n_protein_coding_tx = 0L,
      transcript_length_bp = NA_real_,
      n_exons = NA_integer_,
      GC_5prime = NA_real_,
      GC_3prime = NA_real_
    ))
  }

  exons_by_tx <- GenomicFeatures::exonsBy(
    edb,
    by = "tx",
    filter = AnnotationFilter::GeneIdFilter(gene_id)
  )
  exons_by_tx <- exons_by_tx[names(exons_by_tx) %in% protein_coding_tx]

  if (length(exons_by_tx) == 0) {
    return(tibble(
      gene_id = gene_id,
      selected_tx_id = NA_character_,
      n_protein_coding_tx = 0L,
      transcript_length_bp = NA_real_,
      n_exons = NA_integer_,
      GC_5prime = NA_real_,
      GC_3prime = NA_real_
    ))
  }

  tx_lengths <- vapply(
    exons_by_tx,
    function(exons_gr) sum(width(GenomicRanges::reduce(exons_gr))),
    numeric(1)
  )

  selected_tx_id <- names(tx_lengths)[which.max(tx_lengths)]
  selected_exons <- exons_by_tx[[selected_tx_id]]

  tibble(
    gene_id = gene_id,
    selected_tx_id = selected_tx_id,
    n_protein_coding_tx = length(exons_by_tx),
    transcript_length_bp = unname(tx_lengths[selected_tx_id]),
    n_exons = length(selected_exons),
    GC_5prime = gc_at_transcript_end(selected_exons, genome, "5prime", gc_width),
    GC_3prime = gc_at_transcript_end(selected_exons, genome, "3prime", gc_width)
  )
}

gene_features <- imap_dfr(gene_groups, function(genes, group_name) {
  bind_rows(lapply(genes, summarize_gene_features, edb = edb, genome = genome, gc_width = gc_width)) %>%
    mutate(Group = group_name)
})

write.csv(gene_features, "./gene_sequence_features.csv", row.names = FALSE)

# Paired comparison of 5' and 3' GC content within each gene group.
gc_tests <- lapply(split(gene_features, gene_features$Group), function(x) {
  x <- x[complete.cases(x$GC_5prime, x$GC_3prime), ]

  if (nrow(x) == 0) {
    return(tibble(
      n_genes = 0L,
      median_GC_5prime = NA_real_,
      median_GC_3prime = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  test <- wilcox.test(x$GC_3prime, x$GC_5prime, paired = TRUE, exact = FALSE)
  tibble(
    n_genes = nrow(x),
    median_GC_5prime = median(x$GC_5prime),
    median_GC_3prime = median(x$GC_3prime),
    statistic = unname(test$statistic),
    p_value = test$p.value
  )
}) %>%
  bind_rows(.id = "Group")

write.csv(gc_tests, "./GC_5prime_vs_3prime_tests.csv", row.names = FALSE)

gc_plot_df <- gene_features %>%
  select(Group, gene_id, GC_5prime, GC_3prime) %>%
  pivot_longer(
    cols = c(GC_5prime, GC_3prime),
    names_to = "Region",
    values_to = "GC"
  ) %>%
  mutate(
    Region = recode(Region, GC_5prime = "5'", GC_3prime = "3'"),
    Region = factor(Region, levels = c("5'", "3'"))
  )

p_gc <- ggplot(gc_plot_df, aes(x = Region, y = GC)) +
  geom_boxplot(
    fill = "#9CBFAF",
    color = "black",
    linewidth = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.4) +
  facet_wrap(~Group) +
  labs(x = NULL, y = paste0("GC fraction (terminal ", gc_width, " bp)")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black", face = "bold", size = 11),
    axis.title = element_text(color = "black", face = "bold", size = 13),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_gc, "GC_content_by_gene_set.png", width = 8, height = 5)

p_length <- gene_features %>%
  filter(!is.na(transcript_length_bp)) %>%
  ggplot(aes(x = Group, y = transcript_length_bp)) +
  geom_boxplot(
    fill = "#95AAD3",
    color = "black",
    linewidth = 0.8,
    width = 0.6,
    outlier.size = 1.5,
    outlier.alpha = 0.6
  ) +
  labs(x = NULL, y = "Longest protein-coding transcript length (bp)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black", face = "bold", size = 11),
    axis.title = element_text(color = "black", face = "bold", size = 13),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_length, "longest_protein_coding_transcript_length.png", width = 6.5, height = 6)

### ----- 5] Expression-rank heatmap and non-immune specificity -----

assign_quantiles <- function(expr, n_bins = 4) {
  bins <- rep(0L, length(expr))
  names(bins) <- names(expr)

  keep <- expr > 0 & !is.na(expr)
  if (!any(keep)) {
    return(bins)
  }

  ranks <- rank(expr[keep], ties.method = "average")
  bins[keep] <- ceiling(ranks / max(ranks) * n_bins)
  bins
}

count_geneset_bins <- function(bins, geneset, n_bins = 4) {
  genes <- intersect(names(bins), geneset)
  table(factor(bins[genes], levels = 0:n_bins))
}

get_heatmap_counts <- function(seu, geneset, dataset_name, n_bins = 4) {
  # Only within-sample gene ranks are used, so ranking aggregated counts directly
  # is equivalent to ranking any monotonic library-size-normalized transform.
  aggregate_counts <- AggregateExpression(
    seu,
    group.by = "new_id",
    assays = "RNA",
    return.seurat = FALSE
  )$RNA

  quantile_bins <- apply(aggregate_counts, 2, assign_quantiles, n_bins = n_bins)

  out <- sapply(
    colnames(quantile_bins),
    function(individual) count_geneset_bins(quantile_bins[, individual], geneset, n_bins)
  )

  out <- as.data.frame(t(out))
  colnames(out) <- as.character(0:n_bins)
  rownames(out) <- gsub("_", "-", rownames(out))
  out$ds <- dataset_name
  out
}

get_immune_fraction <- function(seu, dataset_name, celltype_col, non_immune_terms) {
  meta <- seu@meta.data

  # DS2 contains multiple tissues for some individuals. The original analysis
  # represented each individual by the tissue contributing the largest cell count.
  if (dataset_name == "DS2" && "tissue" %in% colnames(meta)) {
    out <- meta %>%
      mutate(new_id = as.factor(new_id)) %>%
      group_by(new_id, tissue) %>%
      summarise(
        immune_pct = mean(!(.data[[celltype_col]] %in% non_immune_terms)),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      group_by(new_id) %>%
      slice_max(order_by = n_cells, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(new_id, immune_pct)
  } else {
    out <- meta %>%
      group_by(new_id) %>%
      summarise(
        immune_pct = mean(!(.data[[celltype_col]] %in% non_immune_terms)),
        .groups = "drop"
      )
  }

  out %>%
    mutate(
      row_id = gsub("_", "-", as.character(new_id)),
      ds = dataset_name
    ) %>%
    select(row_id, ds, immune_pct)
}

classify_nonimmune_expression <- function(
    seu,
    celltype_col,
    non_immune_terms,
    min_cells_per_type = 5,
    min_nonimmune_fraction = 0.30,
    max_immune_fraction = 0.10) {

  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  celltypes <- as.character(seu@meta.data[[celltype_col]])
  names(celltypes) <- colnames(seu)

  celltype_counts <- table(celltypes)
  valid_celltypes <- names(celltype_counts[celltype_counts >= min_cells_per_type])
  nonimmune_celltypes <- intersect(non_immune_terms, valid_celltypes)
  immune_celltypes <- setdiff(valid_celltypes, nonimmune_celltypes)

  high_in_nonimmune <- if (length(nonimmune_celltypes) == 0) {
    character(0)
  } else {
    Reduce(union, lapply(nonimmune_celltypes, function(ct) {
      cells <- names(celltypes)[celltypes == ct]
      detected <- Matrix::rowSums(counts[, cells, drop = FALSE] > 0) >
        ceiling(min_nonimmune_fraction * length(cells))
      rownames(counts)[detected]
    }))
  }

  low_in_immune <- if (length(immune_celltypes) == 0) {
    character(0)
  } else {
    Reduce(intersect, lapply(immune_celltypes, function(ct) {
      cells <- names(celltypes)[celltypes == ct]
      low <- Matrix::rowSums(counts[, cells, drop = FALSE] > 0) <
        ceiling(max_immune_fraction * length(cells))
      rownames(counts)[low]
    }))
  }

  list(
    high_in_nonimmune = strip_ensembl_version(high_in_nonimmune),
    low_in_immune = strip_ensembl_version(low_in_immune)
  )
}

non_immune_terms <- c(
  "endothelial cell",
  "fibroblast",
  "osteoblast",
  "smooth muscle cell",
  "hepatocyte",
  "epithelial cell of nephron",
  "endothelial",
  "lung epithelial cell",
  "malignant cell",
  "RBCs"
)

heatmap_counts_list <- list()
immune_fraction_list <- list()
nonimmune_expression_results <- list()

for (dataset_name in names(data_list)) {
  seu <- data_list[[dataset_name]]
  celltype_col <- celltype_columns[[dataset_name]]

  heatmap_counts_list[[dataset_name]] <- get_heatmap_counts(
    seu = seu,
    geneset = common_ind,
    dataset_name = dataset_name,
    n_bins = 4
  )

  immune_fraction_list[[dataset_name]] <- get_immune_fraction(
    seu = seu,
    dataset_name = dataset_name,
    celltype_col = celltype_col,
    non_immune_terms = non_immune_terms
  )

  # Only datasets containing at least one predefined non-immune population
  # contribute to the non-immune-specific gene analysis.
  if (any(seu@meta.data[[celltype_col]] %in% non_immune_terms)) {
    nonimmune_expression_results[[dataset_name]] <- classify_nonimmune_expression(
      seu = seu,
      celltype_col = celltype_col,
      non_immune_terms = non_immune_terms
    )
  }

  rm(seu)
  gc()
}

# bind_rows() does not preserve row names reliably; rebuild row IDs explicitly.
heatmap_counts <- imap_dfr(heatmap_counts_list, function(x, dataset_name) {
  x$row_id <- rownames(x)
  x
})

immune_fraction_df <- bind_rows(immune_fraction_list)

dataset_order <- names(data_list)

heatmap_counts <- heatmap_counts %>%
  left_join(immune_fraction_df, by = c("row_id", "ds")) %>%
  mutate(ds = factor(ds, levels = dataset_order)) %>%
  arrange(ds, desc(immune_pct))

write.csv(heatmap_counts, "./biased_gene_expression_quantile_counts.csv", row.names = FALSE)

# A gene is considered non-immune specific in this summary when it is observed
# in >30% of cells in at least one non-immune population and in <10% of cells in
# all immune populations, allowing evidence to come from any eligible dataset.
biased_low <- Reduce(
  union,
  lapply(nonimmune_expression_results, function(x) intersect(x$low_in_immune, common_ind))
)
biased_high <- Reduce(
  union,
  lapply(nonimmune_expression_results, function(x) intersect(x$high_in_nonimmune, common_ind))
)
biased_nonimmune_specific <- intersect(biased_low, biased_high)

all_low <- Reduce(
  union,
  lapply(nonimmune_expression_results, function(x) intersect(x$low_in_immune, detected_genes))
)
all_high <- Reduce(
  union,
  lapply(nonimmune_expression_results, function(x) intersect(x$high_in_nonimmune, detected_genes))
)
all_nonimmune_specific <- intersect(all_low, all_high)

nonimmune_summary <- tibble(
  gene_set = c("Biased genes", "Detected genes"),
  total_genes = c(length(common_ind), length(detected_genes)),
  nonimmune_specific_genes = c(length(biased_nonimmune_specific), length(all_nonimmune_specific))
) %>%
  mutate(percent_nonimmune_specific = 100 * nonimmune_specific_genes / total_genes)

write.csv(nonimmune_summary, "./nonimmune_specific_gene_summary.csv", row.names = FALSE)
writeLines(biased_nonimmune_specific, "./biased_nonimmune_specific_genes.txt")

# Heatmap of biased-gene expression ranks across individuals.
to_plot <- as.matrix(heatmap_counts[, as.character(0:4)])
rownames(to_plot) <- heatmap_counts$row_id
colnames(to_plot) <- c("Expr. 0", "Q1", "Q2", "Q3", "Q4")

n_biased <- length(common_ind)
col_fun <- circlize::colorRamp2(
  c(0, n_biased / 2, n_biased),
  c("white", "#E57710", "#AB590D")
)

png("./biased_gene_expression_heatmap.png", width = 3100, height = 2000, res = 300)
ht <- Heatmap(
  to_plot,
  name = "Count",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%d", to_plot[i, j]),
      x,
      y,
      gp = gpar(
        fontsize = 13,
        col = ifelse(to_plot[i, j] > 0.58 * n_biased, "white", "#3A3A3A")
      )
    )
  },
  left_annotation = rowAnnotation(
    Dataset = heatmap_counts$ds,
    `Immune fraction` = anno_barplot(
      heatmap_counts$immune_pct,
      bar_width = 0.8,
      border = FALSE,
      gp = gpar(fill = "#F5F5F5")
    ),
    col = list(
      Dataset = c(
        DS1 = "#23517E",
        DS2 = "#51ABE2",
        DS3 = "#80CED8",
        DS4 = "#F5C248",
        DS5 = "#CE4B4F",
        DS6 = "#801F3D"
      )
    )
  ),
  heatmap_legend_param = list(
    title = "Count",
    at = round(seq(0, n_biased, length.out = 5)),
    legend_height = unit(4, "cm")
  )
)
draw(ht)
dev.off()

plot_nonimmune_fraction <- function(total, specific, label) {
  df <- tibble(
    category = factor(c("Other", "Non-immune specific"), levels = c("Other", "Non-immune specific")),
    n = c(total - specific, specific)
  )

  ggplot(df, aes(x = "", y = n, fill = category)) +
    geom_col(color = "black", width = 0.5) +
    scale_fill_manual(values = c("Other" = "#95AAD3", "Non-immune specific" = "#D95F02")) +
    annotate(
      "text",
      x = 1,
      y = total * 1.05,
      label = paste0(scales::comma(specific), " genes (", round(100 * specific / total, 1), "%)"),
      fontface = "bold",
      size = 5
    ) +
    coord_cartesian(ylim = c(0, total * 1.15)) +
    labs(x = NULL, y = "Number of genes", fill = NULL, title = label) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(color = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}

save_plot(
  plot_nonimmune_fraction(length(common_ind), length(biased_nonimmune_specific), "Biased genes"),
  "nonimmune_specific_common_ind.png",
  width = 5.5,
  height = 6.5
)

save_plot(
  plot_nonimmune_fraction(length(detected_genes), length(all_nonimmune_specific), "All detected genes"),
  "nonimmune_specific_all_detected_genes.png",
  width = 5.5,
  height = 6.5
)

### ----- 6] Overlap of recurrent biased genes with highly variable features -----

# The Metrics script saves the highly variable features as ./hvf.rds for each
# dataset/tissue. Repeat this block using the corresponding hvf.rds file.
hvf <- readRDS("./hvf.rds")
hvf_overlap <- intersect(common_ind, strip_ensembl_version(hvf))
length(hvf_overlap)

### ----- 7] UMAP and mixing metrics after excluding recurrent biased genes -----
# Specify the dataset to reproduce the corresponding UMAP/mixing analysis.
# split RNA layers by individual, normalize and identify variable features, join
# layers, then evaluate each individual separately.

obj <- readRDS("./prep_FINAL.rds") # Specify the dataset
DefaultAssay(obj) <- "RNA"
keep_features <- rownames(obj)[!strip_ensembl_version(rownames(obj)) %in% common_ind]
obj <- subset(obj, features = keep_features)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$new_id)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- JoinLayers(obj)


obj_by_individual <- SplitObject(obj, split.by = "new_id")
umap_plots <- list()
mixing_metrics <- vector("list", length(obj_by_individual))
names(mixing_metrics) <- names(obj_by_individual)

for (individual in names(obj_by_individual)) {
  x <- obj_by_individual[[individual]]
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 40, verbose = FALSE)
  x <- RunUMAP(x, dims = 1:20, reduction = "pca", verbose = FALSE)

  umap_plots[[individual]] <- DimPlot(
    x,
    reduction = "umap",
    group.by = c("assay", "cell_type_ontology_term_id")
  ) +
    ggtitle(paste("Sample", individual))

  mix_assay <- MixingMetric(
    x,
    grouping.var = "assay",
    reduction = "pca",
    dims = 1:20,
    k = 5,
    max.k = 150,
    eps = 0,
    verbose = FALSE
  )

  mix_celltype <- MixingMetric(
    x,
    grouping.var = "cell_type_ontology_term_id",
    reduction = "pca",
    dims = 1:20,
    k = 5,
    max.k = 150,
    eps = 0,
    verbose = FALSE
  )

  mixing_metrics[[individual]] <- tibble(
    individual = individual,
    assay_mixing = median(150 - mix_assay),
    celltype_mixing = median(150 - mix_celltype)
  )

  obj_by_individual[[individual]] <- x
}

pdf("./DS1_UMAP_without_biased_genes.pdf", width = 12, height = 5)
for (p in umap_plots) print(p)
dev.off()

mixing_metrics <- bind_rows(mixing_metrics)
write.csv(mixing_metrics, "./mixing_metrics_without_biased_genes.csv", row.names = FALSE)

### ----- 8] Recurrence and directionality of significant 3' vs 5' differences -----

p_cutoff <- 1e-50

# Load the per-individual 3' vs 5' differential-expression results generated in UseCase.R
ds1 <- readRDS("./DS1_3_vs_5_sign_by_indiv.rds")
ds2bm <- readRDS("./DS2_3_vs_5_sign_by_indiv_bm.rds")
ds2sp <- readRDS("./DS2_3_vs_5_sign_by_indiv_sp.rds")
ds2kd <- readRDS("./DS2_3_vs_5_sign_by_indiv_kd.rds")
ds2lv <- readRDS("./DS2_3_vs_5_sign_by_indiv_lv.rds")
ds2th <- readRDS("./DS2_3_vs_5_sign_by_indiv_th.rds")
ds3 <- readRDS("./DS3_3_vs_5_sign_by_indiv.rds")
ds4 <- readRDS("./DS4_3_vs_5_sign_by_indiv.rds")
ds5 <- readRDS("./DS5_3_vs_5_sign_by_indiv.rds")
ds6 <- readRDS("./DS6_3_vs_5_sign_by_indiv.rds")

significant_results <- c(
  prefix_list_names(ds1, "DS1"),
  prefix_list_names(ds2bm, "DS2_BM"),
  prefix_list_names(ds2sp, "DS2_SP"),
  prefix_list_names(ds2kd, "DS2_KD"),
  prefix_list_names(ds2lv, "DS2_LV"),
  prefix_list_names(ds2th, "DS2_TH"),
  prefix_list_names(ds3, "DS3"),
  prefix_list_names(ds4, "DS4"),
  prefix_list_names(ds5, "DS5"),
  prefix_list_names(ds6, "DS6")
)

get_de_long <- function(result_list) {
  imap_dfr(result_list, function(res, sample_id) {
    parts <- str_split_fixed(sample_id, "::", 2)

    if (!"avg_log2FC" %in% colnames(res)) {
      stop("Expected an 'avg_log2FC' column in every differential-expression result table.")
    }
    if (!"p_val_adj" %in% colnames(res)) {
      stop("Expected a 'p_val_adj' column in every differential-expression result table.")
    }

    tibble(
      gene = strip_ensembl_version(rownames(res)),
      dataset = parts[, 1],
      individual = parts[, 2],
      sample_id = sample_id,
      p_val_adj = res$p_val_adj,
      logFC = res$avg_log2FC
    )
  })
}

de_long <- get_de_long(significant_results) %>%
  mutate(
    significant = !is.na(p_val_adj) & p_val_adj < p_cutoff,
    direction = case_when(
      significant & logFC > 0 ~ "3prime biased",
      significant & logFC < 0 ~ "5prime biased",
      TRUE ~ "Not significant"
    )
  )

# Number of significant genes per individual. Starting from the full sample list
# automatically retains individuals with zero significant genes.
all_samples <- tibble(sample_id = names(significant_results)) %>%
  separate(sample_id, into = c("dataset", "individual"), sep = "::", remove = FALSE)

individual_counts <- de_long %>%
  filter(significant) %>%
  distinct(sample_id, dataset, gene) %>%
  count(sample_id, dataset, name = "n_significant") %>%
  right_join(all_samples %>% select(sample_id, dataset), by = c("sample_id", "dataset")) %>%
  mutate(n_significant = replace_na(n_significant, 0L))

write.csv(individual_counts, "./significant_gene_counts_per_individual.csv", row.names = FALSE)

p_individual_counts <- ggplot(
  individual_counts,
  aes(x = reorder(sample_id, n_significant), y = n_significant, fill = dataset)
) +
  geom_col(color = "black", width = 0.85) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Number of significant genes (" * FDR < 10^-50 * ")"),
    fill = "Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_individual_counts, "significant_genes_per_individual.png", width = 8.7, height = 6.5)

# Number of individuals supporting each recurrent biased gene.
gene_individual_recurrence <- de_long %>%
  filter(significant, gene %in% common_ind) %>%
  distinct(sample_id, gene) %>%
  count(gene, name = "n_individuals")

recurrence_plot_df <- gene_individual_recurrence %>%
  arrange(n_individuals) %>%
  mutate(rank = row_number())

p_gene_recurrence <- ggplot(recurrence_plot_df, aes(x = rank, y = n_individuals)) +
  geom_line(linewidth = 0.8, color = "#2F5597") +
  geom_point(size = 0.6, alpha = 0.5, color = "#2F5597") +
  labs(x = "Genes ranked by recurrence", y = "Number of individuals") +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
  )

save_plot(p_gene_recurrence, "biased_gene_recurrence_by_individual.png", width = 7.7, height = 6)

# Number of dataset/tissue units supporting each gene.
gene_dataset_recurrence <- de_long %>%
  filter(significant) %>%
  distinct(gene, dataset) %>%
  count(gene, name = "n_datasets")

recurrence_biased <- tibble(gene = common_ind) %>%
  left_join(gene_individual_recurrence, by = "gene") %>%
  left_join(gene_dataset_recurrence, by = "gene") %>%
  mutate(
    n_individuals = replace_na(n_individuals, 0L),
    n_datasets = replace_na(n_datasets, 0L)
  )

write.csv(recurrence_biased, "./biased_gene_recurrence_summary.csv", row.names = FALSE)

# Directional consistency among significant observations.
direction_biased <- de_long %>%
  filter(significant, gene %in% common_ind) %>%
  group_by(gene) %>%
  summarise(
    n_individuals = n_distinct(sample_id),
    n_datasets = n_distinct(dataset),
    n_3prime = sum(direction == "3prime biased"),
    n_5prime = sum(direction == "5prime biased"),
    dominant_direction = case_when(
      n_3prime > n_5prime ~ "3prime biased",
      n_5prime > n_3prime ~ "5prime biased",
      TRUE ~ "Tie"
    ),
    direction_consistency = pmax(n_3prime, n_5prime) / (n_3prime + n_5prime),
    .groups = "drop"
  )

write.csv(direction_biased, "./biased_gene_direction_consistency.csv", row.names = FALSE)

p_direction_consistency <- ggplot(direction_biased, aes(x = direction_consistency)) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 1,
    fill = "#9CBFAF",
    color = "black"
  ) +
  labs(
    x = "Proportion of significant individuals supporting dominant direction",
    y = "Number of biased genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  )

save_plot(p_direction_consistency, "biased_gene_direction_consistency.png", width = 7, height = 6)

### ----- 9] Direction of biased-gene fold changes without a significance threshold -----

ds1_all <- readRDS("./DS1_3_vs_5_all_by_indiv.rds")
ds2bm_all <- readRDS("./DS2_3_vs_5_all_by_indiv_bm.rds")
ds2sp_all <- readRDS("./DS2_3_vs_5_all_by_indiv_sp.rds")
ds2kd_all <- readRDS("./DS2_3_vs_5_all_by_indiv_kd.rds")
ds2lv_all <- readRDS("./DS2_3_vs_5_all_by_indiv_lv.rds")
ds2th_all <- readRDS("./DS2_3_vs_5_all_by_indiv_th.rds")
ds3_all <- readRDS("./DS3_3_vs_5_all_by_indiv.rds")
ds4_all <- readRDS("./DS4_3_vs_5_all_by_indiv.rds")
ds5_all <- readRDS("./DS5_3_vs_5_all_by_indiv.rds")
ds6_all <- readRDS("./DS6_3_vs_5_all_by_indiv.rds")

all_results <- c(
  prefix_list_names(ds1_all, "DS1"),
  prefix_list_names(ds2bm_all, "DS2_BM"),
  prefix_list_names(ds2sp_all, "DS2_SP"),
  prefix_list_names(ds2kd_all, "DS2_KD"),
  prefix_list_names(ds2lv_all, "DS2_LV"),
  prefix_list_names(ds2th_all, "DS2_TH"),
  prefix_list_names(ds3_all, "DS3"),
  prefix_list_names(ds4_all, "DS4"),
  prefix_list_names(ds5_all, "DS5"),
  prefix_list_names(ds6_all, "DS6")
)
all_de_long <- get_de_long(all_results) %>%
  filter(gene %in% common_ind)

gene_counts <- all_de_long %>%
  group_by(sample_id) %>%
  summarise(
    higher3 = sum(logFC > 0, na.rm = TRUE),
    higher5 = sum(logFC < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(higher3, higher5),
    names_to = "direction",
    values_to = "count"
  )

write.csv(gene_counts, "./biased_gene_direction_counts_no_significance_threshold.csv", row.names = FALSE)

p_gene_direction <- ggplot(gene_counts, aes(x = sample_id, y = count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.95) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 2.5
  ) +
  scale_fill_manual(
    values = c(higher3 = "#e31a1c", higher5 = "#1f78b4"),
    labels = c(higher3 = "Higher in 3'", higher5 = "Higher in 5'")
  ) +
  labs(
    x = NULL,
    y = "Number of biased genes",
    fill = NULL,
    title = "Direction of biased-gene fold changes per individual"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 12),
    axis.line = element_line(color = "black", linewidth = 1),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.y = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25)
  )

save_plot(p_gene_direction, "biased_gene_direction_no_significance_threshold.png", width = 8, height = 6.3)

higher3_genes <- split(all_de_long$gene[all_de_long$logFC > 0], all_de_long$sample_id[all_de_long$logFC > 0])
higher5_genes <- split(all_de_long$gene[all_de_long$logFC < 0], all_de_long$sample_id[all_de_long$logFC < 0])

# Ensure samples with no genes in one direction are represented by empty vectors.
all_sample_ids <- names(all_results)
get_direction_genes <- function(x, id) {
  if (is.null(x[[id]])) character(0) else unique(x[[id]])
}
higher3_genes <- setNames(lapply(all_sample_ids, function(id) get_direction_genes(higher3_genes, id)), all_sample_ids)
higher5_genes <- setNames(lapply(all_sample_ids, function(id) get_direction_genes(higher5_genes, id)), all_sample_ids)

common_higher3 <- Reduce(intersect, higher3_genes)
common_higher5 <- Reduce(intersect, higher5_genes)

higher3_frequency <- sort(table(unlist(higher3_genes)), decreasing = TRUE)
higher5_frequency <- sort(table(unlist(higher5_genes)), decreasing = TRUE)

# Threshold used for the incomplete-overlap summary. Change if a different
# recurrence threshold is reported in the manuscript/supplement.
min_directional_samples <- 30

higher3_atleast_n <- names(higher3_frequency[higher3_frequency >= min_directional_samples])
higher5_atleast_n <- names(higher5_frequency[higher5_frequency >= min_directional_samples])

directional_overlap_summary <- tibble(
  category = c(
    "Higher in 3' in all samples",
    "Higher in 5' in all samples",
    paste0("Higher in 3' in >=", min_directional_samples, " samples"),
    paste0("Higher in 5' in >=", min_directional_samples, " samples")
  ),
  n_genes = c(
    length(common_higher3),
    length(common_higher5),
    length(higher3_atleast_n),
    length(higher5_atleast_n)
  )
)

write.csv(directional_overlap_summary, "./directional_overlap_summary.csv", row.names = FALSE)
writeLines(common_higher3, "./genes_higher_in_3prime_all_samples.txt")
writeLines(common_higher5, "./genes_higher_in_5prime_all_samples.txt")


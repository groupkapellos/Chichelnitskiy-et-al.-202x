# ============================================================
# Compositional analysis using propeller
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(speckle)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_root <- "D:/Empysema_Fibrosis_Revisions/Compositional_analysis_tables"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

sample_col <- "LuT"
condition_col <- "disease"
cluster_col <- "general_final_annotation"

remove_celltypes <- c(
  "mixed cells", "Mixed cells",
  "Doublet", "Doublets", "doublet"
)

# ------------------------------------------------------------
# Build SingleCellExperiment object for propeller
# ------------------------------------------------------------

build_sce_for_propeller <- function(seu,
                                    sample_col,
                                    condition_col,
                                    cluster_col,
                                    remove_celltypes = NULL) {
  
  metadata <- seu@meta.data %>%
    transmute(
      sample = as.character(.data[[sample_col]]),
      group = as.character(.data[[condition_col]]),
      cluster = as.character(.data[[cluster_col]])
    ) %>%
    mutate(
      group = na_if(group, ""),
      cluster = na_if(cluster, "")
    )
  
  keep_cells <- complete.cases(metadata)
  
  if (!is.null(remove_celltypes)) {
    keep_cells <- keep_cells &
      !tolower(metadata$cluster) %in% tolower(remove_celltypes)
  }
  
  metadata <- metadata[keep_cells, , drop = FALSE]
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = matrix(
        0,
        nrow = 1,
        ncol = nrow(metadata),
        dimnames = list("dummy_gene", rownames(metadata))
      )
    )
  )
  
  colData(sce)$sample <- factor(metadata$sample)
  colData(sce)$group <- factor(metadata$group)
  colData(sce)$cluster <- factor(metadata$cluster)
  
  message("Keeping ", nrow(metadata), " cells.")
  
  sample_summary <- as.data.frame(table(sce$group, sce$sample)) %>%
    filter(Freq > 0) %>%
    group_by(Var1) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    rename(group = Var1)
  
  print(sample_summary)
  
  valid_groups <- sample_summary %>%
    filter(n_samples >= 2) %>%
    pull(group) %>%
    as.character()
  
  if (length(valid_groups) < 2) {
    stop("Fewer than two groups have at least two biological samples.")
  }
  
  sce <- sce[, sce$group %in% valid_groups]
  sce$sample <- droplevels(sce$sample)
  sce$group <- droplevels(sce$group)
  sce$cluster <- droplevels(sce$cluster)
  
  return(sce)
}

# ------------------------------------------------------------
# Run propeller for one pairwise comparison
# ------------------------------------------------------------

run_propeller_pair <- function(sce, group_1, group_2, transform = "asin") {
  
  sce_pair <- sce[, sce$group %in% c(group_1, group_2)]
  sce_pair$sample <- droplevels(sce_pair$sample)
  sce_pair$group <- droplevels(sce_pair$group)
  sce_pair$cluster <- droplevels(sce_pair$cluster)
  
  message("Running propeller: ", group_2, " vs ", group_1)
  
  propeller_result <- speckle::propeller(
    clusters = sce_pair$cluster,
    sample = sce_pair$sample,
    group = sce_pair$group,
    transform = transform
  )
  
  propeller_result %>%
    as_tibble(rownames = "cell_type") %>%
    mutate(
      comparison = paste0(group_2, "_vs_", group_1),
      group_1 = group_1,
      group_2 = group_2
    ) %>%
    relocate(comparison, group_1, group_2, cell_type)
}

# ------------------------------------------------------------
# Run all pairwise comparisons
# ------------------------------------------------------------

run_propeller_all_pairwise <- function(seu,
                                       dataset_name,
                                       out_root,
                                       sample_col,
                                       condition_col,
                                       cluster_col,
                                       remove_celltypes = NULL,
                                       transform = "asin") {
  
  sce <- build_sce_for_propeller(
    seu = seu,
    sample_col = sample_col,
    condition_col = condition_col,
    cluster_col = cluster_col,
    remove_celltypes = remove_celltypes
  )
  
  groups <- levels(sce$group)
  group_pairs <- combn(groups, 2, simplify = FALSE)
  
  results <- map_dfr(
    group_pairs,
    ~ run_propeller_pair(
      sce = sce,
      group_1 = .x[1],
      group_2 = .x[2],
      transform = transform
    )
  ) %>%
    arrange(comparison, FDR)
  
  output_file <- file.path(
    out_root,
    paste0(dataset_name, "_propeller_all_pairwise_table.csv")
  )
  
  write_csv(results, output_file)
  message("Saved: ", output_file)
  
  return(results)
}

# ------------------------------------------------------------
# Analysis 1: Parenchyma
# Tu.free vs Emphysema vs Fibrosis
# ------------------------------------------------------------

combined_parenchyma <- subset(
  combined_list,
  subset = disease %in% c("Tu.free", "Emphysema", "Fibrosis")
)

propeller_parenchyma <- run_propeller_all_pairwise(
  seu = combined_parenchyma,
  dataset_name = "Parenchyma_Tufree_Emphysema_Fibrosis",
  out_root = out_root,
  sample_col = sample_col,
  condition_col = condition_col,
  cluster_col = cluster_col,
  remove_celltypes = remove_celltypes,
  transform = "asin"
)

# ------------------------------------------------------------
# Analysis 2: Lymph node
# Emphysema vs Fibrosis
# ------------------------------------------------------------

lymphnode_subset <- subset(
  lymphnode,
  subset = disease %in% c("Emphysema", "Fibrosis")
)

propeller_lymphnode <- run_propeller_all_pairwise(
  seu = lymphnode_subset,
  dataset_name = "Lymphnode_Emphysema_Fibrosis",
  out_root = out_root,
  sample_col = sample_col,
  condition_col = condition_col,
  cluster_col = cluster_col,
  remove_celltypes = remove_celltypes,
  transform = "asin"
)

# ------------------------------------------------------------
# Save session information
# ------------------------------------------------------------

session_file <- file.path(out_root, "sessionInfo_propeller_analysis.txt")

sink(session_file)
sessionInfo()
sink()

message("Saved session info: ", session_file)
# COPD vs IPF pseudobulk ssGSEA 

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(GSVA)
  library(dplyr)
  library(purrr)
  library(car)
})

COPD_PATH <- "data/COPD_processed.rds"
IPF_PATH <- "data/IPF_processed.rds"
SIGNATURE_PATH <- "data/signature_genes.csv"
OUT_DIR <- "results/ssGSEA"

CELLTYPES <- c(
  "Macrophages", "Monocytes_Neutrophils", "DCs",
  "NK_cells", "T_cells", "B_cells"
)

MIN_TOTAL_COUNTS <- 5
MIN_SIGNATURE_GENES <- 5
MAX_SIGNATURE_GENES <- 600
MIN_SAMPLES <- 3

# Pseudobulk aggregation
get_counts <- function(obj) {
  Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
}

make_pseudobulk <- function(obj, dataset, celltype) {
  counts <- get_counts(obj)
  meta <- obj@meta.data
  cells <- intersect(colnames(counts), rownames(meta))

  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]

  keep <- (
    meta$coarse_celltype == celltype &
      meta$disease_status %in% c("Disease", "Control") &
      !is.na(meta$sample_id)
  )

  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  if (ncol(counts) == 0L) {
    stop("No ", celltype, " cells found in ", dataset, ".")
  }

  meta$dataset_origin <- dataset
  meta$pb_id <- paste(
    dataset, meta$sample_id, meta$disease_status, celltype,
    sep = "__"
  )

  group <- factor(meta$pb_id)
  design <- Matrix::sparse.model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  pb_counts <- counts %*% design

  pb_meta <- meta %>%
    mutate(
      sample_id = as.character(sample_id),
      disease_status = as.character(disease_status),
      coarse_celltype = as.character(coarse_celltype)
    ) %>%
    count(
      pb_id, dataset_origin, sample_id,
      disease_status, coarse_celltype,
      name = "n_cells"
    ) %>%
    as.data.frame()

  pb_meta <- pb_meta[
    match(colnames(pb_counts), pb_meta$pb_id),
    ,
    drop = FALSE
  ]
  rownames(pb_meta) <- pb_meta$pb_id

  stopifnot(identical(colnames(pb_counts), rownames(pb_meta)))
  list(counts = pb_counts, meta = pb_meta)
}


# Normalization and ssGSEA
normalize_counts <- function(counts, metadata) {
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ 1
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) >= MIN_TOTAL_COUNTS, ]
  if (nrow(dds) == 0) stop("No genes remained after filtering.")

  dds <- DESeq2::estimateSizeFactors(dds)
  log2(DESeq2::counts(dds, normalized = TRUE) + 1)
}

match_genes <- function(genes, matrix_genes) {
  lookup <- setNames(matrix_genes, toupper(matrix_genes))
  matched <- unname(lookup[toupper(unique(genes))])
  unique(matched[!is.na(matched)])
}

run_ssgsea <- function(expression, gene_sets) {
  param <- GSVA::ssgseaParam(
    exprData = expression,
    geneSets = gene_sets,
    minSize = MIN_SIGNATURE_GENES,
    maxSize = MAX_SIGNATURE_GENES,
    normalize = TRUE
  )
  GSVA::gsva(param, verbose = FALSE)
}


# Statistical comparison
select_test <- function(values, groups) {
  split_values <- split(values, groups)

  normality <- vapply(
    split_values,
    function(x) {
      if (length(x) >= 3L && stats::sd(x) > 0) {
        stats::shapiro.test(x)$p.value
      } else {
        NA_real_
      }
    },
    numeric(1)
  )

  normal <- all(!is.na(normality)) && all(normality > 0.05)

  levene_p <- tryCatch(
    car::leveneTest(
      values ~ groups,
      data = data.frame(values, groups)
    )$`Pr(>F)`[1],
    error = function(e) NA_real_
  )

  if (normal) {
    test <- stats::t.test(
      values ~ groups,
      var.equal = !is.na(levene_p) && levene_p > 0.05
    )
  } else {
    test <- stats::wilcox.test(values ~ groups, exact = FALSE)
  }

  list(method = test$method, p_value = test$p.value)
}

compare_scores <- function(scores, metadata, signature_info, celltype) {
  disease_meta <- metadata %>%
    filter(disease_status == "Disease") %>%
    mutate(
      group = factor(
        dataset_origin,
        levels = c("COPD", "IPF")
      )
    )

  group_n <- table(disease_meta$group)
  if (
    !all(c("COPD", "IPF") %in% names(group_n)) ||
      any(group_n[c("COPD", "IPF")] < MIN_SAMPLES)
  ) {
    return(NULL)
  }

  purrr::map_dfr(rownames(scores), function(signature_id) {
    values <- as.numeric(
      scores[signature_id, disease_meta$pb_id]
    )

    copd <- values[disease_meta$group == "COPD"]
    ipf <- values[disease_meta$group == "IPF"]

    mean_copd <- mean(copd)
    mean_ipf <- mean(ipf)
    shift <- max(0, -min(values) + PSEUDOCOUNT)
    test <- select_test(values, disease_meta$group)

    data.frame(
      signature_id = signature_id,
      tested_coarse_celltype = celltype,
      n_COPD = length(copd),
      n_IPF = length(ipf),
      mean_score_COPD = mean_copd,
      mean_score_IPF = mean_ipf,
      score_delta_IPF_minus_COPD = mean_ipf - mean_copd,
      plot_log2fc_IPF_over_COPD = log2(
        (mean_ipf + shift + PSEUDOCOUNT) /
          (mean_copd + shift + PSEUDOCOUNT)
      ),
      higher_dataset = case_when(
        mean_ipf > mean_copd ~ "IPF",
        mean_ipf < mean_copd ~ "COPD",
        TRUE ~ "Neutral"
      ),
      test_used = test$method,
      p_value = test$p_value
    )
  }) %>%
    left_join(signature_info, by = "signature_id")
}


# Population-level analysis
run_population <- function(celltype, copd, ipf, signatures) {
  message("Running ", celltype)

  copd_pb <- make_pseudobulk(copd, "COPD", celltype)
  ipf_pb <- make_pseudobulk(ipf, "IPF", celltype)

  common_genes <- intersect(
    rownames(copd_pb$counts),
    rownames(ipf_pb$counts)
  )

  counts <- cbind(
    copd_pb$counts[common_genes, , drop = FALSE],
    ipf_pb$counts[common_genes, , drop = FALSE]
  )

  metadata <- bind_rows(copd_pb$meta, ipf_pb$meta)
  metadata <- metadata[
    match(colnames(counts), metadata$pb_id),
    ,
    drop = FALSE
  ]
  rownames(metadata) <- metadata$pb_id

  expression <- normalize_counts(counts, metadata)

  signature_info <- signatures %>%
    filter(matched_coarse_celltype == celltype) %>%
    group_by(
      signature_id,
      matched_coarse_celltype,
      direction_clean,
      signature_short_label
    ) %>%
    summarise(
      genes = list(unique(Gene)),
      n_genes_raw = n_distinct(Gene),
      .groups = "drop"
    ) %>%
    mutate(
      genes_present = map(
        genes,
        match_genes,
        matrix_genes = rownames(expression)
      ),
      n_genes_present = map_int(genes_present, length)
    ) %>%
    filter(n_genes_present >= MIN_SIGNATURE_GENES)

  if (nrow(signature_info) == 0L) {
    warning("No eligible signatures for ", celltype, ".")
    return(NULL)
  }

  gene_sets <- setNames(
    signature_info$genes_present,
    signature_info$signature_id
  )

  scores <- run_ssgsea(expression, gene_sets)

  compare_scores(
    scores,
    metadata,
    signature_info %>% select(-genes, -genes_present),
    celltype
  )
}


# Execute
copd <- readRDS(COPD_PATH)
ipf <- readRDS(IPF_PATH)
signatures <- read.csv(
  SIGNATURE_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta <- c("sample_id", "disease_status", "coarse_celltype")
required_sig <- c(
  "signature_id", "matched_coarse_celltype",
  "direction_clean", "signature_short_label", "Gene"
)

for (object_name in c("copd", "ipf")) {
  missing <- setdiff(
    required_meta,
    colnames(get(object_name)@meta.data)
  )
  if (length(missing) > 0L) {
    stop(object_name, " is missing: ", paste(missing, collapse = ", "))
  }
}

missing <- setdiff(required_sig, colnames(signatures))
if (length(missing) > 0L) {
  stop("Signature table is missing: ", paste(missing, collapse = ", "))
}

results <- purrr::map_dfr(
  CELLTYPES,
  run_population,
  copd = copd,
  ipf = ipf,
  signatures = signatures
)

if (nrow(results) == 0L) stop("No comparisons passed the filters.")

results <- results %>%
  group_by(direction_clean) %>%
  mutate(
    adj_p_global_direction = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  group_by(tested_coarse_celltype, direction_clean) %>%
  mutate(
    adj_p_population_direction = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup()

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Script: 05_BetaDiversity.R
# Description: Beta diversity analysis using iterative
#              rarefaction and PCoA ordination across
#              Intestine, Sediment, and Water compartments
#              in Litopenaeus vannamei production ponds.
#              Design: 4 treatments (Basal, CFP 5%, CFP 10%,
#              CFP 20%) x 4 ponds, analyzed per compartment.
#              Distance metrics: Aitchison, Jaccard,
#              Weighted UniFrac, Unweighted UniFrac.
# Input: phyloseq object (step4_phyloseq_object.rds)
# Note: Plot aesthetics have been simplified for public release.
#       Final publication figures used customized formatting
#       not included in this repository.
# ============================================================

# =========================================================================
# ENVIRONMENT CLEANUP AND INITIALIZATION
# =========================================================================

cleanup_environment <- function() {
  base_packages <- c("stats", "graphics", "grDevices", "utils", "datasets",
                     "methods", "base")
  loaded_packages    <- search()
  packages_to_detach <- loaded_packages[grepl("^package:", loaded_packages)]
  packages_to_detach <- packages_to_detach[
    !grepl(paste(base_packages, collapse = "|"), packages_to_detach)
  ]
  for (pkg in packages_to_detach) {
    tryCatch(
      detach(pkg, character.only = TRUE, unload = TRUE, force = TRUE),
      error = function(e) invisible(NULL)
    )
  }
  gc(verbose = FALSE)
}

cleanup_environment()

# =========================================================================
# PACKAGE INSTALLATION AND LOADING
# =========================================================================

required_packages <- c(
  "phyloseq", "vegan", "dplyr", "ggplot2", "patchwork", "ggtext",
  "stringr", "parallel", "future.apply", "tidyverse", "tidyr", "janitor",
  "ape", "readr", "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
  library(phyloseq)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
key_packages <- c("phyloseq", "vegan", "ggplot2", "patchwork", "dplyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/16S_V3V4_MicrobiomePipeline/step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Github/Exp2434/BetaDiversity"
results_dir   <- file.path(output_dir, "Statistical_Results")
plots_dir     <- file.path(output_dir, "Plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
n_iterations      <- 100
correction_method <- "fdr"
set.seed(54)
options(contrasts = c("contr.sum", "contr.poly"))

# Type-specific rarefaction depths
# Intestine, Sediment, and Water have different sequencing depth distributions;
# each compartment is rarefied to its own minimum sufficient depth.
rarefaction_depths <- list(
  Intestine = 59815,
  Sediment  = 76637,
  Water     = 85162
)

# Treatment factor order
treatment_order <- c("Basal", "CFP 5%", "CFP 10%", "CFP 20%")

# Colors — Okabe-Ito palette consistent across scripts
type_colors <- c(
  "Intestine" = "#E84646",
  "Sediment"  = "#E69F00",
  "Water"     = "#0072B2"
)
treatment_colors <- c(
  "Basal"   = "#0072B2",
  "CFP 5%"  = "#E84646",
  "CFP 10%" = "#E69F00",
  "CFP 20%" = "#009E73"
)
type_shapes <- c("Intestine" = 21, "Sediment" = 24, "Water" = 22)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

format_pvalue <- function(p) {
  ifelse(is.na(p), "= NA",
         ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
}

clean_sample_metadata <- function(metadata) {
  metadata %>%
    dplyr::filter(trial == "1") %>%
    dplyr::filter(pond %in% c("1", "2", "3", "4")) %>%
    dplyr::filter(is.na(shrimp) | as.numeric(shrimp) <= 5) %>%
    dplyr::mutate(
      type_clean = dplyr::case_when(
        stringr::str_detect(tolower(type), "gut|intestin")  ~ "Intestine",
        stringr::str_detect(tolower(type), "water")         ~ "Water",
        stringr::str_detect(tolower(type), "soil|sediment") ~ "Sediment",
        TRUE ~ as.character(type)
      ),
      trt_clean = dplyr::case_when(
        trt == "1" ~ "Basal",
        trt == "2" ~ "CFP 5%",
        trt == "3" ~ "CFP 10%",
        trt == "4" ~ "CFP 20%",
        TRUE ~ as.character(trt)
      ),
      pond_clean = as.factor(toupper(trimws(as.character(pond))))
    ) %>%
    dplyr::mutate(
      type_clean = factor(type_clean, levels = c("Intestine", "Sediment", "Water")),
      trt_clean  = factor(trt_clean,  levels = treatment_order)
    )
}

create_pvalue_annotation <- function(st, method, permanova_by_type, betadisp_by_type) {
  trt_p  <- tryCatch(
    permanova_by_type[[st]][[method]]$`Pr(>F)`[1],
    error = function(e) NA
  )
  disp_p <- tryCatch(
    betadisp_by_type[[st]][[method]]$anova$`Pr(>F)`[1],
    error = function(e) NA
  )
  stringr::str_c(
    "Diet: P ", format_pvalue(trt_p), "\n",
    "Dispersion: P ", format_pvalue(disp_p)
  )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)
cat(sprintf("Loaded: %d samples, %d ASVs\n",
            phyloseq::nsamples(physeq), phyloseq::ntaxa(physeq)))

# Round to integers — prevents rrarefy warnings from DADA2 float counts
otu_mat_raw <- phyloseq::otu_table(physeq) %>% as("matrix")
if (phyloseq::taxa_are_rows(physeq)) otu_mat_raw <- t(otu_mat_raw)
otu_mat_raw <- round(otu_mat_raw)

sample_metadata <- phyloseq::sample_data(physeq) %>%
  as("data.frame") %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  janitor::clean_names()

sample_metadata_clean <- clean_sample_metadata(sample_metadata)

phyloseq::sample_data(physeq) <- phyloseq::sample_data(
  sample_metadata_clean %>%
    dplyr::select(-sample_id) %>%
    as.data.frame() %>%
    `rownames<-`(sample_metadata_clean$sample_id)
)

# =========================================================================
# STEP 2: TYPE-SPECIFIC RAREFACTION SETUP
# =========================================================================

# Named vector: sample_id -> type label
sample_type_map <- setNames(
  as.character(sample_metadata_clean$type_clean),
  sample_metadata_clean$sample_id
)

# Retain only samples present in rounded OTU matrix after metadata filtering
meta_in_otu <- sample_metadata_clean %>%
  dplyr::filter(sample_id %in% rownames(otu_mat_raw))

# Retain samples with sufficient reads for their type-specific depth
samples_ok <- meta_in_otu %>%
  dplyr::filter(
    mapply(
      function(s, tp) {
        depth <- rarefaction_depths[[tp]]
        !is.null(depth) && rowSums(otu_mat_raw[s, , drop = FALSE]) >= depth
      },
      sample_id, as.character(type_clean)
    )
  ) %>%
  dplyr::pull(sample_id)

otu_mat_filtered     <- otu_mat_raw[samples_ok, , drop = FALSE]
sample_type_filtered <- sample_type_map[samples_ok]
meta_filtered        <- dplyr::filter(sample_metadata_clean, sample_id %in% samples_ok)

cat(sprintf("Samples passing depth filter: %d of %d\n",
            length(samples_ok), nrow(sample_metadata_clean)))
type_counts <- table(sample_type_filtered)
for (tp in names(rarefaction_depths)) {
  n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
  cat(sprintf("  %s: %d samples (depth >= %s)\n",
              tp, n, scales::comma(rarefaction_depths[[tp]])))
}

tree <- phyloseq::phy_tree(physeq, errorIfNULL = FALSE)
if (is.null(tree)) {
  cat("WARNING: No phylogenetic tree found — UniFrac will not be calculated\n")
} else {
  cat("Phylogenetic tree found — UniFrac will be calculated\n")
}

# =========================================================================
# STEP 3: ITERATIVE RAREFACTION
# =========================================================================

run_rarefaction_iteration <- function(i, otu_matrix, sample_type_map,
                                      rarefaction_depths, tree_obj, physeq_obj,
                                      meta_data) {
  if (!requireNamespace("vegan",    quietly = TRUE)) stop("vegan required")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("phyloseq required")
  if (!requireNamespace("ape",      quietly = TRUE)) stop("ape required")
  if (!requireNamespace("dplyr",    quietly = TRUE)) stop("dplyr required")

  set.seed(i * 100)

  # Type-specific rarefaction: each sample is rarefied to its compartment depth
  rarefied           <- matrix(0L, nrow = nrow(otu_matrix), ncol = ncol(otu_matrix))
  rownames(rarefied) <- rownames(otu_matrix)
  colnames(rarefied) <- colnames(otu_matrix)
  for (idx in seq_len(nrow(otu_matrix))) {
    sname         <- rownames(otu_matrix)[idx]
    depth         <- rarefaction_depths[[ sample_type_map[sname] ]]
    rarefied[idx, ] <- vegan::rrarefy(otu_matrix[idx, , drop = FALSE], depth)
  }

  # Aitchison: CLR + Euclidean
  clr_mat <- t(apply(rarefied + 1, 1, function(x) log(x / exp(mean(log(x))))))

  distances <- list(
    jaccard   = vegan::vegdist(rarefied > 0, method = "jaccard"),
    aitchison = vegan::vegdist(clr_mat,      method = "euclidean")
  )

  if (!is.null(tree_obj)) {
    tryCatch({
      present_taxa      <- colnames(rarefied)[colSums(rarefied) > 0]
      tree_tips_to_keep <- intersect(present_taxa, tree_obj$tip.label)
      if (length(tree_tips_to_keep) > 2) {
        tree_pruned <- ape::keep.tip(tree_obj, tree_tips_to_keep)
        temp_phyloseq <- phyloseq::phyloseq(
          phyloseq::otu_table(rarefied[, tree_tips_to_keep], taxa_are_rows = FALSE),
          phyloseq::tax_table(physeq_obj)[tree_tips_to_keep, ],
          phyloseq::sample_data(
            meta_data %>%
              dplyr::select(-sample_id) %>%
              as.data.frame() %>%
              `rownames<-`(meta_data$sample_id)
          ),
          tree_pruned
        )
        distances$wunifrac <- phyloseq::distance(temp_phyloseq, method = "wunifrac")
        distances$uunifrac <- phyloseq::distance(temp_phyloseq, method = "unifrac")
      }
    }, error = function(e) invisible(NULL))
  }

  list(iteration = i, distances = distances, sample_names = rownames(rarefied))
}

n_cores <- max(1L, parallel::detectCores() - 1L)
future::plan(future::multisession, workers = n_cores)
distance_iterations <- future.apply::future_lapply(
  seq_len(n_iterations),
  function(i) tryCatch(
    run_rarefaction_iteration(
      i, otu_mat_filtered, sample_type_filtered,
      rarefaction_depths, tree, physeq, meta_filtered
    ),
    error = function(e) NULL
  ),
  future.seed = TRUE
)
future::plan(future::sequential)

valid_n <- sum(!sapply(distance_iterations, is.null))
cat(sprintf("Completed %d of %d rarefaction iterations\n", valid_n, n_iterations))

# =========================================================================
# STEP 4: AVERAGE DISTANCE MATRICES
# =========================================================================

distance_methods <- c("jaccard", "aitchison", "wunifrac", "uunifrac")
sample_names     <- rownames(otu_mat_filtered)
n_samples        <- length(sample_names)

averaged_distances <- list()
for (method in distance_methods) {
  sum_matrix       <- matrix(0, nrow = n_samples, ncol = n_samples,
                             dimnames = list(sample_names, sample_names))
  valid_iterations <- 0L
  for (iter_result in distance_iterations) {
    if (is.null(iter_result) || !method %in% names(iter_result$distances)) next
    dist_matrix <- as.matrix(iter_result$distances[[method]])
    if (!all(sample_names %in% rownames(dist_matrix))) next
    sum_matrix       <- sum_matrix + dist_matrix[sample_names, sample_names]
    valid_iterations <- valid_iterations + 1L
  }
  if (valid_iterations > 0) {
    averaged_distances[[method]] <- as.dist(sum_matrix / valid_iterations)
    cat(sprintf("  %s: averaged across %d iterations\n", method, valid_iterations))
  }
}
averaged_distances <- averaged_distances[!sapply(averaged_distances, is.null)]

# =========================================================================
# STEP 5: AVERAGE INTESTINE REPLICATES TO POND LEVEL
# =========================================================================
# Intestine has up to 5 shrimp replicates per pond x treatment combination.
# Pond is the true experimental unit; averaging collapses pseudoreplicates
# before ordination and statistical testing.
# =========================================================================

average_intestine <- function(physeq_obj) {
  otu_m   <- as(phyloseq::otu_table(physeq_obj), "matrix")
  if (phyloseq::taxa_are_rows(physeq_obj)) otu_m <- t(otu_m)
  meta_df <- as(phyloseq::sample_data(physeq_obj), "data.frame")

  int_idx  <- meta_df$type_clean == "Intestine"
  oth_idx  <- !int_idx
  int_meta <- meta_df[int_idx, ]
  int_meta$group_id <- paste(int_meta$pond_clean, int_meta$trt_clean, sep = "_")

  avg_otu  <- list()
  avg_meta <- list()
  for (grp in unique(int_meta$group_id)) {
    g <- int_meta$group_id == grp
    new_id              <- paste0("Avg_Int_", grp)
    avg_otu[[new_id]]   <- colMeans(otu_m[int_idx, ][g, , drop = FALSE])
    avg_meta[[new_id]]  <- int_meta[g, ][1, ]
  }

  int_meta_df <- do.call(rbind, avg_meta)
  oth_meta_df <- meta_df[oth_idx, ]
  all_cols    <- union(names(int_meta_df), names(oth_meta_df))
  for (col in all_cols) {
    if (!col %in% names(int_meta_df)) int_meta_df[[col]] <- NA
    if (!col %in% names(oth_meta_df)) oth_meta_df[[col]] <- NA
  }

  final_otu            <- rbind(do.call(rbind, avg_otu), otu_m[oth_idx, , drop = FALSE])
  final_meta           <- rbind(int_meta_df[, all_cols], oth_meta_df[, all_cols])
  rownames(final_meta) <- rownames(final_otu)

  phyloseq::phyloseq(
    phyloseq::otu_table(final_otu,    taxa_are_rows = FALSE),
    phyloseq::sample_data(final_meta),
    phyloseq::tax_table(physeq_obj),
    phyloseq::phy_tree(physeq_obj)
  )
}

temp_physeq <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_mat_filtered, taxa_are_rows = FALSE),
  phyloseq::sample_data(
    meta_filtered %>%
      dplyr::select(-sample_id) %>%
      as.data.frame() %>%
      `rownames<-`(meta_filtered$sample_id)
  ),
  phyloseq::tax_table(physeq)[colnames(otu_mat_filtered), ],
  phyloseq::phy_tree(physeq)
)

physeq_avg <- average_intestine(temp_physeq)
cat(sprintf("Final samples after pond-level averaging: %d\n",
            phyloseq::nsamples(physeq_avg)))

# =========================================================================
# STEP 6: FINAL DISTANCE MATRICES ON AVERAGED DATA
# =========================================================================

otu_avg <- as(phyloseq::otu_table(physeq_avg), "matrix")
if (phyloseq::taxa_are_rows(physeq_avg)) otu_avg <- t(otu_avg)

clr_avg <- t(apply(otu_avg + 1, 1, function(x) log(x / exp(mean(log(x))))))

final_distances <- list(
  jaccard   = vegan::vegdist(otu_avg > 0, method = "jaccard"),
  aitchison = vegan::vegdist(clr_avg,      method = "euclidean")
)

if (!is.null(phyloseq::phy_tree(physeq_avg, errorIfNULL = FALSE))) {
  tryCatch({
    final_distances$wunifrac <- phyloseq::distance(physeq_avg, method = "wunifrac")
    cat("Weighted UniFrac calculated\n")
  }, error = function(e) cat(sprintf("WARNING: wUniFrac failed: %s\n", e$message)))
  tryCatch({
    final_distances$uunifrac <- phyloseq::distance(physeq_avg, method = "unifrac")
    cat("Unweighted UniFrac calculated\n")
  }, error = function(e) cat(sprintf("WARNING: uUniFrac failed: %s\n", e$message)))
}
final_distances <- Filter(Negate(is.null), final_distances)

meta_plot_final <- as(phyloseq::sample_data(physeq_avg), "data.frame") %>%
  dplyr::mutate(
    type_clean = factor(type_clean, levels = c("Intestine", "Sediment", "Water")),
    trt_clean  = factor(trt_clean,  levels = treatment_order)
  )

cat(sprintf("Final distance metrics available: %s\n",
            paste(names(final_distances), collapse = ", ")))

# =========================================================================
# STEP 7: PCoA ORDINATION
# =========================================================================

pcoa_results  <- list()
pcoa_variance <- list()
for (method in names(final_distances)) {
  pcoa_results[[method]]  <- cmdscale(final_distances[[method]], k = 2, eig = TRUE)
  ev                       <- pcoa_results[[method]]$eig
  pcoa_variance[[method]] <- ev / sum(abs(ev))
  cat(sprintf("  %s: PC1=%.1f%%, PC2=%.1f%%\n",
              method,
              pcoa_variance[[method]][1] * 100,
              pcoa_variance[[method]][2] * 100))
}

# =========================================================================
# STEP 8: STATISTICAL ANALYSIS — PERMANOVA + BETADISPER PER COMPARTMENT
#
# adonis2 WITHOUT strata — design rationale:
#   In temporal designs (e.g., Exp2433), strata=pond_clean is valid because
#   each pond has multiple timepoints, enabling within-stratum permutations.
#   In this study, after averaging intestine replicates to pond level, each pond
#   contributes exactly 1 sample per treatment per compartment. Using
#   strata=pond_clean would restrict permutations to within-pond subsets of
#   size 1, making all permutations identical and forcing p = 1.000 regardless
#   of actual treatment signal. Since pseudoreplication is already resolved by
#   pond-level averaging, free permutation across all samples within each
#   compartment (16 averaged samples: 4 ponds x 4 treatments) is the
#   statistically correct approach.
# =========================================================================

permanova_by_type <- list()
betadisp_by_type  <- list()

for (st in levels(meta_plot_final$type_clean)) {
  permanova_by_type[[st]] <- list()
  betadisp_by_type[[st]]  <- list()

  st_samples <- rownames(meta_plot_final)[meta_plot_final$type_clean == st]
  cat(sprintf("  %s: %d samples\n", st, length(st_samples)))

  for (method in names(final_distances)) {
    dist_mat_full <- as.matrix(final_distances[[method]])
    samps_in_dist <- intersect(st_samples, rownames(dist_mat_full))

    if (length(samps_in_dist) < 4) {
      cat(sprintf("    WARNING: %s/%s — insufficient samples, skipping\n", st, method))
      next
    }

    sub_dist <- as.dist(dist_mat_full[samps_in_dist, samps_in_dist])
    sub_meta <- meta_plot_final[samps_in_dist, ]

    if (length(unique(droplevels(sub_meta$trt_clean))) < 2) {
      cat(sprintf("    WARNING: %s/%s — only 1 treatment level, skipping\n", st, method))
      next
    }

    # PERMANOVA: free permutation — strata intentionally omitted (see note above)
    permanova_by_type[[st]][[method]] <- tryCatch(
      vegan::adonis2(sub_dist ~ trt_clean, data = sub_meta,
                     permutations = 9999, by = "terms"),
      error = function(e) {
        cat(sprintf("    WARNING: PERMANOVA failed %s/%s: %s\n", st, method, e$message))
        NULL
      }
    )

    if (!is.null(permanova_by_type[[st]][[method]])) {
      p  <- permanova_by_type[[st]][[method]]$`Pr(>F)`[1]
      r2 <- round(permanova_by_type[[st]][[method]]$R2[1] * 100, 1)
      cat(sprintf("    %s/%s: R2=%.1f%%, P=%.3f\n", st, method, r2, p))
    }

    # BETADISPER: homogeneity of dispersion by treatment
    betadisp_by_type[[st]][[method]] <- tryCatch({
      bd  <- vegan::betadisper(sub_dist, sub_meta$trt_clean)
      anv <- anova(bd)
      list(betadisper = bd, anova = anv)
    }, error = function(e) {
      cat(sprintf("    WARNING: BETADISPER failed %s/%s: %s\n", st, method, e$message))
      NULL
    })
  }
}

# =========================================================================
# STEP 9: SAVE STATISTICAL RESULTS
# =========================================================================

perm_rows <- tibble::tibble()
for (st in names(permanova_by_type)) {
  for (method in names(permanova_by_type[[st]])) {
    res <- permanova_by_type[[st]][[method]]
    if (is.null(res)) next
    perm_rows <- dplyr::bind_rows(perm_rows, tibble::tibble(
      Sample_Type = st,
      Method      = method,
      R2_pct      = round(res$R2[1] * 100, 1),
      p_value     = res$`Pr(>F)`[1]
    ))
  }
}

beta_rows <- tibble::tibble()
for (st in names(betadisp_by_type)) {
  for (method in names(betadisp_by_type[[st]])) {
    res <- betadisp_by_type[[st]][[method]]
    if (is.null(res)) next
    beta_rows <- dplyr::bind_rows(beta_rows, tibble::tibble(
      Sample_Type    = st,
      Method         = method,
      F_value        = round(res$anova$`F value`[1], 3),
      p_value        = res$anova$`Pr(>F)`[1],
      Interpretation = dplyr::if_else(
        res$anova$`Pr(>F)`[1] < 0.05, "Unequal dispersion", "Equal dispersion"
      )
    ))
  }
}

readr::write_csv(perm_rows, file.path(results_dir, "PERMANOVA_Results_by_Type.csv"))
readr::write_csv(beta_rows, file.path(results_dir, "BetaDisper_Results_by_Type.csv"))

readr::write_csv(
  tibble::tibble(
    Method                 = names(final_distances),
    PCoA_PC1_PC2_Variance  = sapply(pcoa_variance,
                                     function(x) round(sum(x[1:2]) * 100, 1)),
    Rarefaction_Iterations = n_iterations,
    FDR_Correction         = correction_method
  ),
  file.path(results_dir, "Ordination_Quality_Summary.csv")
)

readr::write_csv(
  meta_plot_final %>% tibble::rownames_to_column("sample_id"),
  file.path(results_dir, "Final_Metadata.csv")
)

cat("Statistical results saved to:", results_dir, "\n")

# =========================================================================
# STEP 10: CREATE PCoA PLOTS
# =========================================================================

method_titles <- list(
  aitchison = "Aitchison compositional dissimilarity",
  jaccard   = "Jaccard presence-absence dissimilarity",
  wunifrac  = "Weighted UniFrac phylogenetic dissimilarity",
  uunifrac  = "Unweighted UniFrac phylogenetic dissimilarity"
)

# Panel order within each compartment figure: A=Aitchison, B=Jaccard,
# C=Weighted UniFrac, D=Unweighted UniFrac
panel_order <- c("aitchison", "jaccard", "wunifrac", "uunifrac")

# Build per-compartment, per-metric PCoA panels
# PCoA is recomputed on the compartment-specific distance sub-matrix so that
# axis variance explained reflects only the samples in that compartment.
pcoa_panels <- list()

for (st in levels(meta_plot_final$type_clean)) {
  pcoa_panels[[st]] <- list()
  st_samps          <- rownames(meta_plot_final)[meta_plot_final$type_clean == st]
  current_shape     <- type_shapes[st]

  for (method in names(final_distances)) {
    dist_mat_full <- as.matrix(final_distances[[method]])
    samps_avail   <- intersect(st_samps, rownames(dist_mat_full))
    if (length(samps_avail) < 4) next

    sub_dist <- as.dist(dist_mat_full[samps_avail, samps_avail])
    pcoa_sub <- cmdscale(sub_dist, k = 2, eig = TRUE)
    ev       <- pcoa_sub$eig
    var_exp  <- ev / sum(abs(ev))
    pc1_var  <- round(var_exp[1] * 100, 1)
    pc2_var  <- round(var_exp[2] * 100, 1)

    pcoa_df <- as.data.frame(pcoa_sub$points)
    colnames(pcoa_df) <- c("PC1", "PC2")
    pcoa_df <- cbind(pcoa_df, meta_plot_final[samps_avail, ])

    annotation <- create_pvalue_annotation(st, method, permanova_by_type, betadisp_by_type)
    title_str  <- paste0(st, " — ", method_titles[[method]])

    pcoa_panels[[st]][[method]] <- ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = trt_clean, fill = trt_clean),
                 shape = current_shape,
                 size = 4, alpha = 0.6, stroke = 1) +
      stat_ellipse(aes(fill = trt_clean), level = 0.95, type = "norm",
                   alpha = 0.15, linewidth = 0.6, geom = "polygon") +
      annotate("text", x = Inf, y = -Inf, label = annotation,
               hjust = 1.05, vjust = -0.3, size = 3,
               fontface = "italic", lineheight = 0.9) +
      scale_color_manual(values = treatment_colors, name = "Diet") +
      scale_fill_manual(values  = treatment_colors, name = "Diet") +
      labs(
        title = title_str,
        x     = paste0("PC1 (", pc1_var, "%)"),
        y     = paste0("PC2 (", pc2_var, "%)")
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0, size = 10),
        axis.title       = element_text(face = "bold"),
        axis.text        = element_text(color = "black"),
        legend.position  = "right",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }
  cat(sprintf("  %s: %d panels built\n", st, length(pcoa_panels[[st]])))
}

# =========================================================================
# STEP 11: THREE 2x2 COMBINED FIGURES — one per compartment
# =========================================================================

for (st in levels(meta_plot_final$type_clean)) {
  avail_panels <- panel_order[panel_order %in% names(pcoa_panels[[st]])]
  if (length(avail_panels) < 2) {
    cat(sprintf("  WARNING: Skipping %s — fewer than 2 panels available\n", st))
    next
  }

  depth_label <- scales::comma(rarefaction_depths[[st]])

  combined <- patchwork::wrap_plots(
    pcoa_panels[[st]][avail_panels], ncol = 2
  ) +
    patchwork::plot_annotation(
      title = paste0(
        "Compositional dissimilarity of ", tolower(st),
        " microbial communities across Litopenaeus vannamei",
        " fed various levels of corn fermented protein",
        " in semi-intensive, green water production ponds"
      ),
      subtitle = paste0(
        "Based on 16S rRNA amplicon sequencing data rarefied to ",
        depth_label, " reads across ", n_iterations, " iterations. ",
        "Differences in community composition among diet groups were evaluated ",
        "by permutational multivariate analysis of variance (PERMANOVA) at ",
        "9,999 permutations and homogeneity of multivariate dispersion (BETADISPER). ",
        "Panels (A-D): (A) Aitchison compositional dissimilarity, ",
        "(B) Jaccard presence-absence, ",
        "(C) Weighted UniFrac phylogenetic dissimilarity, ",
        "(D) Unweighted UniFrac phylogenetic dissimilarity."
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9,  hjust = 0),
      plot.tag      = element_text(size = 11, face = "bold")
    )

  fname <- paste0("BetaDiversity_", st, "_Combined.png")
  ggsave(
    filename = file.path(plots_dir, fname),
    plot     = combined,
    width    = 16, height = 12, units = "in", dpi = 300
  )
  cat(sprintf("  Saved: %s\n", fname))
}

# =========================================================================
# ANALYSIS SUMMARY
# =========================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("OUTPUTS:\n")
cat(strrep("-", 70), "\n")
cat("TABLES:\n")
cat("  PERMANOVA_Results_by_Type.csv\n")
cat("  BetaDisper_Results_by_Type.csv\n")
cat("  Ordination_Quality_Summary.csv\n")
cat("  Final_Metadata.csv\n")
cat("\nFIGURES:\n")
cat("  BetaDiversity_Intestine_Combined.png\n")
cat("  BetaDiversity_Sediment_Combined.png\n")
cat("  BetaDiversity_Water_Combined.png\n")
cat(sprintf("\nOUTPUT LOCATION:\n  %s\n\n", output_dir))

cat("PERMANOVA summary (by compartment):\n")
print(perm_rows)
cat("\nBETADISPER summary (by compartment):\n")
print(beta_rows)

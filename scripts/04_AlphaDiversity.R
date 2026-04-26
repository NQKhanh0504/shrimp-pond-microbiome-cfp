# ============================================================
# Script: 04_AlphaDiversity.R
# Description: Alpha diversity analysis for CFP
#              treatment comparison experiment.
#              Produces: Observed Richness, Shannon Entropy,
#              Berger-Parker Dominance, and Faith's PD raincloud
#              plots with lm(metric ~ treatment) statistics and
#              Tukey CLD annotations, analyzed separately per
#              sample type (Intestine, Sediment, Water).
# Input: phyloseq object (step4_phyloseq_object.rds)
# Design: 4 treatments (Basal, CFP 5%, CFP 10%, CFP 20%)
#         analyzed by sample compartment using OLS + emmeans.
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

cran_packages <- c(
  "tidyverse", "vegan", "emmeans", "janitor", "picante", "ape",
  "parallel", "future.apply", "multcomp", "multcompView",
  "flextable", "ggpubr", "patchwork", "svglite", "ggtext",
  "ggrepel", "gghalves", "scales", "stringr", "broom"
)

bioc_packages <- c("phyloseq")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
key_packages <- c("phyloseq", "vegan", "emmeans", "multcomp",
                  "flextable", "ggplot2", "dplyr", "tidyr", "picante")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/16S_V3V4_MicrobiomePipeline/step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/AlphaDiversity/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
random_seed  <- 54
n_iterations <- 100
alpha_level  <- 0.05

# Type-specific rarefaction depths — set to actual minimum read counts per compartment
# so no samples are excluded within each type. Adjust if sequencing depths change.
rarefaction_depths <- list(
  Intestine = 59815,
  Sediment  = 76637,
  Water     = 85162
)

# Parallel processing configuration
# Options: "future" (cross-platform, RStudio-safe), "mclapply" (Unix/Mac), "parLapply" (Windows)
parallel_method <- "future"

# Sample filtering
filter_ponds     <- c("C3", "C6", "D1", "D4")
filter_shrimp_max <- 5

# Diversity metrics — all_metrics computed each iteration; primary_metrics plotted
all_metrics <- c(
  "observed_richness", "shannon", "chao1", "simpson", "pielou",
  "ace", "goods_coverage", "faith_pd", "berger_parker", "robbins"
)
primary_metrics <- c("observed_richness", "shannon", "berger_parker", "faith_pd")

display_names <- c(
  "observed_richness" = "Observed Richness",
  "shannon"           = "Shannon Entropy",
  "chao1"             = "Chao1",
  "simpson"           = "Simpson Index",
  "pielou"            = "Pielou's Evenness",
  "ace"               = "ACE",
  "goods_coverage"    = "Good's Coverage",
  "faith_pd"          = "Faith's Phylogenetic Diversity",
  "berger_parker"     = "Berger-Parker Dominance",
  "robbins"           = "Robbins Index"
)

# Visualization parameters — Okabe-Ito palette keyed to treatment label values
color_palette <- c(
  "Basal"   = "#0072B2",
  "CFP 5%"  = "#E69F00",
  "CFP 10%" = "#009E73",
  "CFP 20%" = "#E84646"
)

# Treatment shapes — for scale_shape_manual consistency
type_shapes <- c(
  "Basal"   = 21,
  "CFP 5%"  = 22,
  "CFP 10%" = 23,
  "CFP 20%" = 24
)

# Jitter point shapes keyed to sample type:
#   Intestine = circle (21), Sediment = triangle (24), Water = square (22)
# Diamond (23) is reserved for the mean summary point in stat_summary
sample_type_shapes <- c(
  "Intestine" = 21,
  "Sediment"  = 24,
  "Water"     = 22
)

# Raincloud letter offset as fraction of y-range
letter_offset <- 0.15

set.seed(random_seed)
options(contrasts = c("contr.sum", "contr.poly"))

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

format_pvalue <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
}

# Standardize metadata: rename type labels ("Gut" -> "Intestine", "Soil" -> "Sediment"),
# recode numeric trt codes to descriptive treatment labels used in the palette and plots.
process_metadata <- function(metadata_raw) {
  metadata_raw %>%
    dplyr::mutate(
      type = dplyr::case_when(
        type == "Gut"   ~ "Intestine",
        type == "Soil"  ~ "Sediment",
        type == "Water" ~ "Water",
        TRUE ~ as.character(type)
      ),
      type = factor(type, levels = c("Intestine", "Sediment", "Water")),
      # Recode numeric trt codes to descriptive labels used in palette/plots
      trt = dplyr::case_when(
        trt == "1" ~ "Basal",
        trt == "2" ~ "CFP 5%",
        trt == "3" ~ "CFP 10%",
        trt == "4" ~ "CFP 20%",
        TRUE ~ as.character(trt)
      ),
      trt  = factor(trt, levels = c("Basal", "CFP 5%", "CFP 10%", "CFP 20%")),
      pond = as.factor(pond)
    )
}

filter_samples <- function(metadata_clean, ponds = NULL, shrimp_max = NULL) {
  filtered <- metadata_clean
  if (!is.null(ponds))      filtered <- dplyr::filter(filtered, pond %in% ponds)
  if (!is.null(shrimp_max)) filtered <- dplyr::filter(filtered, is.na(shrimp) | shrimp <= shrimp_max)
  return(filtered)
}

# PSE = pooled standard error = RMSE / sqrt(mean replicate count per treatment).
# Equivalent to SAS LSMEANS stderr; summarizes the precision of treatment mean estimates.
calc_pse <- function(mod, data, metric) {
  e      <- residuals(mod)
  df_e   <- mod$df.residual
  RMSE   <- sqrt(sum(e^2) / df_e)
  n_reps <- mean(table(data$trt))
  RMSE / sqrt(n_reps)
}

# =========================================================================
# ALPHA DIVERSITY CALCULATION (one rarefaction iteration)
# Each sample is rarefied to its compartment-specific depth so that
# depth differences between Intestine, Sediment, and Water do not confound
# richness comparisons within a compartment.
# =========================================================================

calculate_alpha_diversity_one_iter <- function(i, otu_matrix, sample_type_map,
                                               rarefaction_depths, tree) {
  set.seed(i * 100)

  rarefied           <- matrix(0L, nrow = nrow(otu_matrix), ncol = ncol(otu_matrix))
  rownames(rarefied) <- rownames(otu_matrix)
  colnames(rarefied) <- colnames(otu_matrix)

  for (idx in seq_len(nrow(otu_matrix))) {
    sname <- rownames(otu_matrix)[idx]
    depth <- rarefaction_depths[[ sample_type_map[sname] ]]
    rarefied[idx, ] <- vegan::rrarefy(otu_matrix[idx, , drop = FALSE], depth)
  }

  base_metrics <- data.frame(
    sample_id         = rownames(rarefied),
    iteration         = i,
    observed_richness = vegan::specnumber(rarefied),
    shannon           = vegan::diversity(rarefied, index = "shannon"),
    simpson           = vegan::diversity(rarefied, index = "simpson"),
    pielou            = vegan::diversity(rarefied, index = "shannon") /
      log(vegan::specnumber(rarefied)),
    goods_coverage    = sapply(seq_len(nrow(rarefied)), function(x)
      1 - sum(rarefied[x, ] == 1) / sum(rarefied[x, ]))
  )

  richness_estimates <- do.call(rbind, lapply(seq_len(nrow(rarefied)), function(x) {
    est <- vegan::estimateR(rarefied[x, ])
    data.frame(sample_id = rownames(rarefied)[x],
               chao1 = est["S.chao1"], ace = est["S.ACE"])
  }))

  pd_result <- picante::pd(rarefied, tree, include.root = FALSE)
  faith_pd  <- data.frame(sample_id = rownames(pd_result), faith_pd = pd_result$PD)

  additional_metrics <- do.call(rbind, lapply(seq_len(nrow(rarefied)), function(x) {
    row_data    <- rarefied[x, ]
    total_reads <- sum(row_data)
    richness    <- sum(row_data > 0)
    singletons  <- sum(row_data == 1)
    data.frame(
      sample_id     = rownames(rarefied)[x],
      berger_parker = ifelse(total_reads > 0, max(row_data) / total_reads, NA_real_),
      robbins       = ifelse(richness > 0, singletons / (richness + 1), NA_real_)
    )
  }))

  result <- merge(base_metrics, richness_estimates, by = "sample_id")
  result <- merge(result,       faith_pd,           by = "sample_id")
  result <- merge(result,       additional_metrics, by = "sample_id")
  return(result)
}

# =========================================================================
# PLOT FUNCTION
# =========================================================================

# Raincloud plot: one metric x one sample type, x-axis = treatment.
# The global aes() carries NO fill so geom_boxplot after_scale() does not
# override fill in geom_jitter. Each layer maps fill/color explicitly.
create_treatment_box_plot <- function(data, metric, cld_data, treatment_p,
                                      sample_type, treatment_pse = NA,
                                      letter_offset = 0.15) {

  fmt_p <- function(p) {
    if (is.na(p))  return("P-value = NA")
    if (p < 0.001) return("P-value < 0.001")
    paste0("P-value = ", sprintf("%.3f", p))
  }

  show_cld     <- !is.na(treatment_p) && treatment_p < alpha_level
  global_max   <- max(data[[metric]], na.rm = TRUE)
  global_min   <- min(data[[metric]], na.rm = TRUE)
  global_range <- global_max - global_min
  offset_val   <- global_range * letter_offset

  # n = number of ponds per treatment (the true experimental unit after pond-level averaging)
  y_max_group <- data %>%
    dplyr::group_by(trt) %>%
    dplyr::summarise(
      max_val = max(.data[[metric]], na.rm = TRUE),
      min_val = min(.data[[metric]], na.rm = TRUE),
      n       = dplyr::n_distinct(pond),
      .groups = "drop"
    )

  # Resolve treatment palette; fall back gracefully for unregistered levels
  trts    <- levels(data$trt)
  pal     <- color_palette[trts]
  missing <- is.na(pal)
  if (any(missing)) pal[missing] <- scales::hue_pal()(sum(missing))

  # Jitter point shape determined by sample type (circle / triangle / square)
  jitter_shape <- sample_type_shapes[[sample_type]]
  if (is.null(jitter_shape) || is.na(jitter_shape)) jitter_shape <- 21L

  p_label <- paste0(
    fmt_p(treatment_p), "\n",
    "PSE = ", ifelse(is.na(treatment_pse), "NA", sprintf("%.3f", treatment_pse))
  )

  p <- ggplot(data, aes(x = trt, y = .data[[metric]])) +
    gghalves::geom_half_violin(
      aes(fill = trt),
      side = "r", trim = FALSE, alpha = 0.5, color = NA
    ) +
    geom_boxplot(
      aes(color = trt, fill = ggplot2::after_scale(alpha(color, 0.2))),
      width = 0.4, outlier.shape = NA, lwd = 3
    ) +
    stat_summary(
      aes(color = trt),
      fun = mean, geom = "point",
      shape = 23, size = 13, fill = "white", stroke = 4
    ) +
    # fill = trt + color = trt: solid colored filled points with matching stroke.
    # Shape is fixed per sample type (circle/triangle/square); seed for reproducibility.
    geom_jitter(
      aes(color = trt, fill = trt),
      shape    = jitter_shape,
      position = position_jitter(width = 0.25, seed = random_seed),
      size = 13, alpha = 0.6, stroke = 4
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label   = p_label,
      hjust = 1, vjust = 1.1,
      size = 4, fontface = "italic", lineheight = 1.2
    ) +
    scale_fill_manual(values  = pal, name = "Treatment") +
    scale_color_manual(values = pal, name = "Treatment") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(y = display_names[metric], x = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      axis.ticks.x       = element_line(color = "black"),
      axis.text.x        = element_text(face = "bold", color = "black"),
      axis.text.y        = element_text(color = "black"),
      axis.title         = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray92"),
      panel.background   = element_rect(fill = "transparent", color = NA),
      plot.background    = element_rect(fill = "transparent", color = NA),
      legend.position    = "none"
    )

  if (show_cld) {
    cld_plot_data <- cld_data %>%
      dplyr::mutate(Letter = stringr::str_trim(.group)) %>%
      dplyr::select(trt, Letter) %>%
      dplyr::left_join(y_max_group, by = "trt") %>%
      dplyr::mutate(y_pos = max_val + offset_val)

    p <- p + geom_text(
      data = cld_plot_data,
      aes(x = trt, y = y_pos, label = Letter),
      size = 4, fontface = "bold", color = "black"
    )
  }

  n_label_data <- y_max_group %>%
    dplyr::mutate(
      label = paste0("n = ", n),
      y_pos = min_val - global_range * (letter_offset * 0.6)
    )

  p <- p + geom_text(
    data        = n_label_data,
    aes(x = trt, y = y_pos, label = label, color = trt),
    size        = 3,
    fontface    = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  )

  return(p)
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

cat("\n")
cat("STARTING ALPHA DIVERSITY ANALYSIS\n")
cat("\n")
cat("STEP 1: Loading and preparing data\n")

physeq <- readRDS(phyloseq_path)
cat(paste0(
  "Loaded: ", phyloseq::nsamples(physeq), " samples, ",
  phyloseq::ntaxa(physeq), " ASVs"
), "\n")

sample_metadata <- phyloseq::sample_data(physeq) %>%
  tibble::as_tibble(rownames = "sample") %>%
  janitor::clean_names() %>%
  dplyr::filter(!is.na(trial) & trial == "1") %>%
  process_metadata() %>%
  filter_samples(ponds = filter_ponds, shrimp_max = filter_shrimp_max)

cat(paste0(
  "Metadata: ", nrow(sample_metadata), " samples | ",
  "Treatments: ", paste(levels(sample_metadata$trt), collapse = ", "), " | ",
  "Types: ", paste(levels(sample_metadata$type), collapse = ", ")
), "\n")

# Filter to samples with sufficient reads for their type-specific rarefaction depth
otu_mat <- phyloseq::otu_table(physeq) %>% as("matrix")
if (phyloseq::taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
# Round to integers: rrarefy requires count data; non-integer values (e.g. from
# DADA2 sequence table merging) trigger warnings. Rounding preserves library
# size to within 1 read per ASV and does not affect rarefaction validity.
otu_mat <- round(otu_mat)

# Named character vector: sample ID -> type label (no factor encoding needed here)
sample_type_map <- setNames(
  as.character(sample_metadata$type),
  sample_metadata$sample
)

samples_sufficient <- sample_metadata %>%
  dplyr::filter(
    sample %in% rownames(otu_mat),
    mapply(function(s, tp) rowSums(otu_mat[s, , drop = FALSE]) >= rarefaction_depths[[tp]],
           sample, as.character(type))
  ) %>%
  dplyr::pull(sample)

otu_mat_filtered     <- otu_mat[samples_sufficient, ]
# Subset map AFTER building it from the full metadata to preserve names correctly
sample_type_filtered <- sample_type_map[samples_sufficient]
tree                 <- phyloseq::phy_tree(physeq)

cat(paste0(
  "Samples passing depth filter: ",
  length(samples_sufficient), " of ", nrow(sample_metadata)
), "\n")
# Use table() for reliable per-type counts — avoids NA from factor/character mismatch
type_counts <- table(sample_type_filtered)
for (tp in names(rarefaction_depths)) {
  n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
  cat(paste0("  ", tp, ": ", n,
             " samples (depth >= ", rarefaction_depths[[tp]], ")"), "\n")
}

# =========================================================================
# STEP 2: CALCULATE ALPHA DIVERSITY
# =========================================================================

cat(paste0("STEP 2: Iterative rarefaction (", n_iterations, " iterations)\n"))

# Validate with single test iteration before committing to parallel run
test_result <- tryCatch({
  calculate_alpha_diversity_one_iter(1, otu_mat_filtered, sample_type_filtered,
                                     rarefaction_depths, tree)
}, error = function(e) {
  stop(paste0("Test iteration failed: ", e$message,
              "\nParallel processing cannot proceed."))
})
cat(paste0("Test iteration passed — ", nrow(test_result), " samples calculated\n"))

n_cores <- parallel::detectCores() - 1

if (parallel_method == "mclapply") {
  alpha_div_iterations <- parallel::mclapply(seq_len(n_iterations), function(i) {
    tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, sample_type_filtered,
                                         rarefaction_depths, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i))
  }, mc.cores = n_cores, mc.set.seed = TRUE)

} else if (parallel_method == "parLapply") {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("calculate_alpha_diversity_one_iter",
                                "otu_mat_filtered", "sample_type_filtered",
                                "rarefaction_depths", "tree"), envir = environment())
  parallel::clusterEvalQ(cl, { library(vegan); library(picante); library(dplyr) })
  alpha_div_iterations <- parallel::parLapply(cl, seq_len(n_iterations), function(i) {
    tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, sample_type_filtered,
                                         rarefaction_depths, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i))
  })
  parallel::stopCluster(cl)

} else {
  future::plan(future::multisession, workers = n_cores)
  alpha_div_iterations <- future.apply::future_lapply(
    seq_len(n_iterations),
    function(i) tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, sample_type_filtered,
                                         rarefaction_depths, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i)),
    future.seed = TRUE
  )
  future::plan(future::sequential)
}

# Validate iteration results
error_idx <- which(sapply(alpha_div_iterations, function(x)
  is.list(x) && isTRUE(x$error)))
valid_idx  <- setdiff(seq_along(alpha_div_iterations), error_idx)

if (length(error_idx) > 0)
  cat(paste0("Warning: ", length(error_idx), " iterations failed\n"))

if (length(valid_idx) == 0) {
  cat("All parallel iterations failed — falling back to sequential\n")
  alpha_div_iterations <- lapply(seq_len(n_iterations), function(i)
    calculate_alpha_diversity_one_iter(i, otu_mat_filtered, sample_type_filtered,
                                       rarefaction_depths, tree))
  valid_idx <- seq_along(alpha_div_iterations)
}

# Average across valid iterations, then join metadata
alpha_diversity_raw <- do.call(rbind, alpha_div_iterations[valid_idx]) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(dplyr::across(all_of(all_metrics), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop") %>%
  dplyr::left_join(sample_metadata, by = c("sample_id" = "sample")) %>%
  dplyr::filter(complete.cases(trt, type, pond)) %>%
  dplyr::mutate(dplyr::across(c(trt, type, pond), as.factor))

# =========================================================================
# STEP 2B: AVERAGE REPLICATES TO POND LEVEL
# =========================================================================
# Intestine has multiple shrimp replicates per pond x treatment.
# Pond is the true experimental unit; averaging collapses pseudoreplicates
# to give 1 observation per pond x type x treatment before modeling.
# Water and Sediment are unaffected (already 1 sample per pond x treatment).
# =========================================================================

cat("STEP 2B: Averaging replicates to pond level...\n")

alpha_diversity <- alpha_diversity_raw %>%
  dplyr::group_by(pond, trt, type) %>%
  dplyr::summarise(
    dplyr::across(all_of(all_metrics), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dplyr::across(c(trt, type, pond), as.factor))

check <- alpha_diversity %>%
  dplyr::count(pond, trt, type) %>%
  dplyr::filter(n > 1)

if (nrow(check) > 0) {
  cat("WARNING: Duplicate rows remain after averaging\n")
  print(check)
} else {
  cat(paste0(
    "Verified: ", nrow(alpha_diversity), " rows, 1 per pond x trt x type\n"
  ))
}

cat(paste0(
  "Final dataset: ", nrow(alpha_diversity), " observations | ",
  length(unique(alpha_diversity$trt)), " treatments | ",
  length(unique(alpha_diversity$type)), " types | ",
  length(unique(alpha_diversity$pond)), " ponds\n"
))

# =========================================================================
# STEP 3: FIT MODELS AND EXTRACT RESULTS
# =========================================================================
# Model: lm(metric ~ trt), fit separately for each sample type.
# No random effect or AR(1) — pond is not a blocking factor in this design;
# each pond x type combination receives exactly one treatment.
# emmeans + Tukey HSD CLD letters identify which treatments differ significantly.
# =========================================================================

cat("STEP 3: Fitting lm models by sample type\n")

sample_types <- levels(alpha_diversity$type)

# Nested list: models_list[[type]][[metric]]
models_list  <- list()
anova_list   <- list()
emmeans_list <- list()  # CLD tibbles: trt, emmean, SE, df, .group
pse_list     <- list()  # pse_list[[type]][[metric]]

for (tp in sample_types) {
  models_list[[tp]]  <- list()
  anova_list[[tp]]   <- list()
  emmeans_list[[tp]] <- list()
  pse_list[[tp]]     <- list()

  dat_type <- dplyr::filter(alpha_diversity, type == tp)
  cat(paste0("  ", tp, ": ", nrow(dat_type), " samples\n"))

  for (metric in primary_metrics) {
    mod <- lm(as.formula(paste0(metric, " ~ trt")),
              data = dat_type, na.action = na.omit)

    models_list[[tp]][[metric]]  <- mod
    anova_list[[tp]][[metric]]   <- anova(mod)
    pse_list[[tp]][[metric]]     <- calc_pse(mod, dat_type, metric)

    emm <- emmeans::emmeans(mod, ~ trt)
    cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")

    emmeans_list[[tp]][[metric]] <- tibble::as_tibble(cld) %>%
      dplyr::mutate(.group = stringr::str_trim(.group))
  }
}

# =========================================================================
# STEP 4: STATISTICAL SUMMARY TABLES
# =========================================================================

cat("STEP 4: Creating statistical summary tables\n")

# ---- Table 1: Treatment means + CLD + ANOVA p-values, pooled per type ----
# One row per type x metric showing emmean with CLD letter and ANOVA p-value.

table1_rows <- do.call(rbind, lapply(sample_types, function(tp) {
  do.call(rbind, lapply(primary_metrics, function(metric) {
    p_val <- anova_list[[tp]][[metric]]["trt", "Pr(>F)"]
    star  <- dplyr::case_when(
      p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ ""
    )
    emmeans_list[[tp]][[metric]] %>%
      dplyr::mutate(
        `Sample Type` = tp,
        Metric        = display_names[metric],
        `P-value`     = paste0(format_pvalue(p_val), star),
        PSE           = sprintf("%.3f", pse_list[[tp]][[metric]]),
        `Mean (CLD)`  = paste0(sprintf("%.2f", emmean), " ", tolower(.group))
      ) %>%
      dplyr::select(`Sample Type`, Metric, trt, `Mean (CLD)`, PSE, `P-value`)
  }))
})) %>%
  tidyr::pivot_wider(names_from = trt, values_from = `Mean (CLD)`)

ft1 <- flextable::flextable(table1_rows) %>%
  flextable::merge_v(j = c("Sample Type", "Metric")) %>%
  flextable::bold(j = c("Sample Type", "Metric")) %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Sample Type", "Metric"), align = "left", part = "body") %>%
  flextable::add_footer_lines(
    "Marginal means from lm(metric ~ treatment) fit separately per sample type. Different letters indicate significant pairwise differences (Tukey HSD, p < 0.05). *p<0.05, **p<0.01, ***p<0.001"
  ) %>%
  flextable::autofit()

flextable::save_as_docx(ft1, path = file.path(output_dir, "Table1_Treatment_Means_By_Type.docx"))
cat("  Table 1 saved\n")

# ---- Table 2: Pairwise contrasts (all metric x type combinations) ----

table2_rows <- do.call(rbind, lapply(sample_types, function(tp) {
  do.call(rbind, lapply(primary_metrics, function(metric) {
    emm <- emmeans::emmeans(models_list[[tp]][[metric]], ~ trt)
    tibble::as_tibble(emmeans::contrast(emm, method = "pairwise", adjust = "tukey")) %>%
      dplyr::mutate(
        `Sample Type` = tp,
        Metric        = display_names[metric],
        `P-value`     = format_pvalue(p.value),
        Significant   = dplyr::if_else(p.value < alpha_level, "Yes", "No")
      ) %>%
      dplyr::select(`Sample Type`, Metric, contrast, estimate, SE, df,
                    `P-value`, Significant)
  }))
}))

ft2 <- flextable::flextable(table2_rows) %>%
  flextable::merge_v(j = c("Sample Type", "Metric")) %>%
  flextable::bold(j = c("Sample Type", "Metric")) %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Sample Type", "Metric", "contrast"), align = "left", part = "body") %>%
  flextable::colformat_double(j = c("estimate", "SE"), digits = 3) %>%
  flextable::add_footer_lines(
    "Pairwise contrasts from emmeans (Tukey adjustment). Estimate = difference in means."
  ) %>%
  flextable::autofit()

flextable::save_as_docx(ft2, path = file.path(output_dir, "Table2_Pairwise_Contrasts.docx"))
cat("  Table 2 saved\n")

# =========================================================================
# STEP 5: CREATE PLOTS
# =========================================================================

cat("STEP 5: Creating plots — per-type treatment raincloud plots\n")

# Build all 12 individual panels: 3 types x 4 metrics
# all_type_plots[[type]][[metric]]
all_type_plots <- lapply(sample_types, function(tp) {
  dat_type     <- dplyr::filter(alpha_diversity, type == tp)
  metric_plots <- lapply(primary_metrics, function(metric) {
    create_treatment_box_plot(
      data          = dat_type,
      metric        = metric,
      cld_data      = emmeans_list[[tp]][[metric]],
      treatment_p   = anova_list[[tp]][[metric]]["trt", "Pr(>F)"],
      sample_type   = tp,
      treatment_pse = pse_list[[tp]][[metric]],
      letter_offset = letter_offset
    )
  })
  names(metric_plots) <- primary_metrics
  metric_plots
})
names(all_type_plots) <- sample_types

# =========================================================================
# STEP 6: COMBINED FIGURE — 3 rows (types) x 4 columns (metrics)
# Layout: Intestine (A-D), Sediment (E-H), Water (I-L)
# =========================================================================

cat("STEP 6: Assembling combined 3-row figure\n")

# Depth info string for subtitle
depth_str <- paste(
  sapply(names(rarefaction_depths), function(tp)
    paste0(tp, ": ", scales::comma(rarefaction_depths[[tp]]))),
  collapse = ", "
)

row_intestine <- (
  all_type_plots[["Intestine"]][["observed_richness"]] +
    all_type_plots[["Intestine"]][["shannon"]] +
    all_type_plots[["Intestine"]][["berger_parker"]] +
    all_type_plots[["Intestine"]][["faith_pd"]]
) + patchwork::plot_layout(ncol = 4)

row_sediment <- (
  all_type_plots[["Sediment"]][["observed_richness"]] +
    all_type_plots[["Sediment"]][["shannon"]] +
    all_type_plots[["Sediment"]][["berger_parker"]] +
    all_type_plots[["Sediment"]][["faith_pd"]]
) + patchwork::plot_layout(ncol = 4)

row_water <- (
  all_type_plots[["Water"]][["observed_richness"]] +
    all_type_plots[["Water"]][["shannon"]] +
    all_type_plots[["Water"]][["berger_parker"]] +
    all_type_plots[["Water"]][["faith_pd"]]
) + patchwork::plot_layout(ncol = 4)

# Stack the three rows and add shared annotation
combined_plot <- (row_intestine / row_sediment / row_water) +
  patchwork::plot_annotation(
    title    = paste0(
      "Alpha diversity across Litopenaeus vannamei fed various levels of ",
      "corn fermented protein in semi-intensive, green water production ponds"
    ),
    subtitle = paste0(
      "Based on 16S rRNA gene amplicon sequencing data. ",
      "Rarefied to compartment-specific depths (", depth_str, " reads, ",
      n_iterations, " iterations), analyzed using linear models fitted ",
      "separately for each compartment. ",
      "Rows: Intestine (A-D), Sediment (E-H), Water (I-L). ",
      "Columns: (A/E/I) Observed Richness, (B/F/J) Shannon Entropy, ",
      "(C/G/K) Berger-Parker Dominance, (D/H/L) Faith's PD."
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title    = element_text(size = 14, face = "bold",
                                 margin = ggplot2::margin(b = 6)),
    plot.subtitle = element_text(size = 10, lineheight = 1.3),
    plot.tag      = element_text(size = 12, face = "bold"),
    plot.margin   = ggplot2::margin(t = 20, r = 20, b = 20, l = 20),
    legend.position = "none"
  )

ggsave(
  filename  = file.path(output_dir, "AlphaDiversity_Combined.svg"),
  plot      = combined_plot,
  width     = 24, height = 18, units = "in",
  dpi       = 300, bg = "transparent", limitsize = FALSE
)
cat("  Combined figure saved: AlphaDiversity_Combined.svg\n")

ggsave(
  filename = file.path(output_dir, "AlphaDiversity_Combined.png"),
  plot     = combined_plot,
  width    = 24, height = 18, units = "in",
  dpi      = 300
)
cat("  Combined figure saved: AlphaDiversity_Combined.png\n")

# =========================================================================
# STEP 7: SAVE DATA
# =========================================================================

cat("STEP 7: Saving analysis data\n")
readr::write_csv(alpha_diversity, file.path(output_dir, "alpha_diversity_data.csv"))
cat("  alpha_diversity_data.csv saved\n")

# =========================================================================
# ANALYSIS SUMMARY
# =========================================================================

cat("\n")
cat("ANALYSIS COMPLETE — SUMMARY\n")
cat("\n")

cat("\nOUTPUTS:\n")
cat(strrep("-", 70), "\n")
cat("TABLES:\n")
cat("  Table1_Treatment_Means_By_Type.docx\n")
cat("  Table2_Pairwise_Contrasts.docx\n")
cat("\nFIGURES:\n")
cat("  AlphaDiversity_Combined.svg (3 rows x 4 columns)\n")
cat("  AlphaDiversity_Combined.png (3 rows x 4 columns)\n")
cat("\nDATA:\n")
cat("  alpha_diversity_data.csv\n")

cat("\n", strrep("-", 70), "\n")
cat("TREATMENT EFFECT BY SAMPLE TYPE:\n")
for (tp in sample_types) {
  cat(paste0("\n  ", tp, ":\n"))
  for (metric in primary_metrics) {
    p   <- anova_list[[tp]][[metric]]["trt", "Pr(>F)"]
    sig <- dplyr::if_else(p < alpha_level, "SIGNIFICANT *", "Not significant")
    cat(sprintf("    %-35s p = %s  (%s)\n",
                display_names[metric], format_pvalue(p), sig))
  }
}

cat("\n", strrep("-", 70), "\n")
cat("\nOUTPUT LOCATION:\n")
cat(paste0("  ", output_dir, "\n"))
cat("\n")
cat("ALL ANALYSES COMPLETE\n")
cat("\n")

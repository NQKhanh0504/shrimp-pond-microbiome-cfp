# ============================================================
# Script: 02_Rarefaction.R
# Description: Rarefaction curve analysis and sample retention
#              assessment for Litopenaeus vannamei production
#              pond microbiome data.
#              Includes per-sample-type rarefaction depth
#              selection, sample retention threshold analysis,
#              and critical taxa impact assessment.
#              Analysis is performed separately for each sample
#              type (Intestine, Sediment, Water) using
#              type-specific rarefaction depths.
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

cran_packages <- c(
  "vegan", "ggplot2", "dplyr", "patchwork", "RColorBrewer",
  "tidyr", "stringr", "scales", "tidyverse", "ggtext", "readr"
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
key_packages <- c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Github/Exp2434/Rarefaction_Analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters — type-specific rarefaction depths
chosen_depths <- list(
  Intestine = 59815,
  Sediment  = 76637,
  Water     = 85162
)

retention_thresholds        <- seq(0.80, 1.00, by = 0.025)
test_depths                 <- seq(10000, 100000, by = 10000)
use_dense_rarefaction_steps <- TRUE
rarefaction_steps_target    <- 200

# Sample filtering
filter_ponds      <- c("C3", "C6", "D1", "D4")
filter_timepoints <- c("T0", "T1", "T2", "T3", "T4")
filter_shrimp_max <- 5

# Critical taxa patterns
critical_taxa_patterns <- list(
  Genus  = c("Vibrio", "Aeromonas", "Pseudomonas", "Lactobacillus",
             "Bacillus", "Enterococcus", "Shewanella", "Flavobacterium"),
  Family = c("Vibrionaceae", "Lactobacillaceae", "Bacillaceae",
             "Rhodobacteraceae", "Flavobacteriaceae"),
  Order  = c("Vibrionales", "Lactobacillales", "Bacillales",
             "Rhodobacterales", "Pseudomonadales")
)

# Colors
type_colors <- c("Intestine" = "#E84646", "Sediment" = "#E69F00",
                 "Water" = "#0072B2", "Unknown" = "#888888")

# =========================================================================
# ANALYSIS FUNCTIONS
# =========================================================================

analyze_retention_thresholds <- function(sample_sums_df, rarefaction_df,
                                         thresholds = retention_thresholds) {

  do.call(rbind, lapply(thresholds, function(threshold) {
    target_samples  <- ceiling(nrow(sample_sums_df) * threshold)
    sorted_reads    <- sort(sample_sums_df$reads, decreasing = TRUE)
    suggested_depth <- sorted_reads[target_samples]

    closest_data <- rarefaction_df %>%
      dplyr::mutate(depth_diff = abs(subsample_size - suggested_depth)) %>%
      dplyr::filter(depth_diff == min(depth_diff)) %>%
      dplyr::summarise(
        mean_richness   = mean(richness,   na.rm = TRUE),
        median_richness = median(richness, na.rm = TRUE),
        sd_richness     = sd(richness,     na.rm = TRUE)
      )

    data.frame(
      retention_threshold = threshold,
      suggested_depth     = suggested_depth,
      samples_retained    = target_samples,
      samples_excluded    = nrow(sample_sums_df) - target_samples,
      retention_percent   = threshold * 100,
      mean_richness       = closest_data$mean_richness,
      median_richness     = closest_data$median_richness,
      sd_richness         = closest_data$sd_richness
    )
  }))
}

create_raw_averages <- function(rarefaction_df) {
  rarefaction_df %>%
    dplyr::group_by(type, subsample_size) %>%
    dplyr::summarise(
      mean_richness = mean(richness,   na.rm = TRUE),
      sd_richness   = sd(richness,     na.rm = TRUE),
      n_samples     = dplyr::n(),
      .groups       = "drop"
    )
}

analyze_critical_taxa_at_depth <- function(depth, otu_table_mat, taxa_matrix,
                                           critical_patterns, tax_level, sample_metadata) {
  level_col <- if (tax_level %in% colnames(taxa_matrix)) tax_level else tolower(tax_level)
  if (!level_col %in% colnames(taxa_matrix)) return(data.frame())

  do.call(rbind, lapply(critical_patterns, function(critical_taxon) {
    taxa_names <- taxa_matrix[, level_col]
    taxa_names <- gsub("^[a-z]__", "", taxa_names)
    taxa_names[is.na(taxa_names) | taxa_names == ""] <- "Unclassified"

    critical_asv_indices <- grep(paste0("^", critical_taxon), taxa_names)
    if (length(critical_asv_indices) == 0) return(data.frame())

    do.call(rbind, lapply(unique(sample_metadata$type), function(sample_type) {
      type_samples   <- sample_metadata$sample[sample_metadata$type == sample_type]
      type_otu_table <- otu_table_mat[type_samples, critical_asv_indices, drop = FALSE]

      impact_results <- do.call(rbind, lapply(seq_len(nrow(type_otu_table)), function(sample_idx) {
        sample_counts      <- type_otu_table[sample_idx, ]
        sample_total_reads <- sum(otu_table_mat[type_samples[sample_idx], ])

        if (sample_total_reads < depth) {
          lost_asvs  <- sum(sample_counts > 0)
          lost_reads <- sum(sample_counts)
          if (lost_asvs > 0) return(data.frame(asvs_lost = lost_asvs, total_reads_lost = lost_reads, samples_affected = 1))
        } else {
          low_asvs  <- sum(sample_counts > 0 & sample_counts < depth)
          low_reads <- sum(sample_counts[sample_counts > 0 & sample_counts < depth])
          if (low_asvs > 0) return(data.frame(asvs_lost = low_asvs, total_reads_lost = low_reads, samples_affected = 1))
        }
        data.frame(asvs_lost = 0, total_reads_lost = 0, samples_affected = 0)
      }))

      total_asvs  <- sum(impact_results$asvs_lost)
      total_reads <- sum(impact_results$total_reads_lost)
      total_samp  <- sum(impact_results$samples_affected > 0)

      if (total_asvs > 0) {
        data.frame(
          depth                 = depth,
          taxonomic_level       = tax_level,
          critical_taxon        = critical_taxon,
          sample_type           = sample_type,
          asvs_lost             = total_asvs,
          total_reads_lost      = total_reads,
          samples_affected      = total_samp,
          total_samples_of_type = length(type_samples)
        )
      } else {
        data.frame()
      }
    }))
  }))
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)

sample_sums_df <- as_tibble(phyloseq::sample_data(physeq), rownames = "sample") %>%
  dplyr::mutate(
    reads = phyloseq::sample_sums(physeq)[sample],
    timepoint_labeled = factor(
      timepoint,
      levels = c("T0", "T1", "T2", "T3", "T4"),
      labels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")
    ),
    type = dplyr::case_when(
      type == "Gut"  ~ "Intestine",
      type == "Soil" ~ "Sediment",
      TRUE ~ as.character(type)
    )
  ) %>%
  dplyr::filter(
    !is.na(trial) & trial == "1",
    pond %in% filter_ponds,
    timepoint %in% filter_timepoints,
    is.na(shrimp) | shrimp <= filter_shrimp_max
  )

physeq      <- phyloseq::prune_samples(sample_sums_df$sample, physeq)
taxa_matrix <- as.matrix(phyloseq::tax_table(physeq))

sample_types <- unique(sample_sums_df$type)

# =========================================================================
# STEP 2: CALCULATE RAREFACTION CURVES
# =========================================================================

otu_table_mat <- as(phyloseq::otu_table(physeq), "matrix")
if (phyloseq::taxa_are_rows(physeq)) otu_table_mat <- t(otu_table_mat)

max_reads <- max(sample_sums_df$reads)

if (use_dense_rarefaction_steps) {
  rare_steps <- sort(unique(c(
    seq(100,    2000,   by = 100),
    seq(2500,   10000,  by = 250),
    seq(11000,  30000,  by = 500),
    seq(31000,  60000,  by = 1000),
    seq(62000,  100000, by = 2000),
    seq(105000, max_reads, by = 5000)
  )))
  rare_steps <- rare_steps[rare_steps <= max_reads]
  if (length(rare_steps) < rarefaction_steps_target) {
    extra <- seq(min(rare_steps), max(rare_steps), length.out = rarefaction_steps_target)
    rare_steps <- sort(unique(c(rare_steps, extra)))
  }
} else {
  rare_steps <- sort(unique(c(
    seq(100,   1000,   by = 100),
    seq(1500,  10000,  by = 500),
    seq(12000, 30000,  by = 2000),
    seq(35000, max_reads, by = 5000)
  )))
  rare_steps <- rare_steps[rare_steps <= max_reads]
}

n_samples <- nrow(otu_table_mat)

rarefaction_curves <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
  sample_name   <- rownames(otu_table_mat)[i]
  sample_counts <- otu_table_mat[i, ]
  sample_total  <- sum(sample_counts)

  if (sample_total < 100) return(NULL)

  tryCatch({
    valid_steps <- rare_steps[rare_steps <= sample_total]
    if (length(valid_steps) < 2) return(NULL)

    richness_values <- sapply(valid_steps, function(depth) {
      counts_clean <- sample_counts[sample_counts > 0]
      total_clean  <- sum(counts_clean)
      if      (depth >= total_clean) length(counts_clean)
      else if (depth <= 0)           0
      else as.numeric(suppressWarnings(vegan::rarefy(counts_clean, depth)))
    })

    data.frame(
      sample         = sample_name,
      subsample_size = valid_steps,
      richness       = richness_values
    )
  }, error = function(e) {
    NULL
  })
}))

rarefaction_df <- rarefaction_curves %>%
  dplyr::left_join(sample_sums_df, by = "sample") %>%
  dplyr::filter(!is.na(type))

# =========================================================================
# STEP 3: SAMPLE RETENTION ANALYSIS — BY SAMPLE TYPE
# =========================================================================

all_retention_analysis_detailed <- list()

for (sample_type in sample_types) {
  type_sample_sums <- sample_sums_df %>% dplyr::filter(type == sample_type)
  type_rarefaction <- rarefaction_df  %>% dplyr::filter(type == sample_type)

  all_retention_analysis_detailed[[sample_type]] <-
    analyze_retention_thresholds(type_sample_sums, type_rarefaction)
}

# =========================================================================
# STEP 4: RETENTION ACROSS DEPTH RANGE — BY SAMPLE TYPE
# =========================================================================

all_retention_analysis <- do.call(rbind, lapply(sample_types, function(sample_type) {
  type_sample_sums <- sample_sums_df %>% dplyr::filter(type == sample_type)

  do.call(rbind, lapply(test_depths, function(depth) {
    retained <- sum(type_sample_sums$reads >= depth)
    total    <- nrow(type_sample_sums)
    data.frame(
      sample_type       = sample_type,
      depth             = depth,
      samples_retained  = retained,
      total_samples     = total,
      retention_percent = round(100 * retained / total, 1)
    )
  }))
}))

# =========================================================================
# STEP 5: CRITICAL TAXA IMPACT ANALYSIS — BY SAMPLE TYPE
# =========================================================================

all_critical_impacts <- do.call(rbind, lapply(sample_types, function(sample_type) {
  type_sample_sums <- sample_sums_df %>% dplyr::filter(type == sample_type)

  do.call(rbind, lapply(test_depths, function(depth) {
    depth_impacts <- do.call(rbind, lapply(names(critical_taxa_patterns), function(tax_level) {
      analyze_critical_taxa_at_depth(
        depth, otu_table_mat, taxa_matrix,
        critical_taxa_patterns[[tax_level]], tax_level, type_sample_sums
      )
    }))

    if (!is.null(depth_impacts) && nrow(depth_impacts) > 0) depth_impacts else data.frame()
  }))
}))

if (is.null(all_critical_impacts)) all_critical_impacts <- data.frame()

# =========================================================================
# STEP 6: DETAILED ANALYSIS AT CHOSEN DEPTH — BY SAMPLE TYPE
# =========================================================================

all_sample_analysis <- do.call(rbind, lapply(sample_types, function(sample_type) {
  type_sample_sums <- sample_sums_df %>% dplyr::filter(type == sample_type)
  chosen_depth     <- chosen_depths[[sample_type]]

  excluded_samples <- type_sample_sums %>% dplyr::filter(reads < chosen_depth)
  retained_samples <- type_sample_sums %>% dplyr::filter(reads >= chosen_depth)

  retention_pct <- round(100 * nrow(retained_samples) / nrow(type_sample_sums), 1)
  exclusion_pct <- round(100 * nrow(excluded_samples) / nrow(type_sample_sums), 1)

  if (nrow(excluded_samples) > 0) {
    excluded_summary <- excluded_samples %>%
      dplyr::mutate(deficit = chosen_depth - reads) %>%
      dplyr::arrange(reads)

    exclusion_by_type <- excluded_samples %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(n_excluded = dplyr::n(), .groups = "drop") %>%
      dplyr::left_join(
        type_sample_sums %>%
          dplyr::group_by(type) %>%
          dplyr::summarise(total_samples = dplyr::n(), .groups = "drop"),
        by = "type"
      ) %>%
      dplyr::mutate(
        retention_rate = round(100 * (total_samples - n_excluded) / total_samples, 1),
        exclusion_rate = round(100 * n_excluded / total_samples, 1)
      )
  }

  type_sample_sums %>%
    dplyr::mutate(
      retained                 = reads >= chosen_depth,
      chosen_rarefaction_depth = chosen_depth,
      exclusion_deficit        = dplyr::if_else(reads < chosen_depth, chosen_depth - reads, 0L),
      status                   = dplyr::if_else(reads >= chosen_depth, "Retained", "Excluded")
    )
}))

all_chosen_depth_impacts <- data.frame()
if (nrow(all_critical_impacts) > 0) {
  all_chosen_depth_impacts <- do.call(rbind, lapply(sample_types, function(sample_type) {
    chosen_depth <- chosen_depths[[sample_type]]
    all_critical_impacts %>%
      dplyr::filter(sample_type == !!sample_type, depth == chosen_depth)
  }))
}

# =========================================================================
# STEP 7: CREATE RAREFACTION PLOTS
# =========================================================================

raw_avg_curves <- create_raw_averages(rarefaction_df)

has_sample_types <- length(unique(rarefaction_df$type)) > 1 &&
  !"Unknown" %in% unique(rarefaction_df$type)

# Individual rarefaction curves by type — type-specific depth lines
if (has_sample_types) {
  plot_by_type <- ggplot(rarefaction_df,
                         aes(x = subsample_size, y = richness, group = sample)) +
    geom_line(aes(color = type), alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(values = type_colors, name = "Compartment")
} else {
  plot_by_type <- ggplot(rarefaction_df,
                         aes(x = subsample_size, y = richness, group = sample)) +
    geom_line(alpha = 0.8, linewidth = 0.8, color = "#0072B2")
}

for (sample_type in sample_types) {
  plot_by_type <- plot_by_type +
    geom_vline(
      xintercept = chosen_depths[[sample_type]],
      linetype   = "dashed",
      color      = type_colors[[sample_type]],
      linewidth  = 0.8,
      alpha      = 0.8
    )
}

depth_text <- paste(
  sapply(sample_types, function(st) paste0(st, ": ", scales::comma(chosen_depths[[st]]))),
  collapse = ", "
)

plot_by_type <- plot_by_type +
  scale_x_continuous(
    labels = scales::comma_format(),
    breaks = seq(0, 180000, by = 30000),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 8000, 500),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Rarefaction Curves by Compartment",
    x     = "Sequencing Depth (Number of Reads)",
    y     = "Observed Richness (Number of ASVs)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Average curves by type — type-specific depth lines
plot_average_type <- ggplot(raw_avg_curves,
                            aes(x = subsample_size, y = mean_richness, color = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = type_colors, name = "Compartment")

for (sample_type in sample_types) {
  plot_average_type <- plot_average_type +
    geom_vline(
      xintercept = chosen_depths[[sample_type]],
      linetype   = "dashed",
      color      = type_colors[[sample_type]],
      linewidth  = 0.8,
      alpha      = 0.8
    )
}

plot_average_type <- plot_average_type +
  scale_x_continuous(
    labels = scales::comma_format(),
    breaks = seq(0, 180000, by = 30000),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 8000, 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Average Rarefaction Curves by Compartment",
    x     = "Sequencing Depth (Number of Reads)",
    y     = "Mean Observed Richness (Number of ASVs)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Read depth distribution — type-specific depth lines
plot_distribution <- ggplot(sample_sums_df, aes(x = reads, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.7, color = NA) +
  scale_fill_manual(values = type_colors, name = "Compartment")

for (sample_type in sample_types) {
  plot_distribution <- plot_distribution +
    geom_vline(
      xintercept = chosen_depths[[sample_type]],
      linetype   = "dashed",
      color      = type_colors[[sample_type]],
      linewidth  = 0.8,
      alpha      = 0.8
    )
}

plot_distribution <- plot_distribution +
  scale_x_continuous(
    labels = scales::comma_format(),
    breaks = seq(0, 180000, by = 30000),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, 2),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Distribution of Sequencing Depths",
    x     = "Number of Reads per Sample",
    y     = "Number of Samples"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Combined figure
combined_plot <- (plot_by_type) + (plot_distribution) +
  plot_annotation(
    title    = paste0(
      "Rarefaction curves and sequencing depth distribution across compartments ",
      "in Litopenaeus vannamei production ponds"
    ),
    subtitle = paste0(
      "Rarefied to sample type-specific depths (", depth_text, " reads) across ",
      "time points (Week 0-8) within Intestine, Sediment, and Water compartments. ",
      "Colored dashed lines indicate chosen rarefaction depth per compartment. ",
      "Panel A: individual rarefaction curves by compartment. ",
      "Panel B: distribution of sequencing depths per sample."
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9,      hjust = 0),
    plot.tag      = element_text(size = 11,     face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "Rarefaction_Analysis.png"),
  plot     = combined_plot,
  width    = 16, height = 8, units = "in", dpi = 300
)

# =========================================================================
# STEP 8: EXPORT RESULTS — BY SAMPLE TYPE
# =========================================================================

# Per-type comprehensive sample analysis
for (sample_type in sample_types) {
  type_sample_analysis <- all_sample_analysis %>%
    dplyr::filter(type == sample_type)
  readr::write_csv(
    type_sample_analysis,
    file.path(output_dir, paste0("Comprehensive_Sample_Analysis_", sample_type, ".csv"))
  )
}

# Per-type method analysis (retention thresholds + fixed depths)
for (sample_type in sample_types) {
  retention_detailed <- all_retention_analysis_detailed[[sample_type]]

  retention_formatted <- retention_detailed %>%
    dplyr::mutate(
      sample_type = sample_type,
      method      = paste0("Retention ", retention_threshold * 100, "%"),
      method_type = "Sample Retention"
    ) %>%
    dplyr::select(sample_type, method, method_type, suggested_depth, samples_retained,
                  retention_percent, mean_richness, median_richness, sd_richness)

  type_retention_analysis <- all_retention_analysis %>%
    dplyr::filter(sample_type == !!sample_type)

  depth_analysis_formatted <- type_retention_analysis %>%
    dplyr::mutate(
      method          = paste0("Fixed Depth ", scales::comma(depth)),
      method_type     = "Fixed Depth Analysis",
      suggested_depth = depth,
      mean_richness   = NA_real_,
      median_richness = NA_real_,
      sd_richness     = NA_real_
    ) %>%
    dplyr::select(sample_type, method, method_type, suggested_depth, samples_retained,
                  retention_percent, mean_richness, median_richness, sd_richness)

  method_analysis_comprehensive <- dplyr::bind_rows(retention_formatted, depth_analysis_formatted)
  readr::write_csv(
    method_analysis_comprehensive,
    file.path(output_dir, paste0("Comprehensive_Method_Analysis_", sample_type, ".csv"))
  )
}

# Per-type critical taxa impacts
if (nrow(all_critical_impacts) > 0) {
  for (sample_type in sample_types) {
    chosen_depth           <- chosen_depths[[sample_type]]
    type_critical_impacts  <- all_critical_impacts %>%
      dplyr::filter(sample_type == !!sample_type)

    if (nrow(type_critical_impacts) > 0) {
      type_critical_impacts_flagged <- type_critical_impacts %>%
        dplyr::mutate(
          is_chosen_depth   = depth == chosen_depth,
          chosen_depth_flag = dplyr::if_else(depth == chosen_depth, "CHOSEN_DEPTH", "")
        )
      readr::write_csv(
        type_critical_impacts_flagged,
        file.path(output_dir, paste0("Critical_Taxa_Impacts_Analysis_", sample_type, ".csv"))
      )
    }
  }
}

readr::write_csv(rarefaction_df,  file.path(output_dir, "Rarefaction_Curves_Data.csv"))
readr::write_csv(raw_avg_curves,  file.path(output_dir, "Raw_Average_Curves.csv"))

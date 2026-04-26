# ============================================================
# Script: 08_DifferentialAnalysis.R
# Description: ANCOM-BC2 differential abundance analysis for
#              CFP shrimp microbiome data.
#              Treatment comparison design: 4 treatments
#              (Basal, CFP 5%, CFP 10%, CFP 20%), reference = Basal.
#              Analysis run separately within each compartment
#              (Intestine, Sediment, Water).
#              Three comparisons per compartment:
#                CFP 5%  vs Basal (Panel A)
#                CFP 10% vs Basal (Panel B)
#                CFP 20% vs Basal (Panel C)
# Input: step4_phyloseq_object.rds
# Output: Per-compartment volcano figures and CSV result tables.
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

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_packages <- c("phyloseq", "ANCOMBC", "microbiome")
cran_packages <- c("tidyverse", "janitor", "ggrepel", "patchwork",
                   "ggtext", "scales", "parallel", "digest")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
for (pkg in c("phyloseq", "ANCOMBC", "tidyverse", "ggrepel")) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

set.seed(54)

# File paths
phyloseq_path <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/16S_V3V4_MicrobiomePipeline/step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/ANCOMBC_Analysis"
results_dir   <- file.path(output_dir, "Results")
plots_dir     <- file.path(output_dir, "Volcano_Plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# Statistical parameters
alpha_threshold        <- 0.05
prevalence_cutoff_high <- 0.1   # For Phylum and Class levels
prevalence_cutoff_low  <- 0.2   # For Order, Family, Genus levels
taxonomic_levels       <- c("Genus")

# Treatment factor order — Basal is reference (first level)
treatment_order <- c("Basal", "CFP5", "CFP10", "CFP20")

# Rarefaction depths (per compartment, informational)
rarefaction_depths <- c(Intestine = 59815, Sediment = 76637, Water = 85162)

# Three comparisons per compartment (Basal = reference, first factor level):
#   lfc_trt_clean1 = CFP5  vs Basal  (positive = CFP5 enriched)
#   lfc_trt_clean2 = CFP10 vs Basal  (positive = CFP10 enriched)
#   lfc_trt_clean3 = CFP20 vs Basal  (positive = CFP20 enriched)
comparisons_map <- list(
  "CFP5 vs Basal"  = list(
    lfc_col = "lfc_trt_clean1",
    se_col  = "se_trt_clean1",
    p_col   = "p_trt_clean1",
    q_col   = "q_trt_clean1",
    ss_col  = "passed_ss_trt_clean1"
  ),
  "CFP10 vs Basal" = list(
    lfc_col = "lfc_trt_clean2",
    se_col  = "se_trt_clean2",
    p_col   = "p_trt_clean2",
    q_col   = "q_trt_clean2",
    ss_col  = "passed_ss_trt_clean2"
  ),
  "CFP20 vs Basal" = list(
    lfc_col = "lfc_trt_clean3",
    se_col  = "se_trt_clean3",
    p_col   = "p_trt_clean3",
    q_col   = "q_trt_clean3",
    ss_col  = "passed_ss_trt_clean3"
  )
)

# Plotting parameters
max_labels_per_group <- 5

# Colors consistent across scripts
treatment_colors <- c(
  "Basal"  = "#0072B2",
  "CFP5"   = "#E84646",
  "CFP10"  = "#E69F00",
  "CFP20"  = "#009E73"
)
type_shapes <- c("Intestine" = 21, "Sediment" = 24, "Water" = 22)

# Cache/skip flag — set TRUE to delete cache and force a fresh run
force_rerun <- FALSE

# Parallel processing
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy_hierarchical <- function(tax_row) {
  clean_string <- function(s) {
    s <- gsub("^[a-z]__", "", s)
    s <- gsub("_.*", "", s)
    s <- gsub("^(Bacteria|Archaea)_", "", s)
    trimws(s)
  }

  genus  <- clean_string(tax_row["Genus"])
  family <- clean_string(tax_row["Family"])
  order  <- clean_string(tax_row["Order"])
  class  <- clean_string(tax_row["Class"])

  if (!is.na(genus) && genus != "" && !genus %in% c("Bacteria", "Archaea")) {
    return(genus)
  } else if (!is.na(family) && family != "" && !family %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", family))
  } else if (!is.na(order) && order != "" && !order %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", order))
  } else if (!is.na(class) && class != "" && !class %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", class))
  } else {
    return("Unidentified")
  }
}

extract_results <- function(ancombc_result, compartment, alpha = 0.05) {
  if (is.null(ancombc_result) || is.null(ancombc_result$res)) return(NULL)

  res_df <- ancombc_result$res

  # NaN flag: any stat column is NA/NaN/Inf
  stat_cols <- unlist(lapply(comparisons_map, function(x) c(x$lfc_col, x$q_col)))
  stat_cols <- stat_cols[stat_cols %in% colnames(res_df)]
  res_df$NaN_Flag <- apply(
    res_df[, stat_cols, drop = FALSE], 1,
    function(r) any(is.nan(r) | is.infinite(r) | is.na(r))
  )

  results_long <- do.call(rbind, lapply(names(comparisons_map), function(comp) {
    pair <- comparisons_map[[comp]]
    if (!all(c(pair$lfc_col, pair$q_col, pair$ss_col) %in% colnames(res_df))) return(NULL)

    res_df %>%
      dplyr::select(
        taxon, NaN_Flag,
        lfc       = dplyr::all_of(pair$lfc_col),
        q_val     = dplyr::all_of(pair$q_col),
        sens_pass = dplyr::all_of(pair$ss_col)
      ) %>%
      dplyr::mutate(
        comparison  = comp,
        compartment = compartment
      )
  }))

  sig <- dplyr::filter(results_long, q_val < alpha, !NaN_Flag, sens_pass)

  cat(sprintf("  %s: %d significant taxa across all comparisons\n", compartment, nrow(sig)))

  sig_counts <- dplyr::count(sig, comparison)
  for (i in seq_len(nrow(sig_counts))) {
    cat(sprintf("    %s: %d\n", sig_counts$comparison[i], sig_counts$n[i]))
  }

  list(full_results = results_long, significant_results = sig)
}

prepare_volcano_data <- function(full_results, comp,
                                 max_labels = 5, sig_threshold = 0.05) {
  cfp_group <- stringr::str_split(comp, " vs ")[[1]][1]

  comp_data <- dplyr::filter(full_results, comparison == comp)

  # Positive LFC = CFP group enriched; negative LFC = Basal enriched.
  # No sign flip needed — Basal is the reference (first factor level).

  taxa_to_label <- comp_data %>%
    dplyr::filter(
      q_val < sig_threshold, !NaN_Flag,
      is.finite(lfc), is.finite(q_val), sens_pass
    ) %>%
    dplyr::mutate(
      importance = -log10(pmax(q_val, 1e-10)) * abs(lfc),
      group      = dplyr::if_else(lfc > 0, cfp_group, "Basal")
    ) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(importance), .by_group = TRUE) %>%
    dplyr::slice_head(n = max_labels) %>%
    dplyr::ungroup() %>%
    dplyr::pull(taxon)

  comp_data %>%
    dplyr::filter(is.finite(lfc), is.finite(q_val)) %>%
    dplyr::mutate(
      significance = dplyr::if_else(
        q_val < sig_threshold & sens_pass, "Significant", "Not Significant"
      ),
      neg_log_q   = -log10(pmax(q_val, 1e-10)),
      enriched_in = dplyr::if_else(lfc > 0, cfp_group, "Basal"),
      label       = dplyr::if_else(taxon %in% taxa_to_label, taxon, "")
    )
}

make_volcano_panel <- function(plot_data, comp,
                               point_shape = 21, sig_threshold = 0.05) {
  cfp_group  <- stringr::str_split(comp, " vs ")[[1]][1]
  col_values <- treatment_colors[c("Basal", cfp_group)]

  ggplot(plot_data, aes(x = lfc, y = neg_log_q)) +
    geom_hline(
      yintercept = -log10(sig_threshold),
      linetype = "dashed", color = "black", linewidth = 0.8
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed", color = "black", linewidth = 0.8
    ) +
    geom_point(
      aes(
        color = enriched_in,
        fill  = ifelse(significance == "Significant", enriched_in, NA)
      ),
      shape = point_shape, size = 3, alpha = 0.6, stroke = 0.8
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size               = 3,
      max.overlaps       = Inf,
      box.padding        = 0.4,
      point.padding      = 0.2,
      force              = 5,
      force_pull         = 0.1,
      min.segment.length = 0,
      segment.size       = 0.4,
      segment.color      = "black",
      segment.alpha      = 0.6,
      segment.linetype   = 1,
      segment.curvature  = 0
    ) +
    scale_color_manual(values = col_values, name = "Enriched in") +
    scale_fill_manual(
      values   = col_values, name = "Enriched in",
      na.value = NA, guide = "none"
    ) +
    scale_x_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8)) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    labs(
      x = "Log fold change",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

cat("Loading phyloseq object...\n")
physeq <- readRDS(phyloseq_path)
cat(sprintf("Loaded: %d samples, %d ASVs\n", nsamples(physeq), ntaxa(physeq)))

# Apply hierarchical taxonomic naming
cat("Applying hierarchical taxonomic naming...\n")
tax_df             <- as.data.frame(tax_table(physeq))
tax_df$Genus_Clean <- apply(tax_df, 1, clean_taxonomy_hierarchical)

new_tax_table <- tax_table(physeq)
for (i in seq_len(ncol(new_tax_table))) {
  new_tax_table[, i] <- gsub("^[a-z]__", "", new_tax_table[, i])
  new_tax_table[, i] <- gsub("_.*",       "", new_tax_table[, i])
  new_tax_table[, i][is.na(new_tax_table[, i]) | new_tax_table[, i] == ""] <- "Unidentified"
}
new_tax_table[, "Genus"] <- tax_df$Genus_Clean
tax_table(physeq)        <- new_tax_table
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# =========================================================================
# STEP 2: METADATA STANDARDIZATION
# =========================================================================

cat("Standardizing metadata...\n")

meta_df <- as(sample_data(physeq), "data.frame") %>%
  janitor::clean_names() %>%
  dplyr::filter(trial == 1) %>%
  dplyr::mutate(
    type_clean = dplyr::case_when(
      stringr::str_detect(tolower(type), "gut|intestin")  ~ "Intestine",
      stringr::str_detect(tolower(type), "water")         ~ "Water",
      stringr::str_detect(tolower(type), "soil|sediment") ~ "Sediment",
      TRUE ~ as.character(type)
    ),
    trt_clean = dplyr::case_when(
      trt == "1" ~ "Basal",
      trt == "2" ~ "CFP5",
      trt == "3" ~ "CFP10",
      trt == "4" ~ "CFP20",
      TRUE ~ as.character(trt)
    ),
    pond_clean = toupper(trimws(as.character(pond)))
  ) %>%
  dplyr::mutate(
    type_clean = factor(type_clean, levels = c("Intestine", "Sediment", "Water")),
    trt_clean  = factor(trt_clean,  levels = treatment_order)
  )

sample_data(physeq)             <- sample_data(meta_df)
rownames(sample_data(physeq))   <- rownames(meta_df)

physeq_final <- prune_taxa(taxa_sums(physeq) > 0, physeq)
cat(sprintf("Final dataset: %d samples, %d taxa\n",
            nsamples(physeq_final), ntaxa(physeq_final)))

# =========================================================================
# STEP 3: CACHE CHECK
# =========================================================================

cat("Checking for cached results...\n")

results_rds_path <- file.path(output_dir, "ANCOMBC_Results.rds")

current_params <- digest::digest(list(
  alpha      = alpha_threshold,
  prv_high   = prevalence_cutoff_high,
  prv_low    = prevalence_cutoff_low,
  tax_levels = taxonomic_levels,
  design     = "trt_within_compartment_BH"
))

if (force_rerun && file.exists(results_rds_path)) {
  file.remove(results_rds_path)
  cat("Cache deleted — forcing fresh analysis run\n")
}

perform_analysis       <- TRUE
ancombc_by_compartment <- list()
stats_by_compartment   <- list()

if (file.exists(results_rds_path)) {
  cat("Found cached results — loading...\n")
  tryCatch({
    cached <- readRDS(results_rds_path)
    if (!is.null(cached$params_hash) && cached$params_hash == current_params &&
        !is.null(cached$ancombc_by_compartment) &&
        !is.null(cached$stats_by_compartment)) {
      ancombc_by_compartment <- cached$ancombc_by_compartment
      stats_by_compartment   <- cached$stats_by_compartment
      perform_analysis       <- FALSE
      cat("Parameters match — using cached results\n")
    } else {
      cat("Parameters changed — re-running analysis\n")
    }
  }, error = function(e) {
    cat(sprintf("Cache load failed: %s — re-running\n", e$message))
  })
} else {
  cat("No cache found — running ANCOM-BC2 analysis...\n")
}

# =========================================================================
# STEP 4: ANCOM-BC2 ANALYSIS BY COMPARTMENT
# =========================================================================

if (perform_analysis) {
  cat("Running ANCOM-BC2 analysis by compartment...\n")

  for (st in levels(meta_df$type_clean)) {
    cat(sprintf("  Analyzing %s...\n", st))

    physeq_st <- subset_samples(physeq_final, type_clean == st)
    physeq_st <- prune_taxa(taxa_sums(physeq_st) > 0, physeq_st)

    n_st <- nsamples(physeq_st)
    cat(sprintf("    %s: %d samples\n", st, n_st))

    if (n_st < 4) {
      warning(sprintf("    %s: too few samples — skipping", st))
      next
    }

    meta_st    <- as(sample_data(physeq_st), "data.frame")
    trt_levels <- unique(droplevels(meta_st$trt_clean))
    cat(sprintf("    %s: treatments present — %s\n",
                st, paste(trt_levels, collapse = ", ")))

    if (length(trt_levels) < 2) {
      warning(sprintf("    %s: fewer than 2 treatment levels — skipping", st))
      next
    }

    for (tax_level in taxonomic_levels) {
      prv_cutoff <- if (tax_level %in% c("Phylum", "Class")) {
        prevalence_cutoff_high
      } else {
        prevalence_cutoff_low
      }
      key <- paste(st, tax_level, sep = "_")

      ancombc_result <- tryCatch(
        ANCOMBC::ancombc2(
          data         = physeq_st,
          tax_level    = tax_level,
          fix_formula  = "trt_clean",
          rand_formula = NULL,
          p_adj_method = "BH",
          prv_cut      = prv_cutoff,
          lib_cut      = 0,
          s0_perc      = 0.5,
          group        = "trt_clean",
          struc_zero   = TRUE,
          neg_lb       = TRUE,
          alpha        = alpha_threshold,
          n_cl         = n_cores,
          verbose      = FALSE,
          global       = FALSE,
          pairwise     = FALSE,
          dunnet       = FALSE,
          trend        = FALSE,
          pseudo_sens  = TRUE
        ),
        error = function(e) {
          warning(sprintf("    FAILED %s: %s", key, e$message))
          NULL
        }
      )

      if (is.null(ancombc_result)) {
        warning(sprintf("    %s: ancombc2 returned NULL — no results", key))
      } else {
        cat(sprintf("    %s: ancombc2 succeeded — %d taxa in results\n",
                    key, nrow(ancombc_result$res)))
        ancombc_by_compartment[[key]] <- ancombc_result

        results <- extract_results(ancombc_result, st, alpha_threshold)
        if (!is.null(results)) {
          stats_by_compartment[[key]] <- results
          readr::write_csv(
            results$full_results,
            file.path(results_dir,
                      paste0("ANCOMBC_", st, "_", tax_level, "_Full.csv"))
          )
          if (nrow(results$significant_results) > 0) {
            readr::write_csv(
              results$significant_results,
              file.path(results_dir,
                        paste0("ANCOMBC_", st, "_", tax_level, "_Significant.csv"))
            )
          }
        }
      }
      gc()
    }
  }

  saveRDS(
    list(
      ancombc_by_compartment = ancombc_by_compartment,
      stats_by_compartment   = stats_by_compartment,
      params_hash            = current_params,
      params                 = list(
        alpha      = alpha_threshold,
        tax_levels = taxonomic_levels,
        design     = "trt_within_compartment_BH"
      )
    ),
    results_rds_path
  )
  cat("Analysis complete — results cached\n")
}

# =========================================================================
# STEP 5: EXPORT COMPREHENSIVE CSV FILES
# =========================================================================

cat("Exporting per-comparison CSV files...\n")

for (key in names(ancombc_by_compartment)) {
  parts     <- stringr::str_split(key, "_")[[1]]
  st        <- parts[1]
  tax_level <- parts[2]
  res       <- ancombc_by_compartment[[key]]
  if (is.null(res) || is.null(res$res)) next

  res_df <- res$res

  for (comp in names(comparisons_map)) {
    pair <- comparisons_map[[comp]]
    if (!all(c(pair$lfc_col, pair$q_col, pair$ss_col) %in% colnames(res_df))) next

    out <- res_df %>%
      dplyr::select(
        taxon,
        LFC              = dplyr::all_of(pair$lfc_col),
        SE               = dplyr::all_of(pair$se_col),
        P_value          = dplyr::all_of(pair$p_col),
        Adjusted_P_value = dplyr::all_of(pair$q_col),
        Sensitivity_Pass = dplyr::all_of(pair$ss_col)
      ) %>%
      dplyr::mutate(
        Compartment = st,
        Comparison  = comp,
        Significant = Adjusted_P_value < alpha_threshold & Sensitivity_Pass,
        Direction   = dplyr::case_when(
          LFC > 0 ~ stringr::str_split(comp, " vs ")[[1]][1],
          LFC < 0 ~ "Basal",
          TRUE    ~ "No change"
        )
      )

    comp_safe <- stringr::str_replace_all(comp, "[ %]", "_")
    fname     <- paste0("ANCOMBC_", st, "_", tax_level, "_", comp_safe, ".csv")
    readr::write_csv(out, file.path(results_dir, fname))
    cat(sprintf("  Saved: %s (%d significant)\n",
                fname, sum(out$Significant, na.rm = TRUE)))
  }
}

# =========================================================================
# STEP 6: VOLCANO PLOTS
# =========================================================================

cat("Generating volcano plots...\n")

for (st in levels(meta_df$type_clean)) {
  key <- paste0(st, "_Genus")
  if (is.null(stats_by_compartment[[key]])) {
    warning(sprintf("  Skipping %s — no results", st))
    next
  }

  full_res <- stats_by_compartment[[key]]$full_results

  # Build one panel per comparison
  comp_names <- names(comparisons_map)
  panels     <- vector("list", length(comp_names))
  names(panels) <- comp_names

  for (comp in comp_names) {
    plot_data <- prepare_volcano_data(
      full_res, comp,
      max_labels    = max_labels_per_group,
      sig_threshold = alpha_threshold
    )
    if (nrow(plot_data) == 0) {
      panels[[comp]] <- NULL
      next
    }
    panels[[comp]] <- make_volcano_panel(
      plot_data, comp,
      point_shape   = type_shapes[st],
      sig_threshold = alpha_threshold
    )

    # Save individual panel
    comp_safe <- stringr::str_replace_all(comp, "[ %]", "_")
    fname_ind <- paste0("Volcano_", st, "_", comp_safe, ".png")
    ggplot2::ggsave(
      file.path(plots_dir, fname_ind),
      panels[[comp]],
      width = 8, height = 7, dpi = 300
    )
  }

  panels_valid <- Filter(Negate(is.null), panels)

  if (length(panels_valid) == 0) {
    warning(sprintf("  Skipping %s — no plottable panels", st))
    next
  }

  # Combined 1x3 figure: A = CFP5 vs Basal, B = CFP10 vs Basal, C = CFP20 vs Basal
  plot_main_title <- paste0(
    "Volcano plot of genus-level differential abundance across ",
    "*Litopenaeus vannamei* fed various levels of corn fermented protein"
  )
  plot_subtitle <- paste0(
    "ANCOM-BC2 analysis comparing diet groups (CFP 5%, CFP 10%, CFP 20% vs Basal) ",
    "within the ", tolower(st), " compartment. ",
    "FDR-adjusted p-value threshold: 0.05. ",
    "Filled points passed both significance and sensitivity analysis. ",
    "Labels show the ", max_labels_per_group, " most important genera per group. ",
    "Panel (**A**) CFP 5% vs Basal; Panel (**B**) CFP 10% vs Basal; ",
    "Panel (**C**) CFP 20% vs Basal."
  )

  combined <- patchwork::wrap_plots(panels_valid, ncol = 3) +
    patchwork::plot_annotation(
      title      = plot_main_title,
      subtitle   = plot_subtitle,
      tag_levels = "A"
    ) &
    theme(
      plot.title    = ggtext::element_markdown(
        face = "bold", size = 12, hjust = 0,
        margin = ggplot2::margin(t = 10, b = 10)
      ),
      plot.subtitle = ggtext::element_markdown(
        size = 9, hjust = 0, lineheight = 1.3,
        margin = ggplot2::margin(t = 5, b = 10)
      ),
      plot.tag      = element_text(size = 12, face = "bold")
    )

  fname_combined <- paste0("Volcano_Genus_", st, "_Combined.png")
  ggplot2::ggsave(
    file.path(plots_dir, fname_combined),
    combined,
    width = 24, height = 8, dpi = 300
  )
  cat(sprintf("  Saved: %s\n", fname_combined))
}

# =========================================================================
# ANALYSIS SUMMARY
# =========================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("OUTPUTS:\n")
cat(strrep("-", 70), "\n")
cat("TABLES (per compartment x comparison):\n")
cat("  ANCOMBC_{Compartment}_Genus_Full.csv\n")
cat("  ANCOMBC_{Compartment}_Genus_Significant.csv\n")
cat("  ANCOMBC_{Compartment}_Genus_{Comparison}.csv\n")
cat("\nFIGURES:\n")
cat("  Volcano_Genus_Intestine_Combined.png\n")
cat("  Volcano_Genus_Sediment_Combined.png\n")
cat("  Volcano_Genus_Water_Combined.png\n")
cat("  (plus individual panels per compartment x comparison)\n")
cat("\nOUTPUT LOCATION:\n")
cat("  Results:      ", results_dir, "\n")
cat("  Volcano plots:", plots_dir,   "\n\n")

# Significance summary
cat("SIGNIFICANCE SUMMARY:\n")
for (key in names(stats_by_compartment)) {
  st  <- stringr::str_split(key, "_")[[1]][1]
  res <- stats_by_compartment[[key]]
  if (is.null(res)) next
  sig_counts <- dplyr::count(res$significant_results, comparison)
  cat(sprintf("  %s:\n", st))
  for (i in seq_len(nrow(sig_counts))) {
    cat(sprintf("    %s: %d significant genera\n",
                sig_counts$comparison[i], sig_counts$n[i]))
  }
}

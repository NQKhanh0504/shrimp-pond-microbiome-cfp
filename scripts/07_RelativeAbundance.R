# ============================================================
# Script: 07_RelativeAbundance.R
# Description: Stacked bar / alluvial diagram showing genus-level
#              treatment composition across Intestine, Sediment,
#              and Water compartments in Litopenaeus vannamei
#              production ponds fed Basal, CFP 5%, CFP 10%, or
#              CFP 20% diets.
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
  "phyloseq", "vegan", "ggplot2", "patchwork", "ggtext",
  "janitor", "tidyverse", "scales", "parallel", "future",
  "future.apply"
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
key_packages <- c("phyloseq", "vegan", "ggplot2", "patchwork", "dplyr", "tidyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# Rarefaction toggle — TRUE = iterative rarefaction; FALSE = relative abundance
use_rarefaction <- FALSE

# File paths
phyloseq_path <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/16S_V3V4_MicrobiomePipeline/step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/AlluvialPlot_Rarefaction"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cache paths — one per mode so switching modes forces a fresh run
cache_file <- if (use_rarefaction) {
  file.path(output_dir, "alluvial_rarefied_treatment_results.rds")
} else {
  file.path(output_dir, "alluvial_relabund_treatment_results.rds")
}

# Analysis parameters
n_iterations <- 100
random_seed  <- 54
set.seed(random_seed)

# Type-specific rarefaction depths (only used when use_rarefaction = TRUE)
rarefaction_depths <- list(
  Intestine = 59815,
  Sediment  = 76637,
  Water     = 85162
)

# Treatment factor order
treatment_order <- c("Basal", "CFP 5%", "CFP 10%", "CFP 20%")

# Minimum relative abundance thresholds per taxonomic level;
# taxa below threshold are grouped as "Other"
taxa_thresholds <- list(
  phylum = 0.01,
  class  = 0.01,
  order  = 0.01,
  family = 0.02,
  genus  = 0.02
)

# Minimum proportion within a stratum for a taxon label to be shown on plot
label_threshold <- 0.03

# Method label used in filenames
method_label <- if (use_rarefaction) "Rarefied" else "RelAbund"

# Color palette
bright_palette <- c(
  "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FECA57", "#FF9FF3",
  "#54A0FF", "#5F27CD", "#00D2D3", "#FF9F43", "#10AC84", "#EE5A24",
  "#0984E3", "#6C5CE7", "#A29BFE", "#FD79A8", "#FDCB6E", "#6C5CE7",
  "#74B9FF", "#81ECEC", "#FF5252", "#00B894", "#2ED573", "#FF4757",
  "#1E90FF", "#3742FA", "#70A1FF", "#FF6348", "#FF9F1A", "#32FF7E",
  "#7EFFF5", "#18DCFF", "#7D5FFF", "#5352ED", "#C56CF0", "#D6A2E8",
  "#2F3542", "#F368E0", "#00D8D6", "#9B59B6"
)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy <- function(x) {
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("_.*", "", x)
  x[is.na(x) | x == ""] <- "Unidentified"
  return(x)
}

assign_colors <- function(plot_data, taxon_col = "taxon") {
  taxa   <- unique(plot_data[[taxon_col]])
  n      <- length(taxa)
  colors <- bright_palette[seq_len(min(n, length(bright_palette)))]
  if (n > length(bright_palette))
    colors <- rep(bright_palette, ceiling(n / length(bright_palette)))[seq_len(n)]
  if ("Other" %in% taxa)
    colors[which(taxa == "Other")] <- "#808080"
  names(colors) <- taxa
  return(colors)
}

# Recode metadata: Gut -> Intestine, Soil -> Sediment, trt codes -> labels
process_metadata <- function(metadata_raw) {
  metadata_raw %>%
    dplyr::filter(trial == 1) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        type == "Gut"                    ~ "Intestine",
        type == "Water"                  ~ "Water",
        type %in% c("Soil", "Sediment") ~ "Sediment",
        TRUE ~ as.character(type)
      ),
      type = factor(type, levels = c("Intestine", "Sediment", "Water")),
      trt_label = dplyr::case_when(
        trt == "1" ~ "Basal",
        trt == "2" ~ "CFP 5%",
        trt == "3" ~ "CFP 10%",
        trt == "4" ~ "CFP 20%",
        TRUE ~ as.character(trt)
      ),
      trt_label = factor(trt_label, levels = treatment_order),
      pond      = as.factor(pond)
    )
}

# =========================================================================
# ABUNDANCE CALCULATION
# One rarefaction iteration (type-specific depth) or raw counts for rel abund.
# Returns named list of data frames, one per taxonomic level.
# =========================================================================

calculate_abundance_one_iter <- function(i, otu_matrix, sample_type_map,
                                         rarefaction_depths, tax_df,
                                         metadata_df, use_rarefaction) {
  if (use_rarefaction) {
    set.seed(i * 100)
    rarefied           <- matrix(0L, nrow = nrow(otu_matrix), ncol = ncol(otu_matrix))
    rownames(rarefied) <- rownames(otu_matrix)
    colnames(rarefied) <- colnames(otu_matrix)
    for (idx in seq_len(nrow(otu_matrix))) {
      sname            <- rownames(otu_matrix)[idx]
      depth            <- rarefaction_depths[[ sample_type_map[sname] ]]
      rarefied[idx, ]  <- vegan::rrarefy(otu_matrix[idx, , drop = FALSE], depth)
    }
  } else {
    rarefied <- otu_matrix   # use raw counts; relative abundance computed later
  }

  results <- list()
  for (tax_level in c("Phylum", "Class", "Order", "Family", "Genus")) {
    lvl      <- tolower(tax_level)
    taxa_rows <- lapply(seq_len(nrow(rarefied)), function(r) {
      sample_id <- rownames(rarefied)[r]
      nz        <- which(rarefied[r, ] > 0)
      if (length(nz) == 0) return(NULL)
      data.frame(
        sample = sample_id,
        taxon  = tax_df[colnames(rarefied)[nz], tax_level],
        count  = rarefied[r, nz],
        stringsAsFactors = FALSE
      )
    })
    taxa_df_lvl <- dplyr::bind_rows(taxa_rows)
    if (nrow(taxa_df_lvl) == 0) { results[[lvl]] <- NULL; next }

    agg <- aggregate(count ~ sample + taxon, data = taxa_df_lvl, sum)
    res <- dplyr::left_join(agg, metadata_df, by = "sample")
    colnames(res)[colnames(res) == "taxon"] <- lvl
    results[[lvl]] <- res
  }
  return(results)
}

# =========================================================================
# COMBINE ITERATIONS
# Average counts across rarefaction runs and compute relative abundance.
# =========================================================================

combine_iterations <- function(iter_list, use_rarefaction) {
  combined <- list()

  for (lvl in c("phylum", "class", "order", "family", "genus")) {
    cat(sprintf("  Combining %s...\n", lvl))

    level_data <- dplyr::bind_rows(
      Filter(Negate(is.null), lapply(iter_list, `[[`, lvl))
    )
    if (nrow(level_data) == 0) next

    combined[[lvl]] <- level_data %>%
      dplyr::group_by(sample, !!dplyr::sym(lvl), type, trt_label, pond) %>%
      dplyr::summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        total     = sum(mean_count, na.rm = TRUE),
        abundance = dplyr::if_else(total > 0, mean_count / total, 0)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        type      = factor(type,      levels = c("Intestine", "Sediment", "Water")),
        trt_label = factor(trt_label, levels = treatment_order)
      )
  }
  return(combined)
}

# =========================================================================
# STACKED BAR PLOT FUNCTION
# x-axis = trt_label (treatment groups)
# Taxa labels centered inside each bar segment when proportion >= label_threshold
# =========================================================================

create_treatment_plot <- function(data, tax_level, sample_type,
                                  min_abundance = 0.01) {
  tax_col  <- tolower(tax_level)
  filtered <- dplyr::filter(data, type == sample_type, abundance > 0)
  if (nrow(filtered) == 0)
    stop(sprintf("No data for %s in %s", tax_level, sample_type))

  # Group to trt_label x taxon proportions; bin rare taxa as "Other"
  plot_data <- filtered %>%
    dplyr::group_by(trt_label, !!dplyr::sym(tax_col)) %>%
    dplyr::summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(trt_label) %>%
    dplyr::mutate(proportion = total / sum(total, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!dplyr::sym(tax_col)) %>%
    dplyr::mutate(max_prop = max(proportion, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      taxon = dplyr::if_else(max_prop >= min_abundance,
                             !!dplyr::sym(tax_col), "Other"),
      taxon = gsub("_", " ", taxon)
    ) %>%
    dplyr::group_by(trt_label, taxon) %>%
    dplyr::summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(trt_label = factor(trt_label, levels = treatment_order))

  colors <- assign_colors(plot_data)

  ggplot(plot_data, aes(x = trt_label, y = proportion, fill = taxon)) +
    geom_bar(stat = "identity", width = 0.7, position = "stack",
             color = "white", linewidth = 0.1) +
    geom_text(
      aes(label = dplyr::if_else(proportion >= label_threshold, taxon, "")),
      position = position_stack(vjust = 0.5),
      size     = 3.5,
      fontface = "plain",
      color    = "black"
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = c(0.01, 0.01)) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    labs(
      title = paste0(sample_type, " - ", tools::toTitleCase(tax_level)),
      x     = "Treatment",
      y     = "Relative Abundance",
      fill  = tools::toTitleCase(tax_level)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(face = "bold"),
      axis.text.x        = element_text(color = "black", face = "bold"),
      axis.text.y        = element_text(color = "black"),
      legend.title       = element_text(face = "bold"),
      legend.position    = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3),
      plot.background    = element_rect(fill = "transparent", color = NA),
      panel.background   = element_rect(fill = "transparent", color = NA),
      plot.margin        = ggplot2::margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

cat(strrep("=", 70), "\n")
cat("STARTING ALLUVIAL PLOT ANALYSIS\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Mode: %s\n",
            if (use_rarefaction) "Iterative rarefaction" else "Relative abundance"))

perform_processing <- TRUE

if (file.exists(cache_file)) {
  cat(sprintf("Found cached results: %s\n", basename(cache_file)))
  combined_results <- readRDS(cache_file)
  if (is.list(combined_results) && length(combined_results) > 0) {
    cat("Cached results valid - skipping processing\n")
    perform_processing <- FALSE
  } else {
    cat("Cached results invalid - reprocessing\n")
  }
}

if (perform_processing) {

  cat("STEP 1: Loading and preparing data\n")
  physeq <- readRDS(phyloseq_path)
  cat(sprintf("Loaded: %d samples, %d ASVs\n",
              phyloseq::nsamples(physeq), phyloseq::ntaxa(physeq)))

  # Round OTU table to integers (prevents rrarefy warnings from DADA2 floats)
  otu_mat <- phyloseq::otu_table(physeq) %>% as("matrix")
  if (phyloseq::taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  otu_mat <- round(otu_mat)

  # Process metadata
  meta_raw   <- phyloseq::sample_data(physeq) %>%
    tibble::as_tibble(rownames = "sample") %>%
    janitor::clean_names()
  meta_clean <- process_metadata(meta_raw)

  cat(sprintf("Metadata: %d samples | Types: %s | Treatments: %s\n",
              nrow(meta_clean),
              paste(levels(meta_clean$type), collapse = ", "),
              paste(levels(meta_clean$trt_label), collapse = ", ")))

  # Update phyloseq sample_data so taxonomy lookups see group columns
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(
    meta_clean %>% tibble::column_to_rownames("sample")
  )

  # Named vector: sample_id -> type label
  sample_type_map <- setNames(
    as.character(meta_clean$type),
    meta_clean$sample
  )

  # Filter samples with sufficient reads for their type-specific depth
  if (use_rarefaction) {
    samples_ok <- meta_clean %>%
      dplyr::filter(
        sample %in% rownames(otu_mat),
        mapply(
          function(s, tp) rowSums(otu_mat[s, , drop = FALSE]) >= rarefaction_depths[[tp]],
          sample, as.character(type)
        )
      ) %>%
      dplyr::pull(sample)
  } else {
    samples_ok <- intersect(meta_clean$sample, rownames(otu_mat))
  }

  otu_mat_filtered     <- otu_mat[samples_ok, ]
  sample_type_filtered <- sample_type_map[samples_ok]
  meta_filtered        <- dplyr::filter(meta_clean, sample %in% samples_ok)

  cat(sprintf("Samples passing filter: %d of %d\n",
              length(samples_ok), nrow(meta_clean)))

  if (use_rarefaction) {
    type_counts <- table(sample_type_filtered)
    for (tp in names(rarefaction_depths)) {
      n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
      cat(sprintf("  %s: %d (depth >= %s)\n",
                  tp, n, scales::comma(rarefaction_depths[[tp]])))
    }
  }

  # Build taxonomy lookup (rows = ASV IDs)
  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  for (col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
    tax_df[[col]] <- clean_taxonomy(tax_df[[col]])

  # Metadata subset with only columns needed for grouping
  meta_subset <- dplyr::select(meta_filtered, sample, type, trt_label, pond)

  # =========================================================================
  # STEP 2: RUN PROCESSING
  # =========================================================================

  if (use_rarefaction) {
    cat(sprintf("STEP 2: Iterative rarefaction (%d iterations)\n", n_iterations))
    n_cores <- parallel::detectCores() - 1
    future::plan(future::multisession, workers = n_cores)
    cat(sprintf("Running in parallel (%d cores)...\n", n_cores))

    iter_list <- tryCatch(
      future.apply::future_lapply(
        seq_len(n_iterations),
        function(i) tryCatch(
          calculate_abundance_one_iter(
            i, otu_mat_filtered, sample_type_filtered,
            rarefaction_depths, tax_df, meta_subset, TRUE
          ),
          error = function(e) NULL
        ),
        future.seed = TRUE
      ),
      error = function(e) {
        cat(sprintf("Parallel failed (%s), falling back to sequential\n", e$message))
        future::plan(future::sequential)
        lapply(seq_len(n_iterations), function(i)
          tryCatch(
            calculate_abundance_one_iter(
              i, otu_mat_filtered, sample_type_filtered,
              rarefaction_depths, tax_df, meta_subset, TRUE
            ),
            error = function(e2) NULL
          ))
      }
    )
    future::plan(future::sequential)

    valid <- sum(!sapply(iter_list, is.null))
    cat(sprintf("Completed %d of %d iterations\n", valid, n_iterations))

  } else {
    cat("STEP 2: Calculating relative abundance (no rarefaction)\n")
    iter_list <- list(
      calculate_abundance_one_iter(
        1, otu_mat_filtered, sample_type_filtered,
        rarefaction_depths, tax_df, meta_subset, FALSE
      )
    )
  }

  # =========================================================================
  # STEP 3: COMBINE RESULTS
  # =========================================================================

  cat("STEP 3: Combining results\n")
  combined_results <- combine_iterations(
    Filter(Negate(is.null), iter_list),
    use_rarefaction
  )

  saveRDS(combined_results, cache_file)
  cat(sprintf("Results cached: %s\n", basename(cache_file)))

} # end perform_processing

# =========================================================================
# STEP 4: CREATE PLOTS
# =========================================================================

cat("STEP 4: Creating stacked bar plots\n")

treatment_plots <- list()

for (lvl in c("phylum", "class", "order", "family", "genus")) {
  if (is.null(combined_results[[lvl]])) next
  cat(sprintf("  %s\n", lvl))
  treatment_plots[[lvl]] <- list()
  min_abund <- taxa_thresholds[[lvl]]

  for (st in c("Intestine", "Sediment", "Water")) {
    treatment_plots[[lvl]][[st]] <- tryCatch(
      create_treatment_plot(combined_results[[lvl]],
                            tools::toTitleCase(lvl), st, min_abund),
      error = function(e) {
        cat(sprintf("    Error %s: %s\n", st, e$message))
        NULL
      }
    )
  }
}

# =========================================================================
# STEP 5: COMBINED FIGURES — one per taxonomic level
# Three rows: Intestine / Sediment / Water
# =========================================================================

cat("STEP 5: Assembling combined figures\n")

depth_str <- if (use_rarefaction) {
  paste(
    sapply(names(rarefaction_depths), function(tp) {
      sprintf("%s: %s", tp, scales::comma(rarefaction_depths[[tp]]))
    }),
    collapse = "; "
  )
} else NULL

for (lvl in names(treatment_plots)) {
  plots_ok <- !sapply(treatment_plots[[lvl]][c("Intestine", "Sediment", "Water")],
                      is.null)
  if (!all(plots_ok)) {
    cat(sprintf("  Skipping %s combined - missing: %s\n",
                lvl,
                paste(c("Intestine", "Sediment", "Water")[!plots_ok], collapse = ", ")))
    next
  }

  subtitle_text <- if (use_rarefaction) {
    paste0(
      "Based on 16S rRNA gene amplicon sequencing data analyzing ", lvl, "-level ",
      "diversity across Intestine, Sediment, and Water compartments. ",
      "Rarefied to type-specific depths (", depth_str, " reads) with ",
      n_iterations, " iterations. ",
      "Taxa representing less than ", taxa_thresholds[[lvl]] * 100,
      "% relative abundance are grouped as 'Other'."
    )
  } else {
    paste0(
      "Based on 16S rRNA gene amplicon sequencing data analyzing ", lvl, "-level ",
      "diversity across Intestine, Sediment, and Water compartments. ",
      "Relative abundance normalized using total sum scaling. ",
      "Taxa representing less than ", taxa_thresholds[[lvl]] * 100,
      "% relative abundance are grouped as 'Other'."
    )
  }

  combined_plot <- (
    treatment_plots[[lvl]][["Intestine"]] /
    treatment_plots[[lvl]][["Sediment"]]  /
    treatment_plots[[lvl]][["Water"]]
  ) +
    patchwork::plot_annotation(
      title    = paste0(
        "Alluvial plots ", lvl, " composition across Litopenaeus vannamei",
        " fed various levels of corn fermented protein",
        " in semi-intensive, green water production ponds"
      ),
      subtitle   = subtitle_text,
      tag_levels = "A"
    ) &
    theme(
      plot.title    = ggtext::element_markdown(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggtext::element_markdown(size = 10, hjust = 0.5,
                                               lineheight = 1.2),
      plot.tag      = element_text(size = 12, face = "bold"),
      plot.background   = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    )

  fname <- file.path(output_dir, sprintf(
    "Alluvial_%s_Combined_%s.svg",
    tools::toTitleCase(lvl), method_label
  ))
  ggplot2::ggsave(fname, plot = combined_plot,
                  width = 25, height = 35, units = "in",
                  dpi = 300, bg = "transparent")
  cat(sprintf("  Saved: %s\n", basename(fname)))
}

# =========================================================================
# ANALYSIS SUMMARY
# =========================================================================

cat(strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\nOUTPUTS:\n")
cat(strrep("-", 70), "\n")
cat("FIGURES:\n")
for (lvl in names(treatment_plots)) {
  cat(sprintf("  Alluvial_%s_Combined_%s.svg\n",
              tools::toTitleCase(lvl), method_label))
}

cat("\nCACHE:\n")
cat(sprintf("  %s\n", basename(cache_file)))

cat("\nOUTPUT LOCATION:\n")
cat(sprintf("  %s\n", output_dir))

cat("\nMODE:\n")
cat(sprintf("  %s\n",
            if (use_rarefaction) paste("Rarefaction -", n_iterations, "iterations")
            else "Relative abundance (no rarefaction)"))

if (use_rarefaction) {
  cat("\nRAREFACTION DEPTHS:\n")
  for (tp in names(rarefaction_depths))
    cat(sprintf("  %s: %s\n", tp, scales::comma(rarefaction_depths[[tp]])))
}

cat(strrep("=", 70), "\n")
cat("ALL ANALYSES COMPLETE\n")
cat(strrep("=", 70), "\n")

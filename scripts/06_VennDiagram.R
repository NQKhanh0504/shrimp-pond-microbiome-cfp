# ============================================================
# Script: 06_VennDiagram.R
# Description: Venn diagram analysis of shared and unique
#              genera among diet treatments (Basal, CFP 5%,
#              CFP 10%, CFP 20%) within each sample compartment
#              (Intestine, Sediment, Water) in Litopenaeus
#              vannamei production ponds fed corn fermented
#              protein. Analysis uses iterative rarefaction
#              with type-specific depths and consensus threshold.
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
  "ggVennDiagram", "janitor", "future.apply", "parallel", "future", "tidyverse"
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
key_packages <- c("phyloseq", "vegan", "ggplot2", "ggVennDiagram", "dplyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# Caching flag — set TRUE to regenerate plots from saved RDS without re-running rarefaction
skip_analysis <- FALSE

# File paths
phyloseq_path <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Exp2434_CFP_Microbiome/16S_V3V4_MicrobiomePipeline/step4_phyloseq_object.rds"
output_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Github/Exp2434/VennDiagrams"
cache_file    <- file.path(output_dir, "venn_analysis_results.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
n_iterations        <- 100
consensus_threshold <- 0.2   # Taxon retained if present in >= this fraction of iterations
set.seed(54)

# Type-specific rarefaction depths
rarefaction_depths <- list(
  Intestine = 59815,
  Sediment  = 76637,
  Water     = 85162
)

# Treatment order — controls set order in every Venn diagram
treatment_order <- c("Basal", "CFP 5%", "CFP 10%", "CFP 20%")

# Taxonomic level
tax_level <- "Genus"

# Sample filtering
filter_ponds     <- c("C3", "C6", "D1", "D4")
filter_shrimp_max <- 5

# Per-type fill gradient for Venn regions (light to dark of that type's accent color)
type_gradients <- list(
  Intestine = c("#FDEEEE", "#F8BEBE", "#F28E8E", "#EC5E5E", "#E84646"),
  Sediment  = c("#FDF6E6", "#FAE5B3", "#F4CC66", "#F0BB33", "#E69F00"),
  Water     = c("#E6F2FF", "#B3D9FF", "#66B2FF", "#3399FF", "#0072B2")
)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy <- function(tax_names) {
  tax_names <- gsub("^[a-z]__", "", tax_names)
  tax_names <- gsub("_.*", "", tax_names)
  tax_names[is.na(tax_names) | tax_names == ""] <- "Unidentified"
  return(tax_names)
}

process_sample_metadata <- function(physeq) {
  sample_metadata <- tibble::as_tibble(sample_data(physeq), rownames = "sample_id") %>%
    janitor::clean_names() %>%
    dplyr::filter(trial == "1")

  sample_metadata %>%
    dplyr::mutate(
      group_type = dplyr::case_when(
        grepl("gut|intestin", type, ignore.case = TRUE) ~ "Intestine",
        grepl("water",        type, ignore.case = TRUE) ~ "Water",
        grepl("soil|sediment",type, ignore.case = TRUE) ~ "Sediment",
        TRUE ~ as.character(type)
      ),
      group_trt = dplyr::case_when(
        trt == "1" ~ "Basal",
        trt == "2" ~ "CFP 5%",
        trt == "3" ~ "CFP 10%",
        trt == "4" ~ "CFP 20%",
        TRUE ~ as.character(trt)
      ),
      group_trt = factor(group_trt, levels = treatment_order)
    ) %>%
    dplyr::filter(
      pond %in% filter_ponds,
      is.na(shrimp) | shrimp <= filter_shrimp_max
    )
}

get_taxa_for_samples <- function(sample_ids, rarefied, tax_data) {
  if (length(sample_ids) == 0) return(character(0))

  if (length(sample_ids) == 1) {
    present_taxa <- rarefied[sample_ids, ] > 0
  } else {
    present_taxa <- colSums(rarefied[sample_ids, , drop = FALSE] > 0) > 0
  }

  present_taxa_ids <- names(present_taxa)[present_taxa]
  matching_taxa    <- tax_data$taxa_id %in% present_taxa_ids
  tax_names        <- tax_data[[tax_level]][matching_taxa]
  tax_names        <- tax_names[!is.na(tax_names)]
  return(unique(tax_names))
}

# =========================================================================
# RAREFACTION FUNCTION
# =========================================================================

rarefaction_for_venn <- function(iteration, otu_matrix, sample_type_map,
                                 physeq_obj, sample_types, treatments) {
  set.seed(iteration * 100)

  # Type-specific rarefaction: each sample rarefied to its compartment depth
  rarefied           <- matrix(0L, nrow = nrow(otu_matrix), ncol = ncol(otu_matrix))
  rownames(rarefied) <- rownames(otu_matrix)
  colnames(rarefied) <- colnames(otu_matrix)
  for (idx in seq_len(nrow(otu_matrix))) {
    sname <- rownames(otu_matrix)[idx]
    depth <- rarefaction_depths[[ sample_type_map[sname] ]]
    rarefied[idx, ] <- vegan::rrarefy(otu_matrix[idx, , drop = FALSE], depth)
  }

  tax_table_data         <- as.data.frame(tax_table(physeq_obj))
  tax_table_data$taxa_id <- rownames(tax_table_data)
  for (col in setdiff(names(tax_table_data), "taxa_id")) {
    tax_table_data[[col]] <- clean_taxonomy(tax_table_data[[col]])
  }

  sample_data_df           <- as.data.frame(sample_data(physeq_obj))
  sample_data_df$sample_id <- rownames(sample_data_df)

  # Build: for each sample type, list of taxa per treatment
  results <- list()
  for (st in sample_types) {
    taxa_by_trt <- list()
    for (trt in treatments) {
      ids <- sample_data_df$sample_id[
        sample_data_df$group_type == st & as.character(sample_data_df$group_trt) == trt
      ]
      ids <- intersect(ids, rownames(rarefied))
      if (length(ids) > 0) {
        taxa_by_trt[[trt]] <- get_taxa_for_samples(ids, rarefied, tax_table_data)
      }
    }
    taxa_by_trt <- taxa_by_trt[lengths(taxa_by_trt) > 0]
    if (length(taxa_by_trt) >= 2) {
      results[[st]] <- taxa_by_trt[intersect(treatments, names(taxa_by_trt))]
    }
  }

  list(by_sampletype = results)
}

# =========================================================================
# COMBINE RAREFACTION RESULTS
# =========================================================================

combine_rarefaction_results <- function(results_list) {
  all_results  <- lapply(results_list, `[[`, "by_sampletype")
  all_keys     <- unique(unlist(lapply(all_results, names)))

  combined_results <- list()
  for (st in all_keys) {
    st_results       <- Filter(Negate(is.null), lapply(all_results, `[[`, st))
    if (length(st_results) == 0) next
    treatments_found <- unique(unlist(lapply(st_results, names)))
    consensus        <- list()

    for (trt in treatments_found) {
      all_taxa     <- unlist(lapply(st_results, `[[`, trt))
      taxon_counts <- table(all_taxa)
      kept         <- names(taxon_counts)[taxon_counts >= (n_iterations * consensus_threshold)]
      if (length(kept) > 0) consensus[[trt]] <- kept
    }

    if (length(consensus) >= 2) {
      combined_results[[st]] <- consensus[intersect(treatment_order, names(consensus))]
    }
  }

  list(by_sampletype = combined_results)
}

# =========================================================================
# PARALLEL PROCESSING SETUP
# =========================================================================

setup_parallel_processing <- function() {
  tryCatch({
    n_cores <- max(1, parallel::detectCores() - 1)
    future::plan(future::multisession, workers = n_cores)
    return(TRUE)
  }, error = function(e) {
    future::plan(future::sequential)
    return(FALSE)
  })
}

# =========================================================================
# PLOTTING FUNCTIONS
# =========================================================================

get_fill_limits <- function(results_by_type) {
  all_counts <- unlist(lapply(results_by_type, function(grp) {
    if (length(grp) < 2) return(NULL)
    venn <- ggVennDiagram::process_data(ggVennDiagram::Venn(grp))
    ggVennDiagram::venn_region(venn)$count
  }))
  if (length(all_counts) == 0) return(NULL)
  c(min(all_counts, na.rm = TRUE), max(all_counts, na.rm = TRUE))
}

create_venn_plot <- function(data_list, title, colors = NULL, fill_limits = NULL) {
  if (length(data_list) < 2) return(NULL)

  if (is.null(colors)) {
    colors <- c("#EDF8E9", "#BAE4B3", "#74C476", "#238B45")
  }

  p <- ggVennDiagram::ggVennDiagram(
    data_list,
    label               = "both",
    label_alpha         = 1,
    label_color         = "black",
    label_geom          = "text",
    label_percent_digit = 2,
    label_size          = 4,
    set_size            = 4,
    edge_lty            = "solid",
    edge_size           = 0.3,
    show_percentage     = TRUE,
    show_quantity       = TRUE
  ) +
    scale_fill_gradientn(colors = colors, limits = fill_limits) +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.ticks      = element_blank(),
      plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 8),
      axis.title      = element_blank(),
      axis.text       = element_blank()
    )

  return(p)
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

cat("\n")
cat("STARTING VENN DIAGRAM ANALYSIS\n")
cat("\n")

if (!skip_analysis) {

  cat("STEP 1: Loading and preparing data\n")
  physeq <- readRDS(phyloseq_path)
  cat(sprintf("Loaded: %d samples, %d ASVs",
              phyloseq::nsamples(physeq), phyloseq::ntaxa(physeq)), "\n")

  # Round OTU table to integers
  otu_mat_raw <- otu_table(physeq) %>% as("matrix")
  if (taxa_are_rows(physeq)) otu_mat_raw <- t(otu_mat_raw)
  otu_mat_raw <- round(otu_mat_raw)

  # Process metadata and update phyloseq object
  meta_clean <- process_sample_metadata(physeq)
  cat(sprintf(
    "Metadata: %d samples | Types: %s | Treatments: %s",
    nrow(meta_clean),
    paste(sort(unique(meta_clean$group_type)), collapse = ", "),
    paste(levels(meta_clean$group_trt), collapse = ", ")
  ), "\n")

  # Update sample_data on physeq so downstream functions see group_type/group_trt
  sample_data(physeq) <- sample_data(
    meta_clean %>% tibble::column_to_rownames("sample_id")
  )

  # Named character vector: sample_id -> type label (for type-specific rarefaction)
  sample_type_map <- setNames(
    as.character(meta_clean$group_type),
    meta_clean$sample_id
  )

  # Filter samples with sufficient reads for their type-specific depth
  samples_ok <- meta_clean %>%
    dplyr::filter(
      sample_id %in% rownames(otu_mat_raw),
      mapply(
        function(s, tp) rowSums(otu_mat_raw[s, , drop = FALSE]) >= rarefaction_depths[[tp]],
        sample_id, as.character(group_type)
      )
    ) %>%
    dplyr::pull(sample_id)

  otu_mat_filtered     <- otu_mat_raw[samples_ok, ]
  sample_type_filtered <- sample_type_map[samples_ok]

  cat(sprintf("Samples passing depth filter: %d of %d",
              length(samples_ok), nrow(meta_clean)), "\n")
  type_counts <- table(sample_type_filtered)
  for (tp in names(rarefaction_depths)) {
    n <- if (tp %in% names(type_counts)) type_counts[[tp]] else 0L
    cat(sprintf("  %s: %d (depth >= %d)", tp, n, rarefaction_depths[[tp]]), "\n")
  }

  sample_types <- sort(unique(meta_clean$group_type))
  treatments   <- treatment_order

  # =========================================================================
  # STEP 2: RAREFACTION SETUP AND RUN
  # =========================================================================

  cat("STEP 2: Running rarefaction analysis\n")

  parallel_available <- setup_parallel_processing()

  if (parallel_available) {
    rarefaction_results <- tryCatch({
      future.apply::future_lapply(
        seq_len(n_iterations),
        function(i) rarefaction_for_venn(
          i, otu_mat_filtered, sample_type_filtered,
          physeq, sample_types, treatments
        ),
        future.seed = TRUE
      )
    }, error = function(e) {
      cat(sprintf("Parallel failed (%s), falling back to sequential", e$message), "\n")
      future::plan(future::sequential)
      lapply(seq_len(n_iterations), function(i) {
        rarefaction_for_venn(
          i, otu_mat_filtered, sample_type_filtered,
          physeq, sample_types, treatments
        )
      })
    })
  } else {
    rarefaction_results <- lapply(seq_len(n_iterations), function(i) {
      rarefaction_for_venn(
        i, otu_mat_filtered, sample_type_filtered,
        physeq, sample_types, treatments
      )
    })
  }

  future::plan(future::sequential)

  cat(sprintf("Combining %d iterations", n_iterations), "\n")
  combined_results <- combine_rarefaction_results(rarefaction_results)

  # Cache results
  saveRDS(
    list(combined_results = combined_results,
         sample_types     = sample_types,
         treatments       = treatments),
    cache_file
  )
  cat(sprintf("Results cached: %s", cache_file), "\n")

} else {

  cat("STEP 1: Loading cached results\n")
  if (!file.exists(cache_file)) {
    stop("Cache file not found. Run with skip_analysis = FALSE first.")
  }
  cache            <- readRDS(cache_file)
  combined_results <- cache$combined_results
  sample_types     <- cache$sample_types
  treatments       <- cache$treatments
  cat(sprintf("Loaded: %d sample types", length(sample_types)), "\n")

}

# =========================================================================
# STEP 3: CREATE PLOTS
# =========================================================================

cat("STEP 3: Creating Venn diagram figures\n")

results_st  <- combined_results$by_sampletype
fill_limits <- get_fill_limits(results_st)

panels <- list()
for (st in c("Intestine", "Sediment", "Water")) {
  if (!st %in% names(results_st)) next
  grp <- results_st[[st]]
  if (length(grp) < 2) next

  colors <- if (st %in% names(type_gradients)) type_gradients[[st]] else
    c("#D1D5DB", "#9CA3AF", "#6B7280", "#374151")

  panel_title <- sprintf("%s — %s", st, tax_level)

  panels[[st]] <- create_venn_plot(
    grp[intersect(treatment_order, names(grp))],
    title       = panel_title,
    colors      = colors,
    fill_limits = fill_limits
  )
}

# =========================================================================
# STEP 4: COMBINED FIGURE
# =========================================================================

cat("STEP 4: Assembling and saving combined figure\n")

if (length(panels) >= 2) {
  ordered_panels <- panels[intersect(c("Intestine", "Sediment", "Water"), names(panels))]

  trt_str   <- paste(treatments, collapse = ", ")
  depth_str <- paste(
    sapply(names(rarefaction_depths), function(tp) {
      sprintf("%s: %s", tp, scales::comma(rarefaction_depths[[tp]]))
    }),
    collapse = "; "
  )

  fig_combined <- patchwork::wrap_plots(ordered_panels, ncol = length(ordered_panels)) +
    patchwork::plot_annotation(
      title    = sprintf(
        "Shared and unique %s among diet treatments within each compartment",
        tolower(tax_level)
      ),
      subtitle = sprintf(
        paste0(
          "Based on 16S rRNA gene amplicon sequencing data analyzing %s level diversity ",
          "between diet groups (%s) within each compartment. ",
          "Rarefied to compartment-specific depths (%s reads) with %d iterations and ",
          "%g%% consensus threshold. Panels A-C: Intestine, Sediment, Water."
        ),
        tolower(tax_level), trt_str, depth_str,
        n_iterations, consensus_threshold * 100
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5),
      plot.tag      = element_text(size = 10, face = "bold")
    )

  out_fname <- file.path(
    output_dir,
    sprintf("VennDiagram_%s_TreatmentBySampleType.png", tax_level)
  )

  ggsave(
    filename = out_fname,
    plot     = fig_combined,
    width    = 24,
    height   = 8,
    units    = "in",
    dpi      = 300
  )

  cat(sprintf("Saved: %s", basename(out_fname)), "\n")
}

# =========================================================================
# STEP 5: ANALYSIS SUMMARY
# =========================================================================

cat("\n")
cat("ANALYSIS COMPLETE — SUMMARY\n")
cat("\n")

cat("\nOUTPUTS:\n")
cat(strrep("-", 70), "\n")
cat("FIGURES:\n")
cat(sprintf("  VennDiagram_%s_TreatmentBySampleType.png\n", tax_level))
cat("\nCACHE:\n")
cat(sprintf("  %s\n", cache_file))
cat("\nOUTPUT LOCATION:\n")
cat(sprintf("  %s\n", output_dir))

cat("\n")
cat("ALL ANALYSES COMPLETE\n")
cat("\n")

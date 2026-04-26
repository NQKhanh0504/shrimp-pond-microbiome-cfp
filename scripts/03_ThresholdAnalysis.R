# ============================================================
# Script: Threshold_Analysis.R
# Description: Abundance and prevalence filtering threshold
#              analysis for CFP shrimp microbiome data.
#              Tests combined abundance + prevalence filtering
#              and monitors critical taxa impacts across
#              multiple taxonomic levels.
# Input: phyloseq object (seqtab_nochim.rds)
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

required_packages <- c("phyloseq", "stringr", "dplyr", "tidyr", "readr",
                       "parallel")

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
key_packages <- c("phyloseq", "stringr", "dplyr", "readr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

save_detailed_csv <- TRUE

# File paths
base_dir      <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab"
save_dir      <- file.path(base_dir, "Exp2434_CFP_Microbiome/Threshold_Analysis")
gg2_ref       <- file.path(save_dir, "gg2_2024_09_toGenus_trainset.fa.gz")
phyloseq_file <- file.path(save_dir, "seqtab_nochim.rds")

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
csv_dir <- file.path(save_dir, "csv_exports")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# Filtering parameters
threshold_range       <- seq(5, 100, by = 5)
abundance_thresholds  <- threshold_range
prevalence_thresholds <- seq(0, 20, by = 1)

# Analysis scope
analyze_global_prevalence <- FALSE
analyze_group_prevalence  <- TRUE

# Display and export settings
max_displayed_taxa     <- 10
export_sample_specific <- TRUE

# Taxonomic analysis settings
taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
analyze_levels   <- c("Phylum", "Class", "Order", "Family", "Genus")
primary_level    <- "Genus"
secondary_level  <- "Family"

# Phylogenetic annotation handling
aggregate_phylogenetic_variants <- TRUE
phylo_suffix_pattern <- "(_[A-Z0-9]+|_UCG-[0-9]+|_NK[0-9A-Z]+|_unclassified|_incertae_sedis)$"

# Critical taxa for monitoring
critical_taxa <- list(
  Genus = list(
    list(name = "Vibrio",         patterns = c("Vibrio"),         exact_match = FALSE),
    list(name = "Aeromonas",      patterns = c("Aeromonas"),      exact_match = FALSE),
    list(name = "Pseudomonas",    patterns = c("Pseudomonas"),    exact_match = FALSE),
    list(name = "Lactobacillus",  patterns = c("Lactobacillus"),  exact_match = FALSE),
    list(name = "Bacillus",       patterns = c("Bacillus"),       exact_match = FALSE),
    list(name = "Enterococcus",   patterns = c("Enterococcus"),   exact_match = FALSE),
    list(name = "Shewanella",     patterns = c("Shewanella"),     exact_match = FALSE),
    list(name = "Flavobacterium", patterns = c("Flavobacterium"), exact_match = FALSE)
  ),
  Family = list(
    list(name = "Vibrionaceae",      patterns = c("Vibrionaceae"),      exact_match = FALSE),
    list(name = "Lactobacillaceae",  patterns = c("Lactobacillaceae"),  exact_match = FALSE),
    list(name = "Bacillaceae",       patterns = c("Bacillaceae"),       exact_match = FALSE),
    list(name = "Rhodobacteraceae",  patterns = c("Rhodobacteraceae"),  exact_match = FALSE),
    list(name = "Flavobacteriaceae", patterns = c("Flavobacteriaceae"), exact_match = FALSE)
  ),
  Order = list(
    list(name = "Vibrionales",     patterns = c("Vibrionales"),     exact_match = FALSE),
    list(name = "Lactobacillales", patterns = c("Lactobacillales"), exact_match = FALSE),
    list(name = "Bacillales",      patterns = c("Bacillales"),      exact_match = FALSE),
    list(name = "Rhodobacterales", patterns = c("Rhodobacterales"), exact_match = FALSE),
    list(name = "Pseudomonadales", patterns = c("Pseudomonadales"), exact_match = FALSE)
  )
)

# Critical taxa filtering thresholds
critical_abundance_threshold  <- 10
critical_prevalence_threshold <- 0.02
use_critical_filters          <- TRUE

# Taxonomic cleaning settings
remove_prefixes     <- TRUE
handle_unclassified <- TRUE

set.seed(54)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

get_sample_type_column <- function(metadata) {
  potential_cols <- c("type", "sample_type", "sampletype", "group", "category")
  for (col in potential_cols) {
    if (col %in% colnames(metadata)) return(col)
  }
  for (col in colnames(metadata)) {
    if (is.character(metadata[[col]]) || is.factor(metadata[[col]])) {
      n_levels <- length(unique(metadata[[col]]))
      if (n_levels > 1 && n_levels < nrow(metadata) / 2) return(col)
    }
  }
  return(NULL)
}

clean_taxonomic_names <- function(taxa_vector, level_name, preserve_length = FALSE) {
  if (length(taxa_vector) == 0) return(character(0))

  taxa_clean  <- taxa_vector
  na_positions <- is.na(taxa_clean)

  if (remove_prefixes) {
    level_prefixes <- list(
      Kingdom = "^[dk]__", Phylum = "^p__", Class  = "^c__",
      Order   = "^o__",    Family = "^f__", Genus  = "^g__", Species = "^s__"
    )
    level_key <- names(level_prefixes)[
      tolower(names(level_prefixes)) == tolower(level_name)
    ]
    if (length(level_key) > 0) {
      taxa_clean <- stringr::str_remove(taxa_clean, level_prefixes[[level_key[1]]])
    }
  }

  unclassified_pos <- taxa_clean == "" | is.na(taxa_clean) | na_positions
  if (handle_unclassified || preserve_length) {
    taxa_clean[unclassified_pos] <- paste0("Unclassified_", level_name)
  } else if (!preserve_length) {
    taxa_clean <- taxa_clean[!unclassified_pos]
  }

  return(taxa_clean)
}

aggregate_phylogenetic_names <- function(taxa_vector) {
  if (!aggregate_phylogenetic_variants || length(taxa_vector) == 0) return(taxa_vector)
  stringr::str_remove(taxa_vector, phylo_suffix_pattern)
}

match_critical_taxon <- function(taxa_names, pattern, exact_match = FALSE) {
  if (exact_match) {
    which(stringr::str_detect(taxa_names, paste0("^", pattern, "$")))
  } else if (aggregate_phylogenetic_variants) {
    taxa_agg <- aggregate_phylogenetic_names(taxa_names)
    which(stringr::str_detect(taxa_agg, paste0("^", pattern, "$")))
  } else {
    which(stringr::str_detect(taxa_names, paste0("^", pattern, "($|_)")))
  }
}

check_taxonomy_columns <- function(taxa_matrix, levels_to_check = analyze_levels) {
  available <- character()
  for (level in levels_to_check) {
    if      (level          %in% colnames(taxa_matrix)) available <- c(available, level)
    else if (tolower(level) %in% colnames(taxa_matrix)) available <- c(available, tolower(level))
  }
  if (length(available) == 0) stop("No taxonomic levels found. Expected: ",
                                   paste(levels_to_check, collapse = ", "))
  return(available)
}

calculate_threshold_impacts <- function(asv_counts, thresholds) {
  data.frame(
    Threshold                = thresholds,
    ASVs_Remaining           = sapply(thresholds, function(t) sum(asv_counts >= t)),
    Percentage_Remaining     = sapply(thresholds, function(t)
      round(sum(asv_counts >= t) / length(asv_counts) * 100, 2)),
    Reads_Removed            = sapply(thresholds, function(t)
      sum(asv_counts[asv_counts < t])),
    Percentage_Reads_Removed = sapply(thresholds, function(t)
      round(sum(asv_counts[asv_counts < t]) / sum(asv_counts) * 100, 2))
  )
}

apply_abundance_filter <- function(seqtab, threshold) {
  asv_counts <- colSums(seqtab)
  keep_idx   <- which(asv_counts >= threshold)
  list(
    filtered_seqtab = seqtab[, keep_idx, drop = FALSE],
    removed_indices = which(asv_counts < threshold),
    kept_indices    = keep_idx,
    n_removed       = sum(asv_counts < threshold),
    n_kept          = length(keep_idx)
  )
}

apply_prevalence_filter <- function(seqtab, threshold_percent) {
  n_samples       <- nrow(seqtab)
  threshold_count <- ceiling(n_samples * threshold_percent / 100)
  asv_prev        <- colSums(seqtab > 0)
  keep_idx        <- which(asv_prev >= threshold_count)
  list(
    filtered_seqtab = seqtab[, keep_idx, drop = FALSE],
    removed_indices = which(asv_prev < threshold_count),
    kept_indices    = keep_idx,
    n_removed       = sum(asv_prev < threshold_count),
    n_kept          = length(keep_idx),
    threshold_count = threshold_count
  )
}

apply_combined_filter <- function(seqtab, abundance_thresh, prevalence_thresh_percent) {
  abund_result <- apply_abundance_filter(seqtab, abundance_thresh)

  if (prevalence_thresh_percent == 0) {
    return(list(
      filtered_seqtab       = abund_result$filtered_seqtab,
      removed_indices       = abund_result$removed_indices,
      kept_indices          = abund_result$kept_indices,
      n_removed             = abund_result$n_removed,
      n_kept                = abund_result$n_kept,
      removed_by_abundance  = abund_result$n_removed,
      removed_by_prevalence = 0
    ))
  }

  prev_result   <- apply_prevalence_filter(abund_result$filtered_seqtab, prevalence_thresh_percent)
  original_kept <- abund_result$kept_indices
  final_kept    <- original_kept[prev_result$kept_indices]
  final_removed <- setdiff(seq_len(ncol(seqtab)), final_kept)

  list(
    filtered_seqtab       = prev_result$filtered_seqtab,
    removed_indices       = final_removed,
    kept_indices          = final_kept,
    n_removed             = length(final_removed),
    n_kept                = length(final_kept),
    removed_by_abundance  = abund_result$n_removed,
    removed_by_prevalence = prev_result$n_removed
  )
}

analyze_taxonomic_impact <- function(taxa_matrix, removed_indices, taxonomic_level,
                                     critical_taxa_list = NULL, sample_metadata = NULL,
                                     seqtab = NULL, abundance_threshold = NULL,
                                     prevalence_threshold = NULL) {
  empty_result <- list(
    level = taxonomic_level, total_taxa_removed = 0,
    taxa_counts = table(character(0)), critical_impacts = numeric(0),
    sample_type_impacts = data.frame()
  )

  if (length(removed_indices) == 0) return(empty_result)

  level_col <- if (taxonomic_level %in% colnames(taxa_matrix)) taxonomic_level
  else tolower(taxonomic_level)
  if (!level_col %in% colnames(taxa_matrix)) return(empty_result)

  taxa_removed <- taxa_matrix[removed_indices, level_col]
  taxa_clean   <- clean_taxonomic_names(taxa_removed, taxonomic_level, preserve_length = FALSE)
  if (aggregate_phylogenetic_variants && length(taxa_clean) > 0) {
    taxa_clean <- aggregate_phylogenetic_names(taxa_clean)
  }
  if (length(taxa_clean) == 0) return(empty_result)

  taxa_counts         <- table(taxa_clean)
  critical_impacts    <- numeric(0)
  sample_type_impacts <- data.frame()

  if (!is.null(critical_taxa_list) && taxonomic_level %in% names(critical_taxa_list)) {
    critical_list  <- critical_taxa_list[[taxonomic_level]]
    all_taxa_names <- clean_taxonomic_names(taxa_matrix[, level_col], taxonomic_level,
                                            preserve_length = FALSE)

    for (critical_item in critical_list) {
      taxon_name    <- critical_item$name
      all_taxon_idx <- unique(unlist(lapply(critical_item$patterns, function(pat) {
        match_critical_taxon(all_taxa_names, pat, critical_item$exact_match)
      })))

      if (length(all_taxon_idx) == 0) next

      if (use_critical_filters && !is.null(seqtab)) {
        total_abund <- sum(colSums(seqtab[, all_taxon_idx, drop = FALSE]))
        prevalence  <- sum(rowSums(seqtab[, all_taxon_idx, drop = FALSE] > 0) > 0) / nrow(seqtab)
        abund_pass  <- total_abund >= critical_abundance_threshold
        prev_pass   <- prevalence  >= critical_prevalence_threshold
        if (!abund_pass || !prev_pass) next
      }

      removed_count <- length(intersect(all_taxon_idx, removed_indices))
      if (removed_count > 0) critical_impacts[taxon_name] <- removed_count
    }

    if (length(critical_impacts) > 0 && !is.null(sample_metadata) && !is.null(seqtab)) {
      type_column <- get_sample_type_column(sample_metadata)
      if (!is.null(type_column)) {
        for (taxon_name in names(critical_impacts)) {
          critical_item <- NULL
          for (item in critical_list) {
            if (item$name == taxon_name) { critical_item <- item; break }
          }
          if (is.null(critical_item)) next

          all_taxon_idx <- unique(unlist(lapply(critical_item$patterns, function(pat) {
            match_critical_taxon(all_taxa_names, pat, critical_item$exact_match)
          })))

          for (type in unique(sample_metadata[[type_column]])) {
            removed_in_type <- intersect(all_taxon_idx, removed_indices)
            if (length(removed_in_type) > 0) {
              sample_type_impacts <- rbind(sample_type_impacts, data.frame(
                Taxon       = taxon_name,
                Level       = taxonomic_level,
                Sample_Type = type,
                ASVs_Lost   = length(removed_in_type),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }

  list(
    level               = taxonomic_level,
    total_taxa_removed  = length(taxa_counts),
    taxa_counts         = taxa_counts,
    critical_impacts    = critical_impacts,
    sample_type_impacts = sample_type_impacts
  )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

cat("Loading phyloseq object...\n")
ps_object <- readRDS(phyloseq_file)

if (!inherits(ps_object, "phyloseq")) {
  stop("Loaded object is not a phyloseq object. Found: ", class(ps_object)[1])
}

cat(sprintf("Phyloseq object loaded: %d ASVs, %d samples\n",
            phyloseq::ntaxa(ps_object), phyloseq::nsamples(ps_object)))

seqtab.nochim <- as.matrix(phyloseq::otu_table(ps_object))
if (phyloseq::taxa_are_rows(ps_object)) seqtab.nochim <- t(seqtab.nochim)

if (!is.null(phyloseq::sample_data(ps_object, errorIfNULL = FALSE))) {
  sample_metadata <- data.frame(phyloseq::sample_data(ps_object))
} else {
  sample_metadata <- data.frame(
    Sample_ID = rownames(seqtab.nochim),
    group     = "All_Samples",
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- rownames(seqtab.nochim)
}

if (!is.null(phyloseq::tax_table(ps_object, errorIfNULL = FALSE))) {
  taxa_matrix <- as.matrix(phyloseq::tax_table(ps_object))
} else {
  stop("No taxonomy found in phyloseq object.")
}

available_levels <- check_taxonomy_columns(taxa_matrix, analyze_levels)
cat("Available taxonomic levels:", paste(available_levels, collapse = ", "), "\n")

if (!primary_level %in% available_levels && !tolower(primary_level) %in% available_levels) {
  primary_level <- available_levels[1]
}
if (!secondary_level %in% available_levels && !tolower(secondary_level) %in% available_levels) {
  secondary_level <- if (length(available_levels) > 1) available_levels[2] else primary_level
}

type_column <- get_sample_type_column(sample_metadata)
if (is.null(type_column)) {
  sample_metadata$group    <- "All_Samples"
  type_column              <- "group"
  analyze_group_prevalence <- FALSE
}

sample_groups <- unique(sample_metadata[[type_column]])
total_asvs    <- ncol(seqtab.nochim)
asv_counts    <- colSums(seqtab.nochim)
total_reads   <- sum(asv_counts)
singletons    <- sum(asv_counts == 1)
doubletons    <- sum(asv_counts == 2)

cat(sprintf("Total ASVs: %s\n", format(total_asvs, big.mark = ",")))
cat(sprintf("Total reads: %s\n", format(total_reads, big.mark = ",")))
cat(sprintf("Singletons: %d (%.1f%%)\n", singletons, singletons / total_asvs * 100))
cat(sprintf("Doubletons: %d (%.1f%%)\n", doubletons, doubletons / total_asvs * 100))

# =========================================================================
# STEP 2: BASELINE ABUNDANCE AND PREVALENCE ANALYSIS
# =========================================================================

cat("\nAnalyzing baseline abundance and prevalence by sample group...\n")

asv_prevalence <- colSums(seqtab.nochim > 0)
n_samples      <- nrow(seqtab.nochim)

type_asv_lists <- list()
type_summaries <- list()

for (group in sample_groups) {
  group_samples <- sample_metadata[[type_column]] == group
  group_seqtab  <- seqtab.nochim[group_samples, , drop = FALSE]
  group_asv_cnt <- colSums(group_seqtab)
  group_seqtab  <- group_seqtab[, group_asv_cnt > 0, drop = FALSE]

  n_group    <- nrow(group_seqtab)
  group_asvs <- colnames(group_seqtab)
  group_prev <- colSums(group_seqtab > 0)
  group_abund <- colSums(group_seqtab)

  type_asv_lists[[group]] <- group_asvs
  type_summaries[[group]] <- list(
    n_samples         = n_group,
    n_asvs            = length(group_asvs),
    singleton_pct     = round(sum(group_prev == 1) / length(group_asvs) * 100, 1),
    core_pct          = round(sum(group_prev > n_group * 0.5) / length(group_asvs) * 100, 1),
    median_prevalence = median(group_prev),
    median_abundance  = median(group_abund)
  )

  cat(sprintf("  %s: %d samples, %d ASVs\n", group, n_group, length(group_asvs)))
}

comparison_df <- data.frame(
  Sample_Type       = names(type_summaries),
  N_Samples         = sapply(type_summaries, function(x) x$n_samples),
  N_ASVs            = sapply(type_summaries, function(x) x$n_asvs),
  Singleton_Pct     = sapply(type_summaries, function(x) x$singleton_pct),
  Core_Pct          = sapply(type_summaries, function(x) x$core_pct),
  Median_Prevalence = sapply(type_summaries, function(x) x$median_prevalence),
  Median_Abundance  = sapply(type_summaries, function(x) x$median_abundance)
)

readr::write_csv(comparison_df, file.path(csv_dir, "sample_type_comparison.csv"))

# ASV sharing analysis (3 groups)
if (length(sample_groups) == 3) {
  g1 <- sample_groups[1]; g2 <- sample_groups[2]; g3 <- sample_groups[3]

  sharing_df <- data.frame(
    Category = c(
      paste0(g1, "_only"), paste0(g2, "_only"), paste0(g3, "_only"),
      paste0(g1, "_", g2, "_only"), paste0(g1, "_", g3, "_only"),
      paste0(g2, "_", g3, "_only"), "shared_all"
    ),
    N_ASVs = c(
      length(setdiff(type_asv_lists[[g1]], union(type_asv_lists[[g2]], type_asv_lists[[g3]]))),
      length(setdiff(type_asv_lists[[g2]], union(type_asv_lists[[g1]], type_asv_lists[[g3]]))),
      length(setdiff(type_asv_lists[[g3]], union(type_asv_lists[[g1]], type_asv_lists[[g2]]))),
      length(setdiff(intersect(type_asv_lists[[g1]], type_asv_lists[[g2]]), type_asv_lists[[g3]])),
      length(setdiff(intersect(type_asv_lists[[g1]], type_asv_lists[[g3]]), type_asv_lists[[g2]])),
      length(setdiff(intersect(type_asv_lists[[g2]], type_asv_lists[[g3]]), type_asv_lists[[g1]])),
      length(Reduce(intersect, list(type_asv_lists[[g1]], type_asv_lists[[g2]], type_asv_lists[[g3]])))
    )
  )
  sharing_df$Pct_of_Total <- round(sharing_df$N_ASVs / total_asvs * 100, 2)
  readr::write_csv(sharing_df, file.path(csv_dir, "asv_sharing_analysis.csv"))
}

# Overall threshold impact table
overall_results <- calculate_threshold_impacts(asv_counts, threshold_range)
readr::write_csv(overall_results, file.path(csv_dir, "threshold_summary.csv"))

# Recommended thresholds based on read loss
conservative_threshold <- threshold_range[which(overall_results$Percentage_Reads_Removed <= 1)[1]]
moderate_threshold     <- threshold_range[which(overall_results$Percentage_Reads_Removed <= 5)[1]]
aggressive_threshold   <- threshold_range[which(overall_results$Percentage_Reads_Removed <= 10)[1]]

cat(sprintf("\nRecommended thresholds (read loss basis):\n"))
cat(sprintf("  Conservative (<=1%% reads lost): %s\n",
            if (!is.na(conservative_threshold)) conservative_threshold else "none"))
cat(sprintf("  Moderate (<=5%% reads lost): %s\n",
            if (!is.na(moderate_threshold)) moderate_threshold else "none"))
cat(sprintf("  Aggressive (<=10%% reads lost): %s\n",
            if (!is.na(aggressive_threshold)) aggressive_threshold else "none"))

# =========================================================================
# STEP 3: GLOBAL PREVALENCE ANALYSIS
# =========================================================================

if (analyze_global_prevalence) {
  cat("\nRunning global prevalence analysis...\n")
  global_critical_summary <- data.frame()
  filtering_matrix        <- data.frame()

  for (abundance_thresh in abundance_thresholds) {
    for (prevalence_thresh in prevalence_thresholds) {
      filter_result <- apply_combined_filter(seqtab.nochim, abundance_thresh, prevalence_thresh)

      filtering_matrix <- rbind(filtering_matrix, data.frame(
        Abundance_Threshold   = abundance_thresh,
        Prevalence_Threshold  = prevalence_thresh,
        ASVs_Remaining        = filter_result$n_kept,
        ASVs_Removed          = filter_result$n_removed,
        Removed_by_Abundance  = filter_result$removed_by_abundance,
        Removed_by_Prevalence = filter_result$removed_by_prevalence,
        Percent_Remaining     = round(filter_result$n_kept / total_asvs * 100, 2)
      ))

      for (level in available_levels) {
        taxa_impact <- analyze_taxonomic_impact(
          taxa_matrix          = taxa_matrix,
          removed_indices      = filter_result$removed_indices,
          taxonomic_level      = level,
          critical_taxa_list   = critical_taxa,
          sample_metadata      = sample_metadata,
          seqtab               = seqtab.nochim,
          abundance_threshold  = abundance_thresh,
          prevalence_threshold = prevalence_thresh
        )

        if (length(taxa_impact$critical_impacts) > 0) {
          for (taxon in names(taxa_impact$critical_impacts)) {
            global_critical_summary <- rbind(global_critical_summary, data.frame(
              Scope                = "Global",
              Sample_Group         = "All",
              Abundance_Threshold  = abundance_thresh,
              Prevalence_Threshold = prevalence_thresh,
              Taxonomic_Level      = level,
              Taxon                = taxon,
              ASVs_Lost            = taxa_impact$critical_impacts[taxon]
            ))
          }
        }
      }
    }
  }

  readr::write_csv(filtering_matrix,
                   file.path(csv_dir, "filtering_impact_matrix_global.csv"))

  if (nrow(global_critical_summary) > 0) {
    readr::write_csv(global_critical_summary,
                     file.path(csv_dir, "critical_taxa_global.csv"))
  }

  safe_thresholds <- data.frame()
  for (i in seq_len(nrow(filtering_matrix))) {
    row     <- filtering_matrix[i, ]
    abund_t <- row$Abundance_Threshold
    prev_t  <- row$Prevalence_Threshold

    has_critical <- if (nrow(global_critical_summary) > 0) {
      nrow(global_critical_summary[
        global_critical_summary$Abundance_Threshold  == abund_t &
          global_critical_summary$Prevalence_Threshold == prev_t, ]) > 0
    } else FALSE

    if (!has_critical) {
      safe_thresholds <- rbind(safe_thresholds, data.frame(
        Abundance         = abund_t,
        Prevalence        = prev_t,
        ASVs_Remaining    = row$ASVs_Remaining,
        Percent_Remaining = row$Percent_Remaining
      ))
    }
  }

  if (nrow(safe_thresholds) > 0) {
    safe_thresholds <- safe_thresholds[order(safe_thresholds$ASVs_Remaining), ]
    readr::write_csv(safe_thresholds, file.path(csv_dir, "threshold_combinations.csv"))
  }
}

# =========================================================================
# STEP 4: GROUP-SPECIFIC PREVALENCE ANALYSIS
# =========================================================================

if (analyze_group_prevalence) {
  cat("\nRunning group-specific prevalence analysis...\n")
  group_critical_summary <- data.frame()

  for (group in sample_groups) {
    cat(sprintf("  Processing group: %s\n", group))
    group_samples    <- rownames(sample_metadata)[sample_metadata[[type_column]] == group]
    group_seqtab     <- seqtab.nochim[group_samples, , drop = FALSE]
    group_asv_counts <- colSums(group_seqtab)
    group_seqtab     <- group_seqtab[, group_asv_counts > 0, drop = FALSE]
    group_n_asvs     <- ncol(group_seqtab)

    group_matrix <- data.frame()

    for (abundance_thresh in abundance_thresholds) {
      for (prevalence_thresh in prevalence_thresholds) {
        filter_result       <- apply_combined_filter(group_seqtab, abundance_thresh, prevalence_thresh)
        removed_asv_names   <- colnames(group_seqtab)[filter_result$removed_indices]
        removed_global_idx  <- which(colnames(seqtab.nochim) %in% removed_asv_names)

        group_matrix <- rbind(group_matrix, data.frame(
          Sample_Group          = group,
          Abundance_Threshold   = abundance_thresh,
          Prevalence_Threshold  = prevalence_thresh,
          ASVs_Remaining        = filter_result$n_kept,
          ASVs_Removed          = filter_result$n_removed,
          Removed_by_Abundance  = filter_result$removed_by_abundance,
          Removed_by_Prevalence = filter_result$removed_by_prevalence,
          Percent_Remaining     = round(filter_result$n_kept / group_n_asvs * 100, 2)
        ))

        for (level in available_levels) {
          taxa_impact <- analyze_taxonomic_impact(
            taxa_matrix          = taxa_matrix,
            removed_indices      = removed_global_idx,
            taxonomic_level      = level,
            critical_taxa_list   = critical_taxa,
            sample_metadata      = sample_metadata[group_samples, , drop = FALSE],
            seqtab               = seqtab.nochim,
            abundance_threshold  = abundance_thresh,
            prevalence_threshold = prevalence_thresh
          )

          if (length(taxa_impact$critical_impacts) > 0) {
            for (taxon in names(taxa_impact$critical_impacts)) {
              group_critical_summary <- rbind(group_critical_summary, data.frame(
                Scope                = "Group",
                Sample_Group         = group,
                Abundance_Threshold  = abundance_thresh,
                Prevalence_Threshold = prevalence_thresh,
                Taxonomic_Level      = level,
                Taxon                = taxon,
                ASVs_Lost            = taxa_impact$critical_impacts[taxon]
              ))
            }
          }
        }
      }
    }

    clean_group <- stringr::str_replace_all(group, "[^A-Za-z0-9]", "_")
    readr::write_csv(group_matrix,
                     file.path(csv_dir, paste0("filtering_impact_matrix_", clean_group, ".csv")))
  }

  if (nrow(group_critical_summary) > 0) {
    readr::write_csv(group_critical_summary,
                     file.path(csv_dir, "critical_taxa_by_group.csv"))
  }
}

# =========================================================================
# STEP 5: MULTI-LEVEL TAXONOMIC IMPACT ASSESSMENT
# =========================================================================

cat("\nRunning multi-level taxonomic impact assessment across threshold range...\n")

critical_summary <- data.frame(
  Threshold              = character(),
  Taxonomic_Level        = character(),
  Taxon                  = character(),
  Count                  = numeric(),
  Sample_Groups_Affected = character(),
  stringsAsFactors       = FALSE
)

taxa_impact_by_threshold <- list()

for (t in threshold_range) {
  removed_asvs_indices <- which(asv_counts < t)
  threshold_impacts    <- list()

  if (length(removed_asvs_indices) > 0) {
    for (level in available_levels) {
      level_impact <- analyze_taxonomic_impact(
        taxa_matrix          = taxa_matrix,
        removed_indices      = removed_asvs_indices,
        taxonomic_level      = level,
        critical_taxa_list   = critical_taxa,
        sample_metadata      = sample_metadata,
        seqtab               = seqtab.nochim,
        abundance_threshold  = t,
        prevalence_threshold = NULL
      )

      threshold_impacts[[level]] <- level_impact

      if (level_impact$total_taxa_removed > 0) {
        top_n        <- min(max_displayed_taxa, length(level_impact$taxa_counts))
        top_removed  <- sort(level_impact$taxa_counts, decreasing = TRUE)[seq_len(top_n)]
        cat(sprintf("  Threshold %d | %s level - top taxa removed:\n", t, level))
        for (i in seq_along(top_removed)) {
          cat(sprintf("    %s: %d ASVs\n", names(top_removed)[i], top_removed[i]))
        }

        if (length(level_impact$critical_impacts) > 0) {
          cat(sprintf("  CRITICAL %s taxa impacted at threshold %d\n", toupper(level), t))

          for (taxon in names(level_impact$critical_impacts)) {
            count         <- level_impact$critical_impacts[taxon]
            sample_impacts <- level_impact$sample_type_impacts
            taxon_impacts  <- sample_impacts[sample_impacts$Taxon  == taxon &
                                               sample_impacts$Level == level, ]

            affected_groups_summary <- character()
            if (nrow(taxon_impacts) > 0) {
              for (j in seq_len(nrow(taxon_impacts))) {
                affected_groups_summary <- c(
                  affected_groups_summary,
                  paste0(taxon_impacts$Sample_Type[j], "(", taxon_impacts$ASVs_Lost[j], ")")
                )
              }
            }

            critical_summary <- rbind(critical_summary, data.frame(
              Threshold              = as.character(t),
              Taxonomic_Level        = level,
              Taxon                  = taxon,
              Count                  = count,
              Sample_Groups_Affected = paste(affected_groups_summary, collapse = "; "),
              stringsAsFactors       = FALSE
            ))
          }
        }
      }
    }
  } else {
    for (level in available_levels) {
      threshold_impacts[[level]] <- list(
        level = level, total_taxa_removed = 0,
        taxa_counts = table(character(0)), critical_impacts = numeric(0),
        sample_type_impacts = data.frame()
      )
    }
  }

  taxa_impact_by_threshold[[as.character(t)]] <- threshold_impacts
}

# Level-specific summaries
level_summaries <- list()
for (level in available_levels) {
  level_critical <- critical_summary[critical_summary$Taxonomic_Level == level, ]
  level_summaries[[level]] <- list(
    critical_impacts       = level_critical,
    total_critical_events  = nrow(level_critical),
    affected_thresholds    = unique(level_critical$Threshold)
  )
}

taxa_impact_results <- list(
  critical_summary         = critical_summary,
  taxa_impact_by_threshold = taxa_impact_by_threshold,
  level_summaries          = level_summaries,
  analyzed_levels          = available_levels,
  primary_level            = primary_level,
  secondary_level          = secondary_level
)

saveRDS(taxa_impact_results, file.path(save_dir, "step5_taxonomic_impact_multilevel.rds"))

if (nrow(critical_summary) > 0) {
  cat("\nMulti-level critical taxa impact summary:\n")
  for (level in available_levels) {
    level_critical <- critical_summary[critical_summary$Taxonomic_Level == level, ]
    if (nrow(level_critical) > 0) {
      cat(sprintf("  %s level impacts:\n", toupper(level)))
      print(level_critical[, c("Threshold", "Taxon", "Count", "Sample_Groups_Affected")],
            row.names = FALSE)
    } else {
      cat(sprintf("  %s level: no critical taxa impacted\n", toupper(level)))
    }
  }
} else {
  cat("\nNo critical taxa impacted at any level across all tested thresholds\n")
}

# =========================================================================
# STEP 6: DETAILED CSV EXPORT
# =========================================================================

if (save_detailed_csv) {
  cat("\nExporting detailed CSV files...\n")

  # Multi-level critical taxa summary
  if (nrow(critical_summary) > 0) {
    readr::write_csv(critical_summary,
                     file.path(csv_dir, "critical_taxa_impact_multilevel.csv"))

    for (level in available_levels) {
      level_critical <- critical_summary[critical_summary$Taxonomic_Level == level, ]
      if (nrow(level_critical) > 0) {
        readr::write_csv(
          level_critical,
          file.path(csv_dir, paste0("critical_taxa_impact_", tolower(level), ".csv"))
        )
      }
    }
  }

  # Per-threshold detailed removal tables
  results <- lapply(threshold_range, function(t) {
    removed_asvs <- which(asv_counts < t)

    if (length(removed_asvs) > 0) {
      taxa_removed <- taxa_matrix[removed_asvs, , drop = FALSE]

      removed_data <- data.frame(
        ASV_ID     = rownames(taxa_removed),
        Read_Count = asv_counts[removed_asvs],
        stringsAsFactors = FALSE
      )

      for (col in colnames(taxa_removed)) {
        level_name <- if (col %in% available_levels) col else
          available_levels[tolower(available_levels) == col]
        if (length(level_name) > 0) {
          removed_data[[paste0(level_name, "_Raw")]]   <- taxa_removed[, col]
          removed_data[[paste0(level_name, "_Clean")]] <-
            clean_taxonomic_names(taxa_removed[, col], level_name)
        } else {
          removed_data[[col]] <- taxa_removed[, col]
        }
      }

      for (group in unique(sample_metadata[[type_column]])) {
        group_samples    <- rownames(sample_metadata)[sample_metadata[[type_column]] == group]
        group_seqtab     <- seqtab.nochim[group_samples, , drop = FALSE]
        sample_counts    <- numeric(nrow(removed_data))
        sample_percents  <- numeric(nrow(removed_data))
        group_read_cnts  <- numeric(nrow(removed_data))

        for (i in seq_len(nrow(removed_data))) {
          asv_idx <- which(colnames(seqtab.nochim) == removed_data$ASV_ID[i])
          if (length(asv_idx) > 0) {
            asv_in_group       <- group_seqtab[, asv_idx]
            sample_counts[i]   <- sum(asv_in_group > 0)
            sample_percents[i] <- round(sample_counts[i] / nrow(group_seqtab) * 100, 2)
            group_read_cnts[i] <- sum(asv_in_group)
          }
        }

        clean_group <- stringr::str_replace_all(group, "[^A-Za-z0-9]", "_")
        removed_data[[paste0("Present_in_", clean_group, "_samples")]] <- sample_counts
        removed_data[[paste0("Percent_of_", clean_group, "_samples")]] <- sample_percents
        removed_data[[paste0("Total_reads_in_", clean_group)]]         <- group_read_cnts
      }

      removed_data <- removed_data[order(removed_data$Read_Count, decreasing = TRUE), ]
      readr::write_csv(removed_data,
                       file.path(csv_dir, paste0("removed_taxa_threshold_", t, ".csv")))

      if (export_sample_specific) {
        for (group in unique(sample_metadata[[type_column]])) {
          clean_group <- stringr::str_replace_all(group, "[^A-Za-z0-9]", "_")
          group_col   <- paste0("Total_reads_in_", clean_group)

          if (group_col %in% colnames(removed_data)) {
            group_data <- removed_data[removed_data[[group_col]] > 0, ]
            if (nrow(group_data) > 0) {
              group_data_sorted <- group_data[order(group_data[[group_col]],
                                                    decreasing = TRUE), ]
              readr::write_csv(
                group_data_sorted,
                file.path(csv_dir, paste0("removed_taxa_threshold_", t, "_", clean_group, ".csv"))
              )
            }
          }
        }
      }

      return(list(threshold = t, count = nrow(removed_data)))
    } else {
      return(list(threshold = t, count = 0))
    }
  })

  for (res in results) {
    if (res$count > 0) {
      cat(sprintf("  Exported detailed list for threshold %d: %d ASVs\n",
                  res$threshold, res$count))
    }
  }

  # Taxonomic removal summaries per level
  cat("Exporting taxonomic removal summaries for all levels...\n")
  for (level in available_levels) {
    all_taxa <- unique(unlist(lapply(taxa_impact_results$taxa_impact_by_threshold,
                                    function(x) if (level %in% names(x))
                                      names(x[[level]]$taxa_counts) else character(0))))
    all_taxa <- all_taxa[!is.na(all_taxa) & all_taxa != ""]

    if (length(all_taxa) > 0) {
      taxa_removal_summary <- data.frame(Taxon           = sort(all_taxa),
                                         Taxonomic_Level = level,
                                         stringsAsFactors = FALSE)

      for (t in threshold_range) {
        col_name <- paste0("Threshold_", t)
        taxa_removal_summary[[col_name]] <- 0

        t_str <- as.character(t)
        if (t_str %in% names(taxa_impact_results$taxa_impact_by_threshold) &&
            level %in% names(taxa_impact_results$taxa_impact_by_threshold[[t_str]])) {
          taxa_counts_t <- taxa_impact_results$taxa_impact_by_threshold[[t_str]][[level]]$taxa_counts
          for (taxon in names(taxa_counts_t)) {
            idx <- which(taxa_removal_summary$Taxon == taxon)
            if (length(idx) > 0) taxa_removal_summary[idx, col_name] <- taxa_counts_t[taxon]
          }
        }
      }

      readr::write_csv(taxa_removal_summary,
                       file.path(csv_dir, paste0("taxa_removal_summary_", tolower(level), ".csv")))
      cat(sprintf("  Exported %s level removal summary\n", level))
    }
  }

  # Primary level summary for backward compatibility
  if (primary_level %in% available_levels) {
    primary_data     <- taxa_impact_results$taxa_impact_by_threshold
    all_primary_taxa <- unique(unlist(lapply(primary_data,
                                             function(x) if (primary_level %in% names(x))
                                               names(x[[primary_level]]$taxa_counts) else character(0))))

    if (length(all_primary_taxa) > 0) {
      primary_summary <- data.frame(Taxon = sort(all_primary_taxa), stringsAsFactors = FALSE)

      for (t in threshold_range) {
        col_name <- paste0("Threshold_", t)
        primary_summary[[col_name]] <- 0

        t_str <- as.character(t)
        if (t_str %in% names(primary_data) &&
            primary_level %in% names(primary_data[[t_str]])) {
          taxa_counts_t <- primary_data[[t_str]][[primary_level]]$taxa_counts
          for (taxon in names(taxa_counts_t)) {
            idx <- which(primary_summary$Taxon == taxon)
            if (length(idx) > 0) primary_summary[idx, col_name] <- taxa_counts_t[taxon]
          }
        }
      }

      readr::write_csv(primary_summary,
                       file.path(csv_dir, "taxa_removal_summary_primary.csv"))
      cat(sprintf("  Exported primary level (%s) removal summary\n", primary_level))
    }
  }

  cat(sprintf("CSV export complete. Files in: %s\n", csv_dir))
} else {
  cat("Detailed CSV export skipped (disabled in configuration)\n")
}

# =========================================================================
# STEP 7: FINAL SUMMARY AND RECOMMENDATIONS
# =========================================================================

cat("\nFINAL SUMMARY\n")
cat(sprintf("  Total ASVs analyzed: %s\n", format(total_asvs, big.mark = ",")))
cat(sprintf("  Total reads: %s\n", format(total_reads, big.mark = ",")))
cat(sprintf("  Threshold range tested: %d - %d\n", min(threshold_range), max(threshold_range)))
cat(sprintf("  Taxonomic levels analyzed: %s\n", paste(available_levels, collapse = ", ")))

if (nrow(critical_summary) > 0) {
  for (level in available_levels) {
    level_critical <- critical_summary[critical_summary$Taxonomic_Level == level, ]
    if (nrow(level_critical) > 0) {
      affected_thresholds <- as.numeric(level_critical$Threshold)
      safe_thresholds_lvl <- setdiff(threshold_range, affected_thresholds)
      affected_taxa       <- unique(level_critical$Taxon)

      if (length(safe_thresholds_lvl) > 0) {
        cat(sprintf("  %s level - max safe threshold: %d reads\n",
                    level, max(safe_thresholds_lvl)))
      } else {
        cat(sprintf("  %s level - WARNING: critical taxa affected at all thresholds\n", level))
      }
      cat(sprintf("    Critical taxa at risk: %s\n", paste(affected_taxa, collapse = ", ")))
    } else {
      cat(sprintf("  %s level - no critical taxa affected at any threshold\n", level))
    }
  }

  all_affected <- as.numeric(critical_summary$Threshold)
  all_safe     <- setdiff(threshold_range, all_affected)
  if (length(all_safe) > 0) {
    cat(sprintf("  Overall max safe threshold (all levels): %d reads\n", max(all_safe)))
  } else {
    cat("  WARNING: critical taxa affected at all tested thresholds across all levels\n")
  }
} else {
  cat(sprintf("  No critical taxa affected at any level up to threshold %d\n",
              max(threshold_range)))
}

conservative_idx <- which(overall_results$Percentage_Reads_Removed <= 1)
if (length(conservative_idx) > 0) {
  cat(sprintf("  Conservative read retention (<=1%% lost): threshold <= %d\n",
              overall_results$Threshold[max(conservative_idx)]))
  moderate_idx <- which(overall_results$Percentage_Reads_Removed <= 5)
  cat(sprintf("  Moderate read retention (<=5%% lost): threshold <= %d\n",
              overall_results$Threshold[max(moderate_idx)]))
}

# Save final results object
final_results <- list(
  analysis_summary = list(
    total_asvs      = total_asvs,
    total_reads     = total_reads,
    threshold_range = threshold_range
  ),
  abundance_metrics = list(
    singletons        = singletons,
    doubletons        = doubletons,
    singleton_percent = round(singletons / total_asvs * 100, 1),
    median_abundance  = median(asv_counts)
  ),
  overall_threshold_results = overall_results,
  taxonomic_impact          = taxa_impact_results,
  critical_taxa_summary     = critical_summary
)

saveRDS(final_results, file.path(save_dir, "final_threshold_analysis.rds"))
cat(sprintf("\nFinal results saved to: %s\n", file.path(save_dir, "final_threshold_analysis.rds")))
cat(sprintf("Analysis completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

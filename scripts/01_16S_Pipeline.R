# ============================================================
# Script: 16S_V3V4_MicrobiomePipeline.R
# Description: 16S rRNA amplicon sequencing pipeline for
#              CFP shrimp microbiome study. Processes sequence
#              table through filtering, taxonomy assignment,
#              alignment, tree building, and phyloseq object
#              construction.
# Input: seqtab_nochim.rds
# Output: step4_phyloseq_object.rds
# Note: Requires external tools: VSEARCH, VeryFastTree, DECIPHER.
#       Update file paths in CONFIGURATION before running.
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
  "Biostrings", "phyloseq", "dada2", "ShortRead", "ape",
  "phangorn", "parallel", "future.apply", "phytools", "stringr", "DECIPHER"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("phyloseq", "Biostrings", "ShortRead", "DECIPHER", "dada2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
key_packages <- c("phyloseq", "dada2", "DECIPHER", "ape", "phangorn", "Biostrings")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths — update these for your environment
save_dir    <- "/Users/khanhnguyen/Library/CloudStorage/Dropbox/Auburn University/Dr. Davis' Lab/Github/Exp2434/16S_V3V4_MicrobiomePipeline/"
gg2_ref     <- file.path(save_dir, "gg2_2024_09_toGenus_trainset.fa.gz")
seqtab_file <- file.path(save_dir, "seqtab_nochim.rds")

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Parallel processing
cores <- max(1, parallel::detectCores() - 1)

# Data filtering parameters
min_count          <- 5      # Minimum read count for ASVs
low_read_threshold <- 5000   # Minimum reads per sample

# DECIPHER alignment parameters
decipher_iterations   <- 2
decipher_refinements  <- 1
decipher_use_structures <- TRUE
decipher_gap_opening  <- c(-18, -10)
decipher_gap_extension <- -3
decipher_gap_power      <- 1        # Gap cost exponent (default: 1)
decipher_terminal_gap   <- -1       # Terminal gap penalty (default: -1)
decipher_norm_power     <- 1.5      # Normalization power (default: 1.5)
decipher_restrict       <- c(-1000, 2, 10)  # Alignment range restriction
decipher_anchor         <- 0.5      # Structure anchor value (default: 0.5)

# Quality control parameters
col_gap_threshold <- 0.2
seq_gap_threshold <- 0.15
use_chi2_test     <- TRUE
chi2_p_threshold  <- 0.001

# VeryFastTree parameters
bootstrap_replicates <- 10000

# Phylogenetic annotation handling
aggregate_phylogenetic_variants <- TRUE
phylo_suffix_pattern <- "(_[A-Z0-9]+|_UCG-[0-9]+|_NK[0-9A-Z]+|_unclassified|_incertae_sedis)$"

# Critical taxa for protection and monitoring
critical_taxa <- list(
  Genus = list(
    list(name = "Vibrio",         patterns = c("Vibrio"),         exact_match = FALSE),
    list(name = "Photobacterium", patterns = c("Photobacterium"), exact_match = FALSE),
    list(name = "Lactobacillus",  patterns = c("Lactobacillus"),  exact_match = FALSE),
    list(name = "Bacillus",       patterns = c("Bacillus"),       exact_match = FALSE),
    list(name = "Vibrio_A_675473", patterns = c("Vibrio_A_675473"), exact_match = TRUE)
  ),
  Family = list(
    list(name = "Vibrionaceae",    patterns = c("Vibrionaceae"),    exact_match = FALSE),
    list(name = "Lactobacillaceae",patterns = c("Lactobacillaceae"),exact_match = FALSE),
    list(name = "Bacillaceae",     patterns = c("Bacillaceae"),     exact_match = FALSE)
  )
)

# Critical taxa filtering thresholds
critical_abundance_threshold  <- 10
critical_prevalence_threshold <- 0.02
use_critical_filters          <- TRUE

# Group-specific filtering thresholds
# Based on threshold analysis — ASV kept if it passes in ANY sample type
group_abundance_prevalence_filters <- list(
  Gut = list(
    abundance_threshold  = 10,
    prevalence_threshold = 0.02   # ceiling(110 x 0.02) = 3 samples
  ),
  Soil = list(
    abundance_threshold  = 10,
    prevalence_threshold = 0.10   # ceiling(20 x 0.10) = 2 samples
  ),
  Water = list(
    abundance_threshold  = 10,
    prevalence_threshold = 0.10   # ceiling(20 x 0.10) = 2 samples
  )
)

default_group_filter <- list(abundance_threshold = 10, prevalence_threshold = 0.05)

# Taxonomic settings
remove_prefixes     <- TRUE
handle_unclassified <- TRUE
analyze_levels      <- c("Phylum", "Class", "Order", "Family", "Genus")

set.seed(54)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

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

get_sample_type_column <- function(metadata) {
  for (col in c("type", "sample_type", "sampletype", "group", "category")) {
    if (col %in% colnames(metadata)) return(col)
  }
  return(NULL)
}

clean_taxonomic_names <- function(taxa_vector, level_name) {
  if (length(taxa_vector) == 0) return(character(0))
  taxa_clean <- taxa_vector
  if (remove_prefixes) {
    prefix <- switch(tolower(level_name),
                     kingdom = "^[dk]__", phylum = "^p__", class = "^c__",
                     order = "^o__", family = "^f__", genus = "^g__", species = "^s__", NULL)
    if (!is.null(prefix)) taxa_clean <- stringr::str_remove(taxa_clean, prefix)
  }
  taxa_clean[is.na(taxa_clean) | taxa_clean == ""] <- paste0("Unclassified_", level_name)
  return(taxa_clean)
}

check_consistency <- function(seqtab, taxa, alignment = NULL, tree = NULL, step) {
  seqtab_asvs <- colnames(seqtab)
  taxa_asvs   <- rownames(taxa)
  if (length(seqtab_asvs) != length(taxa_asvs) || !all(seqtab_asvs == taxa_asvs)) {
    stop("ASV mismatch at ", step, ": seqtab (", length(seqtab_asvs),
         ") vs taxa (", length(taxa_asvs), ")")
  }
  if (!is.null(alignment)) {
    align_asvs <- names(alignment)
    if (length(align_asvs) != length(seqtab_asvs) || !all(seqtab_asvs == align_asvs)) {
      stop("Alignment mismatch at ", step)
    }
  }
  if (!is.null(tree)) {
    if (!all(seqtab_asvs %in% tree$tip.label)) stop("Tree mismatch at ", step)
  }
}

# =========================================================================
# STEP 1: LOAD AND FILTER SEQUENCE TABLE
# =========================================================================

ps_object <- readRDS(seqtab_file)

if (!inherits(ps_object, "phyloseq")) {
  stop("Input must be a phyloseq object. Found: ", class(ps_object)[1])
}

seqtab.nochim <- as.matrix(phyloseq::otu_table(ps_object))
if (phyloseq::taxa_are_rows(ps_object)) seqtab.nochim <- t(seqtab.nochim)

# Filter low-read samples
sample_depths         <- rowSums(seqtab.nochim)
low_read_samples      <- which(sample_depths < low_read_threshold)
if (length(low_read_samples) > 0) {
  seqtab.nochim <- seqtab.nochim[-low_read_samples, , drop = FALSE]
}

sample_metadata <- if (!is.null(phyloseq::sample_data(ps_object, errorIfNULL = FALSE))) {
  data.frame(phyloseq::sample_data(ps_object))
} else {
  data.frame(Sample_ID = rownames(seqtab.nochim), group = "All_Samples",
             stringsAsFactors = FALSE, row.names = rownames(seqtab.nochim))
}

# =========================================================================
# STEP 2: GROUP-SPECIFIC ABUNDANCE + PREVALENCE FILTERING
# =========================================================================
# Strategy: ASV kept if it passes filters in ANY sample type (logical OR)
# =========================================================================

type_column   <- get_sample_type_column(sample_metadata)
if (is.null(type_column)) {
  sample_metadata$group <- "All_Samples"
  type_column           <- "group"
  group_abundance_prevalence_filters <- list(
    All_Samples = default_group_filter
  )
}

sample_types <- sample_metadata[[type_column]]
unique_types <- unique(sample_types)
total_asvs   <- ncol(seqtab.nochim)

# Add default filter for any unexpected sample types
for (utype in setdiff(unique_types, names(group_abundance_prevalence_filters))) {
  group_abundance_prevalence_filters[[utype]] <- default_group_filter
}

keep_asvs_combined  <- rep(FALSE, total_asvs)
group_filter_summary <- data.frame()

for (type in unique_types) {
  type_samples  <- sample_types == type
  type_matrix   <- seqtab.nochim[type_samples, , drop = FALSE]
  n_type        <- nrow(type_matrix)
  abund_thresh  <- group_abundance_prevalence_filters[[type]]$abundance_threshold
  prev_thresh   <- group_abundance_prevalence_filters[[type]]$prevalence_threshold
  min_samples   <- ceiling(n_type * prev_thresh)

  abund_pass    <- colSums(type_matrix)   >= abund_thresh
  prev_pass     <- colSums(type_matrix > 0) >= min_samples
  keep_in_type  <- abund_pass & prev_pass

  keep_asvs_combined <- keep_asvs_combined | keep_in_type

  group_filter_summary <- rbind(group_filter_summary, data.frame(
    Sample_Type              = type,
    N_Samples                = n_type,
    Abundance_Threshold      = abund_thresh,
    Prevalence_Threshold_Pct = prev_thresh * 100,
    Min_Samples_Required     = min_samples,
    ASVs_Removed_Abundance   = sum(!abund_pass),
    ASVs_Removed_Prevalence  = sum(abund_pass & !prev_pass),
    ASVs_Kept                = sum(keep_in_type),
    stringsAsFactors = FALSE
  ))
}

seqtab.nochim <- seqtab.nochim[, keep_asvs_combined, drop = FALSE]
write.csv(group_filter_summary,
          file.path(save_dir, "step2_group_filtering_summary.csv"),
          row.names = FALSE)

# =========================================================================
# STEP 3: TAXONOMY ASSIGNMENT
# =========================================================================

seqs  <- colnames(seqtab.nochim)
taxa  <- dada2::assignTaxonomy(seqs, gg2_ref, multithread = TRUE, verbose = TRUE)
rownames(taxa) <- seqs

# Report unclassified
unclassified <- sum(apply(taxa, 1, function(x) all(is.na(x))))
if (unclassified / nrow(taxa) > 0.5) {
  warning("Over 50% of ASVs unclassified (", unclassified, "/", nrow(taxa),
          "). Check reference database.")
}

# =========================================================================
# STEP 4: FILTER NON-BACTERIAL SEQUENCES
# =========================================================================

kingdom_clean <- sub("^[dk]__", "", taxa[, "Kingdom"])
phylum_clean  <- sub("^p__",    "", taxa[, "Phylum"])
class_clean   <- sub("^c__",    "", taxa[, "Class"])
order_clean   <- sub("^o__",    "", taxa[, "Order"])
family_clean  <- sub("^f__",    "", taxa[, "Family"])
genus_clean   <- sub("^g__",    "", taxa[, "Genus"])

# Contamination patterns
chloroplast_patterns  <- c("Chloroplast", "Streptophyta", "Bacillariophyta",
                           "Ochrophyta", "Rhodophyta", "Chlorophyceae")
mitochondria_patterns <- c("Mitochondria", "mitochondria", "mitochondrial")
eukaryota_patterns    <- c("Eukaryota", "Fungi", "Metazoa", "Plantae",
                           "Protista", "Viridiplantae", "Opisthokonta")
archaea_patterns      <- c("Archaea", "Crenarchaeota", "Euryarchaeota",
                           "Thaumarchaeota")

detect_any <- function(vec, patterns) {
  grepl(paste(patterns, collapse = "|"), vec, ignore.case = TRUE)
}

is_chloroplast  <- detect_any(kingdom_clean, chloroplast_patterns) |
  detect_any(phylum_clean,  chloroplast_patterns) |
  detect_any(class_clean,   chloroplast_patterns) |
  detect_any(order_clean,   chloroplast_patterns) |
  detect_any(family_clean,  chloroplast_patterns) |
  detect_any(genus_clean,   chloroplast_patterns)

is_mitochondria <- detect_any(kingdom_clean, mitochondria_patterns) |
  detect_any(phylum_clean,  mitochondria_patterns) |
  detect_any(class_clean,   mitochondria_patterns) |
  detect_any(order_clean,   mitochondria_patterns) |
  detect_any(family_clean,  mitochondria_patterns) |
  detect_any(genus_clean,   mitochondria_patterns)

is_eukaryota    <- detect_any(kingdom_clean, eukaryota_patterns) |
  detect_any(phylum_clean,  eukaryota_patterns) |
  detect_any(class_clean,   eukaryota_patterns)

is_archaea      <- detect_any(kingdom_clean, archaea_patterns) |
  detect_any(phylum_clean,  archaea_patterns) |
  detect_any(class_clean,   archaea_patterns)

is_unclassified <- apply(taxa, 1, function(x) all(is.na(x)))

non_bacterial <- is_chloroplast | is_mitochondria | is_eukaryota |
  is_archaea     | is_unclassified
keep_bacterial <- !non_bacterial

# Save contamination report
contam_report <- data.frame(
  Type         = c("Chloroplast", "Mitochondria", "Eukaryota", "Archaea", "Unclassified", "Bacterial"),
  ASV_Count    = c(sum(is_chloroplast), sum(is_mitochondria), sum(is_eukaryota),
                   sum(is_archaea), sum(is_unclassified), sum(keep_bacterial)),
  Percent_ASVs = round(c(sum(is_chloroplast), sum(is_mitochondria), sum(is_eukaryota),
                         sum(is_archaea), sum(is_unclassified), sum(keep_bacterial)) /
                         nrow(taxa) * 100, 2)
)
write.csv(contam_report,
          file.path(save_dir, "step4_contamination_summary.csv"),
          row.names = FALSE)

seqtab.nochim <- seqtab.nochim[, keep_bacterial, drop = FALSE]
taxa          <- taxa[keep_bacterial, , drop = FALSE]

# =========================================================================
# STEP 5: CRITICAL TAXA IMPACT CHECK
# =========================================================================

critical_summary <- data.frame()

for (type in unique_types) {
  type_samples <- rownames(sample_metadata)[sample_metadata[[type_column]] == type]
  type_matrix  <- seqtab.nochim[type_samples, , drop = FALSE]

  for (critical_item in critical_taxa$Genus) {
    taxon_name <- critical_item$name
    all_genera <- clean_taxonomic_names(taxa[, "Genus"], "Genus")

    matched_idx <- unique(unlist(lapply(critical_item$patterns, function(pat) {
      match_critical_taxon(all_genera, pat, critical_item$exact_match)
    })))

    if (length(matched_idx) == 0) next

    type_abund <- sum(colSums(type_matrix[, matched_idx, drop = FALSE]))
    type_prev  <- sum(rowSums(type_matrix[, matched_idx, drop = FALSE] > 0) > 0) /
      nrow(type_matrix)

    filter_status <- if (type_abund >= group_abundance_prevalence_filters[[type]]$abundance_threshold &&
                         type_prev  >= group_abundance_prevalence_filters[[type]]$prevalence_threshold) {
      "PASSED"
    } else "FAILED"

    critical_summary <- rbind(critical_summary, data.frame(
      Sample_Type     = type,
      Taxon           = taxon_name,
      Level           = "Genus",
      Total_ASVs      = length(matched_idx),
      Total_Abundance = type_abund,
      Prevalence_Pct  = round(type_prev * 100, 1),
      Filter_Status   = filter_status,
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(critical_summary) > 0) {
  write.csv(critical_summary,
            file.path(save_dir, "step5_critical_taxa_summary.csv"),
            row.names = FALSE)
}

# =========================================================================
# STEP 6: DECIPHER ALIGNMENT
# =========================================================================

seqs_to_align <- Biostrings::DNAStringSet(colnames(seqtab.nochim))
names(seqs_to_align) <- colnames(seqtab.nochim)

alignment <- DECIPHER::AlignSeqs(
  seqs_to_align,
  guideTree      = NULL,
  iterations     = decipher_iterations,
  refinements    = decipher_refinements,
  useStructures  = decipher_use_structures,
  gapOpening     = decipher_gap_opening,
  gapExtension   = decipher_gap_extension,
  gapPower       = decipher_gap_power,
  terminalGap    = decipher_terminal_gap,
  normPower      = decipher_norm_power,
  restrict       = decipher_restrict,
  anchor         = decipher_anchor,
  processors     = cores,
  verbose        = FALSE
)

names(alignment) <- colnames(seqtab.nochim)

# =========================================================================
# STEP 7: ALIGNMENT QUALITY CONTROL
# =========================================================================

align_matrix    <- as.matrix(alignment)
seq_lengths     <- Biostrings::width(alignment)
gap_freq        <- Biostrings::letterFrequency(alignment, "-")
gap_prop_per_seq <- gap_freq / seq_lengths
ambig_freq      <- Biostrings::letterFrequency(alignment, "NRYKMSWBDHV")
total_ambig     <- rowSums(ambig_freq)
base_freq       <- Biostrings::letterFrequency(alignment, c("A", "C", "G", "T"))
gc_content      <- rowSums(base_freq[, c("G", "C")]) / rowSums(base_freq) * 100

# Chi-squared base composition test
chi2_pvals <- apply(base_freq, 1, function(x) {
  if (sum(x) == 0) return(1)
  expected <- sum(x) * c(0.26, 0.205, 0.342, 0.193)
  chisq.test(x, p = expected / sum(expected))$p.value
})

# Quality filtering — remove sequences with invalid characters or excessive ambiguous bases
valid_chars   <- grepl("^[ACGTN-]+$", as.character(alignment))
excessive_ambig <- total_ambig / seq_lengths > 0.1
remove_seqs   <- (!valid_chars | excessive_ambig)

# Protect critical taxa from removal
all_genera <- clean_taxonomic_names(taxa[, "Genus"], "Genus")
for (critical_item in critical_taxa$Genus) {
  for (pat in critical_item$patterns) {
    protected_idx <- match_critical_taxon(all_genera, pat, critical_item$exact_match)
    remove_seqs[protected_idx] <- FALSE
  }
}

keep_seqs     <- !remove_seqs
alignment     <- alignment[keep_seqs]
seqtab.nochim <- seqtab.nochim[, keep_seqs, drop = FALSE]
taxa          <- taxa[keep_seqs, , drop = FALSE]

# Save QC summary
qc_summary <- data.frame(
  Metric  = c("Total_ASVs", "Alignment_Length", "Gap_Pct", "Sequences_Excessive_Ambig",
              "Sequences_Invalid_Chars", "Chi2_Flagged", "Removed_Total", "Final_ASVs"),
  Value   = c(
    sum(keep_seqs) + sum(!keep_seqs),
    unique(seq_lengths)[1],
    round(sum(gap_freq) / (length(alignment) * unique(seq_lengths)[1]) * 100, 2),
    sum(excessive_ambig),
    sum(!valid_chars),
    sum(chi2_pvals < chi2_p_threshold),
    sum(!keep_seqs),
    ncol(seqtab.nochim)
  )
)
write.csv(qc_summary,
          file.path(save_dir, "step7_alignment_qc_summary.csv"),
          row.names = FALSE)

# =========================================================================
# STEP 8: CHIMERA DETECTION WITH VSEARCH
# =========================================================================

vsearch_fasta   <- file.path(save_dir, "step8_for_chimera_check.fasta")
vsearch_output  <- file.path(save_dir, "step8_chimera_free.fasta")
vsearch_chimeras <- file.path(save_dir, "step8_chimeras.fasta")

ungapped <- Biostrings::DNAStringSet(gsub("-", "", as.character(alignment)))
names(ungapped) <- names(alignment)
Biostrings::writeXStringSet(ungapped, vsearch_fasta)

system2("vsearch", args = c(
  "--uchime_denovo", shQuote(vsearch_fasta),
  "--chimeras",      shQuote(vsearch_chimeras),
  "--nonchimeras",   shQuote(vsearch_output),
  "--mindiffs", 3, "--mindiv", 0.8, "--minh", 0.28, "--threads", 1
), stdout = FALSE, stderr = FALSE)

if (file.exists(vsearch_output)) {
  nonchimeric     <- Biostrings::readDNAStringSet(vsearch_output)
  chimeric_ids    <- if (file.exists(vsearch_chimeras) && file.size(vsearch_chimeras) > 0) {
    names(Biostrings::readDNAStringSet(vsearch_chimeras))
  } else character(0)

  chimera_report <- data.frame(
    Total_ASVs            = length(alignment),
    Chimeras_Detected     = length(chimeric_ids),
    Chimera_Pct           = round(length(chimeric_ids) / length(alignment) * 100, 2),
    Nonchimeric_Sequences = length(nonchimeric)
  )
  write.csv(chimera_report,
            file.path(save_dir, "step8_chimera_summary.csv"),
            row.names = FALSE)

  # Note: chimeras are flagged but not removed here — reported for transparency
}

file.remove(vsearch_fasta)

# =========================================================================
# STEP 9: PHYLOGENETIC TREE WITH VERYFASTTREE
# =========================================================================

align_fasta <- file.path(save_dir, "step9_alignment_for_tree.fasta")
Biostrings::writeXStringSet(alignment, align_fasta)

tree_output <- file.path(save_dir, "step9_veryfasttree.tree")

system(paste(
  "veryfasttree",
  "-nt -gtr -gamma -double-precision",
  "-boot", bootstrap_replicates,
  "-threads", cores,
  "-threads-level 2",
  "-seed 54",
  shQuote(align_fasta),
  ">", shQuote(tree_output)
))

if (!file.exists(tree_output) || file.size(tree_output) == 0) {
  stop("VeryFastTree failed to produce a tree file.")
}

tree <- ape::read.tree(tree_output)

# Process bootstrap support values
if (!is.null(tree$node.label)) {
  support <- as.numeric(tree$node.label)
  support[is.na(support)] <- 0
  if (max(support, na.rm = TRUE) <= 1) support <- support * 100
  tree$node.label <- as.character(round(support))
}

# Root tree using midpoint rooting
if (!ape::is.rooted(tree)) {
  tree <- tryCatch(
    phytools::midpoint.root(tree),
    error = function(e) tree
  )
}

# Save Nexus files
ape::write.nexus(tree, file = file.path(save_dir, "step9_final_tree_simple.nex"))

file.remove(align_fasta)

# =========================================================================
# STEP 10: PHYLOSEQ OBJECT CONSTRUCTION
# =========================================================================

# Align all objects to common ASVs
common_asvs   <- Reduce(intersect, list(
  colnames(seqtab.nochim),
  rownames(taxa),
  tree$tip.label
))

if (length(common_asvs) == 0) stop("No common ASVs across seqtab, taxa, and tree.")

seqtab.nochim <- seqtab.nochim[, common_asvs, drop = FALSE]
taxa          <- taxa[common_asvs, , drop = FALSE]
tree          <- ape::keep.tip(tree, common_asvs)

# Align sample metadata
common_samples <- intersect(rownames(seqtab.nochim), rownames(sample_metadata))
seqtab.nochim  <- seqtab.nochim[common_samples, , drop = FALSE]
sample_metadata <- sample_metadata[common_samples, , drop = FALSE]

# Build phyloseq object
otu_obj    <- phyloseq::otu_table(seqtab.nochim, taxa_are_rows = FALSE)
tax_obj    <- phyloseq::tax_table(taxa)
tree_obj   <- tree
sample_obj <- phyloseq::sample_data(sample_metadata)

physeq <- phyloseq::phyloseq(otu_obj, tax_obj, tree_obj, sample_obj)

# Verify UniFrac
unifrac_test <- tryCatch({
  test_sub <- phyloseq::prune_samples(
    sample(phyloseq::sample_names(physeq), min(5, phyloseq::nsamples(physeq))), physeq
  )
  phyloseq::UniFrac(test_sub, weighted = FALSE)
  "PASSED"
}, error = function(e) paste("FAILED:", e$message))

# Save final summary
final_summary <- data.frame(
  Metric = c("Total_Samples", "Total_ASVs", "Total_Reads", "Min_Reads",
             "Max_Reads", "Median_Reads", "Tree_Tips", "Tree_Rooted",
             "UniFrac_Test"),
  Value  = c(
    phyloseq::nsamples(physeq),
    phyloseq::ntaxa(physeq),
    sum(phyloseq::sample_sums(physeq)),
    min(phyloseq::sample_sums(physeq)),
    max(phyloseq::sample_sums(physeq)),
    median(phyloseq::sample_sums(physeq)),
    ape::Ntip(tree),
    ape::is.rooted(tree),
    unifrac_test
  )
)
write.csv(final_summary,
          file.path(save_dir, "step10_final_phyloseq_summary.csv"),
          row.names = FALSE)

# =========================================================================
# SAVE PHYLOSEQ OBJECT
# =========================================================================

saveRDS(physeq, file.path(save_dir, "step4_phyloseq_object.rds"))

cat("Pipeline complete.\n")
cat("Phyloseq object saved to:", file.path(save_dir, "step4_phyloseq_object.rds"), "\n")
cat("Samples:", phyloseq::nsamples(physeq), "\n")
cat("ASVs:", phyloseq::ntaxa(physeq), "\n")
cat("UniFrac test:", unifrac_test, "\n")

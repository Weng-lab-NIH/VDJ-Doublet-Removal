#!/usr/bin/env Rscript

# Identify VDJ doublets from Cell Ranger filtered contig annotations.
# Doublets are called from productive TRB multiplicity, and when a BCR file is
# provided, from productive IGH multiplicity and TCR/BCR overlap. The script
# removes detected doublets from the provided TCR, BCR, and optional GEX files
# and returns filtered files and a summary of removed cells.

suppressMessages(library(dplyr))
suppressMessages(library(tibble))

args <- commandArgs(trailingOnly = TRUE)

TCR_PATH <- NA_character_
BCR_PATH <- NA_character_
GEX_PATH <- NA_character_
OUT_DIR  <- NA_character_

i <- 1
while (i <= length(args)) {
  a <- args[i]

  if (a == "--tcr") {
    TCR_PATH <- args[i + 1]
    i <- i + 2
  } else if (a == "--bcr") {
    BCR_PATH <- args[i + 1]
    i <- i + 2
  } else if (a == "--gex") {
    GEX_PATH <- args[i + 1]
    i <- i + 2
  } else if (a == "--out") {
    OUT_DIR <- args[i + 1]
    i <- i + 2
  } else {
    stop(paste("Unknown argument:", a))
  }
}

if (is.na(TCR_PATH) || is.na(OUT_DIR)) {
  stop(
    paste(
      "Usage:",
      "Rscript vdj_doublet.R --tcr <filtered_contig_annotations_tcr.csv>",
      "[--bcr <filtered_contig_annotations_bcr.csv>]",
      "[--gex <gex.csv>]",
      "--out <output_directory>"
    )
  )
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Normalize strings for case-insensitive matching.
normalize_string <- function(x) {
  tolower(trimws(ifelse(is.na(x), "", as.character(x))))
}

# Convert common truth-like values to logical TRUE.
is_true_value <- function(x) {
  normalize_string(x) %in% c("true", "t", "1", "yes")
}

# Reduce a barcode to the prefix before the dash suffix.
barcode_prefix <- function(x) {
  sub("-.*$", "", as.character(x))
}

# Read a CSV file if it exists.
safe_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  suppressWarnings(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
}

# Warn when the provided path does not look like a filtered contig annotation file.
warn_if_not_filtered_contig <- function(path, label) {
  if (is.na(path)) {
    return(invisible(NULL))
  }

  file_name <- basename(path)

  if (!grepl("filtered_contig_annotations", file_name, ignore.case = TRUE)) {
    warning(
      paste0(
        label,
        " path does not appear to be a filtered_contig_annotations file: ",
        path
      )
    )
  }
}

# Stop if a required column is missing.
require_columns <- function(df, required_cols, label) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        label,
        " is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

# Return productive barcode prefixes for selected chain names.
get_productive_chain_prefixes <- function(df, chain_values) {
  if (is.null(df)) {
    return(character(0))
  }

  df %>%
    mutate(
      barcode_prefix = barcode_prefix(barcode),
      chain_norm = normalize_string(chain),
      productive_flag = ifelse(is.na(productive), NA, is_true_value(productive))
    ) %>%
    filter(chain_norm %in% normalize_string(chain_values)) %>%
    filter(is.na(productive_flag) | productive_flag) %>%
    filter(!is.na(barcode_prefix), barcode_prefix != "") %>%
    pull(barcode_prefix) %>%
    unique()
}

# Return barcode prefixes with more than one productive contig for a chain class.
get_multi_chain_doublets <- function(df, chain_values, count_name = "n_chain") {
  if (is.null(df)) {
    return(character(0))
  }

  df %>%
    mutate(
      barcode_prefix = barcode_prefix(barcode),
      chain_norm = normalize_string(chain),
      productive_flag = ifelse(is.na(productive), NA, is_true_value(productive))
    ) %>%
    filter(chain_norm %in% normalize_string(chain_values)) %>%
    filter(is.na(productive_flag) | productive_flag) %>%
    filter(!is.na(barcode_prefix), barcode_prefix != "") %>%
    count(barcode_prefix, name = count_name) %>%
    filter(.data[[count_name]] > 1) %>%
    pull(barcode_prefix) %>%
    unique()
}

# Remove rows whose barcode prefix is in the doublet set.
filter_by_doublets <- function(df, barcode_column, doublet_prefixes) {
  if (is.null(df)) {
    return(NULL)
  }

  if (!(barcode_column %in% names(df))) {
    stop(paste("Input file is missing required barcode column:", barcode_column))
  }

  df %>%
    mutate(.barcode_prefix = barcode_prefix(.data[[barcode_column]])) %>%
    filter(!(.barcode_prefix %in% doublet_prefixes)) %>%
    select(-.barcode_prefix)
}

# Build a table listing detected doublet barcodes and their reason.
build_removed_table <- function(doublet_prefixes, tcr_doublets, bcr_doublets, overlap_doublets) {
  if (length(doublet_prefixes) == 0) {
    return(tibble(
      barcode_prefix = character(0),
      reason = character(0)
    ))
  }

  tibble(barcode_prefix = sort(unique(doublet_prefixes))) %>%
    mutate(
      reason = case_when(
        barcode_prefix %in% tcr_doublets &
          barcode_prefix %in% bcr_doublets &
          barcode_prefix %in% overlap_doublets ~ "TRB_MULTI and IGH_MULTI and TCR_BCR_OVERLAP",
        barcode_prefix %in% tcr_doublets &
          barcode_prefix %in% bcr_doublets ~ "TRB_MULTI and IGH_MULTI",
        barcode_prefix %in% tcr_doublets &
          barcode_prefix %in% overlap_doublets ~ "TRB_MULTI and TCR_BCR_OVERLAP",
        barcode_prefix %in% bcr_doublets &
          barcode_prefix %in% overlap_doublets ~ "IGH_MULTI and TCR_BCR_OVERLAP",
        barcode_prefix %in% tcr_doublets ~ "TRB_MULTI",
        barcode_prefix %in% bcr_doublets ~ "IGH_MULTI",
        barcode_prefix %in% overlap_doublets ~ "TCR_BCR_OVERLAP",
        TRUE ~ "UNKNOWN"
      )
    )
}

# Check inputs and load data.
warn_if_not_filtered_contig(TCR_PATH, "TCR")
warn_if_not_filtered_contig(BCR_PATH, "BCR")

tcr_df <- safe_read_csv(TCR_PATH)
if (is.null(tcr_df)) {
  stop("TCR contig annotation file could not be read.")
}

bcr_df <- safe_read_csv(BCR_PATH)
gex_df <- safe_read_csv(GEX_PATH)

require_columns(tcr_df, c("barcode", "chain", "productive"), "TCR contig annotation file")

if (!is.null(bcr_df)) {
  require_columns(bcr_df, c("barcode", "chain", "productive"), "BCR contig annotation file")
}

if (!is.null(gex_df)) {
  require_columns(gex_df, c("cell.barcode"), "GEX file")
}

# Collect productive barcode sets for overlap testing.
tcr_productive_prefixes <- get_productive_chain_prefixes(
  tcr_df,
  c("TRB", "trb")
)
bcr_productive_prefixes <- get_productive_chain_prefixes(
  bcr_df,
  c("IGH")
)

# Call doublets from productive multi-chain barcodes.
tcr_doublet_prefixes <- get_multi_chain_doublets(
  tcr_df,
  c("TRB", "trb"),
  count_name = "n_trb"
)
bcr_doublet_prefixes <- get_multi_chain_doublets(
  bcr_df,
  c("IGH"),
  count_name = "n_igh"
)

# If both TCR and BCR are provided, also call overlap doublets.
overlap_doublet_prefixes <- character(0)
if (!is.null(bcr_df)) {
  overlap_doublet_prefixes <- intersect(
    unique(tcr_productive_prefixes),
    unique(bcr_productive_prefixes)
  )
}

# Use the union of all detected doublet barcode prefixes.
all_doublet_prefixes <- sort(unique(c(
  tcr_doublet_prefixes,
  bcr_doublet_prefixes,
  overlap_doublet_prefixes
)))

# Filter each provided file using the same detected doublet barcodes.
tcr_no_doublets <- filter_by_doublets(tcr_df, "barcode", all_doublet_prefixes)

bcr_no_doublets <- NULL
if (!is.null(bcr_df)) {
  bcr_no_doublets <- filter_by_doublets(bcr_df, "barcode", all_doublet_prefixes)
}

gex_no_doublets <- NULL
if (!is.null(gex_df)) {
  gex_no_doublets <- filter_by_doublets(gex_df, "cell.barcode", all_doublet_prefixes)
}

removed_table <- build_removed_table(
  all_doublet_prefixes,
  tcr_doublet_prefixes,
  bcr_doublet_prefixes,
  overlap_doublet_prefixes
)

# Summarize the input and output counts for each file type.
summary_df <- tibble(
  input_tcr_cells = length(unique(barcode_prefix(tcr_df$barcode))),
  productive_tcr_cells = length(unique(tcr_productive_prefixes)),
  trb_multi_doublets = length(unique(tcr_doublet_prefixes)),
  input_bcr_cells = if (is.null(bcr_df)) NA_integer_ else length(unique(barcode_prefix(bcr_df$barcode))),
  productive_bcr_cells = if (is.null(bcr_df)) NA_integer_ else length(unique(bcr_productive_prefixes)),
  igh_multi_doublets = if (is.null(bcr_df)) NA_integer_ else length(unique(bcr_doublet_prefixes)),
  tcr_bcr_overlap_doublets = if (is.null(bcr_df)) NA_integer_ else length(unique(overlap_doublet_prefixes)),
  total_doublet_barcodes = length(unique(all_doublet_prefixes)),
  tcr_rows_input = nrow(tcr_df),
  tcr_rows_output = nrow(tcr_no_doublets),
  bcr_rows_input = if (is.null(bcr_df)) NA_integer_ else nrow(bcr_df),
  bcr_rows_output = if (is.null(bcr_no_doublets)) NA_integer_ else nrow(bcr_no_doublets),
  gex_rows_input = if (is.null(gex_df)) NA_integer_ else nrow(gex_df),
  gex_rows_output = if (is.null(gex_no_doublets)) NA_integer_ else nrow(gex_no_doublets)
)

# Write output files.
write.csv(
  summary_df,
  file.path(OUT_DIR, "doublet_summary.csv"),
  row.names = FALSE
)

write.csv(
  removed_table,
  file.path(OUT_DIR, "removed_doublets.csv"),
  row.names = FALSE
)

write.csv(
  tcr_no_doublets,
  file.path(OUT_DIR, "tcr_no_doublets.csv"),
  row.names = FALSE
)

if (!is.null(bcr_no_doublets)) {
  write.csv(
    bcr_no_doublets,
    file.path(OUT_DIR, "bcr_no_doublets.csv"),
    row.names = FALSE
  )
}

if (!is.null(gex_no_doublets)) {
  write.csv(
    gex_no_doublets,
    file.path(OUT_DIR, "gex_no_doublets.csv"),
    row.names = FALSE
  )
}
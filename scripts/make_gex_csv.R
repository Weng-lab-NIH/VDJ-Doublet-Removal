#!/usr/bin/env Rscript

# Helper script designed to take a filtered_feature_barcode_matrix file 
# from the output of a cellranger multi run and return a csv of GEX cell
# barcodes for the vdj doublet finder. Optional step for doublet removal

suppressMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)

MATRIX_DIR <- NA_character_
OUT_FILE   <- NA_character_

i <- 1
while (i <= length(args)) {
  a <- args[i]

  if (a == "--matrix") {
    MATRIX_DIR <- args[i + 1]
    i <- i + 2
  } else if (a == "--out") {
    OUT_FILE <- args[i + 1]
    i <- i + 2
  } else {
    stop(paste("Unknown argument:", a))
  }
}

if (is.na(MATRIX_DIR) || is.na(OUT_FILE)) {
  stop(
    paste(
      "Usage:",
      "Rscript make_gex_csv.R",
      "--matrix <filtered_feature_bc_matrix_directory>",
      "--out <output_csv>"
    )
  )
}

if (!dir.exists(MATRIX_DIR)) {
  stop(paste("Matrix directory does not exist:", MATRIX_DIR))
}

# Read the cellranger matrix and select the gene expression assay.
# If cellranger count was used in place of multi, it will not return a list and 
# the matrix file can be used without indexing. Here, an HTO library was part of
# the cellranger multi pipeline so there are two count matrices
gex_all <- Read10X(MATRIX_DIR)

if (is.list(gex_all) && "Gene Expression" %in% names(gex_all)) {
  gex_mat <- gex_all[["Gene Expression"]]
} else if (is.list(gex_all)) {
  gex_mat <- gex_all[[1]]
} else {
  gex_mat <- gex_all
}

# Build a Seurat object to generate standard metadata fields.
gex_seu <- CreateSeuratObject(
  counts = gex_mat,
  min.cells = 3,
  min.features = 200
)

gex_meta <- as.data.frame(gex_seu@meta.data, stringsAsFactors = FALSE)
gex_meta$cell.barcode <- rownames(gex_seu@meta.data)

write.csv(gex_meta, OUT_FILE, row.names = FALSE)
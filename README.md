# VDJ-Doublet-Removal
This repository contains scripts for identifying and removing Scripts for identifying doublets based on TCR CDR3 beta chain and optionally BCR IGH chain. Cell doublets refer to two cells that share the same 10x barcode either through co-encapsulation in one droplet or due to barcode collisions which leads to one barcode capturing transcripts from two cells.

Doublets are identified from:
- more than one productive TCR beta chain (TRB) per barcode
- more than one productive BCR heavy chain (IGH) per barcode when a BCR file is provided
- overlap between cells that have productive TCR and productive BCR barcodes when both are provided

This script detects doublets from filtered_contig_annotations file from the cellranger vdj output. We assume the cellranger framework used is cellranger multi. The same detected doublet barcodes can also be removed from an optional GEX metadata CSV.

# Example pipeline run
The bash wrapper scripts are included in the bash/ directory to create the GEX metadata CSV and run the VDJ doublet script from example data.

Run from the bash root:
1) bash make_gex_csv.bash
2) bash run_doublet_example.bash

## Input files
TCR input:
- Cell Ranger filtered contig annotation found typically at:
  "ExampleID_multi/outs/per_sample_outs/ExampleID_multi/vdj_t/filtered_contig_annotations.csv"
  after a cell ranger run.

Optional BCR input:
- Cell Ranger filtered contig annotation CSV found as the TCR contig above, under vdj_b.

Optional GEX input:
- GEX metadata CSV containing at least: cell.barcode

Preparing the GEX CSV:
- If a GEX CSV is not already available, it can be generated from a Cell Ranger matrix directory:
  "ExampleID_multi/outs/per_sample_outs/ExampleID_multi/count/sample_filtered_feature_bc_matrix/"
  using make_gex_csv.R.

## Output files
The output directory may contain:
- tcr_no_doublets.csv
- bcr_no_doublets.csv
- gex_no_doublets.csv
- removed_doublets.csv
- doublet_summary.csv

### Output File descriptions
The output directory may contain:
- tcr_no_doublets.csv
     - Filtered TCR contig annotation table with detected doublet barcodes removed.
- bcr_no_doublets.csv (If BCR file was provided)
- gex_no_doublets.csv (If GEX metadata was provided)
- removed_doublets.csv
     - A list of removed cells that were identifies as doublets.
- doublet_summary.csv 
     - A summary of total cells, total productive cells, and removed doublet numbers.

## Dependencies
R version 4.0 or higher

Required R packages:
- dplyr
- tibble
- Seurat (required only for make_gex_csv.R)

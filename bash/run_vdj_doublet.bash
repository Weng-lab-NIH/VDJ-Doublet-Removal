#!/bin/bash
  
module load R

Rscript ../scripts/vdj_doublet.R \
	--tcr ../example/filtered_contig_annotations.csv \
	--gex ../example/merged_gex.csv \
	--out ../output/
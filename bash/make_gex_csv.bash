#!/bin/bash
  
module load R

Rscript ../scripts/make_gex_csv.R \
	--matrix ../example/sample_filtered_feature_bc_matrix/ \
	--out ../example/merged_gex.csv
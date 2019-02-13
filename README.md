# OPTM-Analysis

This workflow aims to extract and identify oxidative stress induced post-translational modifications sites, peptides, and proteins as well as quantify their occupancy from high-throughput proteomics raw files. 

## dtaselect_read_step1.R
Run dtaselect_read_step1.R in the directory with the DTA files.
Type:
''' bash
Rscript dtaselect_read_step1.R
'''

## dtaselect_merge_step2.R
Type:
''' bash
Rscript dtaselect_merge_step2.R mass_shift site_mass_shift output_file
'''

## gygi_optm_calculation_step3.R
Upload the peptides to peptidematch, and download the result-file-from-PIR.txt
To calculate O-PTM Occupancy, type:
''' bash
Rscript gygi_optm_calculation_step3.R input_file multiple_peptide_match_data
'''


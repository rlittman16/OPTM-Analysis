# OPTM-Analysis

This workflow aims to quantify oxidative stress induced post-translational modifications sites, peptides, and proteins occupancy from high-throughput proteomics raw files. 

Test data can be found at:
https://drive.google.com/drive/folders/1f_HVhH9rhb9l1ZOesmU8CGO_lod4tQsO?usp=sharing

## dtaselect_read_step1.R
Run dtaselect_read_step1.R in the directory with the DTA files. <br/>
Type:
``` bash
Rscript dtaselect_read_step1.R
```

## dtaselect_merge_step2.R
Type:
``` bash
Rscript dtaselect_merge_step2.R mass_shift site_mass_shift output_file
```

## optm_calculation_step3.R
Upload the peptides to peptidematch, and download the result-file-from-PIR.txt. <br/>
https://research.bioinformatics.udel.edu/peptidematch/batchpeptidematch.jsp <br/>
To calculate O-PTM Occupancy <br/>
Type:
``` bash
Rscript optm_calculation_step3.R input_file multiple_peptide_match_data
```


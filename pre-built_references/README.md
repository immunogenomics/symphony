This directory contains pre-built Symphony reference .rds objects that can be downloaded and used to map new query datasets.

References available for download:

10x PBMCs Atlas (pbmcs_10x_reference.rds)
Pancreatic Islet Cells Atlas (pancreas_plate-based_reference.rds)
Fetal Liver Hematopoiesis Atlas (fetal_liver_reference_3p.rds)
Healthy Fetal Kidney Atlas (kidney_healthy_fetal_reference.rds)
T cell CITE-seq atlas (tbru_ref.rds)
Cross-tissue Fibroblast Atlas (see here)
Cross-tissue Inflammatory Immune Atlas (zhang_reference.rds)
Tabula Muris Senis (FACS) Atlas (TMS_facs_reference.rds)
To read in a reference into R, one may simply execute: reference = readRDS('path/to/reference_name.rds')

Note: To be able to map query datasets into the reference UMAP coordinates, you must also download the corresponding 'uwot_model' file and set the reference$save_uwot_path.

This directory contains pre-built Symphony reference .rds objects that can be downloaded and used to map new query datasets.

The pre-built references available for download include:

- pbmcs_10x_reference.rds: Atlas of PBMCs (20,571 cells) sequenced with three 10x protocols (3'v1, 3'v2, and 5').
- pancreas_plate-based_reference.rds: Atlas of pancreatic islet cells (5,887 cells from 32 donors) from four separate studies.
- fetal_liver_reference_3p.rds: Atlas of fetal liver cells from Popescu et al. (2019) (113,063 cells from 14 donors), sequenced with 10x 3' chemistry.
- (To be released upon publication) Multimodal Memory T cell CITE-seq atlas (500,089 cells from 271 samples)

To read in a reference into R, one may simply execute: reference = readRDS('path/to/reference_name.rds')

Note: In order to map query cells onto the reference UMAP coordinates (e.g. to visualize reference and query cells together), you will need to save the path to the corresponding reference uwot_model file in the reference object's reference$save_uwot_path slot in order to load the uwot model for query mapping. This is due to a technicality of how the uwot package saves and loads UMAP models. If you only wish to map the query cells into the harmonized reference embedding (and compute your own separate UMAP embedding for visualization), you may ignore this step.


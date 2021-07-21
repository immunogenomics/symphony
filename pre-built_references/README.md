Pre-built Symphony reference objects can be used to map new query datasets.

To avoid GitHub memory limits and for clearer versioning, we use Zenodo to house Symphony references.
Download the reference of interest from Zenodo and place in this directory.

Link to Zenodo: [https://zenodo.org/record/5090425](https://zenodo.org/record/5090425)

References available for download: 

| Reference atlas | Filename | Num cells | Description
| --- | ----------- | --------- | ----
| **10x PBMCs Atlas** | pbmcs_10x_reference.rds | 20,571 cells | Healthy human PBMCs sequenced with three 10x protocols (3'v1, 3'v2, and 5')
| **Pancreatic Islet Cells Atlas** | pancreas_plate-based_reference.rds | 5,887 cells from 32 donors | Human pancreatic islet cells from four separate studies (Segerstolpe et al., Lawlor et al., Grun et al., Muraro et al.)
| **Fetal Liver Hematopoeisis Atlas** | fetal_liver_reference_3p.rds | 113,063 cells from 14 donors | Human fetal liver cells (from Popescu et al., 2019), sequenced with 10x (3')
| **Healthy Fetal Kidney Atlas** | kidney_healthy_fetal_reference.rds | 27,203 cells from 6 samples | Healthy human fetal kidney cells (from Stewart et al., 2019).
| **Memory T Cell (CITE-seq) Atlas** | tbru_ref.rds | 500,089 cells from 259 donors | Memory T cells from a tuberculosis cohort assayed with CITE-seq (Nathan et al., 2021)
| **Cross-tissue Fibroblast Atlas** | fibroblast_atlas.rds | 79,148 cells from 74 samples | Human fibroblasts across inflammatory diseases in the lung, gut, synovium, and salivary gland (Korsunsky et al., 2021)
| **Cross-tissue Inflammatory Immune Atlas** | zhang_reference.rds | 307,084 immune cells from 125 donors | Human immune cells across 6 inflammatory diseases (from Zhang et al., 2021)
| **Tabula Muris Senis (FACS) Atlas** | TMS_facs_reference.rds | 110,824 cells from 19 mice | Mouse cells across 23 tissues and organs

To read in a reference into R, simply execute: `reference = readRDS('path/to/reference_name.rds')`

Note: To be able to map query datasets into the reference UMAP coordinates, you must also download the corresponding 'uwot_model' file and set the `reference$save_uwot_path`

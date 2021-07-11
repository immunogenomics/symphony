Pre-built Symphony reference objects that can be downloaded and used to map new query datasets

Link to Zenodo: [https://zenodo.org/record/5090425#.YOqe_hNKhTY](https://zenodo.org/record/5090425#.YOqe_hNKhTY)


References available for download: 

- 10x PBMCs Atlas
- Pancreatic Islet Cells Atlas
- Fetal Liver Hematopoiesis Atlas
- Healthy Fetal Kidney Atlas
- T cell CITE-seq Atlas
- Cross-tissue Fibroblast Atlas
- Cross-tissue Inflammatory Immune Atlas
- Tabula Muris Senis (FACS) Atlas

To read in a reference into R, one may simply execute: reference = readRDS('path/to/reference_name.rds')

Note: To be able to map query datasets into the reference UMAP coordinates, you must also download the corresponding 'uwot_model' file and set the reference$save_uwot_path.

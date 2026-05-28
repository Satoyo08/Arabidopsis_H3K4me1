# Arabidopsis_H3K4me1

Codes used for 'Transcription-coupled and epigenome-encoded mechanisms direct H3K4 methylation (2022)'

for codes for Supplementary Figure.5, please see Figure2 and 3.r

for codes for Supplementary Figure 9,10,11, please see Figure4.r

## Random Forest and SVM analysis pipeline

This repository also includes two machine learning pipelines for genomic feature analysis:

- `RF_functions.R` — custom Random Forest functions for gene group classification
- `RandomForest_arabidopsis.Rmd` — Random Forest workflow used in the paper
- `SVM_training.ipynb` — Support Vector Machine training notebook (downstream visualizaiton: Figure4.r)
- `custom_functions.r` — shared utility functions

This pipeline has been adopted as a core analytical tool in subsequent published work beyond the original 2022 paper. If you use these tools in your research, please cite the original publication (see below).

## Citation

If you use this code, please cite:
Oya, S., Takahashi, M., Takashima, K., Kakutani, T. & Inagaki, S. Transcription-coupled and epigenome-encoded mechanisms direct H3K4 methylation. Nat. Commun. 13, 4521 (2022).
  


## License



## update 2024.01.11

Upon multiple requests, we added the following list of genes in data/list_of_genes folder.

- ATX1/2/R7-marked genes
- ATX3/4/5-marked genes
- ATXR3-marked genes
- ATX1-bound genes
- ATX2-bound genes
- ATXR7-bound genes

As a note for reproducibility, those lists can also be generated (as R object) by running Figure1.R L58- for ATX*-marked genes, and by Supplementary_Figure3.R L8- for ATX*-bound genes.

- **Are you looking for the formatted data table behind scatter plot of Figure1d?** Please use data/Figure1/ChIP1_RPKM.txt and data/Figure1/ChIP1_RPM.txt.
- Are you looking for browser tracks? Sorry, we didn't upload them in GEO nor here — tell us (soya[at]ucdavis.edu) what you need and we'll email you.

Can't find what you are looking for? Pardon my bad code annotation! Please ask via the Issues section or contact soya[at]ucdavis.edu

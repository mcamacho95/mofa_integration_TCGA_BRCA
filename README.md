# mofa_integration_TCGA_BRCA
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![MOFA2](https://img.shields.io/badge/MOFA2-Multi--Omics%20Integration-2E9FDF)](https://biofam.github.io/MOFA2/)

Reproducible pipeline developed for investigating breast cancer (BC) molecular heterogeneity and Triple-Negative Breast Cancer (TNBC) therapeutic targets in **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** cohort.

This framework integrates six molecular modalities—**mRNA, lncRNA, miRNA, DNA methylation, copy-number variation (CNV), and somatic mutations**—via **Multi-Omics Factor Analysis (MOFA+)**, followed by clinical phenotype association testing, differential expression (`limma`), pathway enrichment (`fgsea`), survival modelling (Kaplan–Meier, multivariable Cox, and Elastic Net `glmnet`), and drug–gene interaction profiling (**DGIdb**).

---

## Repository Structure

```text
mofa_integration_TCGA_BRCA/
│
├── README.md                          # Pipeline documentation and setup guide
├── .gitignore                         # Excludes large raw data matrices and HDF5 models
│
├── data/                              # Local directory for input and intermediate datasets
│
├── results/
│   ├── figures/                       # Generated high-resolution figures (UMAP, ROC, KM, Volcano)
│   ├── models/                        # Trained MOFA+ HDF5 models and ELBO summaries
│   └── tables/                        # Statistical summaries, DGE, GSEA, Cox, and DGIdb tables
│
└── src/
    ├── Data_preprocessing.R            # Filtering, normalization (VST/M-values/Z-scores), and alignment
    ├── tnbc_status_annotation.R        # IHC and FISH adjudication for TNBC clinical classification
    ├── Data_visualization_&_Overview.R # Cohort demographics summary and single-omic UMAP projections
    ├── Mofa_train.R                    # Multi-seed MOFA+ model training (10 seeds) & ELBO selection
    ├── Mofa_results.R                  # Variance decomposition, clinical associations, ROC & heatmaps
    ├── DGE_limma.R                     # Differential expression (limma) & Factor 1 weight concordance
    ├── GSEA.R                          # Gene Set Enrichment Analysis (MSigDB Hallmark & GO:BP)
    ├── Survival_analysis.R             # Multivariable/univariate Cox, KM curves & Elastic Net Cox
    └── Drug_repurposing.R              # DGIdb target prioritization & multi-target drug screening

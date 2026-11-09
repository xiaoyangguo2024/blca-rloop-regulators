# R-loop regulator activity in bladder cancer: analysis code

Analysis code for the manuscript

> Xiaoyang Guo, Mingjie Wu, Zhiyue Sun, Zhenghang Liao, Ke Zeng, Bo Wen. **R-Loop Regulator Activity in Bladder Cancer Mainly Reflects Proliferation: An Integrative Single-Cell, Spatial and Bulk Transcriptomic Study.** Manuscript, 2026.
>
> Department of Urology, Shenzhen Hospital of Integrated Traditional Chinese and Western Medicine, Shenzhen, China.

The scripts score 351 evidence-supported R-loop regulators from R-loopBase v2.0 in single-cell, spatial and bulk transcriptomic data from bladder cancer, build a 16-gene R-loop-related risk score (RLRS) and benchmark it against a proliferation score, random gene sets and 97 published bladder cancer signatures.

## Contents

| Path | Purpose |
|---|---|
| `R/00_install.R` | Installs Bioconductor, CRAN and GitHub packages (GitHub packages pinned to fixed commits, see below) |
| `R/00_config.R` | Paths, parameters (random seed 20261107; at most 1,500 cells per single-cell sample) and shared functions; sourced by every script |
| `R/01_geneset.R` | Downloads the R-loopBase regulator lists and builds the R-loop gene sets |
| `R/02_download_bulk.R` | Downloads and harmonises the bulk cohorts: TCGA-BLCA (training); GSE13507, GSE31684, GSE32894 and GSE48075; IMvigor210 and GSE176307 (immunotherapy); GSE169455 (neoadjuvant chemotherapy) |
| `R/03_scRNA.R` | Single-cell analysis of GSE135337 and GSE222315: quality control, Harmony integration, annotation, copyKAT, UCell scoring, pseudobulk differential expression, Hallmark analysis and CellChat |
| `R/03b_sc_cellcycle.R` | R-loop-high versus R-loop-low epithelial cells and the cell cycle (S/G2M fractions; Hallmark comparison within G1-phase cells) |
| `R/03c_copykat_normal_control.R` | copyKAT negative control in epithelial cells from normal samples |
| `R/04_bulk_score.R` | ssGSEA scores in the eight bulk cohorts; associations with clinical variables, molecular subtypes, mutations, APOBEC enrichment and survival |
| `R/04b_proliferation_overlap.R` | RL_core versus proliferation (gene sets without shared genes, gene-level correlations, replication across cohorts); RHI adjusted for EMT |
| `R/05_candidates.R` | Candidate genes: intersection of DESeq2, WGCNA and single-cell differential expression results, then univariable Cox regression |
| `R/06_ml_signature.R` | Machine-learning model combinations (Mime1), selection of RLRS on the two tuning cohorts and comparison with published signatures |
| `R/06b_specificity.R` | Random gene sets of the same size; adjustment for proliferation and stage |
| `R/07_clinical_nomogram.R` | Multivariable Cox models, nomogram, calibration and decision curves |
| `R/07b_nomogram_validation.R` | Bootstrap validation of the nomogram (C-index corrected for optimism) |
| `R/08_tme_therapy.R` | Immune infiltration, immune checkpoint genes, TIDE, immunotherapy and chemotherapy cohorts, oncoPredict |
| `R/08b_response_vs_proliferation.R` | Treatment response in relation to the proliferation score and tumour mutational burden |
| `R/09_hub_gene_spatial.R` | SHAP values, model genes in TCGA-BLCA, single-cell localisation and spatial transcriptomics (GSE319536, 22 Visium sections) |
| `R/09b_spatial_check.R` | Spot-level checks for sequencing depth: partial correlations adjusted for log10 UMI count and proliferation |
| `R/09c_rlrs_gene_check.R` | Cell-cycle gene-set membership of the RLRS genes, their correlation with proliferation and their coverage in each cohort |
| `R/10_tables.R` | Table 1 and the supplementary tables |
| `R/11_figures.R` | Main and supplementary figures |
| `R/compat_rfsrc_var_select.R`, `R/vendor/` | Restore `randomForestSRC::var.select()`, which Mime1 calls but randomForestSRC ≥ 3.4.0 no longer provides |
| `R/tide_run.py` | Runs TIDE locally with tidepy |
| `run_step.sh` | Runs one or more steps and writes a log for each to `logs/` |
| `env/` | conda-forge package list and compiler settings for the R environment |
| `data/geneset/` | R-loop gene sets derived from R-loopBase, and the curated table of published bladder cancer signatures (genes, coefficients and sources) |

## Software environment

The analyses were run with R 4.5.3 from conda-forge on macOS (Apple silicon, 8 GB of memory). Package versions are listed in Table S10 of the manuscript.

```bash
micromamba create -y -p ./conda_env --override-channels -c conda-forge --file env/conda_cran_packages.txt
micromamba install -y -p ./conda_env -c conda-forge r-magick r-locfit
cp env/Makevars.site conda_env/lib/R/etc/Makevars.site   # C17/C23 compiler settings (osx-arm64 compiler names)
micromamba run -p ./conda_env Rscript R/00_install.R
```

`env/Makevars.site` fills in the C17/C23 compiler variables that conda-forge R leaves empty; on other platforms, replace the compiler name with the value of `CC` in `conda_env/lib/R/etc/Makeconf`.

GitHub packages are installed from fixed commits: Mime1 (l-magnificence/Mime, 9a9f6ac), IOBR (IOBR/IOBR, 4a8f4f2), CellChat (jinworks/CellChat, 75253cd), copykat (navinlabcode/copykat, 01924c0), consensusMIBC (cit-bioinfo/consensusMIBC, 8f69c35), CoxBoost (binderh/CoxBoost, 50d5fff), fastAdaboost (souravc83/fastAdaboost, f331ff8) and presto (immunogenomics/presto, b5df6ee). On first use, msigdbr downloads MSigDB 2026.1.Hs from zenodo.org.

TIDE is computed with tidepy, the TIDE authors' Python implementation, in a separate environment:

```bash
micromamba create -y -p ./conda_py --override-channels -c conda-forge python=3.11 pandas numpy pip
./conda_py/bin/python -m pip install tidepy==1.3.9
```

## Running the analysis

Run every script from the repository root, one step at a time, in this order:

01 → 02 → 03 → 03b → 03c → 04 → 04b → 05 → 06 → 06b → 07 → 07b → 08 → 08b → 09 → 09b → 09c → 10 → 11

```bash
bash run_step.sh 01_geneset
bash run_step.sh 02_download_bulk
```

Intermediate files are written to `data/processed/` and `results/<step>/`, so each step can be rerun on its own. The full analysis needs about 30 GB of disk space. With 8 GB of memory, run only one step at a time. Single-cell samples are downsampled to at most 1,500 cells each (`sc_max_cells_per_sample` in `R/00_config.R`); changing this value changes the results.

Before step 08, download `GDSC2_Expr (RMA Normalized and Log Transformed).rds` and `GDSC2_Res.rds` from https://osf.io/c6tfx/ into `data/raw/oncoPredict/`. All other data are downloaded by the scripts: TCGA-BLCA from the Genomic Data Commons (TCGAbiolinks), with survival end points from the TCGA Pan-Cancer Clinical Data Resource; GEO series from the NCBI GEO FTP site; and IMvigor210 from the IMvigor210CoreBiologies package released by the original investigators.

## Notes

- Comments and log messages in the scripts are partly in Chinese. Lines marked `【核对】` ("check") are manual checkpoints, for example the expected number of events in each GEO cohort.
- Mime1's gbm models and RCTD were run on a single core.

## Licence

The code is released under the MIT licence (see `LICENSE`), except `R/vendor/var.select.rfsrc_3.3.1.R`, an unmodified file from randomForestSRC 3.3.1 that remains under the GNU General Public License, version 3 or later (see `R/vendor/README.md`). The curated table of published signatures (`data/geneset/published_BLCA_signatures*.csv`) is released under CC BY 4.0; please also cite the original publications listed in it. The R-loop gene sets are derived from R-loopBase (https://rloopbase.nju.edu.cn); please cite R-loopBase when using them.

## Citation

Please cite the manuscript above; the full reference will be added after publication.

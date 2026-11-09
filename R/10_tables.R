# =============================================================================
# 10_tables.R —— 论文表格：Table 1（队列基线特征）与补充表汇总（一个 Excel 工作簿，每张表一页）
# 输入：02 步各队列 .rds 及 03–09 步的结果表；输出：results/10_tables/
# =============================================================================
source("R/00_config.R")
need(c("survival", "openxlsx"))
`%||%` <- function(a, b) if (is.null(a)) b else a
STEP <- "10_tables"; OUT <- res_dir(STEP)

# -----------------------------------------------------------------------------
# Table 1：只统计进入生存分析（或疗效分析）的患者
# -----------------------------------------------------------------------------
META <- data.frame(
  cohort   = c("TCGA_BLCA", "GSE13507", "GSE31684", "GSE32894", "GSE48075", "IMvigor210", "GSE176307", "GSE169455"),
  label    = c("TCGA-BLCA", "GSE13507", "GSE31684", "GSE32894", "GSE48075", "IMvigor210", "GSE176307", "GSE169455"),
  platform = c("RNA-seq", "Illumina HumanWG-6 v2 (GPL6102)", "Affymetrix HG-U133 Plus 2.0 (GPL570)",
               "Illumina HumanHT-12 v3 (GPL6947)", "Illumina HumanHT-12 v3 (GPL6947)", "RNA-seq", "RNA-seq",
               "Affymetrix HuGene 1.0 ST (GPL6244)"),
  setting  = c("Muscle-invasive, radical cystectomy", "Primary bladder cancer (NMIBC and MIBC)",
               "Radical cystectomy", "Primary bladder cancer (NMIBC and MIBC)", "Primary bladder cancer",
               "Metastatic urothelial carcinoma, atezolizumab", "Metastatic urothelial carcinoma, PD-1/PD-L1 blockade",
               "MIBC, neoadjuvant or induction chemotherapy"),
  role     = c("Training", "Tuning", "Tuning", "Test", "Test", "Test (immunotherapy)", "Test (immunotherapy)",
               "Test (chemotherapy response)"),
  endpoint = c("OS", "OS", "OS", "DSS", "OS", "OS; objective response", "OS; objective response",
               "Pathological response"))

fmt_med <- function(x) { x <- suppressWarnings(as.numeric(x)); x <- x[!is.na(x)]
  if (!length(x)) return("NA"); q <- stats::quantile(x, c(0.5, 0.25, 0.75)); sprintf("%.0f (%.0f-%.0f)", q[1], q[2], q[3]) }
fmt_n <- function(k, n) if (is.na(k) || n == 0) "NA" else sprintf("%d (%.0f%%)", k, 100 * k / n)
t_invasive <- function(cl, n) {
  if (n == "TCGA_BLCA") return(NA)                                   # TCGA-BLCA 全部为肌层浸润（入组标准）
  tt <- if ("T" %in% names(cl)) as.character(cl$T) else if ("cT" %in% names(cl)) as.character(cl$cT) else NULL
  if (is.null(tt)) return(NA)
  ok <- !is.na(tt) & !tt %in% c("X", "x", "")
  if (!any(ok)) return(NA)
  c(k = sum(tt[ok] %in% c("2", "3", "4")), n = sum(ok))
}
rows <- lapply(seq_len(nrow(META)), function(i) {
  n <- META$cohort[i]; cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  if (n == "TCGA_BLCA") cl <- cl[cl$group == "Tumor", ]
  if (n == "GSE169455") cl <- cl[!is.na(cl$pCR) | !is.na(cl$response), ]
  else cl <- cl[!is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0, ]
  N <- nrow(cl)
  male <- if ("gender" %in% names(cl)) { g <- tolower(substr(as.character(cl$gender), 1, 1)); fmt_n(sum(g == "m", na.rm = TRUE), sum(!is.na(g))) } else "NA"
  age <- if ("age" %in% names(cl)) fmt_med(cl$age) else "NA"
  inv <- t_invasive(cl, n)
  inv_txt <- if (n == "TCGA_BLCA") sprintf("%d (100%%)", N) else if (length(inv) == 2) fmt_n(inv[["k"]], inv[["n"]]) else "NA"
  hg <- if ("grade" %in% names(cl)) {
    g <- tolower(as.character(cl$grade)); hi <- g %in% c("high", "g3", "g4"); known <- !is.na(g) & !g %in% c("gx", "")
    fmt_n(sum(hi & known), sum(known)) } else "NA"
  ev <- if (n == "GSE169455") NA else sum(cl$OS == 1)
  fu <- if (n == "GSE169455") NA else {
    km <- survival::survfit(survival::Surv(OS.time, 1 - OS) ~ 1, data = cl)
    as.numeric(summary(km)$table["median"]) / 30.4375 }
  resp <- if ("response" %in% names(cl) && n != "GSE169455") {
    r <- cl$response; fmt_n(sum(r %in% c("CR/PR", "responder", "R", "CR", "PR", 1, "1"), na.rm = TRUE), sum(!is.na(r))) }
    else if (n == "GSE169455") fmt_n(sum(cl$pCR == 1, na.rm = TRUE), sum(!is.na(cl$pCR))) else "NA"
  data.frame(Cohort = META$label[i], Platform = META$platform[i], Setting = META$setting[i], Role = META$role[i],
             Patients = N, `Age, median (IQR)` = age, `Male, n (%)` = male, `Muscle-invasive, n (%)` = inv_txt,
             `High grade, n (%)` = hg, `End point` = META$endpoint[i],
             Events = ifelse(is.na(ev), "NA", as.character(ev)),
             `Median follow-up, months` = ifelse(is.na(fu), "NA", sprintf("%.1f", fu)),
             `Responders, n (%)` = resp, check.names = FALSE)
})
tab1 <- do.call(rbind, rows)
print(tab1, row.names = FALSE)
utils::write.csv(tab1, file.path(OUT, "Table1_cohorts.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 补充表：把各步结果表汇总进一个工作簿（存在才写）
# -----------------------------------------------------------------------------
# 编号按正文首次引用的顺序：S1 基因集 → S2 单细胞 → S3 bulk 评分 → S4 候选基因 → S5 建模 → S6 特异性
#   → S7 临床 → S8 微环境与治疗 → S9 模型基因与空间 → S10 软件版本
R1 <- function(...) file.path(DIR_RES, ...)
SUPP <- list(
  "S1a R-loop regulators"     = R1("01_geneset", "Table_S1_Rloop_regulators_v2.csv"),
  "S1b Gene sets"             = file.path(DIR_GS, "Rloop_genesets_table.csv"),
  "S1c Gene set coverage"     = R1("04_bulk_score", "geneset_coverage_by_cohort.csv"),
  "S2a scRNA samples"         = R1("03_scRNA", "scRNA_samples.csv"),
  "S2b Cluster annotation"    = R1("03_scRNA", "cluster_annotation_auto.csv"),
  "S2c copyKAT reference"     = R1("03_scRNA", "copykat_reference_cells.csv"),
  "S2d copyKAT by sample"     = R1("03_scRNA", "copykat_by_sample.csv"),
  "S2e copyKAT normal control" = R1("03_scRNA", "copykat_normal_control.csv"),
  "S2f Score per sample"      = R1("03_scRNA", "Rloop_score_per_sample_class.csv"),
  "S2g Pseudobulk DE"         = R1("03_scRNA", "pseudobulk_DE_High_vs_Low.csv"),
  "S2h Hallmark paired"       = R1("03_scRNA", "Hallmark_High_vs_Low_paired.csv"),
  "S2i Hallmark G1 only"      = R1("03_scRNA", "Hallmark_all_vs_G1only.csv"),
  "S2j Cell cycle by sample"  = R1("03_scRNA", "cellcycle_by_sample.csv"),
  "S2k CellChat"              = R1("03_scRNA", "CellChat_count_weight_from_malignant.csv"),
  "S3a Spearman TCGA"         = R1("04_bulk_score", "signature_spearman_TCGA.csv"),
  "S3b Prolif by cohort"      = R1("04_bulk_score", "proliferation_correlation_by_cohort.csv"),
  "S3c Prolif overlap"        = R1("04_bulk_score", "proliferation_overlap_summary.csv"),
  "S3d Per-gene rho Prolif"   = R1("04_bulk_score", "proliferation_overlap_per_gene_rho.csv"),
  "S3e Mutation frequency"    = R1("04_bulk_score", "TCGA_mutation_frequency.csv"),
  "S3f Score by mutation"     = R1("04_bulk_score", "Rloop_score_by_mutation.csv"),
  "S3g TCGA Cox"              = R1("04_bulk_score", "TCGA_Cox_Rloop_scores.csv"),
  "S3h External Cox RL_core"  = R1("04_bulk_score", "External_Cox_RL_core.csv"),
  "S3i RHI adjusted for EMT"  = R1("04_bulk_score", "RHI_EMT_adjusted_Cox_by_cohort.csv"),
  "S4 Candidate genes"        = R1("05_candidates", "Table_S3_candidates_unicox.csv"),
  "S5a Model combinations"    = R1("06_ml_signature", "Table_S4_Cindex_101_models.csv"),   # 实际 117 个组合
  "S5b RLRS coefficients"     = R1("09_hub_gene_spatial", "RLRS_coefficients.csv"),        # 线性核 SVM，09 步线性还原
  "S5c RLRS by cohort"        = R1("06_ml_signature", "RLRS_Cox_by_cohort.csv"),
  "S5d RLRS timeROC"          = R1("06_ml_signature", "RLRS_timeROC_AUC.csv"),
  "S5e RLRS genes prolif"     = R1("09_hub_gene_spatial", "RLRS_genes_proliferation.csv"),  # 09c
  "S5f RLRS gene coverage"    = R1("09_hub_gene_spatial", "RLRS_gene_coverage.csv"),        # 09c
  "S6a Random gene sets"      = R1("06b_specificity", "random_signature_benchmark.csv"),
  "S6b Random by cohort"      = R1("06b_specificity", "random_signature_benchmark_by_cohort.csv"),
  "S6c Proliferation-adjusted" = R1("06b_specificity", "Table_S6_proliferation_adjusted.csv"),
  "S6d Published signatures"  = R1("06_ml_signature", "Table_S5_published_models_Cindex.csv"),
  "S6e Signature sources"     = file.path(DIR_GS, "published_BLCA_signatures_sources.csv"),
  "S6f RLRS rank"             = R1("06_ml_signature", "RLRS_rank_among_published.csv"),
  "S7a TCGA Cox"              = R1("07_clinical_nomogram", "TCGA_uni_multi_Cox.csv"),
  "S7b External multivariable" = R1("07_clinical_nomogram", "External_multivariate_Cox.csv"),
  "S7c C-index nomogram"      = R1("07_clinical_nomogram", "Cindex_nomogram_RLRS_stage.csv"),
  "S7d Nomogram validation"   = R1("07_clinical_nomogram", "Nomogram_optimism_corrected.csv"),  # 07b
  "S8a Immune infiltration"   = R1("08_tme_therapy", "TME_High_vs_Low_all_methods.csv"),
  "S8b Checkpoints"           = R1("08_tme_therapy", "Checkpoints_High_vs_Low.csv"),
  "S8c TIDE"                  = R1("08_tme_therapy", "TIDE_High_vs_Low.csv"),
  "S8d Immunotherapy"         = R1("08_tme_therapy", "ICI_response_summary.csv"),
  "S8e Chemotherapy"          = R1("08_tme_therapy", "Chemo_response_GSE169455.csv"),
  "S8f Drug sensitivity"      = R1("08_tme_therapy", "oncoPredict_High_vs_Low.csv"),
  "S8g Response vs prolif"    = R1("08_tme_therapy", "Response_vs_proliferation.csv"),      # 08b
  "S9a SHAP ranking"          = R1("09_hub_gene_spatial", "SHAP_gene_ranking.csv"),
  "S9b Model genes scRNA"     = R1("09_hub_gene_spatial", "hub_scRNA_by_cell_class.csv"),
  "S9c Spatial programmes"    = R1("09_hub_gene_spatial", "Spatial_RLcore_vs_programs_by_slice.csv"),
  "S9d Spatial model genes"   = R1("09_hub_gene_spatial", "Spatial_RLcore_vs_hub_spearman.csv"),
  "S9e Spatial depth check"   = R1("09_hub_gene_spatial", "Spatial_depth_check_by_slice.csv"),
  "S9f Spatial genes adjusted" = R1("09_hub_gene_spatial", "Spatial_hub_partial_by_slice.csv"))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table 1"); openxlsx::writeData(wb, "Table 1", tab1)
for (nm in names(SUPP)) if (file.exists(SUPP[[nm]])) {
  openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, utils::read.csv(SUPP[[nm]], check.names = FALSE))
} else log_msg("（暂缺）", nm, "：", SUPP[[nm]])
# S10：关键 R 包版本（GitHub 包附提交号）
GH_SHA <- c(Mime1 = "l-magnificence/Mime@9a9f6ac", IOBR = "IOBR/IOBR@4a8f4f2", CellChat = "jinworks/CellChat@75253cd",
            copykat = "navinlabcode/copykat@01924c0", consensusMIBC = "cit-bioinfo/consensusMIBC@8f69c35",
            CoxBoost = "binderh/CoxBoost@50d5fff", presto = "immunogenomics/presto@b5df6ee",
            fastAdaboost = "souravc83/fastAdaboost@f331ff8")   # 从固定提交的源码包安装（见 R/00_install.R）
pk <- c("Seurat", "harmony", "copykat", "UCell", "CellChat", "edgeR", "DESeq2", "WGCNA", "GSVA", "msigdbr", "Mime1", "glmnet",
        "randomForestSRC", "survival", "timeROC", "rms", "dcurves", "IOBR", "oncoPredict", "maftools", "consensusMIBC",
        "kernelshap", "TCGAbiolinks", "GEOquery", "HGNChelper")
S10 <- data.frame(package = pk, version = vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                                                        error = function(e) NA_character_), ""),
                 source = vapply(pk, function(p) { d <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
                   if (p %in% names(GH_SHA)) paste0("GitHub ", GH_SHA[[p]])
                   else if (is.null(d)) NA_character_ else if (!is.null(d$RemoteSha)) paste0("GitHub ", d$RemoteRepo %||% "", "@", substr(d$RemoteSha, 1, 7))
                   else if (!is.null(d$biocViews)) "Bioconductor 3.22" else "CRAN" }, ""))
S10 <- rbind(data.frame(package = "R", version = paste(R.version$major, R.version$minor, sep = "."), source = "conda-forge"), S10)
utils::write.csv(S10, file.path(OUT, "Table_S10_package_versions.csv"), row.names = FALSE)
openxlsx::addWorksheet(wb, "S10 Package versions"); openxlsx::writeData(wb, "S10 Package versions", S10)
openxlsx::saveWorkbook(wb, file.path(OUT, "Tables_all.xlsx"), overwrite = TRUE)
save_session(STEP)
log_msg("10 完成")

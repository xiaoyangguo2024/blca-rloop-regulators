# =============================================================================
# 00_install.R —— 安装全部依赖（只需运行一次，约 30–60 分钟）
# 已安装的包会自动跳过；某个包装不上时不会中断，最后统一列出失败的包
# =============================================================================
# 国内网络建议取消下面两行注释，使用清华镜像
# options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(timeout = 3600)
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@")
  options(repos = c(CRAN = "https://cloud.r-project.org"))
# 并行编译；某些沙箱里 detectCores() 返回 NA，此时用 4
local({ nc <- suppressWarnings(parallel::detectCores()); options(Ncpus = if (is.na(nc)) 4L else max(1L, min(4L, nc - 1L))) })

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# 让 install.packages 同时看到 CRAN 和 Bioconductor 仓库（oncoPredict 等 CRAN 包依赖 Bioconductor 包）
options(repos = BiocManager::repositories())

failed <- character(0)
has <- function(p) requireNamespace(p, quietly = TRUE)

inst_cran <- function(pkgs) for (p in pkgs) if (!has(p)) {
  message(">>> CRAN: ", p)
  tryCatch(install.packages(p), error = function(e) NULL)
  if (!has(p)) failed <<- c(failed, p)
}
inst_bioc <- function(pkgs) for (p in pkgs) if (!has(p)) {
  message(">>> Bioconductor: ", p)
  tryCatch(BiocManager::install(p, update = FALSE, ask = FALSE), error = function(e) NULL)
  if (!has(p)) failed <<- c(failed, p)
}
# GitHub 包固定到具体提交（论文方法部分报告这些版本）；用 codeload 下载源码包，不调用 GitHub API，也不需要令牌
inst_gh <- function(repos) for (r in names(repos)) {
  p <- repos[[r]]$pkg
  if (!has(p)) {
    message(">>> GitHub: ", r, "@", substr(repos[[r]]$sha, 1, 7))
    url <- sprintf("https://codeload.github.com/%s/tar.gz/%s", r, repos[[r]]$sha)
    # 只装必需依赖（Depends/Imports/LinkingTo）；IOBR 的 Suggests 很多，用到的那几个在 CRAN 列表里单独装
    tryCatch(remotes::install_url(url, upgrade = "never", dependencies = NA), error = function(e) message(conditionMessage(e)))
    if (!has(p)) failed <<- c(failed, p)
  }
}

# ---- Bioconductor（先装，CRAN 上的 WGCNA 等依赖它们）------------------------
inst_bioc(c(
  "impute", "preprocessCore", "GO.db", "AnnotationDbi", "org.Hs.eg.db",
  "Biobase", "SummarizedExperiment", "GEOquery", "TCGAbiolinks",
  "DESeq2", "edgeR", "limma", "sva",
  "GSVA", "GSEABase", "UCell", "AUCell", "clusterProfiler",
  "maftools", "ComplexHeatmap",
  "cancerclass", "mixOmics", "sparrow",          # Mime1 依赖
  "BiocNeighbors", "glmGamPoi"
))

# ---- CRAN -------------------------------------------------------------------
inst_cran(c(
  # 通用
  "data.table", "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggpubr", "ggrepel",
  "patchwork", "cowplot", "pheatmap", "circlize", "RColorBrewer", "ggsci", "ggbreak", "aplot",
  "Matrix", "jsonlite", "R.utils", "HGNChelper", "msigdbr", "openxlsx", "VennDiagram", "UpSetR",
  # 生存与临床
  "survival", "survminer", "timeROC", "rms", "dcurves", "pROC", "forestplot",
  # 单细胞与空间
  "Seurat", "harmony", "hdf5r", "NMF", "ggalluvial", "scTenifoldKnk",
  # 建模与解释
  "WGCNA", "glmnet", "gbm", "randomForestSRC", "survivalsvm", "plsRcox", "superpc",
  "BART", "miscTools", "compareC", "kernelshap", "shapviz",
  # 数据获取与药敏
  "UCSCXenaTools", "oncoPredict",
  # IOBR 去卷积和 copykat 实际用到的可选依赖
  "e1071", "limSolve", "MASS", "pec", "parallelDist", "mixtools", "MCMCpack", "dlm", "fmsb", "kableExtra"
))

# ---- GitHub -----------------------------------------------------------------
inst_gh(list(
  "binderh/CoxBoost"          = list(pkg = "CoxBoost",      sha = "50d5fff7d7eebc4a1174dfa32adcfb4b4a69ad13"),
  "souravc83/fastAdaboost"    = list(pkg = "fastAdaboost",  sha = "f331ff8ccfe2e7f318448a29bc72af646dfaf2c1"),
  "l-magnificence/Mime"       = list(pkg = "Mime1",         sha = "9a9f6ac89851bf631f9df3868b2fa624bed49df2"),  # 仓库叫 Mime，包名叫 Mime1
  "IOBR/IOBR"                 = list(pkg = "IOBR",          sha = "4a8f4f290b73bdb0d3e6cc6de6843500b014a962"),
  "jinworks/CellChat"         = list(pkg = "CellChat",      sha = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f"),
  "navinlabcode/copykat"      = list(pkg = "copykat",       sha = "01924c0999bb660064da8d3e1561615ac6efb76f"),
  "cit-bioinfo/consensusMIBC" = list(pkg = "consensusMIBC", sha = "8f69c3553fd53ddeaa532489ccf7c169fac5ab6d"),
  "immunogenomics/presto"     = list(pkg = "presto",        sha = "b5df6ee6097eb62522f2e557aa93ac21bab2d05f")  # 可选：加速 Wilcoxon
))

# ---- 可选：较大的包，用到时再装 ---------------------------------------------
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")   # 04 步计算 APOBEC 突变特征（约 800 MB）
# remotes::install_github("dmcable/spacexr")            # 09 步 RCTD 空间反卷积

if (length(failed)) {
  message("\n以下包未能安装，请手动处理（常见原因：网络、缺系统库如 hdf5/gfortran）：\n  ",
          paste(unique(failed), collapse = ", "))
} else {
  message("\n全部依赖已安装。")
}

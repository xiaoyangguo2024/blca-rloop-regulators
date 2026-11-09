# =============================================================================
# 00_config.R —— 全局路径、参数与通用函数
# 其他脚本开头都会 source("R/00_config.R")；必须在项目根目录（BLCA_Rloop/）下运行
# =============================================================================
options(stringsAsFactors = FALSE, timeout = 3600, scipen = 999)

# 核数：某些沙箱/容器里 parallel::detectCores() 返回 NA，会让部分包（future/parallelly）报错
N_CORES <- local({ nc <- suppressWarnings(parallel::detectCores()); if (is.na(nc)) 4L else max(1L, min(4L, nc - 1L)) })
options(mc.cores = N_CORES)
# IOBR 的参考数据（LM22、xCell、EPIC、TIMER、quanTIseq、ESTIMATE）首次使用时从 GitHub 下载；缓存放在项目里，避免每次重下
options(IOBR.cache_dir = file.path(normalizePath("."), "data", "raw", "IOBR_cache"))
# 不能打开本地端口的环境（沙箱、部分服务器）里 BiocParallel 注册默认并行参数会失败，sva::ComBat（oncoPredict 调用）随之报错；
# 显式注册串行参数（只影响未指定 BPPARAM 的调用，结果不变）
if (requireNamespace("BiocParallel", quietly = TRUE))
  try(suppressWarnings(suppressMessages(BiocParallel::register(BiocParallel::SerialParam()))), silent = TRUE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_FALLBACK = N_CORES, OMP_NUM_THREADS = N_CORES)
# 走带认证的 HTTP 代理时（如某些公司网络或沙箱），R 的 curl 包默认先不带凭据试探认证方式，
# 部分代理会直接断开，报 "Proxy CONNECT aborted"。这里让 curl 包新建的连接一律直接发送 Basic 代理认证
# （httr、xml2、GEOquery、UCSCXenaTools 都经由它联网）。没有带认证的代理时这段不起作用
if (grepl("@", Sys.getenv("https_proxy", Sys.getenv("HTTPS_PROXY"))) && requireNamespace("curl", quietly = TRUE)) {
  local({
    nh <- curl::new_handle
    utils::assignInNamespace("new_handle", function(...) { h <- nh(...); curl::handle_setopt(h, proxyauth = 1L); h }, "curl")
  })
  if (requireNamespace("httr", quietly = TRUE)) httr::set_config(httr::config(proxyauth = 1L))
}

if (!file.exists("BLCA_Rloop.Rproj"))
  stop("请先把工作目录设到项目根目录 BLCA_Rloop/（或用 RStudio 打开 BLCA_Rloop.Rproj）")

PROJ    <- normalizePath(".")
DIR_RAW <- file.path(PROJ, "data", "raw")
DIR_PRO <- file.path(PROJ, "data", "processed")
DIR_GS  <- file.path(PROJ, "data", "geneset")
DIR_RES <- file.path(PROJ, "results")
for (d in c(DIR_RAW, DIR_PRO, DIR_GS, DIR_RES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 每个步骤一个结果文件夹
res_dir <- function(step) {
  d <- file.path(DIR_RES, step)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# ----------------------------- 全局参数 --------------------------------------
PARAM <- list(
  seed            = 20261107,
  geneset_tiers   = 1:3,      # 主分析基因集：R-loopBase v2.0 Tier 1–3
  # 单细胞
  sc_min_features = 200,
  sc_max_features = 7000,
  sc_max_mt       = 20,
  sc_max_cells_per_sample = 1500,  # 每个样本随机抽取的细胞上限（读表时即抽样）。8 GB 内存用 1500；32 GB 可设 6000
  sc_npcs         = 30,
  sc_resolution   = 0.8,
  sc_min_cells_pb = 20,       # 伪 bulk：每个"样本×组别"至少多少个细胞才纳入
  sc_deg_lfc      = 0.5,      # 伪 bulk 差异阈值（log2FC）
  sc_deg_fdr      = 0.05,
  copykat_cores   = 2L,       # copyKAT 每次只跑一个样本（并行跑多个样本曾在 8 GB 内存的机器上导致整机卡死），每个样本用几个核
  # bulk
  deg_lfc         = 1,
  deg_padj        = 0.05,
  wgcna_top_mad   = 5000,
  wgcna_mm        = 0.6,
  wgcna_gs        = 0.3,
  min_candidates  = 20,
  unicox_p        = 0.05,
  surv_years      = c(1, 3, 5)
)
set.seed(PARAM$seed)

# 细胞注释用的经典 marker（自动注释后务必看 DotPlot 人工核对）
MARKERS <- list(
  Epithelial  = c("EPCAM", "KRT8", "KRT18", "KRT19", "UPK1B", "UPK2", "KRT5", "GATA3"),
  T_cell      = c("CD3D", "CD3E", "CD2", "CD8A", "CD4", "IL7R"),
  NK          = c("NKG7", "GNLY", "KLRD1", "KLRF1"),
  B_cell      = c("CD79A", "MS4A1", "CD19"),
  Plasma      = c("MZB1", "JCHAIN", "IGHG1"),
  Myeloid     = c("CD68", "CD14", "LYZ", "C1QA", "CD163"),
  DC          = c("CD1C", "CLEC9A", "LAMP3", "FCER1A"),
  Mast        = c("TPSAB1", "CPA3", "KIT"),
  Fibroblast  = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  Myofib_SMC  = c("ACTA2", "RGS5", "MYH11", "TAGLN"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP")
)

# Ayers 2017 T 细胞炎症基因表达谱（GEP18）
GEP18 <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
           "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")

CHECKPOINTS <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
                 "SIGLEC15", "CD276", "VSIR", "BTLA")

COL2 <- c(High = "#D7301F", Low = "#2C7FB8")

# ----------------------------- 通用函数 --------------------------------------
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# 检查并加载 R 包；缺包时提示先运行 00_install.R
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("缺少 R 包：", paste(miss, collapse = ", "), "；请先运行 R/00_install.R")
  invisible(suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE)))
}

save_session <- function(step) {
  writeLines(capture.output(sessionInfo()), file.path(res_dir(step), "sessionInfo.txt"))
}

save_gg <- function(p, file, w = 6, h = 5) ggplot2::ggsave(file, p, width = w, height = h, limitsize = FALSE)

# 打印类图对象（survminer / ComplexHeatmap / base）到 PDF
save_print <- function(obj, file, w = 7, h = 6) {
  grDevices::pdf(file, width = w, height = h, onefile = TRUE)
  on.exit(grDevices::dev.off())
  if (is.list(obj) && !inherits(obj, c("gg", "ggsurvplot"))) for (o in obj) print(o) else print(obj)
}

# 按 GEO2R 的规则判断是否需要 log2；需要时做 log2(x+1)
auto_log2 <- function(m) {
  qx <- as.numeric(stats::quantile(m, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  need_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (need_log) {
    m[m < 0] <- 0
    m <- log2(m + 1)
    log_msg("检测到非 log 数据，已做 log2(x+1)")
  }
  m
}

# 从 GPL 注释表中找基因符号列，返回 data.frame(probe, symbol)
gpl_symbol_map <- function(fd) {
  cand <- c("Symbol", "ILMN_Gene", "Gene Symbol", "GENE_SYMBOL", "Gene symbol", "gene_assignment")
  hit <- cand[cand %in% colnames(fd)][1]
  if (is.na(hit)) stop("GPL 注释中找不到基因符号列，现有列：", paste(colnames(fd), collapse = ", "))
  sym <- as.character(fd[[hit]])
  if (hit == "gene_assignment")   # Affymetrix ST 芯片："NM_xxx // SYMBOL // ..."
    sym <- vapply(strsplit(sym, " // ", fixed = TRUE),
                  function(x) if (length(x) >= 2) trimws(x[2]) else NA_character_, character(1))
  data.frame(probe = as.character(fd$ID), symbol = trimws(sub("\\s*///.*$", "", sym)))
}

# 探针/行 → 基因符号；同一基因多行时保留平均表达最高的一行
collapse_to_symbol <- function(m, symbols) {
  ok <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  m <- m[ok, , drop = FALSE]; symbols <- symbols[ok]
  ord <- order(rowMeans(m, na.rm = TRUE), decreasing = TRUE)
  m <- m[ord, , drop = FALSE]; symbols <- symbols[ord]
  keep <- !duplicated(symbols)
  out <- m[keep, , drop = FALSE]
  rownames(out) <- symbols[keep]
  out
}

probe_to_symbol <- function(expr, p2s) {
  p2s <- p2s[p2s$probe %in% rownames(expr), ]
  collapse_to_symbol(expr[p2s$probe, , drop = FALSE], p2s$symbol)
}

# GEOquery 的 pData 中，"键:ch1" 列已按键名对齐（比 characteristics_ch1.x 可靠）
geo_pheno <- function(eset) {
  pd <- Biobase::pData(eset)
  ch <- grep(":ch1$", colnames(pd), value = TRUE)
  out <- pd[, c("geo_accession", "title", "source_name_ch1", ch), drop = FALSE]
  colnames(out) <- sub(":ch1$", "", colnames(out))
  out
}

MONTH2DAY <- 30.4375

# 每个基因在队列内部做 z-score（跨平台建模必需）
zscore_rows <- function(m) {
  m <- as.matrix(m)
  s <- t(scale(t(m)))
  s[is.na(s)] <- 0
  s
}

# 基因名中的 "-" 等字符在公式里会出问题。Mime1 内部会把 "-" 和 "_" 都改成 "."，
# 这里统一提前改成 "."，保证候选基因名与 Mime1 结果里的基因名一致（如 HLA-DRB1 → HLA.DRB1）
safe_names <- function(x) gsub("[^A-Za-z0-9.]", ".", x)

# 生成 Mime1 输入：ID, OS.time(天), OS(0/1), 各基因（默认队列内 z-score）
# fill_missing = TRUE 时，队列缺失的基因补 0（z-score 时即均值），并给出提示
# scale = FALSE 时保留原始 log2 表达（用于按原文系数计算已发表模型的评分）
build_mime_input <- function(expr, surv, genes, fill_missing = FALSE, scale = TRUE) {
  ids <- intersect(colnames(expr), surv$ID)
  g_have <- intersect(genes, rownames(expr))
  z <- if (scale) zscore_rows(expr[g_have, ids, drop = FALSE]) else as.matrix(expr[g_have, ids, drop = FALSE])
  miss <- setdiff(genes, g_have)
  if (length(miss)) {
    if (!fill_missing) stop("队列缺少基因：", paste(head(miss, 20), collapse = ", "))
    log_msg("提示：缺失 ", length(miss), " 个基因，已补 0：", paste(head(miss, 10), collapse = ", "))
    z <- rbind(z, matrix(0, nrow = length(miss), ncol = ncol(z), dimnames = list(miss, colnames(z))))
  }
  z <- z[genes, , drop = FALSE]
  s <- surv[match(ids, surv$ID), c("ID", "OS.time", "OS")]
  df <- data.frame(s, t(z), check.names = FALSE)
  colnames(df)[-(1:3)] <- safe_names(genes)
  df <- df[!is.na(df$OS.time) & !is.na(df$OS) & df$OS.time > 0, ]
  rownames(df) <- NULL
  df
}

# 从 Mime1 的模型对象中提取最终用到的基因（组合模型 "A + B" 取 B 拟合时用到的变量）
# 取不到时退回 fallback（Mime1 的候选基因），并提示人工核对
model_genes <- function(fit, fallback, quiet = FALSE) {
  g <- tryCatch({
    if (inherits(fit, "rfsrc")) fit$xvar.names
    else if (inherits(fit, "coxph")) names(stats::coef(fit))
    else if (inherits(fit, "cv.glmnet")) {
      cf <- as.matrix(stats::coef(fit, s = fit$lambda.min)); rownames(cf)[cf[, 1] != 0]
    }
    else if (inherits(fit, "CoxBoost")) { cf <- stats::coef(fit); names(cf)[cf != 0] }
    else if (inherits(fit, "plsRcoxmodel")) colnames(fit$dataX)
    else if (inherits(fit, "survivalsvm")) fit$var.names
    else if (is.list(fit) && inherits(fit[[1]], "gbm")) fit[[1]]$var.names
    else if (is.list(fit) && inherits(fit[[1]], "glmnet")) rownames(fit[[1]]$beta)
    else NULL
  }, error = function(e) NULL)
  if (!length(g)) {
    if (!quiet) log_msg("【核对】无法从模型对象（", class(fit)[1], "）自动提取基因，暂用全部候选基因")
    g <- fallback
  }
  unique(g)
}

# ssGSEA（兼容 GSVA ≥1.50 的新接口与旧接口）
run_ssgsea <- function(expr, gsets) {
  expr <- as.matrix(expr)
  gsets <- lapply(gsets, function(g) intersect(g, rownames(expr)))
  gsets <- gsets[lengths(gsets) >= 5]
  if (utils::packageVersion("GSVA") >= "1.50.0") {
    GSVA::gsva(GSVA::ssgseaParam(expr, gsets), verbose = FALSE)
  } else {
    GSVA::gsva(expr, gsets, method = "ssgsea", verbose = FALSE)
  }
}

# msigdbr（兼容 v10 的 collection 参数与旧版的 category 参数）
get_msig <- function(collection, subcollection = NULL) {
  df <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection),
    error = function(e) {
      if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", category = collection)
      else msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
    })
  split(df$gene_symbol, df$gs_name)
}

# 中位数分组（每个队列内部）
median_group <- function(x) factor(ifelse(x >= stats::median(x, na.rm = TRUE), "High", "Low"),
                                   levels = c("High", "Low"))

# KM 曲线（survminer）
km_plot <- function(df, time = "OS.time", status = "OS", group = "group", title = "",
                    xlab = "Time (days)", palette = unname(COL2)) {
  fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, status, group))
  fit <- survminer::surv_fit(fml, data = df)
  survminer::ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                        title = title, xlab = xlab, palette = palette, legend.title = "")
}

# Cox 回归结果整理成表
cox_table <- function(fit) {
  s <- summary(fit)
  data.frame(term = rownames(s$coefficients),
             HR = s$conf.int[, "exp(coef)"], lower = s$conf.int[, "lower .95"],
             upper = s$conf.int[, "upper .95"], p = s$coefficients[, ncol(s$coefficients)],
             row.names = NULL)
}

# 解压 .gz 到指定文件（不依赖 R.utils）
gunzip_to <- function(src, dest) {
  zz <- gzfile(src, "rb"); on.exit(close(zz))
  out <- file(dest, "wb"); on.exit(close(out), add = TRUE)
  repeat {
    b <- readBin(zz, "raw", n = 4e6)
    if (!length(b)) break
    writeBin(b, out)
  }
  invisible(dest)
}

# 下载 GEO 补充文件（直接走 NCBI FTP 的 https 地址，比 getGEOSuppFiles 更可控）
geo_suppl <- function(gse, file, subdir = gse) {
  dest <- file.path(DIR_RAW, "GEO", subdir)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dest, file)
  if (!file.exists(f)) {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
                   sub("\\d{1,3}$", "nnn", gse), gse, file)
    log_msg("下载 ", url)
    utils::download.file(url, f, mode = "wb")
  }
  f
}

# 读取 GEO series matrix（带 GPL 注释）；一个 GSE 含多个平台时按 gpl 取
get_geo_eset <- function(gse, gpl = NULL) {
  dest <- file.path(DIR_RAW, "GEO"); dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  es <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = dest)
  if (length(es) > 1 && !is.null(gpl)) es <- es[grep(gpl, names(es))]
  es[[1]]
}

# 按字段名从 geo_pheno() 的结果中取一列：先按大小写精确匹配，再不区分大小写
# （GSE13507 同时有 "age"（1 个样本的错误注释）和 "AGE"（165 例）两个字段）
pcol <- function(pd, ...) {
  cand <- c(...)
  hit <- cand[cand %in% colnames(pd)][1]
  if (is.na(hit)) {
    hit <- colnames(pd)[match(tolower(cand), tolower(colnames(pd)))]
    hit <- hit[!is.na(hit)][1]
  }
  if (is.na(hit)) {
    log_msg("【核对】找不到字段：", paste(cand, collapse = " / "))
    return(rep(NA_character_, nrow(pd)))
  }
  trimws(as.character(pd[[hit]]))
}

# 行名若是 Ensembl 或 Entrez ID，转换为基因符号
ensure_symbols <- function(m) {
  ids <- sub("\\..*$", "", rownames(m))
  type <- if (mean(grepl("^ENSG", ids)) > 0.5) "ENSEMBL" else if (mean(grepl("^[0-9]+$", ids)) > 0.5) "ENTREZID" else NA
  if (is.na(type)) return(m)
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ids, keytype = type,
                               column = "SYMBOL", multiVals = "first")
  collapse_to_symbol(m, unname(sym))
}

# 下载 GEO 样本（GSM）级补充文件，避免下载整个 RAW.tar
gsm_suppl <- function(gsm, file, subdir) {
  dest <- file.path(DIR_RAW, "GEO", subdir)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dest, file)
  if (!file.exists(f)) {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s",
                   sub("\\d{1,3}$", "nnn", gsm), gsm, file)
    log_msg("下载 ", url)
    utils::download.file(url, f, mode = "wb")
  }
  f
}

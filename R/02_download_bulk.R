# =============================================================================
# 02_download_bulk.R —— 下载并整理 bulk 队列
#   训练：TCGA-BLCA（STAR counts/TPM + TCGA-CDR 生存 + 突变 MAF）
#   验证：GSE13507、GSE32894、GSE31684、GSE48075
#   治疗：IMvigor210、GSE176307（免疫治疗）；GSE169455（新辅助化疗）
# 输出（data/processed/）：每个队列一个 .rds = list(expr = 基因×样本 log2 矩阵, clin = 临床表)
#   clin 统一包含 ID、OS.time（天）、OS（0/1）；GSE32894 的 "OS" 实为 DSS（见 Table S2）
# 字段名已于 2026-11-07 按 GEO series matrix 逐一核对
# =============================================================================
source("R/00_config.R")
need(c("TCGAbiolinks", "SummarizedExperiment", "GEOquery", "Biobase", "UCSCXenaTools", "data.table"))
STEP <- "02_download_bulk"; OUT <- res_dir(STEP)
GDC_DIR <- file.path(DIR_RAW, "GDCdata")

num <- function(x) suppressWarnings(as.numeric(x))

# "T1N0M0" / "cT2N0M0" / "pT3b" / "Ta" / "pTCisN0" → T 分期（"A" = Ta，"IS" = Tis/CIS）
# 已用 GSE169455、GSE13507、GSE48075、GSE32894、GSE31684 的全部实际取值测试
t_stage <- function(x) {
  x <- sub("TCIS", "TIS", toupper(gsub("\\s", "", x)), fixed = TRUE)
  t <- sub("^[CPY]*T(IS|A|X|[0-4]).*$", "\\1", x)
  t[t == x | !nzchar(x)] <- NA
  t
}
n_stage <- function(x) {            # "+" 表示原文只写了 N+
  x <- toupper(gsub("\\s", "", x))
  n <- sub("^.*N([0-3X+]).*$", "\\1", x)
  n[n == x] <- NA
  n
}

COHORTS <- list()   # 汇总信息，最后写成 Table S2
save_cohort <- function(name, expr, clin, role, platform, endpoint) {
  expr <- as.matrix(expr); storage.mode(expr) <- "double"
  ids <- intersect(colnames(expr), clin$ID)
  if (length(ids) < ncol(expr) * 0.8)
    log_msg("【核对】", name, "：表达矩阵与临床表只匹配上 ", length(ids), "/", ncol(expr), " 个样本")
  clin <- clin[match(ids, clin$ID), , drop = FALSE]; rownames(clin) <- NULL
  obj <- list(expr = expr[, ids, drop = FALSE], clin = clin)
  saveRDS(obj, file.path(DIR_PRO, paste0(name, ".rds")))
  has_os <- !is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0
  COHORTS[[name]] <<- data.frame(
    cohort = name, role = role, platform = platform, endpoint = endpoint,
    n_samples = ncol(obj$expr), n_genes = nrow(obj$expr),
    n_with_survival = sum(has_os), n_events = sum(clin$OS[has_os] == 1),
    n_responders = if ("response" %in% colnames(clin)) sum(clin$response %in% c("CR/PR", "pCR"), na.rm = TRUE) else NA)
  log_msg(name, "：", nrow(obj$expr), " 个基因 × ", ncol(obj$expr), " 个样本；有生存 ", sum(has_os),
          " 例，事件 ", sum(clin$OS[has_os] == 1), " 例")
  invisible(obj)
}

geo_expr_symbol <- function(es) {
  m <- auto_log2(Biobase::exprs(es))
  log_msg("样本中位数范围：", paste(round(range(apply(m, 2, stats::median, na.rm = TRUE)), 2), collapse = " ~ "),
          "（差异很大时考虑 limma::normalizeBetweenArrays）")
  probe_to_symbol(m, gpl_symbol_map(Biobase::fData(es)))
}

# =============================================================================
# 1. TCGA-BLCA
# =============================================================================
f_tcga <- file.path(DIR_PRO, "TCGA_BLCA.rds")
if (!file.exists(f_tcga)) {
  log_msg("下载 TCGA-BLCA 表达数据（约 430 个样本，需要较长时间）")
  q <- TCGAbiolinks::GDCquery(project = "TCGA-BLCA", data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  TCGAbiolinks::GDCdownload(q, directory = GDC_DIR, files.per.chunk = 50)
  se <- TCGAbiolinks::GDCprepare(q, directory = GDC_DIR)

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  keep <- rd$gene_type == "protein_coding"
  cnt <- SummarizedExperiment::assay(se, "unstranded")[keep, ]
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")[keep, ]
  sym <- rd$gene_name[keep]
  ord <- order(rowMeans(tpm), decreasing = TRUE)          # 同名基因保留平均表达最高的一行
  pick <- ord[!duplicated(sym[ord])]
  cnt <- cnt[pick, ]; tpm <- tpm[pick, ]
  rownames(cnt) <- rownames(tpm) <- sym[pick]

  # 只保留原发肿瘤（01）和正常组织（11）；同一样本多个 aliquot 时保留排序第一个
  bc <- sort(colnames(cnt))
  bc <- bc[substr(bc, 14, 15) %in% c("01", "11")]
  bc <- bc[!duplicated(substr(bc, 1, 15))]
  cnt <- cnt[, bc]; tpm <- tpm[, bc]
  colnames(cnt) <- colnames(tpm) <- substr(bc, 1, 15)

  # 生存与临床：TCGA-CDR（Liu et al. Cell 2018），UCSC Xena pancanAtlasHub
  cdr <- tryCatch({
    xe <- UCSCXenaTools::XenaGenerate(subset = XenaHostNames == "pancanAtlasHub")
    xe <- UCSCXenaTools::XenaFilter(xe, filterDatasets = "Survival_SupplementalTable_S1_20171025_xena_sp")
    xd <- UCSCXenaTools::XenaDownload(UCSCXenaTools::XenaQuery(xe), destdir = file.path(DIR_RAW, "xena"))
    as.data.frame(UCSCXenaTools::XenaPrepare(xd))
  }, error = function(e) {
    log_msg("UCSCXenaTools 失败（", conditionMessage(e), "），改为直接下载")
    f <- file.path(DIR_RAW, "xena", "Survival_SupplementalTable_S1_20171025_xena_sp.gz")
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    utils::download.file("https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz",
                         f, mode = "wb")
    as.data.frame(data.table::fread(f))
  })
  names(cdr) <- make.names(names(cdr))       # readr/fread 保留列名里的空格（"cancer type abbreviation"）
  cdr <- cdr[cdr$cancer.type.abbreviation == "BLCA", ]
  if (!nrow(cdr)) stop("TCGA-CDR 中没有 BLCA 行，检查列名：", paste(head(names(cdr), 8), collapse = ", "))

  stage <- sub("^Stage ([IV]+)[ABC]?$", "\\1", cdr$ajcc_pathologic_tumor_stage)
  stage[!stage %in% c("I", "II", "III", "IV")] <- NA
  grade <- ifelse(cdr$histological_grade == "High Grade", "High",
                  ifelse(cdr$histological_grade == "Low Grade", "Low", NA))
  cdr_tab <- data.frame(ID = cdr$sample,
                        OS = num(cdr$OS), OS.time = num(cdr$OS.time),
                        DSS = num(cdr$DSS), DSS.time = num(cdr$DSS.time),
                        PFI = num(cdr$PFI), PFI.time = num(cdr$PFI.time),
                        age = num(cdr$age_at_initial_pathologic_diagnosis),
                        gender = tolower(cdr$gender), stage = stage, grade = grade)

  clin <- data.frame(ID = colnames(tpm), patient = substr(colnames(tpm), 1, 12),
                     group = ifelse(substr(colnames(tpm), 14, 15) == "01", "Tumor", "Normal"))
  clin <- merge(clin, cdr_tab, by = "ID", all.x = TRUE, sort = FALSE)
  clin <- clin[match(colnames(tpm), clin$ID), ]
  clin[clin$group == "Normal", c("OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time")] <- NA

  saveRDS(list(counts = cnt, tpm = tpm, expr = log2(tpm + 1), clin = clin), f_tcga)
}
tcga <- readRDS(f_tcga)
print(table(tcga$clin$group))                           # 【核对】预期约 410 例肿瘤 + 19 例正常
tc <- tcga$clin[tcga$clin$group == "Tumor", ]
has_os <- !is.na(tc$OS.time) & tc$OS.time > 0
COHORTS[["TCGA_BLCA"]] <- data.frame(cohort = "TCGA_BLCA", role = "training", platform = "RNA-seq (STAR, log2 TPM+1)",
                                     endpoint = "OS", n_samples = nrow(tc), n_genes = nrow(tcga$expr),
                                     n_with_survival = sum(has_os), n_events = sum(tc$OS[has_os] == 1),
                                     n_responders = NA)

# 体细胞突变（04 步算 TMB、APOBEC，并比较 ARID1A 等突变）
f_maf <- file.path(DIR_PRO, "TCGA_BLCA_maf.rds")
if (!file.exists(f_maf)) {
  qm <- TCGAbiolinks::GDCquery(project = "TCGA-BLCA", data.category = "Simple Nucleotide Variation",
                               access = "open", data.type = "Masked Somatic Mutation",
                               workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
  TCGAbiolinks::GDCdownload(qm, directory = GDC_DIR)
  saveRDS(TCGAbiolinks::GDCprepare(qm, directory = GDC_DIR), f_maf)
}

# =============================================================================
# 2. GEO 预后验证队列
# =============================================================================
# ---- GSE13507（GPL6102）：只用原发膀胱癌 165 例；OS 69 例死亡 ----------------
es <- get_geo_eset("GSE13507"); pd <- geo_pheno(es)
src <- tolower(pcol(pd, "biological source"))
mo  <- num(pcol(pd, "survival month"))
clin <- data.frame(ID = pd$geo_accession,
                   OS.time = mo * MONTH2DAY, OS = as.integer(tolower(pcol(pd, "overall survival")) == "death"),
                   CSS.time = mo * MONTH2DAY, CSS = as.integer(tolower(pcol(pd, "cancer specific survival")) == "death"),
                   age = num(pcol(pd, "AGE")), gender = ifelse(pcol(pd, "SEX") == "M", "male", "female"),
                   T = t_stage(pcol(pd, "stage")), N = n_stage(pcol(pd, "stage")),
                   grade = tolower(pcol(pd, "grade")), invasiveness = pcol(pd, "invasiveness"),
                   systemic_chemo = pcol(pd, "systemic chemo"))
clin <- clin[src %in% "primary bladder cancer", ]
print(table(OS = clin$OS, useNA = "ifany"))           # 【核对】预期 165 例，死亡 69 例
save_cohort("GSE13507", geo_expr_symbol(es), clin, "validation", "GPL6102 Illumina", "OS")

# ---- GSE32894（GPL6947）：DSS，224 例有随访，25 例事件 -------------------------
es <- get_geo_eset("GSE32894"); pd <- geo_pheno(es)
clin <- data.frame(ID = pd$geo_accession,
                   OS.time = num(pcol(pd, "time_to_dod_(months)")) * MONTH2DAY,
                   OS = as.integer(tolower(pcol(pd, "dod_event_(yes/no)")) == "yes"),
                   age = num(pcol(pd, "age")), gender = ifelse(pcol(pd, "gender") == "M", "male", "female"),
                   T = t_stage(pcol(pd, "tumor_stage")), grade = pcol(pd, "tumor_grade"),
                   subtype_Lund = pcol(pd, "molecular_subtype"))
clin$OS[is.na(clin$OS.time)] <- NA
print(table(DSS = clin$OS, useNA = "ifany"))          # 【核对】预期 224 例，事件 25 例
save_cohort("GSE32894", geo_expr_symbol(es), clin, "validation", "GPL6947 Illumina", "DSS")

# ---- GSE31684（GPL570）：根治性膀胱切除 93 例；OS（DOD/DOC 为死亡）-------------
es <- get_geo_eset("GSE31684"); pd <- geo_pheno(es)
st <- toupper(pcol(pd, "last known status"))
clin <- data.frame(ID = pd$geo_accession,
                   OS.time = num(pcol(pd, "survival.months")) * MONTH2DAY,
                   OS = as.integer(st %in% c("DOD", "DOC")),
                   CSS = as.integer(st == "DOD"),
                   RFS.time = num(pcol(pd, "recurrence free survival months (distant and local)")) * MONTH2DAY,
                   RFS = as.integer(tolower(pcol(pd, "recurrence/dod")) == "yes"),
                   age = num(pcol(pd, "age at rc")), gender = tolower(pcol(pd, "gender")),
                   T = t_stage(pcol(pd, "rc_stage")), N = tolower(pcol(pd, "plnd result")),
                   grade = tolower(pcol(pd, "rc grade")))
clin$CSS.time <- clin$OS.time
print(table(status = st, useNA = "ifany"))            # 【核对】预期 DOD 38、DOC 27、NED 28
save_cohort("GSE31684", geo_expr_symbol(es), clin, "validation", "GPL570 Affymetrix", "OS")

# ---- GSE48075（GPL6947）：73 例有随访，OS 45 例事件 ---------------------------
es <- get_geo_eset("GSE48075"); pd <- geo_pheno(es)
cs <- tolower(pcol(pd, "os censor"))
clin <- data.frame(ID = pd$geo_accession,
                   OS.time = num(pcol(pd, "survival (mo)")) * MONTH2DAY,
                   OS = ifelse(cs == "uncensored", 1L, ifelse(cs == "censored", 0L, NA)),
                   age = num(pcol(pd, "age (at specimen collection)")),
                   gender = tolower(pcol(pd, "gender")),
                   T = t_stage(pcol(pd, "cstage")), N = n_stage(pcol(pd, "cstage")))
print(table(OS = clin$OS, useNA = "ifany"))           # 【核对】预期 73 例，事件 45 例
save_cohort("GSE48075", geo_expr_symbol(es), clin, "validation", "GPL6947 Illumina", "OS")

# =============================================================================
# 3. 治疗队列
# =============================================================================
# ---- IMvigor210（阿替利珠单抗）：不安装官方包，直接读压缩包里的 cds.RData ------
imv_dir <- file.path(DIR_RAW, "IMvigor210"); dir.create(imv_dir, showWarnings = FALSE, recursive = TRUE)
imv_tar <- file.path(imv_dir, "IMvigor210CoreBiologies_1.0.1.tar.gz")
if (!file.exists(imv_tar))
  utils::download.file("http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.1.tar.gz",
                       imv_tar, mode = "wb")
imv_rdata <- file.path(imv_dir, "IMvigor210CoreBiologies", "data", "cds.RData")
if (!file.exists(imv_rdata))
  utils::untar(imv_tar, files = "IMvigor210CoreBiologies/data/cds.RData", exdir = imv_dir)
e <- new.env(); load(imv_rdata, envir = e); cds <- e$cds
cnt <- get("counts", envir = attr(cds, "assayData"))                 # 31,286 基因 × 348 例
pd  <- attr(attr(cds, "phenoData"), "data")
fd  <- attr(attr(cds, "featureData"), "data")                        # entrez_id、symbol、length
ok  <- !is.na(fd$symbol) & nzchar(fd$symbol) & !is.na(fd$length) & fd$length > 0
rpk <- cnt[ok, ] / (fd$length[ok] / 1000)
tpm <- t(t(rpk) / colSums(rpk)) * 1e6
imv_expr <- collapse_to_symbol(log2(tpm + 1), fd$symbol[ok])
bor <- pd[["Best Confirmed Overall Response"]]
clin <- data.frame(ID = rownames(pd),
                   OS.time = num(pd$os) * MONTH2DAY, OS = as.integer(pd$censOS),
                   response = as.character(pd$binaryResponse), BOR = as.character(bor),
                   TMB = num(pd[["FMOne mutation burden per MB"]]),
                   neoantigen = num(pd[["Neoantigen burden per MB"]]),
                   immune_phenotype = as.character(pd[["Immune phenotype"]]),
                   IC = as.character(pd[["IC Level"]]), TC = as.character(pd[["TC Level"]]),
                   platinum = as.character(pd[["Received platinum"]]),
                   tissue = as.character(pd$Tissue), ECOG = as.character(pd[["Baseline ECOG Score"]]),
                   Lund2 = as.character(pd$Lund2), TCGA_subtype = as.character(pd[["TCGA Subtype"]]),
                   gender = tolower(as.character(pd$Sex)))
print(table(clin$response, useNA = "ifany"))          # 预期 CR/PR 68、SD/PD 230、NA 50
save_cohort("IMvigor210", imv_expr, clin, "immunotherapy (atezolizumab)", "RNA-seq (counts → log2 TPM+1)", "ORR / OS")

# ---- GSE176307（多种 PD-1/PD-L1 抑制剂）：表达在补充文件里，样本为行 ------------
f <- geo_suppl("GSE176307", "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz")
d <- data.table::fread(f, data.table = FALSE)
m <- t(as.matrix(d[, -1])); colnames(m) <- as.character(d[[1]])
m <- auto_log2(ensure_symbols(m))
es <- get_geo_eset("GSE176307"); pd <- geo_pheno(es)
resp <- toupper(pcol(pd, "io.response"))
clin <- data.frame(ID = sub("^.*(BACI[0-9]+).*$", "\\1", pd$title),
                   OS.time = num(pcol(pd, "overall survival")),     # 单位为天（数值如 727、17、252）
                   OS = as.integer(tolower(pcol(pd, "alive")) == "no"),
                   PFS.time = num(pcol(pd, "pfs")),
                   PFS = as.integer(tolower(pcol(pd, "progressed")) == "yes"),
                   BOR = resp,
                   response = ifelse(resp %in% c("CR", "PR"), "CR/PR", ifelse(resp %in% c("SD", "PD"), "SD/PD", NA)),
                   TMB = num(pcol(pd, "tmb")), therapy = pcol(pd, "io.therapy"),
                   age = num(pcol(pd, "age")), gender = ifelse(pcol(pd, "gender") == "M", "male", "female"))
log_msg("GSE176307 OS.time 范围：", paste(range(clin$OS.time, na.rm = TRUE), collapse = " ~ "),
        "【核对】若最大值 < 150，说明单位是月，需要乘以 30.4375")
print(table(clin$BOR, useNA = "ifany"))               # 预期 CR 7、PR 9、SD 4、PD 69
log_msg("表达矩阵样本名示例：", paste(head(colnames(m), 3), collapse = ", "))
save_cohort("GSE176307", m, clin, "immunotherapy (anti-PD-1/PD-L1)", "RNA-seq (normalized log)", "ORR / OS")

# ---- GSE169455（铂类新辅助/诱导化疗）：按基因汇总的表达在补充文件里 ------------
f <- geo_suppl("GSE169455", "GSE169455_normalized_by_gene.txt.gz")
# 文件开头有 5 行 # 注释和 1 个空行，且注释里也出现 "ID_REF"，所以按"以 ID_REF<Tab> 开头"定位表头
con <- gzfile(f); hdr <- grep("^ID_REF\t", readLines(con, n = 50))[1]; close(con)
d <- data.table::fread(f, skip = hdr - 1, header = TRUE, data.table = FALSE)
m <- auto_log2(collapse_to_symbol(as.matrix(d[, -1]), as.character(d$ID_REF)))
es <- get_geo_eset("GSE169455"); pd <- geo_pheno(es)
colnames(m) <- pd$geo_accession[match(colnames(m), pd$title)]       # 列名是样本标题，如 BLCA_Cx_73
log_msg("GSE169455 表达列匹配到 GSM：", sum(!is.na(colnames(m))), "/", ncol(m))   # 【核对】应为 149/149
m <- m[, !is.na(colnames(m))]
tn <- pcol(pd, "pathologic tn stage rc specimen")
Tp <- t_stage(tn); Np <- n_stage(tn)
clin <- data.frame(ID = pd$geo_accession, OS.time = NA_real_, OS = NA_integer_,
                   ypTN = tn, ypT = Tp, ypN = Np,
                   pCR = as.integer(Tp == "0" & Np == "0"),
                   downstage = as.integer(Tp %in% c("0", "A", "IS", "1") & Np == "0"),
                   setting = pcol(pd, "neoadjuvant vs induction chemotherapy"),
                   consensus = pcol(pd, "consensus classifier subtype"),
                   lundtax = pcol(pd, "lundtax rna subtype"),
                   cT = t_stage(pcol(pd, "clinical tnm-staging")))
clin$response <- ifelse(clin$pCR == 1, "pCR", "non-pCR")
print(table(ypTN = tn, useNA = "ifany"))
print(table(pCR = clin$pCR, downstage = clin$downstage, setting = clin$setting, useNA = "ifany"))   # 【核对】
save_cohort("GSE169455", m, clin, "neoadjuvant/induction chemotherapy", "GPL6244 Affymetrix (gene level)", "pCR / downstaging")

# =============================================================================
# Table S2：队列信息
# =============================================================================
S2 <- do.call(rbind, COHORTS)
print(S2, row.names = FALSE)
utils::write.csv(S2, file.path(OUT, "Table_S2_cohorts.csv"), row.names = FALSE)
save_session(STEP)
log_msg("02 完成")

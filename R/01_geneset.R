# =============================================================================
# 01_geneset.R —— 由 R-loopBase v1/v2.0 构建 5 个人类 R-loop 调控基因集
# 只用 base R；输入文件已放在 data/geneset/，缺失时自动从 R-loopBase 下载
# 输出：data/geneset/Rloop_genesets.rds、.gmt、汇总表、补充表 S1
# =============================================================================
source("R/00_config.R")
STEP <- "01_geneset"; OUT <- res_dir(STEP)

rlb_file <- function(v) {
  f <- file.path(DIR_GS, sprintf("R-loopBase_%s_gene_list_of_R-loop_regulators.tsv.gz", v))
  if (!file.exists(f)) {
    url <- sprintf("https://rloopbase.nju.edu.cn/download/%s/r-loop-regulator/gene_list_of_R-loop_regulators.gz", v)
    log_msg("下载 ", url)
    utils::download.file(url, f, mode = "wb")
  }
  f
}

read_rlb <- function(f) {
  d <- utils::read.delim(gzfile(f), quote = "", check.names = FALSE, comment.char = "",
                         colClasses = "character", na.strings = "")
  d <- d[, nzchar(colnames(d)), drop = FALSE]          # v1 表头末尾多一个制表符
  d <- d[d$Species == "Homo sapiens", , drop = FALSE]
  d$Regulator <- trimws(d$Regulator)
  d
}

v2 <- read_rlb(rlb_file("v2"))
v1 <- read_rlb(rlb_file("v1"))
v2$Tier <- as.integer(v2$Tier)
v2$Direction[is.na(v2$Direction)] <- "."

# ---- 5 个基因集 --------------------------------------------------------------
GS <- list(
  RL_core     = v2$Regulator[v2$Tier %in% PARAM$geneset_tiers],
  RL_all      = v2$Regulator,
  RL_v1       = v1$Regulator,
  RL_promote  = v2$Regulator[v2$Direction == "+"],
  RL_suppress = v2$Regulator[v2$Direction == "-"]
)
GS <- lapply(GS, function(g) sort(unique(g[!is.na(g) & nzchar(g)])))

# 可选：用 HGNChelper 统一基因名。默认关闭：R-loopBase 用的已是当前官方名，
# 而 TCGA/GEO 的注释往往更旧，强行更新反而可能对不上。04 步会报告每个基因集在各队列中的覆盖率，
# 覆盖率明显偏低时再把下面改成 TRUE 重跑
UPDATE_SYMBOLS <- FALSE
if (UPDATE_SYMBOLS && requireNamespace("HGNChelper", quietly = TRUE)) {
  upd <- function(g) {
    chk <- suppressWarnings(HGNChelper::checkGeneSymbols(g, unmapped.as.na = FALSE))
    new <- ifelse(is.na(chk$Suggested.Symbol), g, sub(" ///.*$", "", chk$Suggested.Symbol))
    changed <- g != new
    if (any(changed)) log_msg("HGNChelper 更新了 ", sum(changed), " 个基因名，例如：",
                              paste(head(paste(g[changed], new[changed], sep = "→"), 5), collapse = ", "))
    sort(unique(new))
  }
  GS <- lapply(GS, upd)
}

summ <- data.frame(geneset = names(GS), n_genes = lengths(GS),
                   definition = c("v2.0 Tier 1-3 (low-throughput experimental evidence)",
                                  "v2.0 all human regulators (Tier 1-5)",
                                  "v1 all human regulators",
                                  "v2.0 Direction '+' (promotes R-loop formation)",
                                  "v2.0 Direction '-' (suppresses/resolves R-loops)"))
print(summ, row.names = FALSE)

# 与预期数量核对（R-loopBase 官网更新后数量可能变化，此时只提示，不报错）
expect <- c(RL_core = 351, RL_all = 1709, RL_v1 = 1185, RL_promote = 40, RL_suppress = 250)
diff <- summ$n_genes != expect[summ$geneset]
if (any(diff)) log_msg("提示：以下基因集数量与 2026-11-07 下载版本不同（可能是数据库更新或基因名更新）：",
                       paste(summ$geneset[diff], collapse = ", "))

# 膀胱癌相关的几个代表基因，确认收录情况（方案前言要用）
key <- c("ARID1A", "ZBTB11", "CDK12", "APOBEC3B", "NSUN2", "RNASEH1", "RNASEH2A", "SETX", "AQR",
         "BRCA1", "BRCA2", "TOP1", "SRSF1", "DHX9", "FANCD2", "ERCC2", "KDM6A", "TP53")
key_tab <- data.frame(gene = key,
                      in_v1 = key %in% v1$Regulator,
                      in_v2 = key %in% v2$Regulator,
                      tier_v2 = v2$Tier[match(key, v2$Regulator)],
                      direction_v2 = v2$Direction[match(key, v2$Regulator)])
print(key_tab, row.names = FALSE)

# ---- 保存 --------------------------------------------------------------------
saveRDS(GS, file.path(DIR_GS, "Rloop_genesets.rds"))
writeLines(vapply(names(GS), function(n) paste(c(n, summ$definition[summ$geneset == n], GS[[n]]),
                                                collapse = "\t"), character(1)),
           file.path(DIR_GS, "Rloop_genesets.gmt"))
utils::write.csv(summ, file.path(OUT, "geneset_summary.csv"), row.names = FALSE)
utils::write.csv(key_tab, file.path(OUT, "key_genes_in_RloopBase.csv"), row.names = FALSE)

# 补充表 S1：v2.0 全部人类调控基因的注释（论文补充材料直接用）
S1 <- v2[, intersect(c("Regulator", "Entrez ID", "Description", "Tier", "Direction",
                       "Evidence (PMID)", "Function (PMID)", "References", "Update (PMID)"), colnames(v2))]
S1$in_RL_core <- S1$Tier %in% PARAM$geneset_tiers
S1$in_v1 <- S1$Regulator %in% v1$Regulator
S1 <- S1[order(S1$Tier, S1$Regulator), ]
utils::write.csv(S1, file.path(OUT, "Table_S1_Rloop_regulators_v2.csv"), row.names = FALSE)

# 已发表模型比较用的模板（06 步用；列名为 Mime1 要求的格式）
tpl <- file.path(DIR_GS, "published_BLCA_signatures.csv")
if (!file.exists(tpl))
  writeLines("model,PMID,Cancer,Author,Coef,symbol", tpl)

save_session(STEP)
log_msg("01 完成：基因集保存在 data/geneset/Rloop_genesets.rds")

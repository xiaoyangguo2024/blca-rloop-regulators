# =============================================================================
# 09c_rlrs_gene_check.R —— RLRS 16 个基因的核对（为正文措辞提供依据，不改变任何模型）
#   1. 每个基因是否属于 Hallmark E2F 靶基因 / G2M 检查点、GO 细胞周期（GO:0007049）、GO DNA 复制（GO:0006260）
#   2. TCGA-BLCA 肿瘤中每个基因与增殖评分（Prolif）的 Spearman 相关
#   3. 每个队列中 16 个基因的覆盖情况（缺失基因的 z 值按 0 处理，见 07 步）
# 输入：data/processed/TCGA_BLCA.rds、TCGA_tumor_annot.rds、各队列 rds、results/09_hub_gene_spatial/RLRS_coefficients.csv
# 输出：results/09_hub_gene_spatial/RLRS_genes_proliferation.csv、RLRS_gene_coverage.csv
# 内存：一次只读入一个队列。
# =============================================================================
source("R/00_config.R")
need(c("msigdbr"))
OUT <- res_dir("09_hub_gene_spatial")
genes <- utils::read.csv(file.path(OUT, "RLRS_coefficients.csv"))$term
genes <- setdiff(genes, "(Intercept)")

H <- get_msig("H")
bp <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
go_set <- function(id) unique(bp$gene_symbol[bp$gs_exact_source == id])
cc <- go_set("GO:0007049"); rep_dna <- go_set("GO:0006260")

tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
t1 <- tum[tum$group == "Tumor", ]
ex <- tcga$expr[, t1$ID]
rho <- vapply(genes, function(g) if (g %in% rownames(ex))
  suppressWarnings(stats::cor(as.numeric(ex[g, ]), t1$Prolif, method = "spearman")) else NA_real_, numeric(1))
rm(tcga, ex); invisible(gc())

res <- data.frame(gene = genes,
                  in_HALLMARK_E2F = genes %in% H[["HALLMARK_E2F_TARGETS"]],
                  in_HALLMARK_G2M = genes %in% H[["HALLMARK_G2M_CHECKPOINT"]],
                  in_GO_cell_cycle = genes %in% cc,
                  in_GO_DNA_replication = genes %in% rep_dna,
                  rho_Prolif_TCGA = round(rho, 3))
res <- res[order(-res$rho_Prolif_TCGA), ]
utils::write.csv(res, file.path(OUT, "RLRS_genes_proliferation.csv"), row.names = FALSE)
print(res, row.names = FALSE)
log_msg("Hallmark E2F 或 G2M：", sum(res$in_HALLMARK_E2F | res$in_HALLMARK_G2M), "/", nrow(res),
        "；GO 细胞周期或 DNA 复制：", sum(res$in_GO_cell_cycle | res$in_GO_DNA_replication), "/", nrow(res),
        "；与增殖评分 rho > 0.5：", sum(res$rho_Prolif_TCGA > 0.5, na.rm = TRUE), "/", nrow(res))

# 各队列的基因覆盖
fs <- c(TCGA_BLCA = "TCGA_BLCA", GSE13507 = "GSE13507", GSE32894 = "GSE32894", GSE48075 = "GSE48075", GSE31684 = "GSE31684",
        IMvigor210 = "IMvigor210", GSE176307 = "GSE176307", GSE169455 = "GSE169455")
cov <- do.call(rbind, lapply(names(fs), function(n) {
  x <- readRDS(file.path(DIR_PRO, paste0(fs[[n]], ".rds")))
  rn <- rownames(x$expr); rm(x); invisible(gc())
  data.frame(cohort = n, n_present = sum(genes %in% rn), missing = paste(setdiff(genes, rn), collapse = ", "))
}))
utils::write.csv(cov, file.path(OUT, "RLRS_gene_coverage.csv"), row.names = FALSE)
print(cov, row.names = FALSE)
log_msg("09c 完成")

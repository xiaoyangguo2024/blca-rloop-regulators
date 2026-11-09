# =============================================================================
# 04b_proliferation_overlap.R —— RL_core 评分与增殖高度相关，是不是只因为两组基因有重叠？
#   1. RL_core 与 Hallmark E2F 靶基因 / G2M 检查点基因的重叠
#   2. 去掉重叠基因后重新计算 ssGSEA（RL_core_noProlif），看与增殖评分的相关和预后是否改变
#   3. 每个 RL_core 基因与增殖评分的相关分布
#   4. 8 个队列中 RL_core 与增殖评分的相关（复现性）
# 输入：04 步的 ssgsea_scores_all_cohorts.rds、TCGA_tumor_annot.rds
# 输出：results/04_bulk_score/proliferation_overlap_*.csv
# =============================================================================
source("R/00_config.R")
need(c("GSVA", "msigdbr", "survival"))
OUT <- res_dir("04_bulk_score")
GS <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))
H <- get_msig("H")
prolif_genes <- union(H[["HALLMARK_E2F_TARGETS"]], H[["HALLMARK_G2M_CHECKPOINT"]])
ov <- intersect(GS$RL_core, prolif_genes)
log_msg("RL_core ", length(GS$RL_core), " 个基因，其中 ", length(ov), " 个属于 Hallmark E2F 或 G2M 基因集")

tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
t1 <- tum[tum$group == "Tumor", ]
ex <- tcga$expr[, t1$ID]
s <- run_ssgsea(ex, list(RL_core_noProlif = setdiff(GS$RL_core, prolif_genes)))
t1$RL_core_noProlif <- as.numeric(s["RL_core_noProlif", t1$ID])

ct <- function(a, b) stats::cor.test(a, b, method = "spearman", exact = FALSE)
r1 <- ct(t1$RL_core_noProlif, t1$Prolif); r2 <- ct(t1$RL_core_noProlif, t1$RL_core); r0 <- ct(t1$RL_core, t1$Prolif)
d <- t1[!is.na(t1$OS.time) & t1$OS.time > 0 & !is.na(t1$OS), ]
d$z <- as.numeric(scale(d$RL_core_noProlif))
u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z, data = d))
a <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z + age + stage_bin + Prolif, data = d))
summ <- data.frame(
  item = c("RL_core genes", "RL_core genes in Hallmark E2F/G2M", "RL_core genes left (noProlif set)",
           "rho RL_core vs Prolif", "rho RL_core_noProlif vs Prolif", "rho RL_core_noProlif vs RL_core",
           "OS HR per SD RL_core_noProlif (unadjusted)", "OS HR per SD RL_core_noProlif (adjusted age, stage, Prolif)"),
  value = c(length(GS$RL_core), length(ov), length(setdiff(GS$RL_core, prolif_genes)),
            round(r0$estimate, 2), round(r1$estimate, 2), round(r2$estimate, 2),
            sprintf("%.2f (%.2f-%.2f)", u$HR[1], u$lower[1], u$upper[1]),
            sprintf("%.2f (%.2f-%.2f)", a$HR[a$term == "z"], a$lower[a$term == "z"], a$upper[a$term == "z"])),
  p = c(NA, NA, NA, r0$p.value, r1$p.value, r2$p.value, u$p[1], a$p[a$term == "z"]))
print(summ, row.names = FALSE)
utils::write.csv(summ, file.path(OUT, "proliferation_overlap_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene = sort(ov)), file.path(OUT, "proliferation_overlap_genes.csv"), row.names = FALSE)

# 每个 RL_core 基因与增殖评分的相关
g <- intersect(GS$RL_core, rownames(ex))
rho_g <- vapply(g, function(x) suppressWarnings(stats::cor(as.numeric(ex[x, ]), t1$Prolif, method = "spearman")), numeric(1))
pg <- data.frame(gene = g, rho_Prolif = rho_g, in_E2F_G2M = g %in% prolif_genes)
utils::write.csv(pg[order(-pg$rho_Prolif), ], file.path(OUT, "proliferation_overlap_per_gene_rho.csv"), row.names = FALSE)
log_msg("RL_core 基因与增殖评分的 Spearman rho：中位数 ", round(stats::median(rho_g, na.rm = TRUE), 2),
        "；rho > 0.5 的基因 ", sum(rho_g > 0.5, na.rm = TRUE), " 个（其中不属于 E2F/G2M 的 ",
        sum(rho_g > 0.5 & !pg$in_E2F_G2M, na.rm = TRUE), " 个）")

# 8 个队列中 RL_core 与增殖评分的相关
sc_all <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
rc <- do.call(rbind, lapply(names(sc_all), function(n) {
  x <- sc_all[[n]]
  if (n == "TCGA_BLCA") x <- x[x$ID %in% t1$ID, ]
  r <- ct(x$RL_core, x$Prolif); r2 <- ct(x$RHI, x$Prolif)
  data.frame(cohort = n, n = nrow(x), rho_RLcore_Prolif = unname(r$estimate), p_RLcore = r$p.value,
             rho_RHI_Prolif = unname(r2$estimate), p_RHI = r2$p.value)
}))
print(rc, row.names = FALSE, digits = 3)
utils::write.csv(rc, file.path(OUT, "proliferation_correlation_by_cohort.csv"), row.names = FALSE)

# =============================================================================
# 5. 与之对称的检验：RHI 的预后关联能否由 EMT（基质）信号解释？（RHI 与 EMT 相关 0.61）
#    各有生存数据的队列：RHI 单独、RHI + EMT、EMT 单独（均为每 SD 的 HR；GSE32894 为 DSS）
# =============================================================================
rows <- list()
for (n in c("TCGA_BLCA", "GSE13507", "GSE32894", "GSE31684", "GSE48075", "IMvigor210", "GSE176307")) {
  cl <- if (n == "TCGA_BLCA") t1[, c("ID", "OS.time", "OS")] else readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin[, c("ID", "OS.time", "OS")]
  e <- merge(cl, sc_all[[n]][, c("ID", "RHI", "EMT")], by = "ID")
  e <- e[!is.na(e$OS.time) & e$OS.time > 0 & !is.na(e$OS), ]
  e$RHI_z <- as.numeric(scale(e$RHI)); e$EMT_z <- as.numeric(scale(e$EMT))
  S <- survival::Surv(e$OS.time, e$OS)
  m1 <- cox_table(survival::coxph(S ~ RHI_z, data = e))
  m2 <- cox_table(survival::coxph(S ~ RHI_z + EMT_z, data = e))
  m3 <- cox_table(survival::coxph(S ~ EMT_z, data = e))
  rows[[n]] <- data.frame(cohort = n, n = nrow(e), events = sum(e$OS),
                          rho_RHI_EMT = stats::cor(e$RHI, e$EMT, method = "spearman"),
                          HR_RHI = m1$HR, lower_RHI = m1$lower, upper_RHI = m1$upper, p_RHI = m1$p,
                          HR_RHI_adjEMT = m2$HR[1], lower_RHI_adjEMT = m2$lower[1], upper_RHI_adjEMT = m2$upper[1],
                          p_RHI_adjEMT = m2$p[1], HR_EMT = m3$HR, p_EMT = m3$p)
}
rhi <- do.call(rbind, rows); print(rhi[, c("cohort", "n", "events", "rho_RHI_EMT", "HR_RHI", "p_RHI", "HR_RHI_adjEMT",
                                          "p_RHI_adjEMT", "HR_EMT", "p_EMT")], row.names = FALSE, digits = 3)
utils::write.csv(rhi, file.path(OUT, "RHI_EMT_adjusted_Cox_by_cohort.csv"), row.names = FALSE)
log_msg("04b 完成")

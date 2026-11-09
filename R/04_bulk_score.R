# =============================================================================
# 04_bulk_score.R —— bulk 队列 ssGSEA 评分；TCGA 中与临床、分子亚型、突变、IFN/STING 的关联；生存分析
# 输入：02 步生成的各队列 .rds、01 步的基因集
# 输出：data/processed/ssgsea_scores_all_cohorts.rds、TCGA_tumor_annot.rds；results/04_bulk_score/
# =============================================================================
source("R/00_config.R")
need(c("GSVA", "msigdbr", "survival", "survminer", "ggplot2", "ggpubr"))
STEP <- "04_bulk_score"; OUT <- res_dir(STEP)
GS <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))

COHORTS <- c("TCGA_BLCA", "GSE13507", "GSE32894", "GSE31684", "GSE48075", "IMvigor210", "GSE176307", "GSE169455")
coh <- lapply(stats::setNames(COHORTS, COHORTS), function(n) readRDS(file.path(DIR_PRO, paste0(n, ".rds"))))

# =============================================================================
# 1. 基因集在各队列中的覆盖率（芯片平台基因少，覆盖率低于 80% 要在方法里说明）
# =============================================================================
cov <- do.call(rbind, lapply(COHORTS, function(n) data.frame(
  cohort = n, geneset = names(GS), n_genes = lengths(GS),
  n_found = vapply(GS, function(g) sum(g %in% rownames(coh[[n]]$expr)), integer(1)))))
cov$pct <- round(100 * cov$n_found / cov$n_genes, 1)
utils::write.csv(cov, file.path(OUT, "geneset_coverage_by_cohort.csv"), row.names = FALSE)
print(stats::xtabs(pct ~ cohort + geneset, cov))

# =============================================================================
# 2. ssGSEA：R-loop 基因集 + IFN / cGAS-STING / T 细胞炎症 / DNA 修复 / 增殖
# =============================================================================
H  <- get_msig("H")
C2 <- get_msig("C2")
pick <- function(lst, nm) { if (!nm %in% names(lst)) log_msg("【核对】MSigDB 中找不到 ", nm); lst[[nm]] }
SIG <- c(GS, list(
  IFN_alpha       = pick(H, "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  IFN_gamma       = pick(H, "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  DNA_repair      = pick(H, "HALLMARK_DNA_REPAIR"),
  E2F             = pick(H, "HALLMARK_E2F_TARGETS"),
  G2M             = pick(H, "HALLMARK_G2M_CHECKPOINT"),
  EMT             = pick(H, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  cGAS_STING_KEGG = pick(C2, "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"),
  STING_Reactome  = pick(C2, "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES"),
  GEP18           = GEP18))
SIG <- SIG[lengths(SIG) > 0]

f_sc <- file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds")
if (!file.exists(f_sc)) {
  scores <- lapply(coh, function(o) {
    s <- as.data.frame(t(run_ssgsea(o$expr, SIG)))
    s$RHI    <- as.numeric(scale(s$RL_promote) - scale(s$RL_suppress))   # 探索性指标
    s$Prolif <- as.numeric(scale(s$E2F) + scale(s$G2M)) / 2              # 增殖评分（用于校正）
    s$ID <- rownames(s)
    s
  })
  saveRDS(scores, f_sc)
}
scores <- readRDS(f_sc)

# =============================================================================
# 3. TCGA：肿瘤 vs 正常、分期、分级
# =============================================================================
tc <- merge(coh$TCGA_BLCA$clin, scores$TCGA_BLCA, by = "ID")
tc$group <- factor(tc$group, levels = c("Normal", "Tumor"))
tum <- tc[tc$group == "Tumor", ]
tum$stage_bin <- factor(ifelse(tum$stage %in% c("III", "IV"), "III-IV", ifelse(tum$stage %in% c("I", "II"), "I-II", NA)))

plots <- list()
for (s in c("RL_core", "RHI")) {
  plots[[paste0("TvN_", s)]] <- ggpubr::ggboxplot(tc, "group", s, color = "group", add = "jitter",
                                                  palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(method = "wilcox.test")
  pat <- intersect(tc$patient[tc$group == "Tumor"], tc$patient[tc$group == "Normal"])
  plots[[paste0("paired_", s)]] <- ggpubr::ggpaired(tc[tc$patient %in% pat, ], x = "group", y = s, id = "patient",
                                                   color = "group", line.color = "grey70",
                                                   palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(paired = TRUE, method = "wilcox.test") +
    ggplot2::labs(title = paste0("Paired (n = ", length(pat), ")"))
  plots[[paste0("stage_", s)]] <- ggpubr::ggboxplot(tum[!is.na(tum$stage), ], "stage", s, color = "stage", add = "jitter") +
    ggpubr::stat_compare_means(method = "kruskal.test")
  plots[[paste0("grade_", s)]] <- ggpubr::ggboxplot(tum[!is.na(tum$grade), ], "grade", s, color = "grade", add = "jitter") +
    ggpubr::stat_compare_means(method = "wilcox.test")
}

# =============================================================================
# 4. consensusMIBC 分子亚型
# =============================================================================
if (requireNamespace("consensusMIBC", quietly = TRUE)) {
  cm <- consensusMIBC::getConsensusClass(coh$TCGA_BLCA$expr[, tum$ID], gene_id = "hgnc_symbol")
  tum$consensus <- cm$consensusClass[match(tum$ID, rownames(cm))]
  print(table(tum$consensus, useNA = "ifany"))
  for (s in c("RL_core", "RHI"))
    plots[[paste0("consensus_", s)]] <- ggpubr::ggboxplot(tum[!is.na(tum$consensus), ], "consensus", s,
                                                         color = "consensus", add = "jitter") +
      ggpubr::stat_compare_means(method = "kruskal.test") + ggpubr::rotate_x_text(30)
} else log_msg("未安装 consensusMIBC，跳过分子亚型")

# =============================================================================
# 5. 突变：ARID1A 等膀胱癌常见突变、TMB、APOBEC
# =============================================================================
maf <- readRDS(file.path(DIR_PRO, "TCGA_BLCA_maf.rds"))
maf$sample <- substr(maf$Tumor_Sample_Barcode, 1, 15)
NONSYN <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation",
            "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
ns <- maf[maf$Variant_Classification %in% NONSYN, ]
mut_samples <- unique(maf$sample)
tum$TMB <- ifelse(tum$ID %in% mut_samples,
                  as.numeric(table(factor(ns$sample, levels = mut_samples))[tum$ID]) / 38, NA)   # 外显子按 38 Mb 计
MUT_GENES <- c("ARID1A", "TP53", "FGFR3", "KDM6A", "RB1", "ERCC2", "KMT2D", "PIK3CA", "STAG2", "ELF3", "ATM", "BRCA2")
for (g in MUT_GENES)
  tum[[paste0("mut_", g)]] <- ifelse(tum$ID %in% mut_samples,
                                     ifelse(tum$ID %in% ns$sample[ns$Hugo_Symbol == g], "Mut", "WT"), NA)
mut_freq <- data.frame(gene = MUT_GENES,
                       n_mut = vapply(MUT_GENES, function(g) sum(tum[[paste0("mut_", g)]] == "Mut", na.rm = TRUE), integer(1)),
                       n_tested = sum(tum$ID %in% mut_samples))
mut_freq$freq <- round(mut_freq$n_mut / mut_freq$n_tested, 3)
print(mut_freq)                                       # 前言里 ARID1A 的突变比例以此为准
utils::write.csv(mut_freq, file.path(OUT, "TCGA_mutation_frequency.csv"), row.names = FALSE)

mut_long <- do.call(rbind, lapply(MUT_GENES, function(g)
  data.frame(gene = g, status = tum[[paste0("mut_", g)]], RL_core = tum$RL_core, RHI = tum$RHI)))
mut_long <- mut_long[!is.na(mut_long$status), ]
mut_test <- do.call(rbind, lapply(MUT_GENES, function(g) {
  d <- mut_long[mut_long$gene == g, ]
  wp <- function(f) tryCatch(stats::wilcox.test(f, d)$p.value, error = function(e) NA_real_)
  data.frame(gene = g, n_mut = sum(d$status == "Mut"),
             median_diff_RL_core = stats::median(d$RL_core[d$status == "Mut"]) - stats::median(d$RL_core[d$status == "WT"]),
             p_RL_core = wp(RL_core ~ status), p_RHI = wp(RHI ~ status))
}))
mut_test$FDR_RL_core <- stats::p.adjust(mut_test$p_RL_core, "BH")
mut_test$FDR_RHI <- stats::p.adjust(mut_test$p_RHI, "BH")
utils::write.csv(mut_test, file.path(OUT, "Rloop_score_by_mutation.csv"), row.names = FALSE)
plots$mutation <- ggpubr::ggboxplot(mut_long, "status", "RL_core", color = "status", facet.by = "gene", nrow = 2,
                                    palette = c(Mut = "#D7301F", WT = "#2C7FB8")) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format")
plots$TMB <- ggpubr::ggscatter(tum[!is.na(tum$TMB), ], "RL_core", "TMB", add = "reg.line", cor.coef = TRUE,
                               cor.method = "spearman", size = 1) + ggplot2::scale_y_log10()

if (requireNamespace("maftools", quietly = TRUE) && requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  m2 <- maf; m2$Tumor_Sample_Barcode <- m2$sample
  mo <- maftools::read.maf(m2, verbose = FALSE)
  tnm <- maftools::trinucleotideMatrix(mo, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                                       prefix = if (grepl("^chr", maf$Chromosome[1])) NULL else "chr",
                                       add = TRUE, useSyn = TRUE)
  ap <- as.data.frame(tnm$APOBEC_scores)
  tum$APOBEC_enrich <- ap$APOBEC_Enrichment[match(tum$ID, as.character(ap$Tumor_Sample_Barcode))]
  plots$APOBEC <- ggpubr::ggscatter(tum[!is.na(tum$APOBEC_enrich), ], "RL_core", "APOBEC_enrich", add = "reg.line",
                                    cor.coef = TRUE, cor.method = "spearman", size = 1)
} else log_msg("未安装 maftools 或 BSgenome.Hsapiens.UCSC.hg38，跳过 APOBEC 富集分数")

# =============================================================================
# 6. 与 IFN / cGAS-STING / T 细胞炎症 / DNA 修复 / 增殖的相关性
# =============================================================================
cor_vars <- intersect(c("RL_core", "RL_all", "RL_v1", "RHI", "RL_promote", "RL_suppress", "IFN_alpha", "IFN_gamma",
                        "cGAS_STING_KEGG", "STING_Reactome", "GEP18", "DNA_repair", "E2F", "G2M", "Prolif", "EMT"),
                      colnames(tum))
cm_r <- stats::cor(tum[, cor_vars], method = "spearman", use = "pairwise.complete.obs")
utils::write.csv(cm_r, file.path(OUT, "signature_spearman_TCGA.csv"))
cl <- data.frame(x = rep(colnames(cm_r), each = nrow(cm_r)), y = rep(rownames(cm_r), ncol(cm_r)), r = as.vector(cm_r))
cl$x <- factor(cl$x, levels = cor_vars); cl$y <- factor(cl$y, levels = rev(cor_vars))
plots$cor_heatmap <- ggplot2::ggplot(cl, ggplot2::aes(x, y, fill = r)) + ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 2.5) +
  ggplot2::scale_fill_gradient2(low = "#2C7FB8", high = "#D7301F", limits = c(-1, 1)) +
  ggplot2::theme_minimal() + ggpubr::rotate_x_text(45) + ggplot2::labs(x = NULL, y = NULL)
for (v in intersect(c("IFN_gamma", "cGAS_STING_KEGG", "GEP18", "DNA_repair"), cor_vars))
  plots[[paste0("scatter_", v)]] <- ggpubr::ggscatter(tum, "RL_core", v, add = "reg.line", cor.coef = TRUE,
                                                     cor.method = "spearman", size = 1)

save_print(plots, file.path(OUT, "Fig4_TCGA_associations.pdf"), 7, 6)

# =============================================================================
# 7. 生存：TCGA（OS/DSS/PFI）KM + Cox（每 1 个标准差）；校正年龄、分期和增殖评分
# =============================================================================
EPS <- list(OS = c("OS.time", "OS"), DSS = c("DSS.time", "DSS"), PFI = c("PFI.time", "PFI"))
km_list <- list(); cox_rows <- list()
for (s in c("RL_core", "RHI", "RL_all", "RL_v1")) for (ep in names(EPS)) {
  tt <- EPS[[ep]][1]; ee <- EPS[[ep]][2]
  d <- tum[!is.na(tum[[tt]]) & tum[[tt]] > 0 & !is.na(tum[[ee]]), ]
  d$group <- median_group(d[[s]]); d$z <- as.numeric(scale(d[[s]]))
  km_list[[paste(s, ep)]] <- km_plot(d, tt, ee, "group", title = paste("TCGA-BLCA", s, ep))
  f_uni <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ z", tt, ee))
  f_adj <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ z + age + stage_bin + Prolif", tt, ee))
  u <- cox_table(survival::coxph(f_uni, data = d))
  a <- cox_table(survival::coxph(f_adj, data = d))
  cox_rows[[paste(s, ep)]] <- data.frame(score = s, endpoint = ep, n = nrow(d), events = sum(d[[ee]] == 1),
                                         HR_uni = u$HR[1], p_uni = u$p[1],
                                         HR_adj = a$HR[a$term == "z"], lower_adj = a$lower[a$term == "z"],
                                         upper_adj = a$upper[a$term == "z"], p_adj = a$p[a$term == "z"])
}
cox_tcga <- do.call(rbind, cox_rows)
print(cox_tcga, row.names = FALSE)
utils::write.csv(cox_tcga, file.path(OUT, "TCGA_Cox_Rloop_scores.csv"), row.names = FALSE)
save_print(km_list, file.path(OUT, "Fig4_TCGA_KM.pdf"), 6, 6)

# ---- 外部队列：RL_core 评分本身的预后价值（补充材料）----------------------------
ext <- list(); km_ext <- list()
for (n in setdiff(COHORTS, c("TCGA_BLCA", "GSE169455"))) {
  d <- merge(coh[[n]]$clin, scores[[n]], by = "ID")
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
  if (nrow(d) < 30) next
  d$group <- median_group(d$RL_core); d$z <- as.numeric(scale(d$RL_core))
  km_ext[[n]] <- km_plot(d, title = paste(n, "RL_core"))
  u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ z, data = d))
  ext[[n]] <- data.frame(cohort = n, n = nrow(d), events = sum(d$OS == 1), HR_perSD = u$HR, lower = u$lower,
                         upper = u$upper, p = u$p)
}
utils::write.csv(do.call(rbind, ext), file.path(OUT, "External_Cox_RL_core.csv"), row.names = FALSE)
save_print(km_ext, file.path(OUT, "FigS_external_KM_RL_core.pdf"), 6, 6)

saveRDS(tum, file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
save_session(STEP)
log_msg("04 完成")

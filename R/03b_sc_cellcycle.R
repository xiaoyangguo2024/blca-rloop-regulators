# =============================================================================
# 03b_sc_cellcycle.R —— 单细胞层面检验：恶性细胞的 R-loop 评分是否主要反映细胞周期状态
#   1. 每个样本内：RL_core（UCell）与 S 期、G2/M 期评分的 Spearman 相关
#   2. R-loop 高/低组中处于 S 或 G2/M 期的细胞比例（以样本为单位，配对 Wilcoxon）
#   3. 只用 G1 期恶性细胞：样本内按 RL_core 中位数重新分组，Hallmark 差异是否还在（去掉细胞周期之后还剩什么）
# 输入：03 步的 scRNA_malignant.rds（含 RL_group 与 Hallmark UCell 评分）
# 输出：results/03_scRNA/ 下的 cellcycle_*.csv、Fig3_cellcycle.pdf
# =============================================================================
source("R/00_config.R")
need(c("Seurat", "UCell", "ggplot2", "ggpubr"))
OUT <- res_dir("03_scRNA")
mal <- readRDS(file.path(DIR_PRO, "scRNA_malignant.rds"))
cc <- Seurat::cc.genes.updated.2019
mal <- CellCycleScoring(mal, s.features = intersect(cc$s.genes, rownames(mal)),
                        g2m.features = intersect(cc$g2m.genes, rownames(mal)), set.ident = FALSE)
md <- mal@meta.data
md$cycling <- md$Phase %in% c("S", "G2M")
log_msg("恶性细胞 ", nrow(md), " 个；S/G2M 期比例 ", round(100 * mean(md$cycling), 1), "%")

# ---- 1–2. 每个样本内的相关与周期比例 -----------------------------------------
per_s <- do.call(rbind, lapply(split(md, md$sample), function(x) data.frame(
  sample = x$sample[1], n_cells = nrow(x),
  rho_RLcore_S = stats::cor(x$RL_core_UCell, x$S.Score, method = "spearman"),
  rho_RLcore_G2M = stats::cor(x$RL_core_UCell, x$G2M.Score, method = "spearman"),
  pct_cycling_High = 100 * mean(x$cycling[x$RL_group == "High"]),
  pct_cycling_Low = 100 * mean(x$cycling[x$RL_group == "Low"]))))
utils::write.csv(per_s, file.path(OUT, "cellcycle_by_sample.csv"), row.names = FALSE)
p_cyc <- stats::wilcox.test(per_s$pct_cycling_High, per_s$pct_cycling_Low, paired = TRUE)$p.value
log_msg("样本内 RL_core 与 G2M 评分相关的中位数 rho = ", round(stats::median(per_s$rho_RLcore_G2M), 2),
        "；S/G2M 期比例 高组 vs 低组（中位数 ", round(stats::median(per_s$pct_cycling_High), 1), "% vs ",
        round(stats::median(per_s$pct_cycling_Low), 1), "%），配对 Wilcoxon P = ", signif(p_cyc, 3))

# ---- 3. 只用 G1 期细胞 ---------------------------------------------------------
g1 <- md[md$Phase == "G1", ]
keep <- names(which(table(g1$sample) >= 2 * PARAM$sc_min_cells_pb))
g1 <- g1[g1$sample %in% keep, ]
g1$G1_group <- stats::ave(g1$RL_core_UCell, g1$sample, FUN = function(v) ifelse(v >= stats::median(v), "High", "Low"))
hcols <- grep("^HALLMARK_.*_UCell$", colnames(g1), value = TRUE)
if (!length(hcols)) {                             # 03 步没保存 Hallmark 评分时重新计算
  H <- get_msig("H")
  mal <- UCell::AddModuleScore_UCell(mal, features = H, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
  hcols <- paste0(names(H), "_UCell")
  g1 <- cbind(g1, mal@meta.data[rownames(g1), hcols])
}
hm <- stats::aggregate(g1[, hcols], by = list(sample = g1$sample, group = g1$G1_group), FUN = mean)
hall_g1 <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA))
}))
hall_g1$FDR <- stats::p.adjust(hall_g1$p, "BH")
hall_g1 <- hall_g1[order(hall_g1$mean_diff), ]
utils::write.csv(hall_g1, file.path(OUT, "Hallmark_High_vs_Low_G1only.csv"), row.names = FALSE)
f_all <- file.path(OUT, "Hallmark_High_vs_Low_paired.csv")
if (file.exists(f_all)) {                          # 与全部恶性细胞的结果并排比较
  ha <- utils::read.csv(f_all)
  both <- merge(ha[, c("pathway", "mean_diff", "FDR")], hall_g1[, c("pathway", "mean_diff", "FDR")],
                by = "pathway", suffixes = c("_all", "_G1"))
  utils::write.csv(both, file.path(OUT, "Hallmark_all_vs_G1only.csv"), row.names = FALSE)
  log_msg("G1 期细胞中仍 FDR < 0.05 的 Hallmark：", sum(both$FDR_G1 < 0.05), "（全部恶性细胞中 ", sum(both$FDR_all < 0.05), "）")
}

# ---- 图 --------------------------------------------------------------------------
long <- rbind(data.frame(sample = per_s$sample, group = "High", pct = per_s$pct_cycling_High),
              data.frame(sample = per_s$sample, group = "Low", pct = per_s$pct_cycling_Low))
long$group <- factor(long$group, levels = c("Low", "High"))
p1 <- ggpubr::ggpaired(long, x = "group", y = "pct", id = "sample", color = "group", line.color = "grey70",
                       palette = c(Low = "#2C7FB8", High = "#D7301F")) +
  ggpubr::stat_compare_means(paired = TRUE) + ggplot2::labs(x = "R-loop score group", y = "Cells in S or G2/M phase (%)")
set.seed(PARAM$seed)
sm <- md[sample(nrow(md), min(8000, nrow(md))), ]
p2 <- ggplot2::ggplot(sm, ggplot2::aes(G2M.Score, RL_core_UCell, colour = Phase)) +
  ggplot2::geom_point(size = 0.3, alpha = 0.5) + ggplot2::theme_classic() +
  ggplot2::labs(x = "G2/M score", y = "RL_core score (UCell)") +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1)))
save_print(list(p1, p2), file.path(OUT, "Fig3_cellcycle.pdf"), 5, 4.5)
save_session("03b_sc_cellcycle")
log_msg("03b 完成")

# =============================================================================
# 05_candidates.R —— 候选基因：DESeq2（肿瘤 vs 正常）+ WGCNA（与 RL_core 最相关的模块）+ 单细胞伪 bulk 差异基因
#   候选基因 =（RL_all ∪ WGCNA 关键模块基因）∩ 单细胞高/低差异基因 ∩ 肿瘤 vs 正常差异基因
#   少于 PARAM$min_candidates 个时自动去掉单细胞这一层（论文方法部分需写明）
#   再在 TCGA 中做单因素 Cox（OS）
# 输出：data/processed/candidates.rds；results/05_candidates/（含 Table S3）
# =============================================================================
source("R/00_config.R")
need(c("DESeq2", "WGCNA", "survival", "ggplot2", "ggrepel", "VennDiagram"))
STEP <- "05_candidates"; OUT <- res_dir(STEP)
GS   <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))
tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))

# =============================================================================
# 1. DESeq2：肿瘤 vs 正常
# =============================================================================
f_deg <- file.path(OUT, "TCGA_DESeq2_Tumor_vs_Normal.csv")
if (!file.exists(f_deg)) {
  cl <- tcga$clin
  dds <- DESeq2::DESeqDataSetFromMatrix(round(tcga$counts[, cl$ID]),
                                        colData = data.frame(row.names = cl$ID,
                                                             group = factor(cl$group, levels = c("Normal", "Tumor"))),
                                        design = ~ group)
  dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 10, ]
  dds <- DESeq2::DESeq(dds)
  deg <- as.data.frame(DESeq2::results(dds, contrast = c("group", "Tumor", "Normal")))
  deg$gene <- rownames(deg)
  utils::write.csv(deg, f_deg, row.names = FALSE)
}
deg <- utils::read.csv(f_deg)
deg_genes <- deg$gene[!is.na(deg$padj) & deg$padj < PARAM$deg_padj & abs(deg$log2FoldChange) > PARAM$deg_lfc]
log_msg("TCGA 肿瘤 vs 正常差异基因：", length(deg_genes))

deg$class <- ifelse(deg$gene %in% deg_genes, ifelse(deg$log2FoldChange > 0, "Up", "Down"), "NS")
deg$is_RL <- deg$gene %in% GS$RL_all
lab <- utils::head(deg[deg$is_RL & deg$class != "NS", ][order(deg$padj[deg$is_RL & deg$class != "NS"]), ], 15)
p <- ggplot2::ggplot(deg[!is.na(deg$padj), ], ggplot2::aes(log2FoldChange, -log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = class), size = 0.5, alpha = 0.6) +
  ggplot2::geom_point(data = deg[deg$is_RL & deg$class != "NS", ], shape = 21, color = "black", size = 1.2) +
  ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), size = 3) +
  ggplot2::scale_color_manual(values = c(Up = "#D7301F", Down = "#2C7FB8", NS = "grey80")) +
  ggplot2::theme_classic() + ggplot2::labs(title = "TCGA-BLCA tumor vs normal (circled: R-loop regulators)")
save_gg(p, file.path(OUT, "FigS_volcano_TCGA_DEG.pdf"), 6, 5)

# =============================================================================
# 2. WGCNA：肿瘤样本，MAD 最高的 5000 个基因
# =============================================================================
f_wgcna <- file.path(DIR_PRO, "WGCNA_net.rds")
ex <- tcga$expr[, tum$ID]
datExpr <- t(ex[order(apply(ex, 1, stats::mad), decreasing = TRUE)[seq_len(PARAM$wgcna_top_mad)], ])
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
if (!file.exists(f_wgcna)) {
  WGCNA::allowWGCNAThreads(4)
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)), networkType = "unsigned", verbose = 0)
  power <- sft$powerEstimate
  if (is.na(power)) { power <- 6; log_msg("【核对】无幂次达到无尺度拓扑 R² 0.85，按 WGCNA FAQ 取 6") }
  grDevices::pdf(file.path(OUT, "WGCNA_soft_threshold.pdf"), 9, 4.5); graphics::par(mfrow = c(1, 2))
  fi <- sft$fitIndices
  plot(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, type = "n", xlab = "Soft threshold", ylab = "Scale-free fit (R²)")
  graphics::text(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, labels = fi$Power, col = "red"); graphics::abline(h = 0.85, col = "red")
  plot(fi$Power, fi$mean.k., type = "n", xlab = "Soft threshold", ylab = "Mean connectivity")
  graphics::text(fi$Power, fi$mean.k., labels = fi$Power, col = "red")
  grDevices::dev.off()
  # 其他包（如 Hmisc）可能遮盖 cor()，导致 blockwiseModules 报错，这里临时指定 WGCNA::cor
  cor <- WGCNA::cor
  net <- WGCNA::blockwiseModules(datExpr, power = power, TOMType = "unsigned", minModuleSize = 30,
                                 mergeCutHeight = 0.25, numericLabels = FALSE, maxBlockSize = 6000,
                                 randomSeed = PARAM$seed, verbose = 0)
  cor <- stats::cor
  net$power <- power
  saveRDS(net, f_wgcna)
}
net <- readRDS(f_wgcna)
log_msg("WGCNA 软阈值 = ", net$power, "；模块数 = ", length(unique(net$colors)) - 1)

grDevices::pdf(file.path(OUT, "WGCNA_dendrogram.pdf"), 10, 5)
WGCNA::plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "Module",
                           dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
grDevices::dev.off()

MEs <- WGCNA::orderMEs(net$MEs)
trait <- tum[match(rownames(datExpr), tum$ID), c("RL_core", "RHI", "Prolif")]
mt_cor <- stats::cor(MEs, trait, use = "pairwise.complete.obs")
mt_p <- WGCNA::corPvalueStudent(mt_cor, nrow(datExpr))
grDevices::pdf(file.path(OUT, "Fig5_WGCNA_module_trait.pdf"), 5.5, 9)
graphics::par(mar = c(5, 9, 3, 2))
WGCNA::labeledHeatmap(Matrix = mt_cor, xLabels = colnames(trait), yLabels = colnames(MEs),
                      ySymbols = colnames(MEs), colorLabels = FALSE, colors = WGCNA::blueWhiteRed(50),
                      textMatrix = paste0(signif(mt_cor, 2), "\n(", signif(mt_p, 1), ")"),
                      setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1), main = "Module-trait")
grDevices::dev.off()

me_ok <- setdiff(colnames(MEs), "MEgrey")
best_me <- me_ok[which.max(abs(mt_cor[me_ok, "RL_core"]))]
mod <- sub("^ME", "", best_me)
mod_genes <- colnames(datExpr)[net$colors == mod]
MM  <- as.numeric(stats::cor(datExpr[, mod_genes], MEs[[best_me]], use = "p"))
GSv <- as.numeric(stats::cor(datExpr[, mod_genes], trait$RL_core, use = "p"))
key_mod_genes <- mod_genes[abs(MM) > PARAM$wgcna_mm & abs(GSv) > PARAM$wgcna_gs]
log_msg("关键模块：", mod, "（与 RL_core r = ", round(mt_cor[best_me, "RL_core"], 2), "，", length(mod_genes),
        " 个基因）；|MM|>", PARAM$wgcna_mm, " 且 |GS|>", PARAM$wgcna_gs, "：", length(key_mod_genes), " 个")
grDevices::pdf(file.path(OUT, "Fig5_WGCNA_MM_vs_GS.pdf"), 5, 5)
WGCNA::verboseScatterplot(abs(MM), abs(GSv), xlab = paste("Module membership in", mod),
                          ylab = "Gene significance for RL_core", col = mod, abline = TRUE)
graphics::abline(v = PARAM$wgcna_mm, h = PARAM$wgcna_gs, lty = 2)
grDevices::dev.off()
utils::write.csv(data.frame(gene = mod_genes, MM = MM, GS = GSv, key = mod_genes %in% key_mod_genes),
                 file.path(OUT, paste0("WGCNA_module_", mod, "_genes.csv")), row.names = FALSE)

# =============================================================================
# 3. 取交集
# =============================================================================
f_scdeg <- file.path(DIR_RES, "03_scRNA", "sc_DEG_genes.txt")
sc_deg <- if (file.exists(f_scdeg)) readLines(f_scdeg) else { log_msg("【核对】未找到 03 步的单细胞差异基因"); character(0) }
pool <- union(GS$RL_all, key_mod_genes)
cand <- Reduce(intersect, list(pool, sc_deg, deg_genes))
rule <- "(RL_all ∪ WGCNA key-module genes) ∩ scRNA pseudobulk DEGs ∩ TCGA tumor-vs-normal DEGs"
if (length(cand) < PARAM$min_candidates) {
  log_msg("【核对】三重交集只有 ", length(cand), " 个基因，少于 ", PARAM$min_candidates,
          "；去掉单细胞这一层，论文方法部分需写明")
  cand <- intersect(pool, deg_genes)
  rule <- "(RL_all ∪ WGCNA key-module genes) ∩ TCGA tumor-vs-normal DEGs (scRNA layer dropped)"
}
log_msg("候选基因：", length(cand), " 个；规则：", rule)

venn <- VennDiagram::venn.diagram(list(`RL_all or WGCNA` = pool, `scRNA DEG` = sc_deg, `TCGA DEG` = deg_genes),
                                  filename = NULL, fill = c("#E41A1C", "#377EB8", "#4DAF4A"), alpha = 0.4,
                                  cex = 1.2, cat.cex = 1, disable.logging = TRUE)
grDevices::pdf(file.path(OUT, "Fig5_Venn.pdf"), 5, 5); grid::grid.draw(venn); grDevices::dev.off()

# =============================================================================
# 4. 单因素 Cox（TCGA，OS，每 1 个标准差）
# =============================================================================
d <- tum[!is.na(tum$OS.time) & tum$OS.time > 0 & !is.na(tum$OS), ]
ex_c <- tcga$expr[cand, d$ID, drop = FALSE]
uni <- do.call(rbind, lapply(cand, function(g) {
  z <- as.numeric(scale(ex_c[g, ]))
  if (all(is.na(z))) return(NULL)
  r <- cox_table(survival::coxph(survival::Surv(d$OS.time, d$OS) ~ z))
  data.frame(gene = g, HR = r$HR, lower = r$lower, upper = r$upper, p = r$p)
}))
uni$FDR <- stats::p.adjust(uni$p, "BH")
uni <- uni[order(uni$p), ]
sig <- uni[uni$p < PARAM$unicox_p, ]
log_msg("单因素 Cox P < ", PARAM$unicox_p, "：", nrow(sig), " 个基因")

if (nrow(sig)) {
  fp <- utils::head(sig, 40); fp$gene <- factor(fp$gene, levels = rev(fp$gene))
  p <- ggplot2::ggplot(fp, ggplot2::aes(HR, gene, color = HR > 1)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2C7FB8"), guide = "none") +
    ggplot2::theme_classic() + ggplot2::labs(x = "HR per SD (TCGA-BLCA OS)", y = NULL)
  save_gg(p, file.path(OUT, "Fig5_unicox_forest.pdf"), 5, max(4, nrow(fp) * 0.22))
}

# 候选基因在各验证队列中是否存在（06 步只用所有队列都有的基因）
val <- c("GSE13507", "GSE32894", "GSE31684", "GSE48075")
avail <- do.call(cbind, lapply(stats::setNames(val, val), function(n)
  cand %in% rownames(readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$expr)))
avail <- data.frame(gene = cand, avail, all_cohorts = apply(avail, 1, all), check.names = FALSE)
log_msg("在 4 个验证队列中都有的候选基因：", sum(avail$all_cohorts), "/", length(cand))

# ---- Table S3 -----------------------------------------------------------------
S3 <- merge(data.frame(gene = cand), deg[, c("gene", "log2FoldChange", "padj")], by = "gene", all.x = TRUE)
pbf <- file.path(DIR_RES, "03_scRNA", "pseudobulk_DE_High_vs_Low.csv")
if (file.exists(pbf)) {
  pbd <- utils::read.csv(pbf)
  S3 <- merge(S3, stats::setNames(pbd[, c("gene", "logFC", "FDR")], c("gene", "sc_logFC_High_vs_Low", "sc_FDR")),
              by = "gene", all.x = TRUE)
}
S3$in_RL_core <- S3$gene %in% GS$RL_core
S3$in_RL_all <- S3$gene %in% GS$RL_all
S3$in_WGCNA_key <- S3$gene %in% key_mod_genes
S3 <- merge(S3, uni, by = "gene", all.x = TRUE)
S3 <- merge(S3, avail[, c("gene", "all_cohorts")], by = "gene", all.x = TRUE)
utils::write.csv(S3[order(S3$p), ], file.path(OUT, "Table_S3_candidates_unicox.csv"), row.names = FALSE)

saveRDS(list(candidates = cand, rule = rule, unicox = uni, unicox_sig = sig$gene, key_module = mod,
             key_module_genes = key_mod_genes, availability = avail),
        file.path(DIR_PRO, "candidates.rds"))
save_session(STEP)
log_msg("05 完成")

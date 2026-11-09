# =============================================================================
# 09_hub_gene_spatial.R —— 核心基因：SHAP → TCGA 验证 → 单细胞定位 → 空间转录组 →（可选）RCTD、虚拟敲除
#   空间数据：GSE319536，22 张 MIBC Visium 切片（每张有 h5 矩阵、图像、坐标，2026-11-07 核对）
#   HPA 免疫组化需手动在 proteinatlas.org 截图（见 README）
# 输出：results/09_hub_gene_spatial/
# =============================================================================
source("R/00_config.R")
need(c("shapviz", "survival", "survminer", "ggplot2", "ggpubr", "Seurat", "UCell"))
STEP <- "09_hub_gene_spatial"; OUT <- res_dir(STEP)
N_SHOW      <- 3                                   # 展示 SHAP 排名前几的基因
SHOW_SLICES <- c("BC01", "BC02", "BC03", "BC04")   # 主图展示哪几张切片（看完全部切片后再选）
RUN_RCTD    <- FALSE                               # 需要 spacexr；耗时较长
RUN_KO      <- FALSE                               # 需要 scTenifoldKnk；耗时较长

GS    <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))
model <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
tcga  <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
tum   <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
genes <- model$genes
log_msg("模型基因：", paste(genes, collapse = ", "))

# =============================================================================
# 1. 模型基因对 RLRS 的贡献（SHAP，精确计算）
#   所选模型是线性核的生存 SVM（survivalsvm 默认 kernel = "lin_kernel"），RLRS 是 16 个基因队列内 z 分数的线性组合。
#   这里在训练集上用线性回归还原系数（R² 必须为 1，否则停止），每个基因对每位患者的 SHAP 值即 w_g ×（z_g − 均值），
#   与 kernelshap 对线性模型的结果相同，但是精确值、几秒钟算完。
#   原方案（另拟随机生存森林再用 kernelshap 近似）解释的是替代模型而非 RLRS 本身，且在 8 GB 机器上要十几个小时，已弃用。
# =============================================================================
d <- tum[!is.na(tum$OS.time) & tum$OS.time > 0 & !is.na(tum$OS), c("ID", "OS.time", "OS")]
mi <- build_mime_input(tcga$expr, tcga$clin[tcga$clin$group == "Tumor", ], model$candidate_genes)   # 与 06 步相同的 z 分数
rl_all <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
mi <- merge(rl_all[rl_all$cohort == "TCGA_BLCA", c("ID", "RLRS")], mi, by = "ID")
X <- as.matrix(mi[, safe_names(genes)]); colnames(X) <- genes
lf <- stats::lm(mi$RLRS ~ X)
r2 <- suppressWarnings(summary(lf)$r.squared)
if (r2 < 0.9999) stop("RLRS 不是模型基因 z 分数的线性组合（R² = ", r2, "），不能用线性分解计算 SHAP")
w <- stats::coef(lf)[-1]; names(w) <- genes
log_msg("RLRS = ", round(stats::coef(lf)[1], 3), " + Σ w × z；线性还原 R² = ", format(r2, digits = 8))
utils::write.csv(data.frame(term = c("(Intercept)", genes), coefficient = unname(stats::coef(lf))),
                 file.path(OUT, "RLRS_coefficients.csv"), row.names = FALSE)
S <- sweep(sweep(X, 2, colMeans(X)), 2, w, `*`)
sv <- shapviz::shapviz(S, X = as.data.frame(X), baseline = mean(mi$RLRS))
imp <- sort(colMeans(abs(S)), decreasing = TRUE)
rank_tab <- data.frame(gene = names(imp), coefficient = unname(w[names(imp)]), mean_abs_SHAP = unname(imp),
                       share_pct = unname(100 * imp / sum(imp)))
print(rank_tab, row.names = FALSE)
utils::write.csv(rank_tab, file.path(OUT, "SHAP_gene_ranking.csv"), row.names = FALSE)
hubs <- utils::head(names(imp), N_SHOW)
hub <- hubs[1]
log_msg("核心基因（平均 |SHAP| 最高）：", hub, "；前 ", N_SHOW, " 位：", paste(hubs, collapse = ", "))
save_print(c(list(shapviz::sv_importance(sv, kind = "beeswarm"), shapviz::sv_importance(sv, kind = "bar")),
             lapply(hubs, function(v) shapviz::sv_dependence(sv, v = v))),
           file.path(OUT, "Fig9_SHAP.pdf"), 6, 5)

# =============================================================================
# 2. TCGA：肿瘤 vs 正常（含配对）、KM；外部队列 KM
# =============================================================================
cl <- tcga$clin
pl <- list()
for (g in hubs) {
  dd <- data.frame(cl[, c("ID", "patient", "group")], expr = as.numeric(tcga$expr[g, cl$ID]))
  dd$group <- factor(dd$group, levels = c("Normal", "Tumor"))
  pat <- intersect(dd$patient[dd$group == "Tumor"], dd$patient[dd$group == "Normal"])
  pl[[paste(g, "TvN")]] <- ggpubr::ggboxplot(dd, "group", "expr", color = "group", add = "jitter",
                                             palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(method = "wilcox.test") + ggplot2::labs(title = g, y = "log2(TPM + 1)")
  pl[[paste(g, "paired")]] <- ggpubr::ggpaired(dd[dd$patient %in% pat, ], x = "group", y = "expr", id = "patient",
                                               color = "group", line.color = "grey70",
                                               palette = c("#2C7FB8", "#D7301F")) +
    ggpubr::stat_compare_means(paired = TRUE) + ggplot2::labs(title = paste(g, "paired"), y = "log2(TPM + 1)")
  s <- data.frame(d, expr = as.numeric(tcga$expr[g, d$ID]))
  s$group <- median_group(s$expr)
  pl[[paste(g, "KM")]] <- km_plot(s, title = paste("TCGA-BLCA", g))
  for (n in c("GSE13507", "GSE31684", "GSE48075")) {
    o <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
    if (!g %in% rownames(o$expr)) next
    e <- o$clin[!is.na(o$clin$OS.time) & o$clin$OS.time > 0 & !is.na(o$clin$OS), c("ID", "OS.time", "OS")]
    e$expr <- as.numeric(o$expr[g, e$ID]); e$group <- median_group(e$expr)
    pl[[paste(g, n)]] <- km_plot(e, title = paste(n, g))
  }
}
save_print(pl, file.path(OUT, "Fig9_hub_TCGA_external.pdf"), 6, 6)

# =============================================================================
# 3. 单细胞：核心基因在各细胞类型、R-loop 高/低恶性细胞中的表达
# =============================================================================
f_sc <- file.path(DIR_PRO, "scRNA_integrated.rds")
f_md <- file.path(DIR_RES, "03_scRNA", "cell_metadata.csv")
if (file.exists(f_sc) && file.exists(f_md)) {
  # 内存：只取出核心基因表达和 UMAP 坐标，随即释放整个单细胞对象；恶性细胞的高/低分组从 03 步的 cell_metadata.csv 读取，
  # 不再另外读入 scRNA_malignant.rds（8 GB 内存的机器上避免同时持有两个大对象）。表达值与原 Seurat 对象的 data 层相同。
  sc <- readRDS(f_sc)
  hs <- intersect(hubs, rownames(sc))
  emb <- Embeddings(sc, "umap")
  hx <- as.matrix(LayerData(sc, assay = "RNA", layer = "data")[hs, , drop = FALSE])
  hd <- data.frame(cell = colnames(sc), UMAP1 = emb[colnames(sc), 1], UMAP2 = emb[colnames(sc), 2], t(hx), check.names = FALSE)
  rm(sc, emb, hx); invisible(gc())
  md <- utils::read.csv(f_md, row.names = 1, check.names = FALSE)
  hd <- cbind(hd, md[hd$cell, c("sample", "cell_class", "RL_group")])
  hd$cell_class <- factor(hd$cell_class)
  long <- do.call(rbind, lapply(hs, function(g) data.frame(gene = g, cell_class = hd$cell_class, expr = hd[[g]])))
  p1 <- ggplot2::ggplot(long, ggplot2::aes(cell_class, expr, fill = cell_class)) +
    ggplot2::geom_violin(scale = "width", linewidth = 0.2) + ggplot2::facet_wrap(~ gene, ncol = 1, scales = "free_y") +
    ggplot2::theme_classic() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "Normalised expression")
  p2 <- lapply(hs, function(g) {
    o <- hd[order(hd[[g]]), ]
    ggplot2::ggplot(o, ggplot2::aes(UMAP1, UMAP2, colour = .data[[g]])) + ggplot2::geom_point(size = 0.05) +
      ggplot2::scale_colour_gradient(low = "grey85", high = "#C0392B") + ggplot2::coord_fixed() +
      ggplot2::theme_classic() + ggplot2::labs(title = g, colour = NULL)
  })
  dp <- do.call(rbind, lapply(hs, function(g) do.call(rbind, lapply(split(hd[[g]], hd$cell_class), function(x)
    data.frame(gene = g, avg = mean(x), pct = mean(x > 0))))))
  dp$cell_class <- rep(levels(hd$cell_class), length(hs))
  p3 <- ggplot2::ggplot(dp, ggplot2::aes(gene, cell_class)) + ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg)) +
    ggplot2::scale_colour_gradient(low = "grey85", high = "#C0392B") + ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL, size = "Fraction expressing", colour = "Mean expression")
  # 恶性细胞高/低组：以样本为单位比较（每个样本内两组的平均表达，配对 Wilcoxon）
  hm <- hd[!is.na(hd$RL_group) & hd$cell_class == "Malignant", ]
  agg <- stats::aggregate(hm[, hs, drop = FALSE], by = list(sample = hm$sample, group = hm$RL_group), FUN = mean)
  agg$group <- factor(agg$group, levels = c("Low", "High"))
  agg <- agg[order(agg$sample, agg$group), ]
  p4 <- lapply(hs, function(g) ggpubr::ggpaired(agg, x = "group", y = g, id = "sample", color = "group",
                                               line.color = "grey70", palette = c(Low = "#2C7FB8", High = "#D7301F")) +
                 ggpubr::stat_compare_means(paired = TRUE) + ggplot2::labs(title = paste(g, "(per-sample mean)")))
  utils::write.csv(agg, file.path(OUT, "hub_scRNA_malignant_High_vs_Low_per_sample.csv"), row.names = FALSE)
  utils::write.csv(dp, file.path(OUT, "hub_scRNA_by_cell_class.csv"), row.names = FALSE)
  save_print(c(list(p1), p2, list(p3), p4), file.path(OUT, "Fig9_hub_scRNA.pdf"), 9, 7)
  rm(hd, hm, long); invisible(gc())
} else log_msg("未找到 03 步的单细胞结果，跳过单细胞定位")

# =============================================================================
# 4. 空间转录组：GSE319536（22 张 MIBC 切片）
# =============================================================================
SP <- data.frame(gsm = paste0("GSM", 9519394:9519415), sample = sprintf("BC%02d", 1:22))
sp_root <- file.path(DIR_RAW, "GEO", "GSE319536")
prep_visium <- function(gsm, smp) {              # 整理成 Load10X_Spatial 需要的目录结构
  sdir <- file.path(sp_root, smp); dir.create(file.path(sdir, "spatial"), recursive = TRUE, showWarnings = FALSE)
  h5 <- file.path(sdir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5)) file.copy(gsm_suppl(gsm, paste0(gsm, "_", smp, "_filtered_feature_bc_matrix.h5"), "GSE319536"), h5)
  for (f in c("scalefactors_json.json", "tissue_hires_image.png", "tissue_lowres_image.png", "tissue_positions.csv")) {
    dest <- file.path(sdir, "spatial", f)
    if (!file.exists(dest)) gunzip_to(gsm_suppl(gsm, paste0(gsm, "_", smp, "_", f, ".gz"), "GSE319536"), dest)
  }
  sdir
}

H <- get_msig("H")
SP_SETS <- c(GS[c("RL_core", "RL_promote", "RL_suppress")],
             list(E2F = H[["HALLMARK_E2F_TARGETS"]], G2M = H[["HALLMARK_G2M_CHECKPOINT"]],
                  IFN_gamma = H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], GEP18 = GEP18))
sp_rows <- list(); sp_plots <- list(); sp_prog <- list()
for (i in seq_len(nrow(SP))) {
  smp <- SP$sample[i]
  f_obj <- file.path(DIR_PRO, "spatial", paste0(smp, ".rds"))
  if (!file.exists(f_obj)) {
    dir.create(dirname(f_obj), showWarnings = FALSE, recursive = TRUE)
    sdir <- prep_visium(SP$gsm[i], smp)
    ob <- Load10X_Spatial(data.dir = sdir, filename = "filtered_feature_bc_matrix.h5", slice = smp)
    ob <- subset(ob, subset = nCount_Spatial > 200)
    ob <- NormalizeData(ob)
    ob <- UCell::AddModuleScore_UCell(ob, features = SP_SETS, maxRank = 2500, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
    ob$RHI <- as.numeric(scale(ob$RL_promote_UCell) - scale(ob$RL_suppress_UCell))
    ob$Prolif <- as.numeric(scale(ob$E2F_UCell) + scale(ob$G2M_UCell)) / 2
    saveRDS(ob, f_obj)
  }
  ob <- readRDS(f_obj)
  # 切片内 spot 水平：R-loop 评分与增殖、免疫程序的相关（以切片为单位汇总，不把 spot 当独立样本）
  sp_prog[[smp]] <- data.frame(sample = smp, n_spots = ncol(ob),
    rho_RLcore_Prolif = stats::cor(ob$RL_core_UCell, ob$Prolif, method = "spearman"),
    rho_RLcore_IFNg   = stats::cor(ob$RL_core_UCell, ob$IFN_gamma_UCell, method = "spearman"),
    rho_RLcore_GEP18  = stats::cor(ob$RL_core_UCell, ob$GEP18_UCell, method = "spearman"),
    rho_Prolif_IFNg   = stats::cor(ob$Prolif, ob$IFN_gamma_UCell, method = "spearman"),
    # 去掉增殖后的偏相关：RL_core 与免疫程序（两者分别对增殖评分回归后取残差）
    prho_RLcore_IFNg_given_Prolif = stats::cor(stats::resid(stats::lm(rank(ob$RL_core_UCell) ~ rank(ob$Prolif))),
                                               stats::resid(stats::lm(rank(ob$IFN_gamma_UCell) ~ rank(ob$Prolif)))))
  hs <- intersect(hubs, rownames(ob))
  for (g in hs) {
    x <- as.numeric(LayerData(ob, assay = "Spatial", layer = "data")[g, ])
    ct <- suppressWarnings(stats::cor.test(ob$RL_core_UCell, x, method = "spearman", exact = FALSE))
    sp_rows[[paste(smp, g)]] <- data.frame(sample = smp, gene = g, n_spots = ncol(ob), pct_detected = mean(x > 0),
                                           rho = unname(ct$estimate), p = ct$p.value)
  }
  if (smp %in% SHOW_SLICES)
    sp_plots[[smp]] <- SpatialFeaturePlot(ob, features = c("RL_core_UCell", "Prolif", hs), ncol = length(hs) + 2,
                                          pt.size.factor = 1.6) + patchwork::plot_annotation(title = smp)
}
prog <- do.call(rbind, sp_prog)
utils::write.csv(prog, file.path(OUT, "Spatial_RLcore_vs_programs_by_slice.csv"), row.names = FALSE)
prog_long <- stats::reshape(prog[, c("sample", "rho_RLcore_Prolif", "rho_RLcore_IFNg", "rho_RLcore_GEP18",
                                     "prho_RLcore_IFNg_given_Prolif")], direction = "long", idvar = "sample",
                            varying = list(2:5), v.names = "rho",
                            times = c("Proliferation", "IFN-gamma", "GEP18", "IFN-gamma | proliferation"), timevar = "program")
log_msg("空间：RL_core 与增殖的 spot 水平相关（22 张切片中位数）= ", round(stats::median(prog$rho_RLcore_Prolif), 2))
p_prog <- ggpubr::ggboxplot(prog_long, "program", "rho", add = "jitter", color = "program") +
  ggplot2::geom_hline(yintercept = 0, lty = 2) + ggpubr::rotate_x_text(30) +
  ggplot2::labs(x = NULL, y = "Spot-level Spearman rho with RL_core (per slice)")
save_gg(p_prog, file.path(OUT, "Fig9_spatial_programs.pdf"), 5, 5)
sp_tab <- do.call(rbind, sp_rows)
sp_tab$FDR <- stats::p.adjust(sp_tab$p, "BH")
utils::write.csv(sp_tab, file.path(OUT, "Spatial_RLcore_vs_hub_spearman.csv"), row.names = FALSE)
p <- ggpubr::ggboxplot(sp_tab, "gene", "rho", color = "gene", add = "jitter") + ggplot2::geom_hline(yintercept = 0, lty = 2) +
  ggplot2::labs(y = "Spot-level Spearman rho with R-loop score (22 slices)")
save_print(c(list(p), sp_plots), file.path(OUT, "Fig9_spatial.pdf"), 4 * (N_SHOW + 2), 5)

# ---- 可选：RCTD 反卷积（用 03 步单细胞作参考），看 R-loop 评分是否落在恶性细胞富集区域 ----
if (RUN_RCTD && requireNamespace("spacexr", quietly = TRUE) && file.exists(f_sc)) {
  sc <- readRDS(f_sc)
  keep <- !is.na(sc$cell_class) & !sc$cell_class %in% c("Tumor_epi_undefined")
  set.seed(PARAM$seed)
  cells <- unlist(lapply(split(colnames(sc)[keep], sc$cell_class[keep]), function(x) sample(x, min(length(x), 1000))))
  ref <- spacexr::Reference(LayerData(sc, assay = "RNA", layer = "counts")[, cells],
                            factor(stats::setNames(sc$cell_class[cells], cells)))
  rm(sc); invisible(gc())
  rctd_rows <- list()
  for (smp in SP$sample) {
    ob <- readRDS(file.path(DIR_PRO, "spatial", paste0(smp, ".rds")))
    coords <- GetTissueCoordinates(ob)[, c("x", "y")]
    puck <- spacexr::SpatialRNA(coords, LayerData(ob, assay = "Spatial", layer = "counts"))
    rctd <- spacexr::run.RCTD(spacexr::create.RCTD(puck, ref, max_cores = 1), doublet_mode = "full")
    w <- as.matrix(spacexr::normalize_weights(rctd@results$weights))
    if ("Malignant" %in% colnames(w)) {
      ct <- stats::cor.test(ob$RL_core_UCell[rownames(w)], w[, "Malignant"], method = "spearman", exact = FALSE)
      rctd_rows[[smp]] <- data.frame(sample = smp, rho_RLcore_vs_malignant = unname(ct$estimate), p = ct$p.value)
    }
    saveRDS(w, file.path(DIR_PRO, "spatial", paste0(smp, "_RCTD_weights.rds")))
  }
  utils::write.csv(do.call(rbind, rctd_rows), file.path(OUT, "RCTD_RLcore_vs_malignant.csv"), row.names = FALSE)
}

# ---- 可选：scTenifoldKnk 虚拟敲除核心基因（恶性细胞）→ 受影响基因做 GO 富集 -------
if (RUN_KO && requireNamespace("scTenifoldKnk", quietly = TRUE) && file.exists(f_sc)) {
  mal <- readRDS(file.path(DIR_PRO, "scRNA_malignant.rds"))
  set.seed(PARAM$seed)
  cnt <- LayerData(mal, assay = "RNA", layer = "counts")
  cnt <- cnt[, sample(ncol(cnt), min(1500, ncol(cnt)))]
  cnt <- cnt[Matrix::rowMeans(cnt > 0) >= 0.1 | rownames(cnt) == hub, ]
  ko <- scTenifoldKnk::scTenifoldKnk(countMatrix = as.matrix(cnt), gKO = hub, qc = TRUE, nc_nNet = 10, nc_nCells = 500)
  dr <- ko$diffRegulation
  utils::write.csv(dr, file.path(OUT, paste0("virtualKO_", hub, "_diffRegulation.csv")), row.names = FALSE)
  hit <- setdiff(dr$gene[dr$p.adj < 0.05], hub)
  log_msg("虚拟敲除 ", hub, " 后显著受影响的基因：", length(hit))
  if (length(hit) >= 10 && requireNamespace("clusterProfiler", quietly = TRUE)) {
    ego <- clusterProfiler::enrichGO(hit, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                                     universe = dr$gene, pvalueCutoff = 0.05)
    utils::write.csv(as.data.frame(ego), file.path(OUT, paste0("virtualKO_", hub, "_GO_BP.csv")), row.names = FALSE)
    if (nrow(as.data.frame(ego))) save_gg(enrichplot::dotplot(ego, showCategory = 15),
                                          file.path(OUT, paste0("Fig9_virtualKO_", hub, "_GO.pdf")), 7, 6)
  }
}

log_msg("HPA：请在 https://www.proteinatlas.org 搜索 ", paste(hubs, collapse = "、"),
        "，截取 Tissue（urinary bladder）与 Pathology（urothelial cancer）的免疫组化图，图注写明抗体编号")
save_session(STEP)
log_msg("09 完成")

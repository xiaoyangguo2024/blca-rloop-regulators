# =============================================================================
# 11_figures.R —— 论文主图：统一主题重画（不拼接各步的 PDF）
#   每张图输出 600 dpi PNG（results/figures/FigureN.png）和单页 PDF；面板标 A/B/C，图上不写 "Figure"
#   图中标注的统计量同时写入 results/figures/FigureN_stats.csv，正文数字以此为准
# 用法：Rscript R/11_figures.R            # 全部
#       Rscript R/11_figures.R 4 6        # 只画 Figure 4 和 6
# 图 1（研究设计）是手写 SVG：results/figures/Figure1.svg，这里用 rsvg-convert 转成 PNG
# =============================================================================
source("R/00_config.R")
need(c("ggplot2", "patchwork", "ragg", "survival", "scales"))
FIG <- res_dir("figures")
args <- commandArgs(trailingOnly = TRUE)
WHICH <- if (length(args)) args else c(1:9, paste0("S", 1:3))   # 例：Rscript R/11_figures.R 3 S2

FONT <- "Arial"
BASE <- 7
theme_pub <- function(base = BASE) {
  ggplot2::theme_classic(base_size = base, base_family = FONT) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = base - 0.5),
                   axis.line = ggplot2::element_line(linewidth = 0.3), axis.ticks = ggplot2::element_line(linewidth = 0.3),
                   plot.title = ggplot2::element_text(size = base + 0.5, hjust = 0.5),
                   plot.tag = ggplot2::element_text(size = base + 4, face = "bold", family = FONT),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = base, face = "bold"),
                   legend.key.size = ggplot2::unit(3, "mm"), legend.text = ggplot2::element_text(size = base - 0.5),
                   legend.title = ggplot2::element_text(size = base - 0.5))
}
ggplot2::theme_set(theme_pub())
C_HI <- "#C0392B"; C_LO <- "#2E6DA4"; C_GREY <- "grey55"
PAL2 <- c(High = C_HI, Low = C_LO)

# 图上：P < 0.0001 时只写 "< 0.0001"；统计表（*_stats.csv）里用 exact_p() 保留精确值
fmt_p <- function(p, digits = 2) {
  ifelse(is.na(p), "NA", ifelse(p < 1e-4, "<0.0001", formatC(signif(p, digits), format = "fg", digits = digits)))
}
lab_p <- function(p) ifelse(!is.na(p) & p < 1e-4, "P < 0.0001", paste0("P = ", fmt_p(p)))
exact_p <- function(p) paste0("P = ", signif(p, 2))
lab_q <- function(q) ifelse(!is.na(q) & q < 1e-4, "FDR < 0.0001", paste0("FDR = ", fmt_p(q)))

STATS <- list()
add_stat <- function(fig, panel, what, value, detail = "") {
  STATS[[length(STATS) + 1]] <<- data.frame(figure = fig, panel = panel, statistic = what,
                                            value = as.character(value), detail = detail)
}

save_fig <- function(p, n, w = 180, h = 200) {      # 宽高单位 mm；J Cancer 双栏宽约 17–18 cm
  f <- file.path(FIG, paste0("Figure", n))
  tmp <- tempfile(fileext = ".png")        # ragg 不能直接写含中文的路径：先写到临时目录再复制
  ragg::agg_png(tmp, width = w, height = h, units = "mm", res = 600, background = "white")
  print(p); grDevices::dev.off()
  stopifnot(file.copy(tmp, paste0(f, ".png"), overwrite = TRUE)); unlink(tmp)
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = w / 25.4, height = h / 25.4, family = FONT)
  print(p); grDevices::dev.off()
  st <- do.call(rbind, Filter(function(x) x$figure == n, STATS))
  if (!is.null(st)) utils::write.csv(st, paste0(f, "_stats.csv"), row.names = FALSE)
  log_msg("Figure ", n, " 已保存（", w, " × ", h, " mm）")
}

# Kaplan–Meier（ggplot 重画）：组内中位数分组，注明 log-rank P 和 HR（高 vs 低）
km_gg <- function(time, status, group, title = NULL, xlab = "Years", ylab = "Survival probability",
                  cols = PAL2, years = TRUE) {
  t <- if (years) time / 365.25 else time
  d <- data.frame(t = t, s = status, g = factor(group, levels = names(cols)))
  d <- d[stats::complete.cases(d), ]
  fit <- survival::survfit(survival::Surv(t, s) ~ g, data = d)
  st <- summary(fit, censored = TRUE)
  grp <- sub("^g=", "", names(fit$strata))
  curve <- do.call(rbind, lapply(seq_along(grp), function(i) {
    ss <- fit[i]
    data.frame(g = grp[i], t = c(0, ss$time), surv = c(1, ss$surv), cens = c(0, ss$n.censor))
  }))
  lr <- survival::survdiff(survival::Surv(t, s) ~ g, data = d)
  p <- stats::pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  cx <- summary(survival::coxph(survival::Surv(t, s) ~ relevel(g, ref = names(cols)[2]), data = d))$conf.int
  ntab <- table(d$g)
  lab <- sprintf("HR %.2f (%.2f-%.2f)\nLog-rank %s", cx[1, 1], cx[1, 3], cx[1, 4], lab_p(p))
  curve$g <- factor(curve$g, levels = names(cols),
                    labels = sprintf("%s (n = %d)", names(cols), as.integer(ntab[names(cols)])))
  cl <- stats::setNames(cols, levels(curve$g))
  # 图例放左下（曲线从左上角出发，左下角总是空的）；HR 放右侧：后三分之一时间内曲线都低于 75% 时放右上，否则放右下
  late <- curve$surv[curve$t >= 2 / 3 * max(curve$t)]
  hr_top <- length(late) > 0 && max(late) < 0.75
  gg <- ggplot2::ggplot(curve, ggplot2::aes(t, surv, colour = g)) +
    ggplot2::geom_step(linewidth = 0.45) +
    ggplot2::geom_point(data = curve[curve$cens > 0, ], shape = 3, size = 0.6, stroke = 0.3) +
    ggplot2::annotate("text", x = Inf, y = if (hr_top) 1 else 0.02, hjust = 1.03, vjust = if (hr_top) 1 else 0, label = lab,
                      size = (BASE - 1) / ggplot2::.pt, family = FONT, lineheight = 0.9) +
    ggplot2::scale_colour_manual(values = cl, name = NULL) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.01, 0.01), legend.justification = c(0, 0),
                   legend.background = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = BASE - 1.5),
                   legend.key.height = ggplot2::unit(2.5, "mm"))
  attr(gg, "stats") <- c(HR = cx[1, 1], lower = cx[1, 3], upper = cx[1, 4], p_logrank = p)
  gg
}

# 森林图（HR 与 95% CI，对数刻度）
forest_gg <- function(df, label = "label", group = NULL, xlab = "Hazard ratio per SD (95% CI)", cols = NULL) {
  df$label <- factor(df[[label]], levels = rev(unique(df[[label]])))
  aes_c <- if (is.null(group)) ggplot2::aes(HR, label) else ggplot2::aes(HR, label, colour = .data[[group]])
  g <- ggplot2::ggplot(df, aes_c) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.4, orientation = "y",
                           position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 1.3, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_x_log10() + ggplot2::labs(x = xlab, y = NULL)
  if (!is.null(cols)) g <- g + ggplot2::scale_colour_manual(values = cols, name = NULL)
  g
}

# =============================================================================
# Figure 1：研究设计（ggplot 画流程图，与其他图同一字体和导出方式）
#   坐标单位 mm，画布 180 × 121；方框内容只写设计，不写结果
#   行高按正文行数计算：标题 + n 行正文 ≈ 5 + 2.7n + 1.5 mm
# =============================================================================
fig1 <- function() {
  W <- 180; H <- 121
  FILL <- c(sc = "#F7E1DD", bulk = "#DDE9F5", sp = "#E0F0E1", gs = "#EDEDED", mod = "#FBF0D5", ev = "#FFFFFF")
  box <- function(id, x0, x1, y1, h, fill, title, body)
    data.frame(id = id, x0 = x0, x1 = x1, y1 = y1, y0 = y1 - h, fill = fill, title = title, body = body)
  # 行的上沿（自上而下）：A 120, B 105.5, C 85.5, D 63, E 45.5, F 16.5
  bx <- rbind(
    box("gs", 2, 178, 120, 9.5, "gs", "R-loop regulator gene sets (R-loopBase v2.0)",
        paste0("RL_core: 351 tier 1–3 regulators  ·  RL_all: 1,709  ·  RL_v1: 1,185  ·  ",
               "R-loop homeostasis index (RHI) = promoting (40) − suppressing (250)")),
    box("sc", 2, 58, 105.5, 15, "sc", "Single-cell RNA-seq", "GSE135337, GSE222315\n16 tumours, 5 normal samples"),
    box("bulk", 62, 118, 105.5, 15, "bulk", "Bulk transcriptomes",
        "TCGA-BLCA (406 tumours, 19 normal)\n+ 7 GEO and treatment cohorts\n1,545 patients with outcome data"),
    box("sp", 122, 178, 105.5, 15, "sp", "Spatial transcriptomics",
        "GSE319536\n22 muscle-invasive tumour sections\n(10x Genomics Visium)"),
    box("sc2", 2, 58, 85.5, 17.5, "sc", "Tumour epithelial cells",
        "Harmony integration, annotation\ncopyKAT with normal-epithelium control\nUCell scores; R-loop-high vs -low cells\n(pseudobulk DE, cell cycle, CellChat)"),
    box("bulk2", 62, 118, 85.5, 17.5, "bulk", "R-loop scores in bulk tumours",
        "ssGSEA scores vs proliferation,\nmutations, consensus classes\nand survival in 8 cohorts"),
    box("sp2", 122, 178, 85.5, 17.5, "sp", "Spot-level programmes",
        "R-loop score vs proliferation\nand immune programmes\n(sections as units of analysis)"),
    box("cand", 2, 118, 63, 12.5, "mod", "Candidate genes",
        "(RL_all or WGCNA key-module genes) also differentially expressed in scRNA-seq\nand in TCGA tumour vs normal tissue; univariable Cox regression in TCGA"),
    box("mime", 2, 170, 45.5, 24, "mod", "Mime: 10 survival algorithms, 117 combinations",
        "Within-cohort z-scores; fixed-direction C-index → R-loop-related risk score (RLRS)"),
    box("train", 4, 56.67, 36, 12.5, "ev", "Training", "TCGA-BLCA (OS)"),
    box("tune", 59.67, 112.33, 36, 12.5, "ev", "Tuning (model selection only)", "GSE13507, GSE31684 (OS)"),
    box("test", 115.33, 168, 36, 12.5, "ev", "Independent testing",
        "GSE32894 (DSS), GSE48075 (OS)\nIMvigor210, GSE176307, GSE169455"),
    box("spec", 2, 41.75, 16.5, 15, "ev", "Specificity",
        "1,000 random gene sets × 2\nproliferation adjustment\n97 published signatures"),
    box("clin", 44.75, 84.5, 16.5, 15, "ev", "Clinical model", "Multivariable Cox\nnomogram, calibration\ndecision curves"),
    box("ther", 87.5, 127.25, 16.5, 15, "ev", "Treatment response",
        "ICI: IMvigor210, GSE176307\nchemotherapy: GSE169455\nIOBR, TIDE, oncoPredict"),
    box("genes", 130.25, 170, 16.5, 15, "mod", "Model genes",
        "SHAP (linear model)\ntumour vs normal tissue\nsingle-cell and spatial"))
  bx$xm <- (bx$x0 + bx$x1) / 2
  fx <- bx$xm[match(c("spec", "clin", "ther", "genes"), bx$id)]
  ar <- data.frame(x    = c(30, 90, 150, 30, 90, 150, 30, 90, 60, fx),
                   y    = c(110.5, 110.5, 110.5, 90.5, 90.5, 90.5, 68, 68, 50.5, rep(21.5, 4)),
                   xend = c(30, 90, 150, 30, 90, 150, 30, 90, 60, fx),
                   yend = c(105.5, 105.5, 105.5, 85.5, 85.5, 85.5, 63, 63, 45.5, rep(16.5, 4)))
  # 空间转录组 → 模型基因：沿右侧空白走虚线
  sp_line <- data.frame(x = c(175, 175), y = c(68, 9))
  sp_arr  <- data.frame(x = 175, y = 9, xend = 170, yend = 9)
  pt <- ggplot2::.pt
  arw <- ggplot2::arrow(length = ggplot2::unit(1.4, "mm"), type = "closed")
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bx, ggplot2::aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill),
                       colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_text(data = bx, ggplot2::aes(x = xm, y = y1 - 1.6, label = title), vjust = 1, fontface = "bold",
                       size = BASE / pt, family = FONT) +
    ggplot2::geom_text(data = bx, ggplot2::aes(x = xm, y = y1 - 5, label = body), vjust = 1, lineheight = 0.95,
                       size = (BASE - 0.5) / pt, family = FONT) +
    ggplot2::geom_segment(data = ar, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35,
                          colour = "grey30", arrow = arw) +
    ggplot2::geom_path(data = sp_line, ggplot2::aes(x, y), linewidth = 0.35, colour = "grey30", linetype = "22") +
    ggplot2::geom_segment(data = sp_arr, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.35,
                          colour = "grey30", linetype = "22", arrow = arw) +
    ggplot2::scale_fill_manual(values = FILL, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
    ggplot2::coord_fixed() + ggplot2::theme_void()
  save_fig(p, 1, W, H)
}

# =============================================================================
# Figure 4：bulk 队列中的 R-loop 评分
# =============================================================================
fig4 <- function() {
  tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
  sc_all <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
  tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
  tc <- merge(tcga$clin, sc_all$TCGA_BLCA, by = "ID")

  # A：配对的肿瘤与癌旁
  pat <- intersect(tc$patient[tc$group == "Tumor"], tc$patient[tc$group == "Normal"])
  pr <- tc[tc$patient %in% pat, c("patient", "group", "RL_core", "RHI")]
  pr <- pr[!duplicated(pr[, c("patient", "group")]), ]
  long <- rbind(data.frame(pr[, c("patient", "group")], score = "RL_core", value = pr$RL_core),
                data.frame(pr[, c("patient", "group")], score = "RHI", value = pr$RHI))
  long$group <- factor(long$group, levels = c("Normal", "Tumor"), labels = c("Normal", "Tumour"))
  long$score <- factor(long$score, levels = c("RL_core", "RHI"))
  pv <- vapply(levels(long$score), function(s) {
    w <- stats::reshape(long[long$score == s, c("patient", "group", "value")], idvar = "patient", timevar = "group",
                        direction = "wide")
    stats::wilcox.test(w$value.Tumour, w$value.Normal, paired = TRUE)$p.value
  }, numeric(1))
  unp <- vapply(c("RL_core", "RHI"), function(s) stats::wilcox.test(tc[[s]] ~ tc$group)$p.value, numeric(1))
  for (s in names(pv)) {
    add_stat(4, "A", paste(s, "paired tumour vs normal, Wilcoxon signed-rank P"), signif(pv[[s]], 3), paste(length(pat), "pairs"))
    add_stat(4, "A", paste(s, "all tumour vs normal, Wilcoxon P"), signif(unp[[s]], 3),
             sprintf("%d tumours, %d normals", sum(tc$group == "Tumor"), sum(tc$group == "Normal")))
  }
  pA <- ggplot2::ggplot(long, ggplot2::aes(group, value)) +
    ggplot2::geom_line(ggplot2::aes(group = patient), colour = "grey75", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 0.9) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = group), fill = NA, width = 0.35, outlier.shape = NA, linewidth = 0.35) +
    ggplot2::facet_wrap(~ score, scales = "free_y", nrow = 1) +
    ggplot2::geom_text(data = data.frame(score = factor(names(pv), levels = levels(long$score)), lab = lab_p(pv)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = (BASE - 0.5) / ggplot2::.pt,
                       inherit.aes = FALSE, family = FONT) +
    ggplot2::scale_colour_manual(values = c(Normal = C_LO, Tumour = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::labs(x = NULL, y = "ssGSEA score", title = sprintf("TCGA paired samples (n = %d)", length(pat)))

  # B：相关热图
  vars <- c(RL_core = "RL_core", RL_all = "RL_all", RL_v1 = "RL_v1", RHI = "RHI", E2F = "E2F targets",
            G2M = "G2M checkpoint", Prolif = "Proliferation", DNA_repair = "DNA repair", IFN_gamma = "IFN-γ response",
            IFN_alpha = "IFN-α response", cGAS_STING_KEGG = "Cytosolic DNA sensing", GEP18 = "T cell-inflamed GEP",
            EMT = "EMT")
  t1 <- tum[tum$group == "Tumor", ]
  cm <- stats::cor(t1[, names(vars)], method = "spearman", use = "pairwise.complete.obs")
  rows <- c("RL_core", "RL_all", "RL_v1", "RHI")
  hm <- data.frame(x = rep(colnames(cm), each = length(rows)), y = rep(rows, ncol(cm)),
                   r = as.vector(cm[rows, ]))
  hm$x <- factor(vars[hm$x], levels = vars); hm$y <- factor(hm$y, levels = rev(rows))
  for (v in c("Prolif", "E2F", "G2M", "EMT", "IFN_gamma", "GEP18", "DNA_repair"))
    for (s in c("RL_core", "RHI")) add_stat(4, "C", sprintf("Spearman rho %s vs %s (TCGA tumours)", s, v), round(cm[s, v], 2), paste(nrow(t1), "tumours"))
  pB <- ggplot2::ggplot(hm, ggplot2::aes(x, y, fill = r)) + ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_fill_gradient2(low = C_LO, mid = "white", high = C_HI, limits = c(-1, 1), name = "Spearman ρ") +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::coord_fixed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), legend.position = "right")

  # 只进统计表、不画在主图上的关联（补充图 S3 用）：分期、分级、consensus 分子亚型、RHI 与突变
  mt <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "Rloop_score_by_mutation.csv"))
  for (s in c("RL_core", "RHI")) {
    st <- t1[!is.na(t1$stage_bin), ]
    add_stat(4, "suppl", paste(s, "stage III-IV vs I-II, Wilcoxon P"), signif(stats::wilcox.test(st[[s]] ~ st$stage_bin)$p.value, 3),
             sprintf("median %.3f vs %.3f", stats::median(st[[s]][st$stage_bin == "III-IV"]), stats::median(st[[s]][st$stage_bin == "I-II"])))
    gr <- t1[t1$grade %in% c("High", "Low"), ]
    add_stat(4, "suppl", paste(s, "high vs low grade, Wilcoxon P"), signif(stats::wilcox.test(gr[[s]] ~ gr$grade)$p.value, 3),
             sprintf("n = %d high, %d low", sum(gr$grade == "High"), sum(gr$grade == "Low")))
    cs <- t1[!is.na(t1$consensus), ]
    md_c <- tapply(cs[[s]], cs$consensus, stats::median)
    add_stat(4, "suppl", paste(s, "by consensus class, Kruskal-Wallis P"), signif(stats::kruskal.test(cs[[s]] ~ cs$consensus)$p.value, 3),
             paste(sprintf("%s %.3f", names(md_c), md_c), collapse = "; "))
  }
  for (g in c("FGFR3", "KDM6A", "ELF3", "TP53", "RB1")) {
    v <- t1[[paste0("mut_", g)]]; ok <- !is.na(v)
    add_stat(4, "suppl", paste("RHI median, mutant vs wild-type:", g),
             sprintf("%.2f vs %.2f", stats::median(t1$RHI[ok & v == "Mut"]), stats::median(t1$RHI[ok & v == "WT"])),
             paste0("FDR = ", signif(mt$FDR_RHI[mt$gene == g], 2)))
  }

  # C：RL_core 与增殖评分
  rho <- stats::cor.test(t1$RL_core, t1$Prolif, method = "spearman", exact = FALSE)
  add_stat(4, "B", "Spearman rho RL_core vs proliferation score", round(rho$estimate, 2),
           paste0(exact_p(rho$p.value), "; n = ", nrow(t1)))
  pC <- ggplot2::ggplot(t1, ggplot2::aes(Prolif, RL_core)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6, colour = "grey30") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = C_HI, linewidth = 0.5, se = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, family = FONT, size = (BASE - 0.5) / ggplot2::.pt,
                      label = sprintf("ρ = %.2f, %s", rho$estimate, lab_p(rho$p.value))) +
    ggplot2::labs(x = "Proliferation score (E2F + G2M)", y = "RL_core score", title = "TCGA-BLCA tumours")

  # D：突变（TP53、RB1、ERCC2）
  mt <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "Rloop_score_by_mutation.csv"))
  genes_m <- c("TP53", "RB1", "ERCC2", "ARID1A")
  ml <- do.call(rbind, lapply(genes_m, function(g) data.frame(gene = g, status = t1[[paste0("mut_", g)]], RL_core = t1$RL_core)))
  ml <- ml[!is.na(ml$status), ]
  ml$status <- factor(ml$status, levels = c("WT", "Mut")); ml$gene <- factor(ml$gene, levels = genes_m)
  labd <- data.frame(gene = factor(genes_m, levels = genes_m),
                     lab = lab_q(mt$FDR_RL_core[match(genes_m, mt$gene)]))
  for (g in genes_m) add_stat(4, "D", paste("RL_core mutant vs wild-type, Wilcoxon FDR:", g),
                              signif(mt$FDR_RL_core[mt$gene == g], 3),
                              sprintf("%d mutant / %d tested", mt$n_mut[mt$gene == g], sum(!is.na(t1[[paste0("mut_", g)]]))))
  pD <- ggplot2::ggplot(ml, ggplot2::aes(status, RL_core, colour = status)) +
    ggplot2::geom_jitter(width = 0.2, size = 0.35, alpha = 0.6) +
    ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.35, colour = "black") +
    ggplot2::facet_wrap(~ gene, nrow = 1) +
    ggplot2::geom_text(data = labd, ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::scale_colour_manual(values = c(WT = C_LO, Mut = C_HI), guide = "none") +
    ggplot2::labs(x = NULL, y = "RL_core score")

  # E：APOBEC 富集与 TMB
  ap <- t1[!is.na(t1$APOBEC_enrich), ]
  r_ap <- stats::cor.test(ap$RL_core, ap$APOBEC_enrich, method = "spearman", exact = FALSE)
  tm <- t1[!is.na(t1$TMB), ]
  r_tmb <- stats::cor.test(tm$RL_core, tm$TMB, method = "spearman", exact = FALSE)
  r_ap_rhi <- stats::cor.test(ap$RHI, ap$APOBEC_enrich, method = "spearman", exact = FALSE)
  add_stat(4, "E", "Spearman rho RL_core vs APOBEC enrichment", round(r_ap$estimate, 2), paste0(exact_p(r_ap$p.value), "; n = ", nrow(ap)))
  add_stat(4, "E", "Spearman rho RHI vs APOBEC enrichment", round(r_ap_rhi$estimate, 2), paste0(exact_p(r_ap_rhi$p.value), "; n = ", nrow(ap)))
  add_stat(4, "E", "Spearman rho RL_core vs TMB", round(r_tmb$estimate, 2), paste0(exact_p(r_tmb$p.value), "; n = ", nrow(tm)))
  pE <- ggplot2::ggplot(ap, ggplot2::aes(APOBEC_enrich, RL_core)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6, colour = "grey30") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = C_HI, linewidth = 0.5, se = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, family = FONT, size = (BASE - 0.5) / ggplot2::.pt,
                      label = sprintf("ρ = %.2f, %s", r_ap$estimate, lab_p(r_ap$p.value))) +
    ggplot2::labs(x = "APOBEC enrichment score", y = "RL_core score")

  # F：TCGA 中各 R-loop 评分的 Cox（每 SD），未校正与校正（年龄、分期、增殖评分）
  cx <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "TCGA_Cox_Rloop_scores.csv"))
  cx <- cx[cx$endpoint == "OS", ]
  # 未校正 CI 需要重算（04 步只存了 HR 与 P）
  d_os <- t1[!is.na(t1$OS.time) & t1$OS.time > 0 & !is.na(t1$OS), ]
  fr <- do.call(rbind, lapply(cx$score, function(s) {
    u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ scale(d_os[[s]]), data = d_os))
    rbind(data.frame(score = s, model = "Unadjusted", HR = u$HR, lower = u$lower, upper = u$upper, p = u$p),
          data.frame(score = s, model = "Adjusted", HR = cx$HR_adj[cx$score == s], lower = cx$lower_adj[cx$score == s],
                     upper = cx$upper_adj[cx$score == s], p = cx$p_adj[cx$score == s]))
  }))
  for (i in seq_len(nrow(fr))) add_stat(4, "F", sprintf("TCGA OS HR per SD, %s, %s", fr$score[i], fr$model[i]),
                                        sprintf("%.2f (%.2f-%.2f)", fr$HR[i], fr$lower[i], fr$upper[i]), exact_p(fr$p[i]))
  fr$label <- fr$score
  fr$model <- factor(fr$model, levels = c("Adjusted", "Unadjusted"))     # 躲避位置：未校正在上
  pF <- forest_gg(fr, group = "model") +
    ggplot2::geom_text(ggplot2::aes(x = max(fr$upper) * 1.25, label = fmt_p(p)), hjust = 0, show.legend = FALSE,
                       size = (BASE - 1) / ggplot2::.pt, family = FONT, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = c(Unadjusted = "grey35", Adjusted = C_HI), name = NULL,
                                 breaks = c("Unadjusted", "Adjusted"),
                                 labels = c("Unadjusted", "Adjusted for age, stage, proliferation")) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme(plot.margin = ggplot2::margin(5, 22, 5, 5)) +
    ggplot2::labs(title = "TCGA-BLCA, overall survival") +
    ggplot2::theme(legend.position = "bottom")

  # G：外部队列中 RL_core 的 Cox
  ex <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "External_Cox_RL_core.csv"))
  ex <- rbind(data.frame(cohort = "TCGA-BLCA", n = nrow(d_os), events = sum(d_os$OS),
                         fr[fr$score == "RL_core" & fr$model == "Unadjusted", c("HR", "lower", "upper", "p")]),
              data.frame(ex[, c("cohort", "n", "events")], HR = ex$HR_perSD, ex[, c("lower", "upper", "p")]))
  ex$label <- sprintf("%s (%d/%d)", ex$cohort, ex$events, ex$n)
  ex$label[ex$cohort == "GSE32894"] <- paste0(ex$label[ex$cohort == "GSE32894"], " DSS")
  for (i in seq_len(nrow(ex))) add_stat(4, "G", paste("RL_core HR per SD,", ex$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", ex$HR[i], ex$lower[i], ex$upper[i]), exact_p(ex$p[i]))
  pG <- forest_gg(ex) + ggplot2::geom_point(size = 1.3, colour = C_HI) +
    ggplot2::geom_text(ggplot2::aes(x = max(ex$upper) * 1.25, label = fmt_p(p)), hjust = 0, size = (BASE - 1) / ggplot2::.pt,
                       family = FONT) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme(plot.margin = ggplot2::margin(5, 22, 5, 5)) +
    ggplot2::labs(title = "RL_core score across cohorts (events/patients)")

  top <- (pA | pC) + patchwork::plot_layout(widths = c(1.15, 1))
  mid <- pB
  low <- (pD | pE) + patchwork::plot_layout(widths = c(2.2, 1))
  bot <- (pF | pG) + patchwork::plot_layout(widths = c(1, 1.1))
  p <- (top / mid / low / bot) + patchwork::plot_layout(heights = c(1, 0.55, 0.9, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 4, 180, 230)
}

# =============================================================================
# 单细胞数据：只读一次，取出作图需要的部分后释放（整个对象约 1.7 GB）
# =============================================================================
CT_NAMES <- c(Epithelial = "Epithelial", T_cell = "T cell", NK = "Cytotoxic T/NK", B_cell = "B cell", Plasma = "Plasma cell",
              Myeloid = "Myeloid", DC = "Dendritic cell", Mast = "Mast cell", Fibroblast = "Fibroblast",
              Myofib_SMC = "Myofibroblast/SMC", Endothelial = "Endothelial")
CT_COLS <- c(Epithelial = "#C0392B", `T cell` = "#2E6DA4", `Cytotoxic T/NK` = "#5DADE2", `B cell` = "#8E44AD",
             `Plasma cell` = "#C39BD3", Myeloid = "#E67E22", `Dendritic cell` = "#F5B041", `Mast cell` = "#A04000",
             Fibroblast = "#27AE60", `Myofibroblast/SMC` = "#82E0AA", Endothelial = "#B7950B")
MK_SHOW <- list(Epithelial = c("EPCAM", "KRT8", "KRT19", "UPK1B"), `T cell` = c("CD3D", "CD3E"),
                `Cytotoxic T/NK` = c("CD8A", "NKG7", "GNLY"), `B cell` = c("CD79A", "MS4A1"), `Plasma cell` = c("MZB1", "JCHAIN"),
                Myeloid = c("CD68", "CD14", "LYZ"), Fibroblast = c("COL1A1", "DCN"), `Myofibroblast/SMC` = c("ACTA2", "RGS5"),
                Endothelial = c("PECAM1", "VWF"))
SC <- NULL
load_sc <- function() {
  if (is.null(SC)) {
    need(c("Seurat", "Matrix"))
    sc <- readRDS(file.path(DIR_PRO, "scRNA_integrated.rds"))
    md <- sc@meta.data
    emb <- SeuratObject::Embeddings(sc, "umap")
    md$UMAP1 <- emb[rownames(md), 1]; md$UMAP2 <- emb[rownames(md), 2]
    md$ct <- factor(CT_NAMES[md$celltype], levels = intersect(CT_NAMES, CT_NAMES[unique(md$celltype)]))
    g <- intersect(unique(unlist(MK_SHOW)), rownames(sc))
    ex <- SeuratObject::LayerData(sc, assay = "RNA", layer = "data")[g, rownames(md), drop = FALSE]
    rm(sc); invisible(gc())
    SC <<- list(md = md, mk = ex)
  }
  SC
}

# =============================================================================
# Figure 2：单细胞全景、copyKAT（含正常上皮阴性对照）、R-loop 评分
# =============================================================================
fig2 <- function() {
  s <- load_sc(); md <- s$md
  set.seed(PARAM$seed); md <- md[sample(nrow(md)), ]           # 打乱绘制顺序，避免某一类总压在最上层
  cen <- stats::aggregate(cbind(UMAP1, UMAP2) ~ ct, data = md, FUN = stats::median)
  add_stat(2, "A", "cells after QC", nrow(md), sprintf("%d samples (%d tumour, %d normal)", length(unique(md$sample)),
                                                      length(unique(md$sample[md$tissue == "Tumor"])),
                                                      length(unique(md$sample[md$tissue == "Normal"]))))
  for (k in levels(md$ct)) add_stat(2, "A", paste("cells:", k), sum(md$ct == k), sprintf("%.1f%%", 100 * mean(md$ct == k)))
  pA <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = ct)) + ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggrepel::geom_text_repel(data = cen, ggplot2::aes(label = ct), colour = "black", size = BASE / ggplot2::.pt, family = FONT,
                             seed = 1, box.padding = 0.5, min.segment.length = 0, segment.size = 0.2,
                             bg.color = "white", bg.r = 0.12) +
    ggplot2::scale_colour_manual(values = CT_COLS, guide = "none") + ggplot2::coord_fixed() +
    ggplot2::labs(title = sprintf("%s cells, %d samples", format(nrow(md), big.mark = ","), length(unique(md$sample)))) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  md$tissue_lab <- factor(ifelse(md$tissue == "Tumor", "Tumour", "Normal/adjacent"), levels = c("Tumour", "Normal/adjacent"))
  pB <- ggplot2::ggplot(md[order(md$tissue_lab == "Normal/adjacent"), ], ggplot2::aes(UMAP1, UMAP2, colour = tissue_lab)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +      # 正常/癌旁细胞较少，画在上层
    ggplot2::scale_colour_manual(values = c(Tumour = C_HI, `Normal/adjacent` = C_LO), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) + ggplot2::coord_fixed() +
    ggplot2::labs(title = "Tissue of origin") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0))

  # C：marker 点图（平均表达按基因缩放到 0–1，点大小为表达比例）
  ex <- s$mk; ct <- s$md$ct
  mk <- unlist(MK_SHOW); mk <- mk[mk %in% rownames(ex)]
  dp <- do.call(rbind, lapply(levels(ct), function(k) {
    m <- ex[mk, ct == k, drop = FALSE]
    data.frame(ct = k, gene = mk, avg = Matrix::rowMeans(m), pct = Matrix::rowMeans(m > 0))
  }))
  dp$avg_s <- stats::ave(dp$avg, dp$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  dp$gene <- factor(dp$gene, levels = unique(mk)); dp$ct <- factor(dp$ct, levels = rev(levels(ct)))
  pC <- ggplot2::ggplot(dp, ggplot2::aes(gene, ct)) + ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 2.6, name = "Fraction\nexpressing", breaks = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   legend.box = "horizontal")

  # D：copyKAT 判定（肿瘤样本上皮细胞）+ 阴性对照（正常样本上皮细胞，03c）
  te <- md[md$celltype == "Epithelial" & md$tissue == "Tumor", ]
  te$call <- ifelse(te$copykat %in% "aneuploid", "Aneuploid", ifelse(te$copykat %in% "diploid", "Diploid", "Not analysed/undefined"))
  cb <- as.data.frame(table(sample = te$sample, call = te$call))
  cb$set <- "Tumour"
  f_ctl <- file.path(DIR_RES, "03_scRNA", "copykat_normal_control.csv")
  if (file.exists(f_ctl)) {
    ctl <- utils::read.csv(f_ctl)
    cb <- rbind(cb, data.frame(sample = rep(ctl$sample, 3), call = rep(c("Aneuploid", "Diploid", "Not analysed/undefined"), each = nrow(ctl)),
                               Freq = c(ctl$n_aneuploid, ctl$n_diploid, ctl$n_undefined), set = "Normal"))
    for (i in seq_len(nrow(ctl))) add_stat(2, "D", paste("copyKAT negative control, aneuploid fraction:", ctl$sample[i]),
                                           ctl$aneuploid_frac[i], paste(ctl$n_epithelial[i], "normal epithelial cells"))
  }
  tot <- stats::aggregate(Freq ~ sample + set, data = cb, FUN = sum)
  cb$frac <- cb$Freq / tot$Freq[match(paste(cb$sample, cb$set), paste(tot$sample, tot$set))]
  an <- cb[cb$call == "Aneuploid" & cb$set == "Tumour", ]
  for (i in seq_len(nrow(an))) add_stat(2, "D", paste("aneuploid fraction of tumour epithelial cells:", an$sample[i]),
                                        round(an$frac[i], 3), paste(tot$Freq[tot$sample == an$sample[i] & tot$set == "Tumour"], "epithelial cells"))
  add_stat(2, "D", "tumour epithelial cells called aneuploid (all samples)", sum(te$call == "Aneuploid"),
           sprintf("of %d (%.1f%%)", nrow(te), 100 * mean(te$call == "Aneuploid")))
  ord <- tot$sample[order(tot$set != "Tumour", tot$sample)]
  cb$sample <- factor(cb$sample, levels = unique(ord))
  cb$set <- factor(cb$set, levels = c("Tumour", "Normal"))
  pD <- ggplot2::ggplot(cb, ggplot2::aes(sample, frac, fill = call)) + ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(~ set, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c(Aneuploid = C_HI, Diploid = C_LO, `Not analysed/undefined` = "grey80"), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Epithelial cells (copyKAT)", title = "copyKAT calls") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5),
                   legend.position = "top", strip.clip = "off")

  # E：RL_core UCell 评分
  pE <- ggplot2::ggplot(md[order(md$RL_core_UCell), ], ggplot2::aes(UMAP1, UMAP2, colour = RL_core_UCell)) +
    ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggplot2::scale_colour_gradientn(colours = c("grey88", "#F5B7B1", C_HI, "#641E16"), name = "RL_core\n(UCell)") +
    ggplot2::coord_fixed() + ggplot2::labs(title = "R-loop regulator activity") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  # F：每个样本的平均 RL_core：肿瘤上皮细胞（肿瘤样本的全部上皮细胞）vs 正常上皮（样本为单位）
  #    cell_class "Malignant" 在 MALIGNANT_RULE = "all_tumor_epi" 下即肿瘤样本的全部上皮细胞；P2–P5 有配对的癌旁样本，连线
  cls <- c(Malignant = "Tumour", Normal_epi = "Normal")
  sub_md <- md[md$cell_class %in% names(cls), ]
  agg <- stats::aggregate(RL_core_UCell ~ sample + cell_class, data = sub_md, FUN = mean)
  nn <- stats::aggregate(RL_core_UCell ~ sample + cell_class, data = sub_md, FUN = length)
  agg <- agg[nn$RL_core_UCell[match(paste(agg$sample, agg$cell_class), paste(nn$sample, nn$cell_class))] >= 20, ]   # 每组至少 20 个细胞
  agg$cls <- factor(cls[agg$cell_class], levels = cls)
  agg$patient <- ifelse(grepl("-(BCa|NAT)$", agg$sample), sub("-(BCa|NAT)$", "", agg$sample), agg$sample)
  p_mn <- stats::wilcox.test(agg$RL_core_UCell[agg$cell_class == "Malignant"], agg$RL_core_UCell[agg$cell_class == "Normal_epi"])$p.value
  add_stat(2, "F", "per-sample mean RL_core, tumour vs normal epithelium, Wilcoxon P", signif(p_mn, 3),
           sprintf("%d tumour samples vs %d normal samples", sum(agg$cell_class == "Malignant"), sum(agg$cell_class == "Normal_epi")))
  add_stat(2, "F", "median of per-sample mean RL_core: tumour / normal epithelium",
           sprintf("%.3f / %.3f", stats::median(agg$RL_core_UCell[agg$cell_class == "Malignant"]),
                   stats::median(agg$RL_core_UCell[agg$cell_class == "Normal_epi"])))
  pr <- intersect(agg$patient[agg$cell_class == "Malignant"], agg$patient[agg$cell_class == "Normal_epi"])
  pair <- agg[agg$patient %in% pr, ]
  for (pp in pr) add_stat(2, "F", paste("paired patient, RL_core tumour / normal epithelium:", pp),
                          sprintf("%.3f / %.3f", pair$RL_core_UCell[pair$patient == pp & pair$cell_class == "Malignant"],
                                  pair$RL_core_UCell[pair$patient == pp & pair$cell_class == "Normal_epi"]))
  pF <- ggplot2::ggplot(agg, ggplot2::aes(cls, RL_core_UCell)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    ggplot2::geom_line(data = pair, ggplot2::aes(group = patient), colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = cls), size = 1, position = ggplot2::position_jitter(width = 0.06, seed = 1)) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 1.3, label = lab_p(p_mn), size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(C_HI, C_LO), cls), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = "Epithelial cells", y = "RL_core (per-sample mean)", title = "Samples as the unit")

  # 手动加标签、逐排包装：避免下排 y 轴标题跟 C 的长标签对齐而留出空白
  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), list(pA, pB, pC, pD, pE, pF), LETTERS[1:6])
  p <- patchwork::wrap_plots(patchwork::wrap_elements(full = pl[[1]] | pl[[2]]), patchwork::wrap_elements(full = pl[[3]]),
                             patchwork::wrap_elements(full = (pl[[4]] | pl[[5]] | pl[[6]]) + patchwork::plot_layout(widths = c(1.35, 1, 0.7))),
                             ncol = 1, heights = c(1.05, 0.62, 1))
  save_fig(p, 2, 180, 215)
}

# =============================================================================
# Figure 3：R-loop 高/低肿瘤上皮细胞（伪 bulk 差异、Hallmark、细胞周期、细胞通讯）
# =============================================================================
fig3 <- function(d3 = file.path(DIR_RES, "03_scRNA"), id = 3, cells = "Tumour epithelial cells",
                 all_lab = "All tumour epithelial cells",
                 lab_genes = c("FEN1", "CDCA7", "FANCI", "CHEK1", "CDC6", "TYMS", "HELLS", "BRCA1", "GINS2", "MCM10", "POLQ",
                               "PIGR", "LCN2", "SPRR1A", "UPK1A", "UPK2", "CXCL17")) {
  GS <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))
  de <- utils::read.csv(file.path(d3, "pseudobulk_DE_High_vs_Low.csv"))
  de$sig <- ifelse(de$FDR < PARAM$sc_deg_fdr & de$logFC > PARAM$sc_deg_lfc, "Up in R-loop-high",
                   ifelse(de$FDR < PARAM$sc_deg_fdr & de$logFC < -PARAM$sc_deg_lfc, "Up in R-loop-low", "NS"))
  add_stat(id, "A", "pseudobulk DEGs up in R-loop-high", sum(de$sig == "Up in R-loop-high"), "FDR < 0.05, log2FC > 0.5")
  add_stat(id, "A", "pseudobulk DEGs up in R-loop-low", sum(de$sig == "Up in R-loop-low"), "FDR < 0.05, log2FC < -0.5")
  add_stat(id, "A", "DEGs that are R-loop regulators (RL_all)", sum(de$sig != "NS" & de$gene %in% GS$RL_all),
           sprintf("of %d DEGs; RL_core: %d", sum(de$sig != "NS"), sum(de$sig != "NS" & de$gene %in% GS$RL_core)))
  lab <- de[de$sig != "NS" & de$gene %in% lab_genes, ]                      # 正文提到的基因
  if (nrow(lab) < 8) lab <- rbind(lab, utils::head(de[de$sig != "NS" & !de$gene %in% lab$gene, ][order(de$FDR[de$sig != "NS" & !de$gene %in% lab$gene]), ], 8 - nrow(lab)))
  pA <- ggplot2::ggplot(de, ggplot2::aes(logFC, -log10(FDR), colour = sig)) + ggplot2::geom_point(size = 0.35) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * PARAM$sc_deg_lfc, linetype = 2, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(PARAM$sc_deg_fdr), linetype = 2, linewidth = 0.25, colour = "grey60") +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), size = (BASE - 1.5) / ggplot2::.pt, family = FONT,
                             max.overlaps = Inf, segment.size = 0.2, fontface = "italic", show.legend = FALSE,
                             box.padding = 0.35, point.padding = 0.15, min.segment.length = 0, force = 3, seed = 1) +
    ggplot2::scale_colour_manual(values = c(`Up in R-loop-high` = C_HI, `Up in R-loop-low` = C_LO, NS = "grey80"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.2))) +      # 上方留空给基因标签
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.14))) +
    ggplot2::labs(x = "log2 fold change (R-loop-high vs -low)", y = "-log10 FDR", title = paste0(cells, ", pseudobulk")) +
    ggplot2::theme(legend.position = "bottom") + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))

  # B：Hallmark，全部肿瘤上皮细胞 vs 仅 G1 期细胞
  hb <- utils::read.csv(file.path(d3, "Hallmark_all_vs_G1only.csv"))
  # 先纳入任一分析中 FDR < 0.05 的通路，再按差异绝对值补足到 18 条（正文提到的显著通路都要在图上）
  ord <- hb[order(-abs(hb$mean_diff_all)), ]
  keep <- unique(c(ord$pathway[ord$FDR_all < 0.05 | ord$FDR_G1 < 0.05], ord$pathway))
  top <- hb[match(keep[seq_len(min(max(18, sum(ord$FDR_all < 0.05 | ord$FDR_G1 < 0.05)), length(keep)))], hb$pathway), ]
  hl <- rbind(data.frame(pathway = top$pathway, set = all_lab, diff = top$mean_diff_all, FDR = top$FDR_all),
              data.frame(pathway = top$pathway, set = "G1-phase cells only", diff = top$mean_diff_G1, FDR = top$FDR_G1))
  hl$pathway <- factor(gsub("_", " ", hl$pathway), levels = gsub("_", " ", top$pathway[order(top$mean_diff_all)]))
  hl$sig <- ifelse(hl$FDR < 0.05, "FDR < 0.05", "NS")
  for (i in seq_len(nrow(hb))) if (hb$FDR_all[i] < 0.05 || hb$FDR_G1[i] < 0.05)
    add_stat(id, "B", paste("Hallmark High-Low difference (all / G1):", hb$pathway[i]),
             sprintf("%.4f / %.4f", hb$mean_diff_all[i], hb$mean_diff_G1[i]),
             sprintf("FDR %s / %s", signif(hb$FDR_all[i], 2), signif(hb$FDR_G1[i], 2)))
  add_stat(id, "B", "Hallmark pathways with FDR < 0.05 (all tumour epithelial / G1 only)",
           sprintf("%d / %d", sum(hb$FDR_all < 0.05), sum(hb$FDR_G1 < 0.05)))
  pB <- ggplot2::ggplot(hl, ggplot2::aes(diff, pathway, fill = set, alpha = sig)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = stats::setNames(c(C_HI, "#7F8C8D"), c(all_lab, "G1-phase cells only")), name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`FDR < 0.05` = 1, NS = 0.35), name = NULL) +
    ggplot2::labs(x = "UCell difference, R-loop-high minus -low\n(per-sample mean)", y = NULL) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical", axis.text.y = ggplot2::element_text(size = BASE - 1.5))

  # C：细胞周期
  cc <- utils::read.csv(file.path(d3, "cellcycle_by_sample.csv"))
  pc <- stats::wilcox.test(cc$pct_cycling_High, cc$pct_cycling_Low, paired = TRUE)$p.value
  add_stat(id, "C", "S/G2M cells, R-loop-high vs -low, median %",
           sprintf("%.1f vs %.1f", stats::median(cc$pct_cycling_High), stats::median(cc$pct_cycling_Low)),
           paste0(exact_p(pc), " (paired Wilcoxon, ", nrow(cc), " samples)"))
  add_stat(id, "C", "within-sample Spearman rho RL_core vs G2M score, median (range)",
           sprintf("%.2f (%.2f to %.2f)", stats::median(cc$rho_RLcore_G2M), min(cc$rho_RLcore_G2M), max(cc$rho_RLcore_G2M)))
  add_stat(id, "C", "within-sample Spearman rho RL_core vs S score, median (range)",
           sprintf("%.2f (%.2f to %.2f)", stats::median(cc$rho_RLcore_S), min(cc$rho_RLcore_S), max(cc$rho_RLcore_S)))
  cl <- rbind(data.frame(sample = cc$sample, grp = "Low", pct = cc$pct_cycling_Low),
              data.frame(sample = cc$sample, grp = "High", pct = cc$pct_cycling_High))
  cl$grp <- factor(cl$grp, levels = c("Low", "High"))
  pC <- ggplot2::ggplot(cl, ggplot2::aes(grp, pct)) +
    ggplot2::geom_line(ggplot2::aes(group = sample), colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = grp), size = 1.1) +
    ggplot2::annotate("text", x = 1.5, y = Inf, vjust = 1.3, label = lab_p(pc), size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_colour_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = "R-loop group", y = "Cells in S or G2/M phase (%)")
  rl <- rbind(data.frame(sample = cc$sample, score = "S phase", rho = cc$rho_RLcore_S),
              data.frame(sample = cc$sample, score = "G2/M phase", rho = cc$rho_RLcore_G2M))
  pD <- ggplot2::ggplot(rl, ggplot2::aes(score, rho)) + ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.35) +
    ggplot2::geom_jitter(width = 0.1, size = 1, colour = C_HI) +
    ggplot2::labs(x = NULL, y = "Spearman ρ with RL_core\n(within sample)")

  # E：CellChat，肿瘤上皮细胞 → 其他细胞
  f_cc <- file.path(d3, "CellChat_count_weight_from_malignant.csv")
  pE <- NULL
  if (file.exists(f_cc)) {
    cw <- utils::read.csv(f_cc)
    cw <- cw[!cw$target %in% c("Malig_High", "Malig_Low"), ]
    cw$grp <- factor(ifelse(cw$source == "Malig_High", "High", "Low"), levels = c("Low", "High"))
    cw$target <- factor(ifelse(cw$target %in% names(CT_NAMES), CT_NAMES[cw$target], cw$target))
    tot <- stats::aggregate(cbind(count, weight) ~ grp, data = cw, FUN = sum)
    add_stat(id, "E", "CellChat interactions sent by tumour epithelial cells to non-epithelial cells (Low / High)",
             sprintf("%d / %d", tot$count[tot$grp == "Low"], tot$count[tot$grp == "High"]),
             sprintf("total weight %.3f / %.3f", tot$weight[tot$grp == "Low"], tot$weight[tot$grp == "High"]))
    pE <- ggplot2::ggplot(cw, ggplot2::aes(count, target, fill = grp)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::scale_fill_manual(values = PAL2, name = "R-loop group") +
      ggplot2::labs(x = "Inferred interactions (CellChat)", y = NULL, title = paste(cells, "as senders")) +
      ggplot2::theme(legend.position = "bottom")
  }
  # 手动加标签，再把每一排包成一个整体：避免下排的坐标轴跟上排火山图对齐而留出大块空白
  pl <- Filter(Negate(is.null), list(pA, pB, pC, pD, pE))
  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), pl, LETTERS[seq_along(pl)])
  top <- patchwork::wrap_elements(full = (pl[[1]] | pl[[2]]) + patchwork::plot_layout(widths = c(1, 1.05)))
  mid <- patchwork::wrap_elements(full = pl[[3]] | pl[[4]])
  p <- if (length(pl) == 5) patchwork::wrap_plots(top, mid, patchwork::wrap_elements(full = pl[[5]]), ncol = 1, heights = c(1.15, 0.6, 0.55)) else
    patchwork::wrap_plots(top, mid, ncol = 1, heights = c(1.15, 0.6))
  save_fig(p, id, 180, 200)
}
figS2 <- function() fig3(file.path(DIR_RES, "03_scRNA", "sensitivity_copykat_rule"), "S2",
                         cells = "Aneuploid epithelial cells", all_lab = "All aneuploid epithelial cells")

# =============================================================================
# 共用：各队列的 RLRS 与临床数据（06 步产出）
# =============================================================================
COH_LAB  <- c(TCGA_BLCA = "TCGA-BLCA", GSE13507 = "GSE13507", GSE31684 = "GSE31684", GSE32894 = "GSE32894",
              GSE48075 = "GSE48075", IMvigor210 = "IMvigor210", GSE176307 = "GSE176307", GSE169455 = "GSE169455")
COH_ROLE <- c(TCGA_BLCA = "training", GSE13507 = "tuning", GSE31684 = "tuning", GSE32894 = "test", GSE48075 = "test",
              IMvigor210 = "test, ICI", GSE176307 = "test, ICI", GSE169455 = "test, chemotherapy")
COH_COLS <- c(`TCGA-BLCA` = "grey30", GSE13507 = "#2E6DA4", GSE31684 = "#5DADE2", GSE32894 = "#C0392B",
              GSE48075 = "#E59866", IMvigor210 = "#8E44AD", GSE176307 = "#C39BD3")
SURV_COH <- c("TCGA_BLCA", "GSE13507", "GSE31684", "GSE32894", "GSE48075")
coh_clin <- function(n) {
  if (n == "TCGA_BLCA") readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds")) else readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
}
rlrs_all <- function() readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
rlrs_cohort <- function(n) {
  r <- rlrs_all(); r <- r[r$cohort == n, c("ID", "RLRS", "RLRS_z", "RLRS_group")]
  merge(r, coh_clin(n), by = "ID")
}
txt <- function(...) ggplot2::annotate("text", size = (BASE - 0.5) / ggplot2::.pt, family = FONT, ...)

# =============================================================================
# Figure 5：候选基因与 RLRS 的构建、验证
# =============================================================================
fig5 <- function() {
  GS <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))
  cands <- readRDS(file.path(DIR_PRO, "candidates.rds"))
  model <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
  d6 <- file.path(DIR_RES, "06_ml_signature")
  sc_deg <- readLines(file.path(DIR_RES, "03_scRNA", "sc_DEG_genes.txt"))
  uc <- utils::read.csv(file.path(DIR_RES, "05_candidates", "Table_S3_candidates_unicox.csv"))

  # A：候选基因筛选（各步剩余基因数）
  pool <- union(GS$RL_all, cands$key_module_genes)
  n_steps <- c(length(pool), length(intersect(pool, sc_deg)), length(cands$candidates), sum(uc$all_cohorts),
               sum(uc$all_cohorts & uc$p < PARAM$unicox_p), length(model$genes))
  lab_steps <- c("RL_all or WGCNA key-module genes", "+ DEG, R-loop-high vs -low cells", "+ DEG, tumour vs normal (TCGA)",
                 "+ measured in all validation cohorts", "+ univariate Cox P < 0.05 (TCGA)", "Genes in the selected model")
  for (i in seq_along(n_steps)) add_stat(5, "A", paste("candidate selection:", lab_steps[i]), n_steps[i])
  st <- data.frame(step = factor(lab_steps, levels = rev(lab_steps)), n = n_steps)
  pA <- ggplot2::ggplot(st, ggplot2::aes(y = step)) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = n, yend = step), linewidth = 3.2, colour = "#D98880") +
    ggplot2::geom_text(ggplot2::aes(x = n, label = format(n, big.mark = ",")), hjust = -0.2,
                       size = (BASE - 0.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(x = "Genes (log scale)", y = NULL, title = "Candidate genes")

  # B：16 个模型基因的单因素 Cox（TCGA）
  ug <- uc[uc$gene %in% model$genes, ]
  ug <- ug[order(-ug$HR), ]
  ug$dir <- ifelse(ug$HR > 1, "Higher risk", "Lower risk")
  add_stat(5, "B", "model genes with HR > 1 in TCGA (univariate, per SD)", sum(ug$HR > 1), paste("of", nrow(ug)))
  pB <- forest_gg(ug, label = "gene", group = "dir", xlab = "HR per SD (95% CI)",
                  cols = c(`Higher risk` = C_HI, `Lower risk` = C_LO)) +
    ggplot2::labs(title = "Model genes, TCGA-BLCA OS") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "none")

  # C：全部算法组合：调参队列 vs 检验队列的平均 C-index
  ci <- utils::read.csv(file.path(d6, "Table_S4_Cindex_101_models.csv"), check.names = FALSE)
  ci$sel <- ci$Model == model$model
  prof <- apply(round(ci[, SURV_COH], 6), 1, paste, collapse = "_")
  rk <- rank(-ci$mean_test)[ci$sel]
  rt <- stats::cor.test(ci$mean_tune, ci$mean_test, method = "spearman", exact = FALSE)
  add_stat(5, "C", "algorithm combinations evaluated", nrow(ci), sprintf("%d distinct C-index profiles", length(unique(prof))))
  add_stat(5, "C", "selected model", model$model,
           sprintf("mean C-index tuning %.3f, test %.3f; %d genes", ci$mean_tune[ci$sel], ci$mean_test[ci$sel], ci$n_genes[ci$sel]))
  add_stat(5, "C", "rank of the selected model by mean test C-index", rk,
           sprintf("of %d; mean test C-index range %.3f-%.3f", nrow(ci), min(ci$mean_test), max(ci$mean_test)))
  add_stat(5, "C", "mean tuning C-index range across combinations", sprintf("%.3f-%.3f", min(ci$mean_tune), max(ci$mean_tune)))
  add_stat(5, "C", "Spearman rho, mean tuning vs mean test C-index", sprintf("%.2f", rt$estimate), exact_p(rt$p.value))
  pC <- ggplot2::ggplot(ci, ggplot2::aes(mean_tune, mean_test)) +
    ggplot2::geom_point(data = ci[!ci$sel, ], colour = "grey60", size = 0.8, alpha = 0.8) +
    ggplot2::geom_point(data = ci[ci$sel, ], colour = C_HI, size = 1.8) +
    txt(x = ci$mean_tune[ci$sel], y = ci$mean_test[ci$sel], label = "Selected", hjust = 0.75, vjust = 2.1, colour = C_HI) +
    txt(x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, label = sprintf("Spearman ρ = %.2f", rt$estimate)) +
    ggplot2::labs(x = "Mean C-index, tuning cohorts", y = "Mean C-index, test cohorts",
                  title = sprintf("%d algorithm combinations", nrow(ci)))

  # D：所选模型在各队列的 C-index（方向固定；95% CI）
  cx <- utils::read.csv(file.path(d6, "RLRS_Cox_by_cohort.csv"))
  cd <- data.frame(cohort = cx$cohort, C = cx$Cindex, lo = cx$Cindex - 1.96 * cx$Cindex_se, hi = cx$Cindex + 1.96 * cx$Cindex_se)
  for (n in c("IMvigor210", "GSE176307")) {
    d <- rlrs_cohort(n); d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
    cc <- survival::concordance(survival::Surv(OS.time, OS) ~ RLRS, data = d, reverse = TRUE)
    cd <- rbind(cd, data.frame(cohort = n, C = cc$concordance, lo = cc$concordance - 1.96 * sqrt(cc$var),
                               hi = cc$concordance + 1.96 * sqrt(cc$var)))
  }
  for (i in seq_len(nrow(cd))) add_stat(5, "D", paste("C-index (fixed direction):", cd$cohort[i]),
                                        sprintf("%.3f (%.3f-%.3f)", cd$C[i], cd$lo[i], cd$hi[i]), COH_ROLE[cd$cohort[i]])
  cd$lab <- factor(COH_LAB[cd$cohort], levels = rev(COH_LAB[cd$cohort]))
  cd$role <- factor(sub(",.*", "", COH_ROLE[cd$cohort]), levels = c("training", "tuning", "test"))
  pD <- ggplot2::ggplot(cd, ggplot2::aes(C, lab, colour = role)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.4, orientation = "y") +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = c(training = "grey30", tuning = C_LO, test = C_HI), name = NULL) +
    ggplot2::labs(x = "C-index (95% CI)", y = NULL, title = "Selected model") +
    ggplot2::theme(legend.position = "bottom")

  # J：时间依赖 AUC
  au <- utils::read.csv(file.path(d6, "RLRS_timeROC_AUC.csv"))
  for (i in seq_len(nrow(au))) add_stat(5, "J", sprintf("time-dependent AUC, %d year(s): %s", au$year[i], au$cohort[i]), sprintf("%.3f", au$AUC[i]))
  au$lab <- factor(COH_LAB[au$cohort], levels = COH_LAB[SURV_COH])
  pE <- ggplot2::ggplot(au, ggplot2::aes(factor(year), AUC, colour = lab, group = lab)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = COH_COLS, name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0.45, 0.9)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = "Years", y = "Time-dependent AUC", title = "Discrimination over time") +
    ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(2.5, "mm"),
                   legend.text = ggplot2::element_text(size = BASE - 1.5), legend.margin = ggplot2::margin(0, 0, 0, 0))

  # E–I：KM（队列内中位数分组）
  kms <- lapply(SURV_COH, function(n) {
    d <- rlrs_cohort(n); d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
    g <- km_gg(d$OS.time, d$OS, d$RLRS_group, title = sprintf("%s (%s)", COH_LAB[n], COH_ROLE[n]),
               ylab = if (n == "GSE32894") "Disease-specific survival" else "Overall survival")
    s <- attr(g, "stats")
    add_stat(5, "E-I", paste("KM RLRS high vs low:", n), sprintf("HR %.2f (%.2f-%.2f)", s["HR"], s["lower"], s["upper"]),
             paste0("log-rank ", exact_p(s["p_logrank"]), "; n = ", nrow(d), ", events = ", sum(d$OS)))
    g
  })
  # 版式：A|B；C|D；KM 两行三列（第六格放时间依赖 AUC）。下半部分单独排版后整体嵌入（wrap_elements），
  # 避免上半部分的长标签影响下半部分的对齐；因此面板字母手动标注
  tagged <- Map(function(g, t) g + ggplot2::labs(tag = t), c(list(pA, pB, pC, pD), kms, list(pE)), LETTERS[1:10])
  top <- (tagged[[1]] | tagged[[2]]) / (patchwork::free(tagged[[3]]) | tagged[[4]])
  bottom <- patchwork::wrap_elements(full = patchwork::wrap_plots(tagged[5:10], ncol = 3))
  p <- patchwork::wrap_plots(list(patchwork::wrap_elements(full = top), bottom), ncol = 1, heights = c(1.05, 1))
  save_fig(p, 5, 180, 225)
}

# =============================================================================
# Figure 6：特异性：随机基因集、增殖、已发表模型
# =============================================================================
adv_stage <- function(n, c0) {     # 与 06b 相同：TCGA 为 III–IV 期；混合队列为 T2–4；MIBC 队列为 T3–4
  if (n == "TCGA_BLCA") return(ifelse(c0$stage %in% c("III", "IV"), 1L, ifelse(c0$stage %in% c("I", "II"), 0L, NA)))
  if (!"T" %in% colnames(c0)) return(rep(NA_integer_, nrow(c0)))
  hi <- if (n %in% c("GSE13507", "GSE32894")) c("2", "3", "4") else c("3", "4")
  ifelse(c0$T %in% hi, 1L, ifelse(is.na(c0$T), NA, 0L))
}
fig6 <- function() {
  d6b <- file.path(DIR_RES, "06b_specificity"); d6 <- file.path(DIR_RES, "06_ml_signature")
  # A：随机基因集零分布（16 个基因，与 RLRS 基因用同样的 ridge Cox 拟合）
  nl <- readRDS(file.path(DIR_PRO, "random_signature_null.rds"))
  bm <- utils::read.csv(file.path(d6b, "random_signature_benchmark.csv"))
  obs <- mean(nl$obs)
  nA <- rowMeans(nl$nullA); nB <- rowMeans(nl$nullB)
  pv <- c((1 + sum(nA >= obs, na.rm = TRUE)) / (1 + sum(!is.na(nA))), (1 + sum(nB >= obs, na.rm = TRUE)) / (1 + sum(!is.na(nB))))
  if (any(abs(pv - bm$p_RLRS_ridge) > 1e-6)) log_msg("【核对】Figure 6A 的经验 P 与 06b 不一致")
  nm <- c("Random expressed genes", "Random prognostic genes")
  dn <- rbind(data.frame(null = nm[1], C = nA), data.frame(null = nm[2], C = nB))
  lb <- data.frame(null = nm, p = pv)
  for (i in 1:2) add_stat(6, "A", paste("random-signature null:", nm[i]),
                          sprintf("mean %.3f, 95th percentile %.3f", mean(c(nA, nB)[dn$null == nm[i]]),
                                  stats::quantile(dn$C[dn$null == nm[i]], 0.95)),
                          sprintf("RLRS genes (ridge) %.3f, empirical P = %.3g; selected model %.3f, P = %.3g",
                                  obs, pv[i], bm$RLRS_best_mean[i], bm$p_RLRS_best[i]))
  pA <- ggplot2::ggplot(dn, ggplot2::aes(C)) + ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = NA) +
    ggplot2::geom_vline(xintercept = obs, colour = C_HI, linewidth = 0.5) +
    ggplot2::geom_text(data = lb, ggplot2::aes(x = obs, y = Inf, label = paste0("RLRS genes\nempirical ", lab_p(p))),
                       hjust = 1.05, vjust = 1.3, size = (BASE - 1) / ggplot2::.pt, family = FONT, lineheight = 0.9) +
    ggplot2::facet_wrap(~ null, ncol = 1) +
    ggplot2::labs(x = "Mean C-index, four validation cohorts", y = "Random 16-gene sets", title = "Random-signature benchmark")

  # B/C：RLRS 与增殖评分
  sc <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
  r <- rlrs_all()
  t1 <- merge(r[r$cohort == "TCGA_BLCA", ], sc$TCGA_BLCA[, c("ID", "Prolif")], by = "ID")
  rr <- stats::cor.test(t1$RLRS, t1$Prolif, method = "spearman", exact = FALSE)
  pB <- ggplot2::ggplot(t1, ggplot2::aes(Prolif, RLRS_z)) + ggplot2::geom_point(size = 0.45, colour = "grey45", alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = C_HI, linewidth = 0.5) +
    txt(x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, label = sprintf("ρ = %.2f", rr$estimate)) +
    ggplot2::labs(x = "Proliferation score", y = "RLRS (z-score)", title = "TCGA-BLCA")
  adj <- utils::read.csv(file.path(d6b, "Table_S6_proliferation_adjusted.csv"))
  for (i in seq_len(nrow(adj))) add_stat(6, "C", paste("Spearman rho RLRS vs proliferation / RL_core:", adj$cohort[i]),
                                         sprintf("%.2f / %.2f", adj$rho_RLRS_Prolif[i], adj$rho_RLRS_RLcore[i]))
  rl <- rbind(data.frame(cohort = adj$cohort, what = "Proliferation", rho = adj$rho_RLRS_Prolif),
              data.frame(cohort = adj$cohort, what = "RL_core", rho = adj$rho_RLRS_RLcore))
  rl$lab <- factor(COH_LAB[rl$cohort], levels = rev(COH_LAB[adj$cohort]))
  pC <- ggplot2::ggplot(rl, ggplot2::aes(rho, lab, fill = what)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::scale_fill_manual(values = c(Proliferation = C_HI, RL_core = "grey60"), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = "Spearman ρ with RLRS", y = NULL, title = "Correlation with RLRS") +
    ggplot2::theme(legend.position = "bottom")

  # D：C-index：增殖评分、RLRS、两者合并（各队列内拟合的 Cox 模型）
  cc <- rbind(data.frame(cohort = adj$cohort, model = "Proliferation", C = adj$C_Prolif),
              data.frame(cohort = adj$cohort, model = "RLRS", C = adj$C_RLRS),
              data.frame(cohort = adj$cohort, model = "Proliferation + RLRS", C = adj$C_Prolif_RLRS))
  for (i in seq_len(nrow(adj))) add_stat(6, "D", paste("Cox C-index proliferation / RLRS / both:", adj$cohort[i]),
                                         sprintf("%.3f / %.3f / %.3f", adj$C_Prolif[i], adj$C_RLRS[i], adj$C_Prolif_RLRS[i]),
                                         paste("LRT for adding RLRS", exact_p(adj$LRT_p_add_RLRS[i])))
  cc$lab <- factor(COH_LAB[cc$cohort], levels = rev(COH_LAB[adj$cohort]))
  cc$model <- factor(cc$model, levels = c("Proliferation", "RLRS", "Proliferation + RLRS"))
  pD <- ggplot2::ggplot(cc, ggplot2::aes(C, lab, colour = model, shape = model)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(size = 1.5, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_colour_manual(values = c(Proliferation = "grey50", RLRS = C_HI, `Proliferation + RLRS` = "#641E16"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Proliferation = 16, RLRS = 17, `Proliferation + RLRS` = 15), name = NULL) +
    ggplot2::labs(x = "C-index", y = NULL, title = "Added value over proliferation") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")

  # E：RLRS 的 HR：未校正、校正增殖、校正增殖 + 分期（复算以取得 95% CI；HR 与 06b 表核对）
  fr <- list()
  for (n in adj$cohort) {
    c0 <- coh_clin(n)
    d <- merge(data.frame(ID = c0$ID, OS.time = c0$OS.time, OS = c0$OS, adv = adv_stage(n, c0)),
               r[r$cohort == n, c("ID", "RLRS")], by = "ID")
    d <- merge(d, sc[[n]][, c("ID", "Prolif")], by = "ID")
    d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
    d$RLRS_z <- as.numeric(scale(d$RLRS)); d$Prolif_z <- as.numeric(scale(d$Prolif))
    one <- function(f, what) {
      m <- survival::coxph(f, data = d); s <- summary(m)
      data.frame(cohort = n, model = what, HR = s$conf.int["RLRS_z", 1], lower = s$conf.int["RLRS_z", 3],
                 upper = s$conf.int["RLRS_z", 4], p = s$coefficients["RLRS_z", 5])
    }
    o <- rbind(one(survival::Surv(OS.time, OS) ~ RLRS_z, "Unadjusted"),
               one(survival::Surv(OS.time, OS) ~ RLRS_z + Prolif_z, "Adjusted for proliferation"))
    if (sum(!is.na(d$adv)) > 30 && length(unique(stats::na.omit(d$adv))) == 2)
      o <- rbind(o, one(survival::Surv(OS.time, OS) ~ RLRS_z + Prolif_z + adv, "Adjusted for proliferation and stage"))
    if (abs(o$HR[2] - adj$HR_RLRS_adjProlif[adj$cohort == n]) > 1e-6) log_msg("【核对】Figure 6E 与 06b 的 HR 不一致：", n)
    fr[[n]] <- o
  }
  fr <- do.call(rbind, fr)
  for (i in seq_len(nrow(fr))) add_stat(6, "E", paste("RLRS HR per SD,", fr$model[i], ":", fr$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", fr$HR[i], fr$lower[i], fr$upper[i]), exact_p(fr$p[i]))
  fr$lab <- COH_LAB[fr$cohort]
  fr$model <- factor(fr$model, levels = c("Unadjusted", "Adjusted for proliferation", "Adjusted for proliferation and stage"))
  pE <- forest_gg(fr, label = "lab", group = "model", xlab = "HR per SD of RLRS (95% CI)",
                  cols = c(Unadjusted = "grey55", `Adjusted for proliferation` = C_HI, `Adjusted for proliferation and stage` = "#641E16")) +
    ggplot2::labs(title = "Survival models") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")

  # F：与已发表的预后模型比较（各外部队列）
  pub <- utils::read.csv(file.path(d6, "Table_S5_published_models_Cindex.csv"))
  rk <- utils::read.csv(file.path(d6, "RLRS_rank_among_published.csv"))
  mw <- utils::read.csv(file.path(d6, "published_models_mean_external_Cindex.csv"), check.names = FALSE)
  ext <- c("GSE13507", "GSE32894", "GSE31684", "GSE48075", "IMvigor210", "GSE176307")
  pub <- pub[pub$cohort %in% ext, ]
  rk <- rk[match(ext, rk$cohort), ]
  ylab_rk <- stats::setNames(sprintf("%s (%d/%d)", COH_LAB[ext], rk$RLRS_rank, rk$n_models), ext)   # 括号内：RLRS 名次/模型数
  pub$lab <- factor(ylab_rk[pub$cohort], levels = rev(ylab_rk))
  add_stat(6, "F", "published signatures evaluated", length(unique(pub$model[pub$model != "RLRS (this study)"])),
           "in at least one external cohort")
  for (i in seq_len(nrow(rk))) add_stat(6, "F", paste("RLRS rank among models by C-index:", rk$cohort[i]),
                                        sprintf("%d of %d", rk$RLRS_rank[i], rk$n_models[i]),
                                        sprintf("RLRS %.3f; published median %.3f, max %.3f", rk$RLRS_Cindex[i], rk$published_median[i], rk$published_max[i]))
  mw <- mw[order(-mw$mean_external), ]
  add_stat(6, "F", "rank of RLRS by mean C-index across the six external cohorts", which(mw$model == "RLRS (this study)"),
           sprintf("of %d models computable in all six cohorts (RLRS %.3f; best %.3f)", nrow(mw),
                   mw$mean_external[mw$model == "RLRS (this study)"], max(mw$mean_external)))
  pF <- ggplot2::ggplot(pub[pub$model != "RLRS (this study)", ], ggplot2::aes(Cindex, lab)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, colour = "grey45") +
    ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.15, seed = 1), size = 0.25, colour = "grey70") +
    ggplot2::geom_point(data = pub[pub$model == "RLRS (this study)", ], colour = C_HI, size = 2.2, shape = 18) +
    ggplot2::labs(x = "C-index", y = NULL, title = "RLRS (red) vs published signatures")

  row1 <- (pA | pB | pC) + patchwork::plot_layout(widths = c(1.05, 0.9, 1))
  row2 <- (pD | pE | pF) + patchwork::plot_layout(widths = c(0.9, 1.05, 1.1))
  p <- (row1 / row2) + patchwork::plot_layout(heights = c(1, 1.1)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 6, 180, 165)
}

# =============================================================================
# Figure 7：临床意义（多因素 Cox、列线图、校准、时间 ROC、DCA）
# =============================================================================
fig7 <- function() {
  need(c("rms", "gridGraphics"))
  d7 <- file.path(DIR_RES, "07_clinical_nomogram")
  um <- utils::read.csv(file.path(d7, "TCGA_uni_multi_Cox.csv"))
  um$model <- factor(ifelse(um$model == "Univariate", "Univariable", "Multivariable"), levels = c("Univariable", "Multivariable"))
  um$variable <- sub("^Gender", "Sex", um$variable)
  um$variable <- sub("^Gender", "Sex", um$variable)
  for (i in seq_len(nrow(um))) add_stat(7, "A", paste(um$model[i], "Cox, TCGA OS:", um$variable[i]),
                                        sprintf("%.2f (%.2f-%.2f)", um$HR[i], um$lower[i], um$upper[i]), exact_p(um$p[i]))
  pA <- forest_gg(um, label = "variable", group = "model", xlab = "HR (95% CI)",
                  cols = c(Univariable = "grey55", Multivariable = C_HI)) +
    ggplot2::labs(title = "TCGA-BLCA, overall survival") + ggplot2::theme(legend.position = "bottom")

  ex <- utils::read.csv(file.path(d7, "External_multivariate_Cox.csv"))
  ex <- ex[ex$term == "RLRS_z", ]
  cv <- c(`age10+gender+MIBC` = "age, sex, MIBC", `age10+MIBC` = "age, MIBC")
  tc <- um[um$model == "Multivariable" & um$variable == "RLRS (per SD)", ]
  eb <- rbind(data.frame(lab = "TCGA-BLCA\n(age, sex, stage, grade)", HR = tc$HR, lower = tc$lower, upper = tc$upper, p = tc$p, LRT = NA),
              data.frame(lab = paste0(COH_LAB[ex$cohort], "\n(", cv[ex$covariates], ")"), HR = ex$HR, lower = ex$lower,
                         upper = ex$upper, p = ex$p, LRT = ex$LRT_p_add_RLRS))
  for (i in seq_len(nrow(ex))) add_stat(7, "B", paste("RLRS HR per SD adjusted for clinical factors:", ex$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", ex$HR[i], ex$lower[i], ex$upper[i]),
                                        sprintf("%s; LRT %s; C-index clinical %.3f vs clinical + RLRS %.3f", exact_p(ex$p[i]),
                                                exact_p(ex$LRT_p_add_RLRS[i]), ex$C_clinical[i], ex$C_clinical_RLRS[i]))
  eb$plab <- lab_p(eb$p)
  pB <- forest_gg(eb, label = "lab", xlab = "HR per SD of RLRS (95% CI)") +
    ggplot2::geom_text(ggplot2::aes(x = upper, label = plab), hjust = -0.15, size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.45))) +
    ggplot2::labs(title = "Multivariable models")

  nom <- readRDS(file.path(DIR_PRO, "nomogram_object.rds"))
  names(nom)[match(c("age", "stage_bin"), names(nom))] <- c("Age (years)", "Stage")
  pC <- patchwork::wrap_elements(full = ~ {
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5), family = FONT, cex = 0.56)
    plot(nom, xfrac = 0.28, lmgp = 0.18, cex.axis = 0.85, cex.var = 1, col.grid = grDevices::gray(c(0.85, 0.95)))
  })

  cal <- utils::read.csv(file.path(d7, "calibration_bootstrap.csv"))
  cal$yr <- factor(paste0(cal$year, "-year"), levels = c("1-year", "3-year", "5-year"))
  yr_cols <- c(`1-year` = "#E59866", `3-year` = C_HI, `5-year` = "#641E16")
  pD <- ggplot2::ggplot(cal, ggplot2::aes(mean.predicted, colour = yr)) +
    ggplot2::geom_abline(linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(ggplot2::aes(y = KM.corrected), linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(y = KM.corrected), size = 1) +
    ggplot2::scale_colour_manual(values = yr_cols, name = NULL) +
    ggplot2::coord_equal(xlim = c(0.2, 1), ylim = c(0.2, 1)) +
    ggplot2::labs(x = "Predicted survival (nomogram)", y = "Observed survival\n(bias-corrected)", title = "Calibration") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = ggplot2::element_blank())

  cn <- utils::read.csv(file.path(d7, "Cindex_nomogram_RLRS_stage.csv"))
  for (i in seq_len(nrow(cn))) add_stat(7, "E", paste("C-index, TCGA:", cn$model[i]), sprintf("%.3f", cn$Cindex[i]))
  ta <- utils::read.csv(file.path(d7, "timeROC_nomogram_RLRS_stage.csv"))
  tl <- stats::reshape(ta, direction = "long", varying = c("Nomogram", "RLRS", "Stage"), v.names = "AUC",
                       timevar = "model", times = c("Nomogram", "RLRS", "Stage"), idvar = "year")
  for (i in seq_len(nrow(tl))) add_stat(7, "E", sprintf("time-dependent AUC, %d year(s): %s", tl$year[i], tl$model[i]), sprintf("%.3f", tl$AUC[i]))
  m_cols <- c(Nomogram = "#641E16", RLRS = C_HI, Stage = "grey50")
  pE <- ggplot2::ggplot(tl, ggplot2::aes(factor(year), AUC, colour = model, group = model)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = m_cols, name = NULL) + ggplot2::coord_cartesian(ylim = c(0.45, 0.8)) +
    ggplot2::labs(x = "Years", y = "Time-dependent AUC", title = "TCGA-BLCA") +
    ggplot2::theme(legend.position = "bottom")

  dc <- utils::read.csv(file.path(d7, "DCA_net_benefit.csv"))
  dc <- dc[dc$year == 3 & dc$threshold <= 0.6, ]
  dc$label <- factor(dc$label, levels = c("Treat All", "Treat None", "Nomogram", "RLRS", "Stage"))
  pF <- ggplot2::ggplot(dc, ggplot2::aes(threshold, net_benefit, colour = label, linetype = label)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = c(`Treat All` = "grey65", `Treat None` = "black", m_cols), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(`Treat All` = 2, `Treat None` = 1, Nomogram = 1, RLRS = 1, Stage = 1), name = NULL) +
    ggplot2::coord_cartesian(ylim = c(-0.05, max(dc$net_benefit[dc$label != "Treat All"], na.rm = TRUE) * 1.15)) +
    ggplot2::labs(x = "Threshold probability (death within 3 years)", y = "Net benefit", title = "Decision curves, 3 years") +
    ggplot2::theme(legend.position = "bottom")

  p <- ((pA | pB) + patchwork::plot_layout(widths = c(1, 1))) /
    ((pC | pD) + patchwork::plot_layout(widths = c(1.45, 1))) /
    ((patchwork::free(pE) | pF) + patchwork::plot_layout(widths = c(1, 1.3))) +
    patchwork::plot_layout(heights = c(1, 1.15, 1)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 7, 180, 205)
}

# =============================================================================
# Figure 8：免疫微环境与治疗反应
# =============================================================================
TME_MAP <- data.frame(variable = c(
  "T_cells_CD8_CIBERSORT", "CD8_Tcells_EPIC", "CD8_T_cells_MCPcounter", "T_cells_CD8_quantiseq", "T_cell_CD8_TIMER", "CD8+_T-cells_xCell",
  "T_cells_CD4_memory_activated_CIBERSORT", "CD4_Tcells_EPIC", "T_cells_CD4_quantiseq", "T_cell_CD4_TIMER", "CD4+_T-cells_xCell",
  "T_cells_regulatory_(Tregs)_CIBERSORT", "Tregs_quantiseq", "Tregs_xCell",
  "NK_cells_activated_CIBERSORT", "NKcells_EPIC", "NK_cells_MCPcounter", "NK_cells_quantiseq", "NK_cells_xCell",
  "B_cells_naive_CIBERSORT", "Bcells_EPIC", "B_lineage_MCPcounter", "B_cells_quantiseq", "B_cell_TIMER", "B-cells_xCell",
  "Macrophages_M1_CIBERSORT", "Macrophages_M1_quantiseq", "Macrophages_M1_xCell",
  "Macrophages_M2_CIBERSORT", "Macrophages_M2_quantiseq", "Macrophages_M2_xCell",
  "Macrophages_EPIC", "Monocytic_lineage_MCPcounter", "Macrophage_TIMER", "Macrophages_xCell",
  "Dendritic_cells_activated_CIBERSORT", "Myeloid_dendritic_cells_MCPcounter", "Dendritic_cells_quantiseq", "DC_TIMER", "DC_xCell",
  "Neutrophils_CIBERSORT", "Neutrophils_MCPcounter", "Neutrophils_quantiseq", "Neutrophil_TIMER", "Neutrophils_xCell",
  "CAFs_EPIC", "Fibroblasts_MCPcounter", "Fibroblasts_xCell",
  "Endothelial_EPIC", "Endothelial_cells_MCPcounter", "Endothelial_cells_xCell",
  "ImmuneScore_estimate", "StromalScore_estimate", "TumorPurity_estimate"),
  cell = c(rep("CD8 T cell", 6), rep("CD4 T cell", 5), rep("Treg", 3), rep("NK cell", 5), rep("B cell", 6),
           rep("M1 macrophage", 3), rep("M2 macrophage", 3), rep("Macrophage/monocyte", 4), rep("Dendritic cell", 5),
           rep("Neutrophil", 5), rep("Fibroblast", 3), rep("Endothelial cell", 3), "Immune score", "Stromal score", "Tumour purity"),
  method = c("CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "EPIC", "quanTIseq", "TIMER", "xCell",
             "CIBERSORT", "quanTIseq", "xCell", "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "xCell",
             "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "quanTIseq", "xCell",
             "CIBERSORT", "quanTIseq", "xCell", "EPIC", "MCPcounter", "TIMER", "xCell",
             "CIBERSORT", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CIBERSORT", "MCPcounter", "quanTIseq", "TIMER", "xCell",
             "EPIC", "MCPcounter", "xCell", "EPIC", "MCPcounter", "xCell", "ESTIMATE", "ESTIMATE", "ESTIMATE"),
  check.names = FALSE)
TME_NONCELL <- c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT", "Other_quantiseq", "otherCells_EPIC")

fig8 <- function() {
  need("pROC")
  d8 <- file.path(DIR_RES, "08_tme_therapy")
  # A：免疫浸润（7 种算法；Spearman ρ 与 FDR）
  tme <- utils::read.csv(file.path(d8, "TME_High_vs_Low_all_methods.csv"), check.names = FALSE)
  cellv <- tme[!tme$variable %in% TME_NONCELL & !grepl("_estimate$", tme$variable) &
                 !tme$variable %in% c("ImmuneScore_xCell", "StromaScore_xCell", "MicroenvironmentScore_xCell"), ]
  add_stat(8, "A", "cell-type estimates (7 methods) with FDR < 0.05: higher / lower in RLRS-high",
           sprintf("%d / %d", sum(cellv$FDR < 0.05 & cellv$median_High > cellv$median_Low),
                   sum(cellv$FDR < 0.05 & cellv$median_High < cellv$median_Low)), paste("of", nrow(cellv), "estimates"))
  ta <- merge(TME_MAP, tme, by = "variable")
  if (nrow(ta) < nrow(TME_MAP)) log_msg("【核对】TME_MAP 有变量没匹配上：", paste(setdiff(TME_MAP$variable, ta$variable), collapse = ", "))
  for (i in seq_len(nrow(ta))) if (ta$method[i] == "ESTIMATE" || ta$FDR[i] < 0.05)
    add_stat(8, "A", paste("rho with RLRS:", ta$cell[i], "/", ta$method[i]), sprintf("%.2f", ta$rho_RLRS[i]), lab_q(ta$FDR[i]))
  ta$cell <- factor(ta$cell, levels = rev(unique(TME_MAP$cell)))
  ta$method <- factor(ta$method, levels = c("CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "ESTIMATE"))
  ta$sig <- ifelse(ta$FDR < 0.05, "FDR < 0.05", "NS")
  pA <- ggplot2::ggplot(ta, ggplot2::aes(method, cell)) +
    ggplot2::geom_point(ggplot2::aes(size = pmin(-log10(FDR), 10), fill = rho_RLRS, colour = sig), shape = 21, stroke = 0.35) +
    ggplot2::scale_fill_gradient2(low = C_LO, mid = "white", high = C_HI, limits = c(-0.7, 0.7), oob = scales::squish,
                                  name = "Spearman ρ\nwith RLRS") +
    ggplot2::scale_size_area(max_size = 3.2, name = "-log10 FDR", breaks = c(2, 5, 10), labels = c("2", "5", "≥10")) +
    ggplot2::scale_colour_manual(values = c(`FDR < 0.05` = "black", NS = "grey80"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Immune and stromal cells, TCGA-BLCA") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.box = "vertical")

  # B：免疫检查点
  ck <- utils::read.csv(file.path(d8, "Checkpoints_High_vs_Low.csv"))
  for (i in seq_len(nrow(ck))) add_stat(8, "B", paste("checkpoint gene rho with RLRS:", ck$variable[i]),
                                        sprintf("%.2f", ck$rho_RLRS[i]), lab_q(ck$FDR[i]))
  ck <- ck[order(ck$rho_RLRS), ]; ck$variable <- factor(ck$variable, levels = ck$variable)
  ck$sig <- ifelse(ck$FDR < 0.05, ifelse(ck$rho_RLRS > 0, "Positive", "Negative"), "NS")
  pB <- ggplot2::ggplot(ck, ggplot2::aes(rho_RLRS, variable, fill = sig)) + ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(Positive = C_HI, Negative = C_LO, NS = "grey75"), guide = "none") +
    ggplot2::labs(x = "Spearman ρ with RLRS", y = NULL, title = "Immune checkpoints") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))

  # C：TIDE（tidepy）
  r <- rlrs_all(); rt <- r[r$cohort == "TCGA_BLCA", ]
  td <- utils::read.csv(file.path(d8, "TIDE_output.csv"), check.names = FALSE)
  td <- merge(rt, td, by.x = "ID", by.y = "Patient")
  tl <- rbind(data.frame(score = "TIDE", v = td$TIDE, g = td$RLRS_group), data.frame(score = "Dysfunction", v = td$Dysfunction, g = td$RLRS_group),
              data.frame(score = "Exclusion", v = td$Exclusion, g = td$RLRS_group))
  tl$score <- factor(tl$score, levels = c("TIDE", "Dysfunction", "Exclusion")); tl$g <- factor(tl$g, levels = c("Low", "High"))
  tp <- do.call(rbind, lapply(levels(tl$score), function(s) data.frame(score = s, p = stats::wilcox.test(v ~ g, data = tl[tl$score == s, ])$p.value)))
  for (i in seq_len(nrow(tp))) add_stat(8, "C", paste("TIDE component, RLRS high vs low, Wilcoxon:", tp$score[i]), exact_p(tp$p[i]),
                                        sprintf("median high %.2f vs low %.2f; n = %d", stats::median(tl$v[tl$score == tp$score[i] & tl$g == "High"]),
                                                stats::median(tl$v[tl$score == tp$score[i] & tl$g == "Low"]), nrow(td)))
  tp$score <- factor(tp$score, levels = levels(tl$score))
  pC <- ggplot2::ggplot(tl, ggplot2::aes(g, v, fill = g)) +
    ggplot2::geom_boxplot(outlier.size = 0.2, width = 0.6, linewidth = 0.3) +
    ggplot2::geom_text(data = tp, ggplot2::aes(x = 1.5, y = Inf, label = lab_p(p)), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ score, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(x = NULL, y = "Score", title = "TIDE, TCGA-BLCA")

  # D/E：免疫治疗队列（客观缓解；AUC 方向固定为"RLRS 越高越可能缓解"）
  ici <- utils::read.csv(file.path(d8, "ICI_response_summary.csv"))
  for (i in seq_len(nrow(ici))) add_stat(8, "D-E", paste("ICI response:", ici$cohort[i]),
                                         sprintf("ORR high %.1f%% vs low %.1f%%, Fisher %s", 100 * ici$ORR_High[i], 100 * ici$ORR_Low[i], exact_p(ici$p_fisher[i])),
                                         sprintf("AUC %.3f (%.3f-%.3f), Wilcoxon %s; %d evaluable, %d responders", ici$AUC[i], ici$AUC_lower[i],
                                                 ici$AUC_upper[i], exact_p(ici$p_wilcox[i]), ici$n_evaluable[i], ici$n_responders[i]))
  orr <- rbind(data.frame(cohort = ici$cohort, g = "High", orr = ici$ORR_High), data.frame(cohort = ici$cohort, g = "Low", orr = ici$ORR_Low))
  orr$g <- factor(orr$g, levels = c("Low", "High")); orr$cohort <- factor(orr$cohort, levels = ici$cohort)
  ici$cohort <- factor(ici$cohort, levels = ici$cohort)
  pD <- ggplot2::ggplot(orr, ggplot2::aes(g, orr, fill = g)) + ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * orr)), vjust = -0.4, size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::geom_text(data = ici, ggplot2::aes(x = 1.5, y = Inf, label = paste("Fisher", lab_p(p_fisher))), vjust = 1.3,
                       inherit.aes = FALSE, size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ cohort, nrow = 1) +
    ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 0.45), expand = c(0, 0)) +
    ggplot2::labs(x = "RLRS group", y = "Objective response rate", title = "Anti-PD-L1/PD-1 therapy")
  roc_df <- do.call(rbind, lapply(c("IMvigor210", "GSE176307"), function(n) {
    e <- rlrs_cohort(n); e <- e[e$response %in% c("CR/PR", "SD/PD"), ]
    ro <- pROC::roc(factor(e$response, levels = c("SD/PD", "CR/PR")), e$RLRS, levels = c("SD/PD", "CR/PR"), direction = "<", quiet = TRUE)
    data.frame(cohort = n, fpr = 1 - ro$specificities, tpr = ro$sensitivities, auc = as.numeric(ro$auc))
  }))
  roc_df <- roc_df[order(roc_df$cohort, roc_df$fpr, roc_df$tpr), ]
  al <- stats::aggregate(auc ~ cohort, data = roc_df, FUN = function(x) x[1])
  al$lab <- sprintf("%s  AUC %.2f", al$cohort, al$auc)
  roc_df$cohort <- factor(roc_df$cohort, levels = c("IMvigor210", "GSE176307"), labels = al$lab[match(c("IMvigor210", "GSE176307"), al$cohort)])
  pE <- ggplot2::ggplot(roc_df, ggplot2::aes(fpr, tpr, colour = cohort)) +
    ggplot2::geom_abline(linetype = 2, linewidth = 0.3, colour = "grey50") + ggplot2::geom_path(linewidth = 0.45) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#8E44AD", "#C39BD3"), levels(roc_df$cohort)), name = NULL) +
    ggplot2::coord_equal() + ggplot2::labs(x = "1 - specificity", y = "Sensitivity", title = "Response (CR/PR)") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0), legend.background = ggplot2::element_blank())

  # F：IMvigor210 OS
  e <- rlrs_cohort("IMvigor210"); e <- e[!is.na(e$OS.time) & e$OS.time > 0 & !is.na(e$OS), ]
  tunit <- if (max(e$OS.time, na.rm = TRUE) > 100) TRUE else FALSE     # IMvigor210 生存时间单位：天 → 年；若为月，另行换算
  pF <- km_gg(if (tunit) e$OS.time else e$OS.time * 30.44, e$OS, e$RLRS_group, title = "IMvigor210", ylab = "Overall survival")
  s <- attr(pF, "stats")
  add_stat(8, "F", "IMvigor210 OS, RLRS high vs low", sprintf("HR %.2f (%.2f-%.2f)", s["HR"], s["lower"], s["upper"]),
           paste0("log-rank ", exact_p(s["p_logrank"]), "; n = ", nrow(e)))

  # G：GSE169455 新辅助化疗 pCR（7/16 个模型基因在该数据中缺失）
  ch <- utils::read.csv(file.path(d8, "Chemo_response_GSE169455.csv"))
  for (i in seq_len(nrow(ch))) add_stat(8, "G", paste("GSE169455", ch$set[i], ch$endpoint[i]),
                                        sprintf("AUC %.3f (%.3f-%.3f); Wilcoxon %s", ch$AUC[i], ch$AUC_lower[i], ch$AUC_upper[i], exact_p(ch$p_wilcox[i])),
                                        sprintf("OR per SD %.2f (%s); adjusted for consensus class %.2f (%s); n = %d, events = %d",
                                                ch$OR_perSD[i], exact_p(ch$p_logistic[i]), ch$OR_perSD_adj_subtype[i], exact_p(ch$p_adj_subtype[i]), ch$n[i], ch$n_yes[i]))
  e <- rlrs_cohort("GSE169455"); e <- e[e$setting == "neoadjuvant", ]
  gl <- rbind(data.frame(ep = "pCR", g = e$RLRS_group, y = e$pCR), data.frame(ep = "Downstaging", g = e$RLRS_group, y = e$downstage))
  gl <- gl[!is.na(gl$y), ]
  ga <- stats::aggregate(y ~ ep + g, data = gl, FUN = mean)
  ga$g <- factor(ga$g, levels = c("Low", "High")); ga$ep <- factor(ga$ep, levels = c("pCR", "Downstaging"))
  gp <- data.frame(ep = factor(c("pCR", "Downstaging"), levels = c("pCR", "Downstaging")),
                   p = c(ch$p_wilcox[ch$set == "neoadjuvant" & ch$endpoint == "pCR"], ch$p_wilcox[ch$set == "neoadjuvant" & ch$endpoint == "downstage"]),
                   auc = c(ch$AUC[ch$set == "neoadjuvant" & ch$endpoint == "pCR"], ch$AUC[ch$set == "neoadjuvant" & ch$endpoint == "downstage"]))
  pG <- ggplot2::ggplot(ga, ggplot2::aes(g, y, fill = g)) + ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * y)), vjust = -0.4, size = (BASE - 1) / ggplot2::.pt, family = FONT) +
    ggplot2::geom_text(data = gp, ggplot2::aes(x = 1.5, y = Inf, label = sprintf("AUC %.2f\n%s", auc, lab_p(p))), vjust = 1.2,
                       inherit.aes = FALSE, size = (BASE - 1.5) / ggplot2::.pt, family = FONT, lineheight = 0.9) +
    ggplot2::facet_wrap(~ ep, nrow = 1) + ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 0.85), expand = c(0, 0)) +
    ggplot2::labs(x = "RLRS group", y = "Patients with the outcome", title = "Neoadjuvant chemotherapy (GSE169455)")

  # H：oncoPredict 预测的 IC50（如有）
  f_op <- file.path(d8, "oncoPredict_High_vs_Low.csv")
  pH <- NULL
  if (file.exists(f_op)) {
    op <- utils::read.csv(f_op, check.names = FALSE)
    pred <- utils::read.csv(file.path(d8, "oncoPredict", "calcPhenotype_Output", "DrugPredictions.csv"), row.names = 1, check.names = FALSE)
    add_stat(8, "H", "drugs with predicted IC50 differing between RLRS groups (FDR < 0.05): lower / higher in RLRS-high",
             sprintf("%d / %d", sum(op$FDR < 0.05 & op$median_High < op$median_Low), sum(op$FDR < 0.05 & op$median_High > op$median_Low)),
             paste("of", nrow(op), "drugs"))
    focus <- grep("^(Cisplatin|Gemcitabine|Docetaxel|Olaparib|AZD6738|Erlotinib)_", colnames(pred), value = TRUE)
    focus <- focus[!duplicated(sub("_.*", "", focus))]
    for (v in focus) add_stat(8, "H", paste("predicted IC50, RLRS high vs low:", v),
                              sprintf("median %.3g vs %.3g", op$median_High[op$variable == v], op$median_Low[op$variable == v]),
                              lab_q(op$FDR[op$variable == v]))
    dl <- do.call(rbind, lapply(focus, function(v) data.frame(drug = sub("_.*", "", v), id = rownames(pred), v = log10(pred[[v]]))))
    dl$g <- factor(rt$RLRS_group[match(dl$id, rt$ID)], levels = c("Low", "High"))
    dl <- dl[!is.na(dl$g), ]
    qq <- data.frame(drug = sub("_.*", "", focus), q = op$FDR[match(focus, op$variable)])
    pH <- ggplot2::ggplot(dl, ggplot2::aes(g, v, fill = g)) + ggplot2::geom_boxplot(outlier.size = 0.2, width = 0.6, linewidth = 0.3) +
      ggplot2::geom_text(data = qq, ggplot2::aes(x = 1.5, y = Inf, label = lab_q(q)), vjust = 1.3, inherit.aes = FALSE,
                         size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
      ggplot2::facet_wrap(~ drug, ncol = 3, scales = "free_y") + ggplot2::scale_fill_manual(values = PAL2, guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
      ggplot2::labs(x = "RLRS group", y = "Predicted log10 IC50", title = "Predicted drug sensitivity (GDSC2, oncoPredict)")
  }
  # 手动加标签、每排单独排版：避免各排坐标轴互相对齐造成的大块空白
  pl <- Filter(Negate(is.null), list(pA, pB, pC, pD, pE, pF, pG, pH))
  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), pl, LETTERS[seq_along(pl)])
  top <- patchwork::wrap_elements(full = (pl[[1]] | (pl[[2]] / pl[[3]] + patchwork::plot_layout(heights = c(1, 0.8)))) +
                                    patchwork::plot_layout(widths = c(1.15, 1)))
  mid <- patchwork::wrap_elements(full = (pl[[4]] | pl[[5]] | pl[[6]]) + patchwork::plot_layout(widths = c(1.1, 0.9, 1)))
  bot <- patchwork::wrap_elements(full = if (length(pl) == 8) (pl[[7]] | pl[[8]]) + patchwork::plot_layout(widths = c(0.75, 1.25)) else pl[[7]])
  p <- patchwork::wrap_plots(top, mid, bot, ncol = 1, heights = c(1.3, 0.85, 1.0))
  save_fig(p, 8, 180, 245)
}

# =============================================================================
# 补充图 S1：单细胞聚类、数据集、各簇 marker、各样本细胞组成与质控指标
# =============================================================================
figS1 <- function() {
  s <- load_sc(); md0 <- s$md
  lev <- as.character(sort(unique(as.integer(as.character(md0$seurat_clusters)))))
  md0$cluster <- factor(as.character(md0$seurat_clusters), levels = lev)
  set.seed(PARAM$seed); md <- md0[sample(nrow(md0)), ]
  cen <- stats::aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = md, FUN = stats::median)
  add_stat("S1", "A", "clusters", nlevels(md$cluster), "Harmony integration, resolution 0.8")
  noax <- ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  pA <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = cluster)) + ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggrepel::geom_text_repel(data = cen, ggplot2::aes(label = cluster), colour = "black", size = BASE / ggplot2::.pt, family = FONT,
                             seed = 1, box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.2,
                             bg.color = "white", bg.r = 0.12) +
    ggplot2::scale_colour_manual(values = grDevices::hcl.colors(nlevels(md$cluster), "Dark 3"), guide = "none") +
    ggplot2::coord_fixed() + ggplot2::labs(title = sprintf("%d clusters", nlevels(md$cluster))) + noax
  pB <- ggplot2::ggplot(md, ggplot2::aes(UMAP1, UMAP2, colour = dataset)) + ggplot2::geom_point(size = 0.25, stroke = 0) +
    ggplot2::scale_colour_manual(values = c(GSE135337 = "#E67E22", GSE222315 = "#2E6DA4"), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) + ggplot2::coord_fixed() +
    ggplot2::labs(title = "Dataset") + noax +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0))

  # C：各簇的 marker 点图（y 轴：簇号与注释）
  ex <- s$mk; mk <- unlist(MK_SHOW); mk <- mk[mk %in% rownames(ex)]
  ann <- tapply(as.character(md0$ct), md0$cluster, function(v) names(sort(table(v), decreasing = TRUE))[1])
  dp <- do.call(rbind, lapply(lev, function(k) {
    m <- ex[mk, md0$cluster == k, drop = FALSE]
    data.frame(cluster = k, gene = mk, avg = Matrix::rowMeans(m), pct = Matrix::rowMeans(m > 0))
  }))
  dp$avg_s <- stats::ave(dp$avg, dp$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  ylab <- paste0(lev, "  ", ann[lev])
  dp$y <- factor(paste0(dp$cluster, "  ", ann[dp$cluster]), levels = rev(ylab))
  dp$gene <- factor(dp$gene, levels = unique(mk))
  for (k in lev) add_stat("S1", "C", paste("cluster annotation:", k), ann[[k]], paste(sum(md0$cluster == k), "cells"))
  pC <- ggplot2::ggplot(dp, ggplot2::aes(gene, y)) + ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 2.4, name = "Fraction\nexpressing", breaks = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Marker genes by cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   axis.text.y = ggplot2::element_text(size = BASE - 1.5), legend.box = "horizontal")

  # D：各样本的细胞类型组成
  cp <- as.data.frame(table(sample = md0$sample, ct = md0$ct))
  tot <- tapply(cp$Freq, cp$sample, sum); cp$frac <- cp$Freq / tot[as.character(cp$sample)]
  tis <- tapply(as.character(md0$tissue), md0$sample, function(v) v[1])
  cp$tissue <- factor(ifelse(tis[as.character(cp$sample)] == "Tumor", "Tumour", "Normal"),
                      levels = c("Tumour", "Normal"))
  for (k in names(tot)) add_stat("S1", "D", paste("cells after QC:", k), tot[[k]],
                                 sprintf("epithelial %.1f%%", 100 * mean(md0$ct[md0$sample == k] == "Epithelial")))
  pD <- ggplot2::ggplot(cp, ggplot2::aes(sample, frac, fill = ct)) + ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_grid(~ tissue, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = CT_COLS, name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Cells", title = "Cell-type composition") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5), strip.clip = "off",
                   legend.key.size = ggplot2::unit(2.4, "mm"), legend.text = ggplot2::element_text(size = BASE - 1.5))

  # E：质控指标（过滤后）
  mtc <- grep("percent.?mt|pct.?mt|mito", names(md0), value = TRUE, ignore.case = TRUE)[1]
  qv <- c(nFeature_RNA = "Detected genes", stats::setNames("Mitochondrial reads (%)", mtc))
  qv <- qv[names(qv) %in% names(md0)]
  ql <- do.call(rbind, lapply(names(qv), function(v) data.frame(sample = md0$sample, metric = qv[[v]], value = md0[[v]])))
  ql$sample <- factor(ql$sample, levels = levels(factor(cp$sample)))
  pE <- ggplot2::ggplot(ql, ggplot2::aes(sample, value)) +
    ggplot2::geom_violin(scale = "width", fill = "grey80", colour = "grey40", linewidth = 0.2) +
    ggplot2::facet_wrap(~ metric, ncol = 1, scales = "free_y") + ggplot2::labs(x = NULL, y = NULL, title = "Quality control") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 1.5))

  pl <- Map(function(g, t) g + ggplot2::labs(tag = t), list(pA, pB, pC, pD, pE), LETTERS[1:5])
  p <- patchwork::wrap_plots(patchwork::wrap_elements(full = pl[[1]] | pl[[2]]), patchwork::wrap_elements(full = pl[[3]]),
                             patchwork::wrap_elements(full = (pl[[4]] | pl[[5]]) + patchwork::plot_layout(widths = c(1.25, 1))),
                             ncol = 1, heights = c(1, 1.1, 0.95))
  save_fig(p, "S1", 180, 235)
}

# =============================================================================
# 补充图 S3：bulk 中 RL_core 与 RHI 按分子亚型、分期、分级和突变的分布（TCGA-BLCA）
# =============================================================================
figS3 <- function() {
  tum <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
  t1 <- tum[tum$group == "Tumor", ]
  mt <- utils::read.csv(file.path(DIR_RES, "04_bulk_score", "Rloop_score_by_mutation.csv"))
  sc2 <- c(RL_core = "RL_core", RHI = "RHI")
  box_by <- function(d, var, lev, test = c("kw", "wilcox"), title) {
    test <- match.arg(test)
    d <- d[!is.na(d[[var]]) & d[[var]] %in% lev, ]
    l <- do.call(rbind, lapply(names(sc2), function(s) data.frame(score = sc2[[s]], grp = factor(d[[var]], levels = lev), v = d[[s]])))
    l$score <- factor(l$score, levels = sc2)
    pv <- vapply(sc2, function(s) {
      x <- l[l$score == s, ]
      if (test == "kw") stats::kruskal.test(v ~ grp, data = x)$p.value else stats::wilcox.test(v ~ grp, data = x)$p.value
    }, numeric(1))
    for (s in sc2) add_stat("S3", title, paste(s, "by", var, if (test == "kw") "Kruskal-Wallis P" else "Wilcoxon P"),
                            signif(pv[[s]], 3), paste(sprintf("%s n = %d", lev, as.integer(table(d[[var]])[lev])), collapse = "; "))
    ggplot2::ggplot(l, ggplot2::aes(grp, v)) +
      ggplot2::geom_boxplot(ggplot2::aes(colour = grp), outlier.shape = NA, width = 0.6, linewidth = 0.35) +
      ggplot2::geom_jitter(ggplot2::aes(colour = grp), width = 0.15, size = 0.25, alpha = 0.6) +
      ggplot2::geom_text(data = data.frame(score = factor(sc2, levels = sc2), lab = lab_p(pv)),
                         ggplot2::aes(x = (length(lev) + 1) / 2, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                         size = (BASE - 1) / ggplot2::.pt, family = FONT) +
      ggplot2::facet_wrap(~ score, scales = "free_y", nrow = 1) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
      ggplot2::labs(x = NULL, y = "ssGSEA score") + ggplot2::theme(legend.position = "none")
  }
  cs_lev <- c("LumP", "LumNS", "LumU", "Stroma-rich", "Ba/Sq", "NE-like")
  cs_col <- c(LumP = "#2E6DA4", LumNS = "#5DADE2", LumU = "#8E44AD", `Stroma-rich` = "#27AE60", `Ba/Sq` = C_HI, `NE-like` = "#E67E22")
  pA <- box_by(t1, "consensus", cs_lev, "kw", "A") + ggplot2::scale_colour_manual(values = cs_col) +
    ggplot2::labs(title = "Consensus molecular classes") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  pB <- box_by(t1, "stage_bin", c("I-II", "III-IV"), "wilcox", "B") +
    ggplot2::scale_colour_manual(values = c(`I-II` = C_LO, `III-IV` = C_HI)) + ggplot2::labs(title = "AJCC stage")
  pC <- box_by(t1, "grade", c("Low", "High"), "wilcox", "C") +
    ggplot2::scale_colour_manual(values = c(Low = C_LO, High = C_HI)) + ggplot2::labs(title = "Histological grade")
  # D：RHI 与常见突变（FDR 来自 04 步，12 个基因校正）
  gm <- c("FGFR3", "KDM6A", "ELF3", "TP53", "RB1")
  ml <- do.call(rbind, lapply(gm, function(g) data.frame(gene = g, status = t1[[paste0("mut_", g)]], RHI = t1$RHI)))
  ml <- ml[!is.na(ml$status), ]; ml$status <- factor(ml$status, levels = c("WT", "Mut")); ml$gene <- factor(ml$gene, levels = gm)
  for (g in gm) add_stat("S3", "D", paste("RHI mutant vs wild-type, Wilcoxon FDR:", g), signif(mt$FDR_RHI[mt$gene == g], 3),
                         sprintf("%d mutant / %d tested", sum(ml$gene == g & ml$status == "Mut"), sum(ml$gene == g)))
  pD <- ggplot2::ggplot(ml, ggplot2::aes(status, RHI)) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = status), outlier.shape = NA, width = 0.6, linewidth = 0.35) +
    ggplot2::geom_jitter(ggplot2::aes(colour = status), width = 0.15, size = 0.25, alpha = 0.6) +
    ggplot2::geom_text(data = data.frame(gene = factor(gm, levels = gm), lab = lab_q(mt$FDR_RHI[match(gm, mt$gene)])),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, nrow = 1) + ggplot2::scale_colour_manual(values = c(WT = "grey50", Mut = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(x = NULL, y = "RHI (ssGSEA)", title = "RHI by somatic mutation")
  p <- (pA / (pB | pC) / pD) + patchwork::plot_layout(heights = c(1.1, 1, 0.9)) + patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, "S3", 180, 190)
}

# =============================================================================
# Figure 9：模型基因（SHAP）→ bulk → 单细胞 → 空间转录组
# =============================================================================
fig9 <- function() {
  d9 <- file.path(DIR_RES, "09_hub_gene_spatial")
  rk <- utils::read.csv(file.path(d9, "SHAP_gene_ranking.csv"))
  cf <- utils::read.csv(file.path(d9, "RLRS_coefficients.csv"))
  genes <- rk$gene; hubs <- utils::head(genes, 3)
  HUB_COLS <- stats::setNames(c(C_HI, "#E67E22", "#2E6DA4"), hubs)

  # A：精确 SHAP（线性模型：w ×（z − 均值）），TCGA-BLCA
  tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
  mi <- build_mime_input(tcga$expr, tcga$clin[tcga$clin$group == "Tumor", ], genes)
  r0 <- rlrs_all(); mi <- merge(r0[r0$cohort == "TCGA_BLCA", c("ID", "RLRS")], mi, by = "ID")
  X <- as.matrix(mi[, safe_names(genes)]); colnames(X) <- genes
  w <- stats::setNames(cf$coefficient[match(genes, cf$term)], genes)
  S <- sweep(sweep(X, 2, colMeans(X)), 2, w, `*`)
  chk <- max(abs(colMeans(abs(S)) - rk$mean_abs_SHAP))
  if (chk > 1e-6) log_msg("【核对】Figure 9A 的 SHAP 与 09 步不一致（最大差 ", signif(chk, 3), "）")
  fit_chk <- max(abs(drop(cf$coefficient[1] + X %*% w) - mi$RLRS))
  if (fit_chk > 1e-3) log_msg("【核对】由系数还原的 RLRS 与 06 步不一致（最大差 ", signif(fit_chk, 3), "）")
  for (i in seq_len(nrow(rk))) add_stat(9, "A", paste("mean |SHAP|, share of total:", rk$gene[i]),
                                        sprintf("%.2f (%.1f%%)", rk$mean_abs_SHAP[i], rk$share_pct[i]),
                                        sprintf("coefficient %.2f", rk$coefficient[i]))
  ylab <- stats::setNames(sprintf("%s (%.1f%%)", rk$gene, rk$share_pct), rk$gene)
  sl <- data.frame(gene = factor(rep(genes, each = nrow(S)), levels = rev(genes)), shap = as.vector(S),
                   z = pmax(pmin(as.vector(X), 2.5), -2.5))
  pA <- ggplot2::ggplot(sl, ggplot2::aes(shap, gene, colour = z)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.28, width = 0, seed = 1), size = 0.2, stroke = 0) +
    ggplot2::scale_colour_gradient2(low = C_LO, mid = "grey85", high = C_HI, name = "Expression\n(z-score)", limits = c(-2.5, 2.5)) +
    ggplot2::scale_y_discrete(labels = ylab) +
    ggplot2::labs(x = "SHAP value (contribution to RLRS)", y = NULL, title = "Gene contributions, TCGA-BLCA") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom",
                   legend.key.width = ggplot2::unit(5, "mm"), legend.key.height = ggplot2::unit(2, "mm"))

  # B：肿瘤 vs 正常（19 对配对样本）
  cl <- tcga$clin
  pat <- intersect(cl$patient[cl$group == "Tumor"], cl$patient[cl$group == "Normal"])
  pr <- cl[cl$patient %in% pat, ]; pr <- pr[!duplicated(pr[, c("patient", "group")]), ]
  bl <- do.call(rbind, lapply(hubs, function(g) data.frame(gene = g, patient = pr$patient, group = pr$group,
                                                          v = as.numeric(tcga$expr[g, pr$ID]))))
  bl$group <- factor(bl$group, levels = c("Normal", "Tumor"), labels = c("Normal", "Tumour")); bl$gene <- factor(bl$gene, levels = hubs)
  pvB <- vapply(hubs, function(g) {
    x <- stats::reshape(bl[bl$gene == g, c("patient", "group", "v")], idvar = "patient", timevar = "group", direction = "wide")
    stats::wilcox.test(x$v.Tumour, x$v.Normal, paired = TRUE)$p.value
  }, numeric(1))
  for (g in hubs) {
    xa <- as.numeric(tcga$expr[g, cl$ID[cl$group == "Tumor"]]); xb <- as.numeric(tcga$expr[g, cl$ID[cl$group == "Normal"]])
    add_stat(9, "B", paste("tumour vs normal, paired Wilcoxon P (all-sample Wilcoxon P):", g), signif(pvB[[g]], 3),
             sprintf("%d pairs; all samples P = %s; median tumour %.2f vs normal %.2f", length(pat),
                     signif(stats::wilcox.test(xa, xb)$p.value, 3), stats::median(xa), stats::median(xb)))
  }
  pB <- ggplot2::ggplot(bl, ggplot2::aes(group, v)) +
    ggplot2::geom_line(ggplot2::aes(group = patient), colour = "grey75", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 0.7) +
    ggplot2::geom_text(data = data.frame(gene = factor(hubs, levels = hubs), lab = lab_p(pvB)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = c(Normal = C_LO, Tumour = C_HI), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggplot2::labs(x = NULL, y = "log2(TPM + 1)", title = sprintf("TCGA-BLCA paired samples (n = %d)", length(pat))) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  # C：各队列中单个基因的预后价值（每 SD 的 HR）
  hr <- list()
  for (n in SURV_COH) {
    o <- if (n == "TCGA_BLCA") tcga else readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
    c0 <- coh_clin(n); c0 <- c0[!is.na(c0$OS.time) & c0$OS.time > 0 & !is.na(c0$OS) & c0$ID %in% colnames(o$expr), ]
    for (g in intersect(hubs, rownames(o$expr))) {
      z <- as.numeric(scale(as.numeric(o$expr[g, c0$ID])))
      cs <- summary(survival::coxph(survival::Surv(c0$OS.time, c0$OS) ~ z))
      hr[[paste(n, g)]] <- data.frame(cohort = n, gene = g, HR = cs$conf.int[1, 1], lower = cs$conf.int[1, 3],
                                      upper = cs$conf.int[1, 4], p = cs$coefficients[1, 5], n = nrow(c0), events = sum(c0$OS))
    }
  }
  hr <- do.call(rbind, hr)
  for (i in seq_len(nrow(hr))) add_stat(9, "C", paste("HR per SD:", hr$gene[i], hr$cohort[i]),
                                        sprintf("%.2f (%.2f-%.2f)", hr$HR[i], hr$lower[i], hr$upper[i]), exact_p(hr$p[i]))
  hr$lab <- factor(COH_LAB[hr$cohort], levels = COH_LAB[SURV_COH]); hr$gene <- factor(hr$gene, levels = hubs)
  pC <- forest_gg(hr, label = "lab", group = "gene", cols = HUB_COLS) +
    ggplot2::labs(title = "Overall survival (DSS in GSE32894)") +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(face = "italic"))

  # D：单细胞，各细胞类型中的表达（03 步的细胞分类）
  cc <- utils::read.csv(file.path(d9, "hub_scRNA_by_cell_class.csv"))
  CLS <- c(Malignant = "Tumour epithelial", Normal_epi = "Normal epithelial", T_cell = "T cell", NK = "Cytotoxic T/NK",
           B_cell = "B cell", Plasma = "Plasma cell", Myeloid = "Myeloid", DC = "Dendritic cell", Mast = "Mast cell",
           Fibroblast = "Fibroblast", Myofib_SMC = "Myofibroblast/SMC", Endothelial = "Endothelial")
  cc$cls <- factor(ifelse(cc$cell_class %in% names(CLS), CLS[cc$cell_class], cc$cell_class), levels = rev(CLS[CLS %in% CLS[cc$cell_class]]))
  cc$gene <- factor(cc$gene, levels = hubs)
  cc$avg_s <- stats::ave(cc$avg, cc$gene, FUN = function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
  for (i in seq_len(nrow(cc))) add_stat(9, "D", paste("fraction expressing / mean expression:", cc$gene[i], cc$cell_class[i]),
                                        sprintf("%.3f / %.3f", cc$pct[i], cc$avg[i]))
  pD <- ggplot2::ggplot(cc, ggplot2::aes(gene, cls)) + ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_s)) +
    ggplot2::scale_size_area(max_size = 3.2, name = "Fraction\nexpressing", breaks = c(0.05, 0.25, 0.5, 0.75)) +
    ggplot2::scale_colour_gradient(low = "grey90", high = C_HI, name = "Scaled\nexpression", breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Single-cell expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"), legend.box = "horizontal",
                   legend.position = "right")

  # E：R-loop 高/低肿瘤上皮细胞（样本内平均，配对）
  hv <- utils::read.csv(file.path(d9, "hub_scRNA_malignant_High_vs_Low_per_sample.csv"), check.names = FALSE)
  el <- do.call(rbind, lapply(hubs, function(g) data.frame(gene = g, sample = hv$sample, grp = hv$group, v = hv[[g]])))
  el$grp <- factor(el$grp, levels = c("Low", "High")); el$gene <- factor(el$gene, levels = hubs)
  pvE <- vapply(hubs, function(g) {
    x <- stats::reshape(el[el$gene == g, c("sample", "grp", "v")], idvar = "sample", timevar = "grp", direction = "wide")
    suppressWarnings(stats::wilcox.test(x$v.High, x$v.Low, paired = TRUE, exact = FALSE)$p.value)
  }, numeric(1))
  for (g in hubs) {
    x <- stats::reshape(el[el$gene == g, c("sample", "grp", "v")], idvar = "sample", timevar = "grp", direction = "wide")
    add_stat(9, "E", paste("R-loop-high vs -low tumour epithelial cells, paired Wilcoxon P:", g), signif(pvE[[g]], 3),
             sprintf("higher in R-loop-high in %d of %d samples", sum(x$v.High > x$v.Low), nrow(x)))
  }
  pE <- ggplot2::ggplot(el, ggplot2::aes(grp, v)) +
    ggplot2::geom_line(ggplot2::aes(group = sample), colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.9) +
    ggplot2::geom_text(data = data.frame(gene = factor(hubs, levels = hubs), lab = lab_p(pvE)),
                       ggplot2::aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, inherit.aes = FALSE,
                       size = (BASE - 1.5) / ggplot2::.pt, family = FONT) +
    ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 1) + ggplot2::scale_colour_manual(values = PAL2, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggplot2::labs(x = "R-loop group", y = "Mean expression per sample", title = "Tumour epithelial cells, samples as the unit") +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"))

  # F：空间分布（两张示例切片；各特征在切片内按百分位着色）
  slim <- readRDS(file.path(DIR_PRO, "spatial_spots_slim.rds"))
  show <- SP_SHOW[SP_SHOW %in% slim$sample][1]
  feats <- c(RL_core = "RL_core", Prolif = "Proliferation", UMI = "UMI count", stats::setNames(hubs[1], hubs[1]))
  fl <- do.call(rbind, lapply(show, function(smp) {
    d <- slim[slim$sample == smp, ]
    do.call(rbind, lapply(names(feats), function(f) data.frame(sample = smp, feature = feats[[f]], x = d$x, y = d$y,
                                                               v = rank(d[[f]], ties.method = "average") / nrow(d))))
  }))
  fl$feature <- factor(fl$feature, levels = feats); fl$sample <- factor(fl$sample, levels = show)
  sp_rho <- function(a, b) stats::cor(slim[[a]][slim$sample == show], slim[[b]][slim$sample == show], method = "spearman")
  add_stat(9, "F", paste("example section", show, "spot-level rho RL_core vs proliferation / UMI / ANXA1"),
           sprintf("%.2f / %.2f / %.2f", sp_rho("RL_core", "Prolif"), sp_rho("RL_core", "UMI"), sp_rho("RL_core", hubs[1])),
           paste(sum(slim$sample == show), "spots"))
  pF <- ggplot2::ggplot(fl, ggplot2::aes(y, -x, colour = v)) + ggplot2::geom_point(size = 0.38, stroke = 0, shape = 16) +
    ggplot2::facet_wrap(~ feature, nrow = 1) + ggplot2::coord_fixed() +
    ggplot2::scale_colour_viridis_c(option = "inferno", name = "Percentile\nwithin section",
                                    breaks = c(0, 0.5, 1), labels = c("0", "50", "100")) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Visium section %s (%s spots)", show,
                                                       format(sum(slim$sample == show), big.mark = ","))) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())

  # G：各切片内 spot 水平的相关（22 张切片），以及校正测序深度（和增殖）后的偏相关
  dc <- utils::read.csv(file.path(d9, "Spatial_depth_check_by_slice.csv"))
  gv <- c(rho_UMI_RLcore = "Sequencing depth (UMI)",
          rho_RLcore_Prolif = "Proliferation", prho_RLcore_Prolif_given_UMI = "Proliferation | UMI",
          rho_RLcore_IFNg = "IFN-γ response", prho_RLcore_IFNg_given_Prolif_UMI = "IFN-γ | proliferation, UMI",
          rho_RLcore_GEP18 = "T cell-inflamed GEP", prho_RLcore_GEP18_given_Prolif_UMI = "GEP | proliferation, UMI")
  gl <- do.call(rbind, lapply(names(gv), function(v) data.frame(sample = dc$sample, what = gv[[v]], rho = dc[[v]],
                                                                adj = grepl("^prho", v))))
  gl$what <- factor(gl$what, levels = rev(gv)); gl$adj[gl$what == gv[["rho_UMI_RLcore"]]] <- NA
  mr <- function(x) sprintf("%.2f (%.2f to %.2f)", stats::median(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  for (v in names(gv)) add_stat(9, "G", paste("spot-level correlation with RL_core, median (range) over sections:", gv[[v]]),
                                mr(dc[[v]]), sprintf("%d sections", sum(!is.na(dc[[v]]))))
  for (v in c("rho_UMI_RLcore", "rho_UMI_Prolif", "rho_UMI_IFNg"))
    add_stat(9, "G", paste("spot-level Spearman rho with log10 UMI, median (range):", v), mr(dc[[v]]))
  for (v in c("prho_RLcore_IFNg_given_UMI", "rho_Prolif_IFNg", "prho_Prolif_IFNg_given_UMI"))
    add_stat(9, "G", paste("spot-level correlation, median (range):", v), mr(dc[[v]]))
  pG <- ggplot2::ggplot(gl, ggplot2::aes(rho, what)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = adj), outlier.shape = NA, width = 0.6, linewidth = 0.3) +
    ggplot2::geom_jitter(height = 0.12, width = 0, size = 0.5, colour = "grey25") +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#F5B7B1"), na.value = "grey55", guide = "none") +
    ggplot2::labs(x = "Spot-level correlation with RL_core\n(one point per section)", y = NULL, title = "Programmes")

  # H：核心基因与 RL_core 的 spot 水平相关（未校正 vs 校正增殖和深度）
  hp <- utils::read.csv(file.path(d9, "Spatial_hub_partial_by_slice.csv"))
  hl <- rbind(data.frame(gene = hp$gene, what = "Unadjusted", rho = hp$rho),
              data.frame(gene = hp$gene, what = "Adjusted for proliferation and UMI", rho = hp$prho_given_Prolif_UMI))
  hl <- hl[!is.na(hl$rho), ]
  hl$gene <- factor(hl$gene, levels = rev(hubs)); hl$what <- factor(hl$what, levels = c("Unadjusted", "Adjusted for proliferation and UMI"))
  for (g in hubs) {
    h <- hp[hp$gene == g, ]
    add_stat(9, "H", paste("spot-level rho with RL_core, median (range): unadjusted / adjusted:", g),
             paste(mr(h$rho), "/", mr(h$prho_given_Prolif_UMI)),
             sprintf("detected in median %.1f%% of spots; %d sections with detection; adjusted rho > 0 in %d/%d; UMI-only adjusted %s",
                     100 * stats::median(h$pct_detected), sum(!is.na(h$rho)), sum(h$prho_given_Prolif_UMI > 0, na.rm = TRUE),
                     sum(!is.na(h$prho_given_Prolif_UMI)), mr(h$prho_given_UMI)))
  }
  pH <- ggplot2::ggplot(hl, ggplot2::aes(rho, gene, fill = what)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_boxplot(outlier.size = 0.3, width = 0.65, linewidth = 0.3, position = ggplot2::position_dodge(width = 0.75)) +
    ggplot2::scale_fill_manual(values = c(Unadjusted = "grey85", `Adjusted for proliferation and UMI` = "#F5B7B1"), name = NULL) +
    ggplot2::labs(x = "Spot-level correlation with RL_core\n(22 sections)", y = NULL, title = "Model genes") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom", legend.direction = "vertical")

  r1 <- (pA | (pB / pC + patchwork::plot_layout(heights = c(0.8, 1)))) + patchwork::plot_layout(widths = c(1, 1.1))
  r2 <- (pD | pE) + patchwork::plot_layout(widths = c(1, 1.15))
  r4 <- (pG | pH) + patchwork::plot_layout(widths = c(1.2, 1))
  p <- (r1 / r2 / pF / r4) + patchwork::plot_layout(heights = c(1.45, 0.8, 0.62, 0.8)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_fig(p, 9, 180, 240)
}
SP_SHOW <- c("BC01", "BC13")    # Figure 9F 的示例切片（看过 09 步全部切片后选定）

for (n in WHICH) {
  f <- get0(paste0("fig", n))
  if (is.null(f)) { log_msg("Figure ", n, " 的绘图函数还没写，跳过"); next }
  f()
}
log_msg("11 完成")

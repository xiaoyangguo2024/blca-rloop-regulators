# =============================================================================
# 08_tme_therapy.R —— 免疫微环境与治疗反应
#   TCGA：IOBR 7 种免疫浸润算法、免疫检查点、TIDE（网页版）、oncoPredict（GDSC2）
#   免疫治疗：IMvigor210（阿替利珠单抗）、GSE176307（多种 PD-1/PD-L1 抑制剂）
#   化疗：GSE169455（铂类新辅助 125 例为主分析，加诱导 24 例为敏感性分析）
# 手动准备：TIDE 网页结果；oncoPredict 的 GDSC2 训练文件（见 README）
# 输出：results/08_tme_therapy/
# =============================================================================
source("R/00_config.R")
need(c("survival", "survminer", "pROC", "ggplot2", "ggpubr"))
STEP <- "08_tme_therapy"; OUT <- res_dir(STEP)
tcga <- readRDS(file.path(DIR_PRO, "TCGA_BLCA.rds"))
tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
rlrs <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
d <- merge(tum, rlrs[rlrs$cohort == "TCGA_BLCA", c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
d$group <- factor(d$RLRS_group, levels = c("High", "Low"))
PAL <- c(High = "#D7301F", Low = "#2C7FB8")

wilcox_table <- function(df, vars, group = "group") {
  out <- do.call(rbind, lapply(vars, function(v) {
    x <- df[[v]]; g <- df[[group]]
    data.frame(variable = v,
               median_High = stats::median(x[g == "High"], na.rm = TRUE),
               median_Low = stats::median(x[g == "Low"], na.rm = TRUE),
               rho_RLRS = suppressWarnings(stats::cor(x, df$RLRS, method = "spearman", use = "complete.obs")),
               p = tryCatch(stats::wilcox.test(x ~ g)$p.value, error = function(e) NA_real_))
  }))
  out$FDR <- stats::p.adjust(out$p, "BH")
  out
}

# =============================================================================
# 1. 免疫浸润（IOBR，输入非 log 的 TPM）
# =============================================================================
f_tme <- file.path(DIR_PRO, "TCGA_IOBR_tme.rds")
if (!file.exists(f_tme) && requireNamespace("IOBR", quietly = TRUE)) {
  tpm <- tcga$tpm[, d$ID]
  meths <- c("cibersort", "mcpcounter", "xcell", "epic", "timer", "quantiseq", "estimate")
  tme <- lapply(meths, function(m) tryCatch({
    args <- list(eset = tpm, method = m, arrays = FALSE)
    if (m == "cibersort") args$perm <- 200
    if (m == "timer") args$group_list <- rep("blca", ncol(tpm))
    if (m == "quantiseq") { args$tumor <- TRUE; args$scale_mrna <- TRUE }
    log_msg("IOBR：", m)
    do.call(IOBR::deconvo_tme, args)
  }, error = function(e) { log_msg("IOBR ", m, " 失败：", conditionMessage(e)); NULL }))
  tme <- Reduce(function(a, b) merge(a, b, by = "ID"), Filter(Negate(is.null), tme))
  saveRDS(tme, f_tme)
}
if (file.exists(f_tme)) {
  tme <- readRDS(f_tme)
  dt <- merge(d[, c("ID", "RLRS", "group")], tme, by = "ID")
  cells <- setdiff(colnames(tme), "ID")
  cells <- cells[vapply(dt[, cells], function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  wt <- wilcox_table(dt, cells)
  utils::write.csv(wt[order(wt$p), ], file.path(OUT, "TME_High_vs_Low_all_methods.csv"), row.names = FALSE)

  cib <- grep("_CIBERSORT$", cells, value = TRUE)
  long <- do.call(rbind, lapply(cib, function(v) data.frame(cell = sub("_CIBERSORT$", "", v), value = dt[[v]], group = dt$group)))
  p1 <- ggpubr::ggboxplot(long, "cell", "value", color = "group", palette = PAL, outlier.shape = NA) +
    ggpubr::stat_compare_means(ggplot2::aes(group = group), method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
    ggpubr::rotate_x_text(60) + ggplot2::labs(x = NULL, y = "CIBERSORT fraction")
  top <- utils::head(wt[order(-abs(wt$rho_RLRS)), ], 30)
  top$variable <- factor(top$variable, levels = top$variable[order(top$rho_RLRS)])
  p2 <- ggplot2::ggplot(top, ggplot2::aes(rho_RLRS, variable, color = FDR < 0.05)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = rho_RLRS, yend = variable)) + ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "grey60")) +
    ggplot2::theme_classic() + ggplot2::labs(x = "Spearman rho with RLRS", y = NULL)
  est <- grep("_estimate$", cells, value = TRUE)
  p3 <- lapply(est, function(v) ggpubr::ggboxplot(dt, "group", v, color = "group", palette = PAL, add = "jitter") +
                 ggpubr::stat_compare_means(method = "wilcox.test"))
  save_print(c(list(p1, p2), p3), file.path(OUT, "Fig8_TME.pdf"), 11, 6)
} else log_msg("未安装 IOBR 或反卷积失败，跳过免疫浸润")

# =============================================================================
# 2. 免疫检查点
# =============================================================================
cp <- intersect(CHECKPOINTS, rownames(tcga$expr))
dc <- data.frame(d[, c("ID", "RLRS", "group")], t(tcga$expr[cp, d$ID]), check.names = FALSE)
utils::write.csv(wilcox_table(dc, cp), file.path(OUT, "Checkpoints_High_vs_Low.csv"), row.names = FALSE)
long <- do.call(rbind, lapply(cp, function(g) data.frame(gene = g, value = dc[[g]], group = dc$group)))
p <- ggpubr::ggboxplot(long, "gene", "value", color = "group", palette = PAL, outlier.shape = NA) +
  ggpubr::stat_compare_means(ggplot2::aes(group = group), method = "wilcox.test", label = "p.signif") +
  ggpubr::rotate_x_text(45) + ggplot2::labs(x = NULL, y = "log2(TPM + 1)")
save_gg(p, file.path(OUT, "Fig8_checkpoints.pdf"), 9, 4.5)

# =============================================================================
# 3. TIDE：优先用本地 tidepy（TIDE 作者的官方 Python 实现，环境见 env/README.md）；
#    没有本地环境时，把 TIDE_input.txt 上传到 tide.dfci.harvard.edu，结果存为 TIDE_output.csv 再重跑
# =============================================================================
tide_in <- tcga$expr[, d$ID]
tide_in <- tide_in - rowMeans(tide_in)                 # 按基因中心化（TIDE 要求）
utils::write.table(round(tide_in, 4), file.path(OUT, "TIDE_input.txt"), sep = "\t", quote = FALSE, col.names = NA)
f_tide <- file.path(OUT, "TIDE_output.csv")
py <- file.path(PROJ, "conda_py", "bin", "python")
if (!file.exists(f_tide) && file.exists(py)) {
  log_msg("用本地 tidepy 计算 TIDE")
  system2(py, c("-I", shQuote(file.path(PROJ, "R", "tide_run.py")), shQuote(file.path(OUT, "TIDE_input.txt")), shQuote(f_tide)))
}
if (file.exists(f_tide)) {
  ti <- utils::read.csv(f_tide, check.names = FALSE)
  dt <- merge(d[, c("ID", "RLRS", "group")], ti, by.x = "ID", by.y = "Patient")
  vars <- intersect(c("TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2", "IFNG", "CD8"), colnames(dt))
  utils::write.csv(wilcox_table(dt, vars), file.path(OUT, "TIDE_High_vs_Low.csv"), row.names = FALSE)
  pl <- lapply(vars, function(v) ggpubr::ggviolin(dt, "group", v, fill = "group", palette = PAL, add = "boxplot") +
                 ggpubr::stat_compare_means(method = "wilcox.test"))
  if ("Responder" %in% colnames(dt)) {
    tb <- table(dt$group, dt$Responder); print(tb)
    pl$resp <- ggplot2::ggplot(dt, ggplot2::aes(group, fill = Responder)) + ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(y = "Proportion", title = sprintf("Chi-square P = %.3g", stats::chisq.test(tb)$p.value)) +
      ggplot2::theme_classic()
  }
  save_print(pl, file.path(OUT, "Fig8_TIDE.pdf"), 4.5, 4.5)
} else log_msg("尚无 TIDE_output.csv：请上传 ", file.path(OUT, "TIDE_input.txt"), " 到 TIDE 网站（癌种选 Other）")

# =============================================================================
# 4. 免疫治疗队列：IMvigor210、GSE176307
# =============================================================================
ici_plots <- list(); ici_rows <- list()
for (n in c("IMvigor210", "GSE176307")) {
  cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  e <- merge(cl, rlrs[rlrs$cohort == n, c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
  e$group <- factor(e$RLRS_group, levels = c("High", "Low"))
  r <- e[e$response %in% c("CR/PR", "SD/PD"), ]
  r$response <- factor(r$response, levels = c("SD/PD", "CR/PR"))
  # AUC 方向固定（RLRS 越高 → 越可能缓解）；不用 direction = "auto"，否则 AUC 永远 ≥ 0.5、偏乐观。AUC < 0.5 表示缓解者 RLRS 更低
  roc <- pROC::roc(r$response, r$RLRS, levels = c("SD/PD", "CR/PR"), direction = "<", quiet = TRUE)
  ci  <- as.numeric(pROC::ci.auc(roc))
  tb  <- table(r$group, r$response)
  ici_rows[[n]] <- data.frame(cohort = n, n_evaluable = nrow(r), n_responders = sum(r$response == "CR/PR"),
                              AUC = ci[2], AUC_lower = ci[1], AUC_upper = ci[3],
                              p_wilcox = stats::wilcox.test(RLRS ~ response, r)$p.value,
                              ORR_High = mean(r$response[r$group == "High"] == "CR/PR"),
                              ORR_Low = mean(r$response[r$group == "Low"] == "CR/PR"),
                              p_fisher = stats::fisher.test(tb)$p.value)
  ici_plots[[paste(n, "box")]] <- ggpubr::ggboxplot(r, "response", "RLRS", color = "response", add = "jitter") +
    ggpubr::stat_compare_means(method = "wilcox.test") + ggplot2::ggtitle(n)
  ici_plots[[paste(n, "bar")]] <- ggplot2::ggplot(r, ggplot2::aes(group, fill = response)) +
    ggplot2::geom_bar(position = "fill") + ggplot2::theme_classic() +
    ggplot2::labs(y = "Proportion", title = sprintf("%s  Fisher P = %.3g", n, stats::fisher.test(tb)$p.value))
  ici_plots[[paste(n, "roc")]] <- pROC::ggroc(roc, legacy.axes = TRUE) + ggplot2::geom_abline(linetype = 2) +
    ggplot2::theme_classic() + ggplot2::ggtitle(sprintf("%s  AUC = %.3f (%.3f-%.3f)", n, ci[2], ci[1], ci[3]))
  s <- e[!is.na(e$OS.time) & e$OS.time > 0 & !is.na(e$OS), ]
  ici_plots[[paste(n, "KM")]] <- km_plot(s, title = paste(n, "OS"))
  if ("TMB" %in% colnames(s) && sum(!is.na(s$TMB)) > 40) {
    t2 <- s[!is.na(s$TMB), ]
    t2$TMB_group <- ifelse(t2$TMB >= stats::median(t2$TMB), "TMB-H", "TMB-L")
    t2$combo <- factor(paste(t2$TMB_group, paste0("RLRS-", t2$RLRS_group), sep = " + "))
    ici_plots[[paste(n, "KM_TMB")]] <- km_plot(t2, group = "combo", title = paste(n, "RLRS x TMB"),
                                               palette = c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A"))
  }
  if (n == "IMvigor210") {
    ip <- e[!is.na(e$immune_phenotype), ]
    ici_plots$IMv_phenotype_box <- ggpubr::ggboxplot(ip, "immune_phenotype", "RLRS", color = "immune_phenotype",
                                                     add = "jitter") + ggpubr::stat_compare_means(method = "kruskal.test")
    ici_plots$IMv_phenotype_bar <- ggplot2::ggplot(ip, ggplot2::aes(group, fill = immune_phenotype)) +
      ggplot2::geom_bar(position = "fill") + ggplot2::theme_classic() +
      ggplot2::labs(y = "Proportion", title = sprintf("Chi-square P = %.3g",
                                                      stats::chisq.test(table(ip$group, ip$immune_phenotype))$p.value))
  }
}
ici_tab <- do.call(rbind, ici_rows); print(ici_tab, row.names = FALSE)
utils::write.csv(ici_tab, file.path(OUT, "ICI_response_summary.csv"), row.names = FALSE)
save_print(ici_plots, file.path(OUT, "Fig8_immunotherapy.pdf"), 6, 6)

# =============================================================================
# 5. 新辅助化疗：GSE169455（主分析：新辅助 125 例；敏感性分析：全部 149 例）
# =============================================================================
cl <- readRDS(file.path(DIR_PRO, "GSE169455.rds"))$clin
e <- merge(cl, rlrs[rlrs$cohort == "GSE169455", c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
chemo_rows <- list(); chemo_plots <- list()
for (set in c("neoadjuvant", "all")) {
  s <- if (set == "all") e else e[e$setting == "neoadjuvant", ]
  for (ep in c("pCR", "downstage")) {
    x <- s[!is.na(s[[ep]]), ]
    x$outcome <- factor(ifelse(x[[ep]] == 1, "Yes", "No"), levels = c("No", "Yes"))
    roc <- pROC::roc(x$outcome, x$RLRS, levels = c("No", "Yes"), direction = "<", quiet = TRUE)   # 方向固定，同上
    ci <- as.numeric(pROC::ci.auc(roc))
    lg <- summary(stats::glm(outcome ~ RLRS_z, data = x, family = stats::binomial()))$coefficients
    lg_adj <- tryCatch(summary(stats::glm(outcome ~ RLRS_z + consensus, data = x, family = stats::binomial()))$coefficients,
                       error = function(err) NULL)
    chemo_rows[[paste(set, ep)]] <- data.frame(
      set = set, endpoint = ep, n = nrow(x), n_yes = sum(x$outcome == "Yes"), AUC = ci[2], AUC_lower = ci[1],
      AUC_upper = ci[3], p_wilcox = stats::wilcox.test(RLRS ~ outcome, x)$p.value,
      OR_perSD = exp(lg["RLRS_z", 1]), p_logistic = lg["RLRS_z", 4],
      OR_perSD_adj_subtype = if (!is.null(lg_adj)) exp(lg_adj["RLRS_z", 1]) else NA,
      p_adj_subtype = if (!is.null(lg_adj)) lg_adj["RLRS_z", 4] else NA)
    if (set == "neoadjuvant") {
      chemo_plots[[paste(ep, "box")]] <- ggpubr::ggboxplot(x, "outcome", "RLRS", color = "outcome", add = "jitter") +
        ggpubr::stat_compare_means(method = "wilcox.test") + ggplot2::labs(x = ep, title = "GSE169455 neoadjuvant")
      chemo_plots[[paste(ep, "roc")]] <- pROC::ggroc(roc, legacy.axes = TRUE) + ggplot2::geom_abline(linetype = 2) +
        ggplot2::theme_classic() + ggplot2::ggtitle(sprintf("%s  AUC = %.3f (%.3f-%.3f)", ep, ci[2], ci[1], ci[3]))
    }
  }
}
chemo_plots$subtype <- ggpubr::ggboxplot(e[!is.na(e$consensus), ], "consensus", "RLRS", color = "consensus",
                                         add = "jitter") + ggpubr::stat_compare_means(method = "kruskal.test")
chemo_tab <- do.call(rbind, chemo_rows); print(chemo_tab, row.names = FALSE)
utils::write.csv(chemo_tab, file.path(OUT, "Chemo_response_GSE169455.csv"), row.names = FALSE)
save_print(chemo_plots, file.path(OUT, "Fig8_chemotherapy.pdf"), 5, 5)

# =============================================================================
# 6. 药敏预测：oncoPredict（GDSC2）
# =============================================================================
op_dir <- file.path(DIR_RAW, "oncoPredict")
f_gx <- file.path(op_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")
f_gr <- file.path(op_dir, "GDSC2_Res.rds")
f_pred <- file.path(OUT, "oncoPredict", "calcPhenotype_Output", "DrugPredictions.csv")
if (!file.exists(f_pred) && file.exists(f_gx) && file.exists(f_gr) && requireNamespace("oncoPredict", quietly = TRUE)) {
  dir.create(file.path(OUT, "oncoPredict"), showWarnings = FALSE)
  old <- setwd(file.path(OUT, "oncoPredict"))                 # calcPhenotype 把结果写到当前目录
  oncoPredict::calcPhenotype(trainingExprData = readRDS(f_gx), trainingPtype = exp(readRDS(f_gr)),  # GDSC2_Res 为 ln(IC50)
                             testExprData = as.matrix(tcga$expr[, d$ID]), batchCorrect = "eb",
                             powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2, minNumSamples = 10,
                             printOutput = TRUE, removeLowVaringGenesFrom = "rawData",
                             folder = TRUE)   # oncoPredict ≥ 1.3 默认 folder = FALSE，不写 DrugPredictions.csv
  setwd(old)
}
if (file.exists(f_pred)) {
  pred <- utils::read.csv(f_pred, row.names = 1, check.names = FALSE)
  dp <- data.frame(d[match(rownames(pred), d$ID), c("ID", "RLRS", "group")], pred, check.names = FALSE)
  drugs <- colnames(pred)
  wt <- wilcox_table(dp, drugs)
  utils::write.csv(wt[order(wt$p), ], file.path(OUT, "oncoPredict_High_vs_Low.csv"), row.names = FALSE)
  focus <- grep("^(Cisplatin|Gemcitabine|Docetaxel|Paclitaxel|Vinblastine|Vinorelbine|Mitomycin|Doxorubicin|Epirubicin|Methotrexate|Olaparib|Talazoparib|Niraparib|AZD6738|VE-822|Erdafitinib)",
                drugs, value = TRUE)
  # 药名含连字符等字符（如 VE-822_1613），ggpubr 会把列名当表达式解析而报错：改用 ggplot2 的 .data[[v]]
  pl <- lapply(focus, function(v) ggplot2::ggplot(dp, ggplot2::aes(group, .data[[v]], colour = group)) +
                 ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, size = 0.5) +
                 ggplot2::scale_colour_manual(values = PAL) + ggpubr::stat_compare_means(method = "wilcox.test") +
                 ggplot2::theme_classic() + ggplot2::labs(x = NULL, y = paste(v, "predicted IC50")))
  save_print(pl, file.path(OUT, "Fig8_drug_sensitivity.pdf"), 4, 4.5)
} else log_msg("未运行 oncoPredict：请按 README 下载 GDSC2 训练文件到 ", op_dir)

save_session(STEP)
log_msg("08 完成")

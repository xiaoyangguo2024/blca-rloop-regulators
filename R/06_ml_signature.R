# =============================================================================
# 06_ml_signature.R —— Mime1：10 种算法、101 种组合构建 R-loop 相关风险评分（RLRS）
#   训练：TCGA-BLCA；外部队列：GSE13507、GSE32894（DSS）、GSE31684、GSE48075
#   每个队列内部按基因 z-score（Mime1 不会替你做，已查源码确认）
#   只用 2 个"调参队列"（TUNE）的平均 C-index 选模型，相差 < 0.005 时选基因少的；
#   另 2 个外部队列（TEST）和 3 个治疗队列不参与选模型，作为真正独立的检验（避免"用验证集挑模型"的乐观偏倚）
# 输出：data/processed/Mime_res.rds、RLRS_model.rds、RLRS_all_cohorts.rds；results/06_ml_signature/
# 耗时：16 核约数小时（结果缓存，重跑直接读取）
# =============================================================================
source("R/00_config.R")
need(c("Mime1", "survival", "survminer", "timeROC", "ggplot2"))
# Mime1 的 NAMESPACE 没有 import，内部直接调用 rfsrc()、CoxBoost()、cv.glmnet()、superpc.train()、survivalsvm()、
# gbm()、plsRcox()、rownames_to_column()、%>%、insert_right() 等，需要这些包在搜索路径上（与 Mime1 教程的用法一致）
for (p in c("randomForestSRC", "glmnet", "CoxBoost", "plsRcox", "superpc", "gbm", "survivalsvm", "dplyr", "tibble", "aplot"))
  suppressPackageStartupMessages(library(p, character.only = TRUE))
source("R/compat_rfsrc_var_select.R")   # randomForestSRC ≥ 3.4.0 删除了 var.select()，补回 3.3.1 的原版
STEP <- "06_ml_signature"; OUT <- res_dir(STEP)
TRAIN <- "TCGA_BLCA"
VALID <- c("GSE13507", "GSE32894", "GSE31684", "GSE48075")
TUNE  <- c("GSE13507", "GSE31684")     # 用于挑模型
TEST  <- setdiff(VALID, TUNE)          # 不参与挑模型
THER  <- c("IMvigor210", "GSE176307", "GSE169455")
BEST_MODEL <- NULL     # 想手动指定模型时填模型名，例如 "StepCox[both] + Ridge"

cands <- readRDS(file.path(DIR_PRO, "candidates.rds"))
genes <- cands$availability$gene[cands$availability$all_cohorts]   # 只用 4 个验证队列都有的候选基因
log_msg("进入 Mime1 的候选基因：", length(genes), " 个（Mime1 内部会再按训练集单因素 Cox P < ",
        PARAM$unicox_p, " 筛选）")
if (length(genes) < 10) log_msg("【核对】候选基因少于 10 个，建模意义有限，考虑放宽 05 步阈值")

load_cohort <- function(n) {
  o <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
  if (n == "TCGA_BLCA") o$clin <- o$clin[o$clin$group == "Tumor", ]
  o
}
coh <- lapply(stats::setNames(c(TRAIN, VALID, THER), c(TRAIN, VALID, THER)), load_cohort)
inp <- lapply(coh[c(TRAIN, VALID)], function(o) build_mime_input(o$expr, o$clin, genes))
print(vapply(inp, function(x) c(n = nrow(x), events = sum(x$OS)), numeric(2)))

# =============================================================================
# 1. 101 种组合
# =============================================================================
f_res <- file.path(DIR_PRO, "Mime_res.rds")
if (!file.exists(f_res)) {
  res <- Mime1::ML.Dev.Prog.Sig(train_data = inp[[TRAIN]], list_train_vali_Data = inp,
                                unicox.filter.for.candi = TRUE, unicox_p_cutoff = PARAM$unicox_p,
                                candidate_genes = safe_names(genes), mode = "all", nodesize = 5,
                                seed = PARAM$seed, cores_for_parallel = 1)
  # cores_for_parallel 只用于 gbm 的交叉验证；gbm 在 n.cores > 1 时建 PSOCK 集群（需要本地端口，沙箱里会失败）。
  # n.cores = 1 时 gbm 逐折串行，每折仍用 gbm 内部固定的种子，模型设置不变。
  # CoxBoost 的 parallel = TRUE 依赖 snowfall（未安装），自动退回串行
  saveRDS(res, f_res)
}
res <- readRDS(f_res)

# =============================================================================
# 2. 选模型：验证集平均 C-index（不看训练集），接近时选基因少的
# =============================================================================
# Mime1 的 Cindex.res 在每个队列重新拟合 coxph(Surv ~ RS)，方向不固定：外部队列里关联反向的模型也会得到 > 0.5。
# 这里改用固定方向的 C-index：方向由训练集（TCGA）的 Cox 系数确定，外部队列不再调整
ci_mime <- res$Cindex.res
dir_c <- function(rs, d) survival::concordance(survival::Surv(d$OS.time, d$OS) ~ rs, reverse = TRUE)$concordance
if (!is.null(res$riskscore)) {
  ci <- do.call(rbind, lapply(names(res$riskscore), function(m) {
    rsl <- res$riskscore[[m]]
    s <- tryCatch(sign(stats::coef(survival::coxph(survival::Surv(OS.time, OS) ~ RS, data = rsl[[TRAIN]]))),
                  error = function(e) 1)
    if (!is.finite(s) || s == 0) s <- 1
    data.frame(ID = names(rsl), Cindex = vapply(names(rsl), function(n) dir_c(s * rsl[[n]]$RS, rsl[[n]]), numeric(1)),
               Model = m)
  }))
  res$Cindex.res <- ci          # 之后 Mime1 的作图函数显示的也是固定方向的 C-index
  chg <- merge(ci, ci_mime, by = c("Model", "ID"), suffixes = c("_dir", "_mime"))
  log_msg("固定方向后 C-index 改变 > 0.001 的 模型×队列：", sum(abs(chg$Cindex_dir - chg$Cindex_mime) > 0.001),
          " / ", nrow(chg), "（其中外部队列反向关联 C < 0.5：", sum(chg$Cindex_dir < 0.5), "）")
} else {
  log_msg("【核对】Mime1 结果里没有 riskscore，只能用 Mime1 原始 C-index 选模型")
  ci <- ci_mime
}
ci_w <- stats::reshape(ci, idvar = "Model", timevar = "ID", direction = "wide")
colnames(ci_w) <- sub("^Cindex\\.", "", colnames(ci_w))
ci_w$mean_tune  <- rowMeans(ci_w[, TUNE, drop = FALSE])
ci_w$mean_test  <- rowMeans(ci_w[, TEST, drop = FALSE])
ci_w$mean_valid <- rowMeans(ci_w[, VALID, drop = FALSE])
ci_w$mean_all <- rowMeans(ci_w[, c(TRAIN, VALID), drop = FALSE])
ci_w$n_genes <- vapply(ci_w$Model, function(m) length(model_genes(res$ml.res[[m]], res$Sig.genes, quiet = TRUE)),
                       integer(1))
ci_w <- ci_w[order(-ci_w$mean_tune), ]
near <- ci_w[ci_w$mean_tune >= max(ci_w$mean_tune) - 0.005, ]
best <- if (!is.null(BEST_MODEL)) BEST_MODEL else near$Model[order(near$n_genes, -near$mean_tune)][1]
log_msg("最优模型：", best, "（调参队列平均 C-index = ", round(ci_w$mean_tune[ci_w$Model == best], 3),
        "；独立检验队列平均 C-index = ", round(ci_w$mean_test[ci_w$Model == best], 3),
        "，", ci_w$n_genes[ci_w$Model == best], " 个基因）【核对】")
log_msg("该模型在 101 个组合中按检验队列平均 C-index 排第 ", rank(-ci_w$mean_test)[ci_w$Model == best], " 位")
print(utils::head(ci_w[, c("Model", TRAIN, VALID, "mean_tune", "mean_test", "n_genes")], 10), row.names = FALSE)
utils::write.csv(ci_w, file.path(OUT, "Table_S4_Cindex_101_models.csv"), row.names = FALSE)

save_print(Mime1::cindex_dis_all(res, validate_set = VALID, order = c(TRAIN, VALID), width = 0.35),
           file.path(OUT, "Fig6_Cindex_101_models.pdf"), 10, 16)
save_print(Mime1::cindex_dis_select(res, model = best, order = c(TRAIN, VALID)),
           file.path(OUT, "Fig6_Cindex_best_model.pdf"), 6, 4)

# 模型基因（安全名 → 原始基因名）
sig_safe <- model_genes(res$ml.res[[best]], res$Sig.genes)
sig_genes <- genes[match(sig_safe, safe_names(genes))]
log_msg("模型基因（", length(sig_genes), "）：", paste(sig_genes, collapse = ", "))
fit_best <- res$ml.res[[best]]
coef_tab <- tryCatch({
  if (inherits(fit_best, "coxph")) data.frame(gene = names(stats::coef(fit_best)), coef = stats::coef(fit_best))
  else if (inherits(fit_best, "cv.glmnet")) {
    cf <- as.matrix(stats::coef(fit_best, s = fit_best$lambda.min)); data.frame(gene = rownames(cf), coef = cf[, 1])
  } else NULL
}, error = function(e) NULL)
if (!is.null(coef_tab)) {
  coef_tab <- coef_tab[coef_tab$coef != 0, ]
  coef_tab$gene <- genes[match(coef_tab$gene, safe_names(genes))]
  utils::write.csv(coef_tab, file.path(OUT, "RLRS_coefficients.csv"), row.names = FALSE)
} else log_msg("模型对象（", class(fit_best)[1], "）不能直接取系数；survivalsvm 默认线性核时 RLRS 仍是基因 z 分数的线性组合，系数在 09 步用线性回归还原（RLRS_coefficients.csv）")

# =============================================================================
# 3. RLRS：方向统一为"越高风险越大"
#   Mime1 用 coxph(Surv ~ RS) 的 concordance 计算 C-index，不区分方向；这里按训练集的 Cox 系数统一方向
# =============================================================================
rs_train <- res$riskscore[[best]][[TRAIN]]
sgn <- sign(stats::coef(survival::coxph(survival::Surv(OS.time, OS) ~ RS, data = rs_train)))
if (sgn < 0) log_msg("RS 与风险方向相反，已乘以 -1")

# 治疗队列：没有生存信息的样本放占位值（只为计算 RS，生存分析在 08 步用真实临床数据）
inp_ther <- lapply(coh[THER], function(o) {
  s <- o$clin[, c("ID", "OS.time", "OS")]
  s$OS.time[is.na(s$OS.time) | s$OS.time <= 0] <- 1
  s$OS[is.na(s$OS)] <- 0
  build_mime_input(o$expr, s, genes, fill_missing = TRUE)
})
rs_ther <- Mime1::cal_RS_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = inp[[TRAIN]],
                                inputmatrix.list = inp_ther, mode = "all")

rlrs <- do.call(rbind, c(
  lapply(c(TRAIN, VALID), function(n) data.frame(cohort = n, ID = res$riskscore[[best]][[n]]$ID,
                                                 RLRS = sgn * res$riskscore[[best]][[n]]$RS)),
  lapply(THER, function(n) data.frame(cohort = n, ID = rs_ther[[best]][[n]]$ID, RLRS = sgn * rs_ther[[best]][[n]]$RS))))
rlrs$RLRS_z <- stats::ave(rlrs$RLRS, rlrs$cohort, FUN = function(x) as.numeric(scale(x)))
rlrs$RLRS_group <- stats::ave(rlrs$RLRS, rlrs$cohort, FUN = function(x) ifelse(x >= stats::median(x), "High", "Low"))
saveRDS(rlrs, file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
saveRDS(list(model = best, genes = sig_genes, genes_safe = sig_safe, sign = sgn, candidate_genes = genes,
             coef = coef_tab, cindex = ci_w), file.path(DIR_PRO, "RLRS_model.rds"))

# =============================================================================
# 4. 各队列 KM、Cox（每 1 个标准差）、时间依赖 ROC
# =============================================================================
km_list <- list(); cox_rows <- list(); auc_rows <- list()
grDevices::pdf(file.path(OUT, "Fig6_timeROC.pdf"), 5, 5)
for (n in c(TRAIN, VALID)) {
  d <- merge(inp[[n]][, c("ID", "OS.time", "OS")], rlrs[rlrs$cohort == n, ], by = "ID")
  d$group <- factor(d$RLRS_group, levels = c("High", "Low"))
  ep <- if (n == "GSE32894") "DSS" else "OS"
  km_list[[n]] <- km_plot(d, title = paste0(n, " (", ep, ")"))
  u <- cox_table(survival::coxph(survival::Surv(OS.time, OS) ~ RLRS_z, data = d))
  # Mime1 的 C-index 在每个队列重新拟合 coxph(Surv ~ RS)，方向不固定（反向关联也会 > 0.5）；
  # 论文报告固定方向（RLRS 越高风险越高）的 C-index
  cc <- survival::concordance(survival::Surv(OS.time, OS) ~ RLRS, data = d, reverse = TRUE)
  cox_rows[[n]] <- data.frame(cohort = n, role = if (n == TRAIN) "training" else if (n %in% TUNE) "tuning" else "test",
                              endpoint = ep, n = nrow(d), events = sum(d$OS), HR_perSD = u$HR,
                              lower = u$lower, upper = u$upper, p = u$p,
                              Cindex = cc$concordance, Cindex_se = sqrt(cc$var),
                              Cindex_mime = ci_mime$Cindex[ci_mime$Model == best & ci_mime$ID == n])
  yrs <- PARAM$surv_years[PARAM$surv_years * 365.25 < max(d$OS.time[d$OS == 1])]
  if (length(yrs)) {
    roc <- timeROC::timeROC(T = d$OS.time, delta = d$OS, marker = d$RLRS, cause = 1,
                            times = yrs * 365.25, iid = FALSE)
    auc_rows[[n]] <- data.frame(cohort = n, year = yrs, AUC = as.numeric(roc$AUC))
    cols <- c("#E41A1C", "#377EB8", "#4DAF4A")[seq_along(yrs)]
    for (k in seq_along(yrs)) plot(roc, time = yrs[k] * 365.25, col = cols[k], add = k > 1, title = FALSE)
    graphics::legend("bottomright", bty = "n", col = cols, lwd = 2,
                     legend = sprintf("%d-year AUC = %.3f", yrs, roc$AUC))
    graphics::title(paste0(n, " (", ep, ")"))
  }
}
grDevices::dev.off()
save_print(km_list, file.path(OUT, "Fig6_KM_RLRS.pdf"), 6, 6)
cox_out <- do.call(rbind, cox_rows); print(cox_out, row.names = FALSE)
utils::write.csv(cox_out, file.path(OUT, "RLRS_Cox_by_cohort.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, auc_rows), file.path(OUT, "RLRS_timeROC_AUC.csv"), row.names = FALSE)

# =============================================================================
# 5. 与已发表的膀胱癌预后模型比较
#    模型来自 Europe PMC 开放全文中印出的风险评分公式（逐条核对，见 published_BLCA_signatures_sources.csv）
#    按原文系数 × 各队列 log2 表达计算（不做 z-score）；基因名先用 HGNC 现行名，队列里没有时用原文写法
#    某队列中模型基因覆盖 < 80% 时，该模型在该队列不计算
#    C-index 固定方向（原文系数下分数越高风险越高），与 RLRS 的固定方向 C-index 可直接比较
# =============================================================================
f_pub <- file.path(DIR_GS, "published_BLCA_signatures.csv")
pub <- if (file.exists(f_pub)) utils::read.csv(f_pub) else data.frame()
if (nrow(pub) > 0) {
  PUBCOH <- c(TRAIN, VALID, "IMvigor210", "GSE176307")       # 有 OS（GSE32894 为 DSS）的队列
  surv_of <- function(n) {
    cl <- coh[[n]]$clin; cl <- cl[!is.na(cl$OS.time) & !is.na(cl$OS) & cl$OS.time > 0, c("ID", "OS.time", "OS")]
    cl[cl$ID %in% colnames(coh[[n]]$expr), ]
  }
  tab <- do.call(rbind, lapply(PUBCOH, function(n) {
    ex <- coh[[n]]$expr; sv <- surv_of(n)
    do.call(rbind, lapply(split(pub, pub$model), function(ct) {
      use <- ifelse(ct$symbol %in% rownames(ex), ct$symbol,
                    ifelse(ct$symbol_published %in% rownames(ex), ct$symbol_published, NA))
      ok <- !is.na(use)
      cover <- mean(ok)
      ci <- if (cover >= 0.8 && sum(ok) >= 2) {
        rs <- as.numeric(crossprod(as.matrix(ex[use[ok], sv$ID, drop = FALSE]), ct$Coef[ok]))
        if (stats::sd(rs) > 0) dir_c(rs, sv) else NA_real_
      } else NA_real_
      data.frame(model = ct$model[1], PMID = ct$PMID[1], cohort = n, n_genes = nrow(ct), n_present = sum(ok),
                 coverage = round(cover, 2), Cindex = ci)
    }))
  }))
  rl <- do.call(rbind, lapply(PUBCOH, function(n) {
    d <- merge(surv_of(n), rlrs[rlrs$cohort == n, c("ID", "RLRS")], by = "ID")
    data.frame(model = "RLRS (this study)", PMID = NA, cohort = n, n_genes = length(sig_genes),
               n_present = length(sig_genes), coverage = 1, Cindex = dir_c(d$RLRS, d))
  }))
  tab <- rbind(rl, tab)
  tab$role <- ifelse(tab$cohort == TRAIN, "training", ifelse(tab$cohort %in% TUNE, "tuning", "test"))
  utils::write.csv(tab, file.path(OUT, "Table_S5_published_models_Cindex.csv"), row.names = FALSE)

  ev <- tab[!is.na(tab$Cindex), ]
  rk <- do.call(rbind, lapply(split(ev, ev$cohort), function(x) data.frame(
    cohort = x$cohort[1], role = x$role[1], n_models = nrow(x),
    RLRS_Cindex = x$Cindex[x$model == "RLRS (this study)"],
    RLRS_rank = rank(-x$Cindex, ties.method = "min")[x$model == "RLRS (this study)"],
    published_median = stats::median(x$Cindex[x$model != "RLRS (this study)"]),
    published_max = max(x$Cindex[x$model != "RLRS (this study)"]))))
  rk <- rk[match(PUBCOH, rk$cohort), ]
  print(rk, row.names = FALSE)
  utils::write.csv(rk, file.path(OUT, "RLRS_rank_among_published.csv"), row.names = FALSE)
  # 外部队列平均（不含训练集）：要求模型在全部外部队列都可计算
  ext <- setdiff(PUBCOH, TRAIN)
  mw <- stats::reshape(ev[ev$cohort %in% ext, c("model", "cohort", "Cindex")], idvar = "model", timevar = "cohort",
                       direction = "wide")
  mw <- mw[stats::complete.cases(mw), ]
  mw$mean_external <- rowMeans(mw[, -1, drop = FALSE])
  mw <- mw[order(-mw$mean_external), ]
  log_msg("在全部 ", length(ext), " 个外部队列都可计算的模型：", nrow(mw) - 1, " 个已发表模型 + RLRS；RLRS 外部平均 C-index 排第 ",
          which(mw$model == "RLRS (this study)"), " 位")
  utils::write.csv(mw, file.path(OUT, "published_models_mean_external_Cindex.csv"), row.names = FALSE)

  ev$is_rlrs <- ev$model == "RLRS (this study)"
  ev$cohort <- factor(ev$cohort, levels = PUBCOH)
  lab <- data.frame(cohort = factor(rk$cohort, levels = PUBCOH), y = max(ev$Cindex) + 0.03,
                    txt = sprintf("%d/%d", rk$RLRS_rank, rk$n_models))
  p <- ggplot2::ggplot(ev[!ev$is_rlrs, ], ggplot2::aes(cohort, Cindex)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey45", fill = "grey92") +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.8, colour = "grey40", alpha = 0.6) +
    ggplot2::geom_point(data = ev[ev$is_rlrs, ], colour = "#C0392B", shape = 18, size = 4) +
    ggplot2::geom_text(data = lab, ggplot2::aes(y = y, label = txt), size = 3.2) +
    ggplot2::labs(x = NULL, y = "C-index (fixed direction)",
                  subtitle = "Red diamond: RLRS; grey: published BLCA signatures; numbers: RLRS rank") +
    ggplot2::theme_bw(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  save_gg(p, file.path(OUT, "Fig6_Cindex_vs_published.pdf"), 7, 4.5)
  log_msg("已比较 ", length(unique(pub$model)), " 个已发表模型")
} else log_msg("published_BLCA_signatures.csv 为空，跳过与已发表模型的比较")

save_session(STEP)
log_msg("06 完成")

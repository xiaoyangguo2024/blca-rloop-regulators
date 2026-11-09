# =============================================================================
# 07_clinical_nomogram.R —— RLRS 的临床价值
#   单/多因素 Cox（TCGA + 外部队列）、列线图（rms）、1/3/5 年校准曲线（1000 次 bootstrap）、
#   DCA（dcurves）、时间依赖 ROC 与 C-index 比较、亚组 KM
# 输出：results/07_clinical_nomogram/
# =============================================================================
source("R/00_config.R")
need(c("survival", "survminer", "rms", "dcurves", "timeROC", "ggplot2"))
STEP <- "07_clinical_nomogram"; OUT <- res_dir(STEP)
NOMO_VARS <- c("RLRS", "age", "stage_bin")      # 列线图变量；多因素 Cox 结果出来后可按需调整

tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
rlrs <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
d <- merge(tum, rlrs[rlrs$cohort == "TCGA_BLCA", c("ID", "RLRS", "RLRS_z", "RLRS_group")], by = "ID")
d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
d$time_y    <- d$OS.time / 365.25
d$age10     <- d$age / 10
d$age_group <- factor(ifelse(d$age > 65, ">65", "<=65"), levels = c("<=65", ">65"))
d$gender    <- factor(d$gender, levels = c("female", "male"))
d$stage_bin <- factor(d$stage_bin, levels = c("I-II", "III-IV"))
d$grade     <- factor(d$grade, levels = c("Low", "High"))
d$group     <- factor(d$RLRS_group, levels = c("High", "Low"))
log_msg("TCGA 纳入 ", nrow(d), " 例（事件 ", sum(d$OS), "）；低级别仅 ", sum(d$grade == "Low", na.rm = TRUE),
        " 例，分级的 HR 不稳定，解读时注意")

# =============================================================================
# 1. 单因素与多因素 Cox（TCGA）
# =============================================================================
vars <- c(RLRS_z = "RLRS (per SD)", age10 = "Age (per 10 years)", gender = "Gender (male vs female)",
          stage_bin = "Stage (III-IV vs I-II)", grade = "Grade (high vs low)")
uni <- do.call(rbind, lapply(names(vars), function(v) {
  r <- cox_table(survival::coxph(stats::as.formula(paste("survival::Surv(OS.time, OS) ~", v)), data = d))
  data.frame(variable = vars[[v]], r[, -1])
}))
fit_multi <- survival::coxph(survival::Surv(OS.time, OS) ~ RLRS_z + age10 + gender + stage_bin + grade, data = d)
multi <- data.frame(variable = unname(vars), cox_table(fit_multi)[, -1])
cox_all <- rbind(data.frame(model = "Univariate", uni), data.frame(model = "Multivariate", multi))
print(cox_all, row.names = FALSE)
utils::write.csv(cox_all, file.path(OUT, "TCGA_uni_multi_Cox.csv"), row.names = FALSE)

forest <- function(tab, title) {
  tab$variable <- factor(tab$variable, levels = rev(tab$variable))
  tab$label <- sprintf("%.2f (%.2f-%.2f), P=%s", tab$HR, tab$lower, tab$upper, format.pval(tab$p, digits = 2, eps = 0.001))
  ggplot2::ggplot(tab, ggplot2::aes(HR, variable)) + ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), color = "#D7301F") +
    ggplot2::geom_text(ggplot2::aes(x = max(tab$upper) * 1.3, label = label), hjust = 0, size = 3) +
    ggplot2::scale_x_log10() + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic() + ggplot2::theme(plot.margin = ggplot2::margin(5, 150, 5, 5)) +
    ggplot2::labs(title = title, x = "Hazard ratio (log scale)", y = NULL)
}
save_print(list(forest(uni, "Univariate Cox (TCGA-BLCA, OS)"), forest(multi, "Multivariate Cox (TCGA-BLCA, OS)")),
           file.path(OUT, "Fig7_forest_uni_multi.pdf"), 8, 3.5)

# ---- 外部队列的多因素 Cox（有哪些临床变量就用哪些）------------------------------
ext_rows <- list()
for (n in c("GSE13507", "GSE32894", "GSE31684", "GSE48075")) {
  cl <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin
  e <- merge(cl, rlrs[rlrs$cohort == n, c("ID", "RLRS_z")], by = "ID")
  e <- e[!is.na(e$OS.time) & e$OS.time > 0 & !is.na(e$OS), ]
  if ("T" %in% colnames(e)) e$MIBC <- factor(ifelse(e$T %in% c("2", "3", "4"), "T2-4", ifelse(is.na(e$T), NA, "Ta-T1")),
                                             levels = c("Ta-T1", "T2-4"))
  e$age10 <- e$age / 10
  cov <- intersect(c("age10", "gender", "MIBC"), colnames(e))
  cov <- cov[vapply(cov, function(v) length(unique(stats::na.omit(e[[v]]))) > 1, logical(1))]
  e <- e[stats::complete.cases(e[, c("OS.time", "OS", "RLRS_z", cov)]), ]
  f <- stats::as.formula(paste("survival::Surv(OS.time, OS) ~", paste(c("RLRS_z", cov), collapse = " + ")))
  m_full <- survival::coxph(f, data = e)
  r <- cox_table(m_full)
  # RLRS 在临床变量之外是否增加信息：似然比检验（同一批完整病例）；C-index 为队列内拟合（两模型同等条件）
  m_clin <- if (length(cov)) survival::coxph(stats::as.formula(paste("survival::Surv(OS.time, OS) ~",
                                                                     paste(cov, collapse = " + "))), data = e) else NULL
  ext_rows[[n]] <- data.frame(cohort = n, n = nrow(e), events = sum(e$OS), covariates = paste(cov, collapse = "+"), r,
                              LRT_p_add_RLRS = if (is.null(m_clin)) NA else stats::anova(m_clin, m_full)[2, "Pr(>|Chi|)"],
                              C_clinical = if (is.null(m_clin)) NA else survival::concordance(m_clin)$concordance,
                              C_clinical_RLRS = survival::concordance(m_full)$concordance)
}
ext_tab <- do.call(rbind, ext_rows); print(ext_tab, row.names = FALSE, digits = 3)
utils::write.csv(ext_tab, file.path(OUT, "External_multivariate_Cox.csv"), row.names = FALSE)

# =============================================================================
# 2. 列线图（rms）
# =============================================================================
dn <- d[stats::complete.cases(d[, c("time_y", "OS", NOMO_VARS)]), ]
dd <- rms::datadist(dn); options(datadist = "dd")
# rms 解析公式时用未加前缀的 Surv（survival 已由 need() 加载）
f_nomo <- stats::as.formula(paste("Surv(time_y, OS) ~", paste(NOMO_VARS, collapse = " + ")))
fit_cph <- rms::cph(f_nomo, data = dn, x = TRUE, y = TRUE, surv = TRUE)
S <- rms::Survival(fit_cph)
nom <- rms::nomogram(fit_cph, lp = FALSE,
                     fun = list(function(x) S(1, x), function(x) S(3, x), function(x) S(5, x)),
                     funlabel = c("1-year OS", "3-year OS", "5-year OS"),
                     fun.at = c(0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1))
grDevices::pdf(file.path(OUT, "Fig7_nomogram.pdf"), 9, 6); plot(nom); grDevices::dev.off()
saveRDS(nom, file.path(DIR_PRO, "nomogram_object.rds"))          # 11 步重画主图用
f_nomo_cox <- stats::as.formula(paste("survival::Surv(time_y, OS) ~", paste(NOMO_VARS, collapse = " + ")))
utils::write.csv(cox_table(survival::coxph(f_nomo_cox, data = dn)), file.path(OUT, "nomogram_model_Cox.csv"), row.names = FALSE)

# ---- 校准曲线：1/3/5 年，1000 次 bootstrap ------------------------------------
grDevices::pdf(file.path(OUT, "Fig7_calibration.pdf"), 5, 5)
cols <- c("#E41A1C", "#377EB8", "#4DAF4A")
cal_rows <- list()
for (k in seq_along(PARAM$surv_years)) {
  t <- PARAM$surv_years[k]
  fit_t <- rms::cph(f_nomo, data = dn, x = TRUE, y = TRUE, surv = TRUE, time.inc = t)
  set.seed(PARAM$seed)
  cal <- rms::calibrate(fit_t, cmethod = "KM", method = "boot", u = t, m = floor(nrow(dn) / 4), B = 1000)
  cal_rows[[k]] <- data.frame(year = t, as.data.frame(unclass(cal)[, , drop = FALSE]), check.names = FALSE)
  plot(cal, lwd = 2, lty = 1, errbar.col = cols[k], col = cols[k], add = k > 1,
       xlim = c(0, 1), ylim = c(0, 1), subtitles = FALSE,
       xlab = "Nomogram-predicted OS", ylab = "Observed OS (Kaplan-Meier)")
}
graphics::abline(0, 1, lty = 2, col = "grey50")
graphics::legend("bottomright", legend = paste0(PARAM$surv_years, "-year"), col = cols, lwd = 2, bty = "n")
grDevices::dev.off()
utils::write.csv(do.call(rbind, cal_rows), file.path(OUT, "calibration_bootstrap.csv"), row.names = FALSE)

# =============================================================================
# 3. C-index、时间依赖 ROC、DCA：列线图 vs RLRS vs 分期
# =============================================================================
fit_rlrs  <- survival::coxph(survival::Surv(time_y, OS) ~ RLRS, data = dn)
fit_stage <- survival::coxph(survival::Surv(time_y, OS) ~ stage_bin, data = dn)
fit_nomo  <- survival::coxph(f_nomo, data = dn)
lp <- data.frame(Nomogram = stats::predict(fit_nomo, type = "lp"), RLRS = stats::predict(fit_rlrs, type = "lp"),
                 Stage = stats::predict(fit_stage, type = "lp"))
cidx <- vapply(lp, function(x) survival::concordance(survival::Surv(dn$time_y, dn$OS) ~ x, reverse = TRUE)$concordance,
               numeric(1))
print(round(cidx, 3))
utils::write.csv(data.frame(model = names(cidx), Cindex = cidx), file.path(OUT, "Cindex_nomogram_RLRS_stage.csv"),
                 row.names = FALSE)

grDevices::pdf(file.path(OUT, "Fig7_timeROC_compare.pdf"), 5, 5)
yrs <- PARAM$surv_years[PARAM$surv_years < max(dn$time_y[dn$OS == 1])]
aucs <- sapply(lp, function(x) timeROC::timeROC(T = dn$time_y, delta = dn$OS, marker = x, cause = 1,
                                                times = yrs, iid = FALSE)$AUC)
graphics::matplot(yrs, aucs, type = "b", pch = 16, lty = 1, col = c("#D7301F", "#2C7FB8", "grey40"),
                  ylim = c(0.5, 0.9), xlab = "Years", ylab = "Time-dependent AUC")
graphics::legend("topright", legend = colnames(aucs), col = c("#D7301F", "#2C7FB8", "grey40"), lwd = 2, bty = "n")
grDevices::dev.off()
utils::write.csv(data.frame(year = yrs, matrix(aucs, nrow = length(yrs), dimnames = list(NULL, colnames(lp)))),
                 file.path(OUT, "timeROC_nomogram_RLRS_stage.csv"), row.names = FALSE)

risk_at <- function(fit, t) {          # 每个患者 t 年时的预测死亡风险
  sf <- survival::survfit(fit, newdata = dn)
  1 - as.numeric(summary(sf, times = t, extend = TRUE)$surv)
}
dca_obj <- lapply(c(1, 3, 5)[c(1, 3, 5) < max(dn$time_y)], function(t) {
  dd_t <- data.frame(time_y = dn$time_y, OS = dn$OS, Nomogram = risk_at(fit_nomo, t),
                     RLRS = risk_at(fit_rlrs, t), Stage = risk_at(fit_stage, t))
  list(t = t, dca = dcurves::dca(Surv(time_y, OS) ~ Nomogram + RLRS + Stage, data = dd_t, time = t,
                                 thresholds = seq(0, 0.7, 0.01)))
})
dca_list <- lapply(dca_obj, function(o) plot(o$dca, smooth = TRUE) + ggplot2::ggtitle(paste0(o$t, "-year DCA")))
save_print(dca_list, file.path(OUT, "Fig7_DCA.pdf"), 6, 4.5)
utils::write.csv(do.call(rbind, lapply(dca_obj, function(o) data.frame(year = o$t, as.data.frame(o$dca$dca)))),
                 file.path(OUT, "DCA_net_benefit.csv"), row.names = FALSE)

# =============================================================================
# 4. 亚组 KM（RLRS 高/低，按 TCGA 全队列中位数分组）
# =============================================================================
sub_km <- list()
for (v in c("stage_bin", "age_group", "gender", "grade")) for (lv in levels(d[[v]])) {
  s <- d[!is.na(d[[v]]) & d[[v]] == lv, ]
  if (nrow(s) < 30 || length(unique(s$group)) < 2) next
  sub_km[[paste(v, lv)]] <- km_plot(s, title = paste0(v, ": ", lv, " (n = ", nrow(s), ")"))
}
save_print(sub_km, file.path(OUT, "Fig7_subgroup_KM.pdf"), 6, 6)

save_session(STEP)
log_msg("07 完成")

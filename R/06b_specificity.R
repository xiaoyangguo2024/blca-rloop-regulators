# =============================================================================
# 06b_specificity.R —— RLRS 是 R-loop 特异的信号，还是只是增殖/分期的另一种标签？
#   1. 随机基因集基准（Venet et al. 2011 的思路）：基因数与 RLRS 相同，用同一种简单模型
#      （岭回归 Cox：TCGA 训练、4 个外部队列检验）比较 RLRS 基因与 1000 组随机基因的外部 C-index
#      零分布 A：全部表达基因中随机抽取
#      零分布 B：从 TCGA 单因素 Cox P < 0.05 的基因中随机抽取（与 RLRS 基因经历的预筛选条件对等）
#   2. 与增殖的关系：RLRS 与增殖评分（E2F + G2M 的 ssGSEA）的相关；
#      同时纳入增殖评分（及分期）的 Cox；似然比检验 RLRS 在增殖评分之外是否还有信息
# 输出：results/06b_specificity/；data/processed/random_signature_null.rds
# 耗时：约 15–30 分钟（结果缓存）
# =============================================================================
source("R/00_config.R")
need(c("survival", "glmnet", "ggplot2"))
STEP <- "06b_specificity"; OUT <- res_dir(STEP)
TRAIN <- "TCGA_BLCA"; VALID <- c("GSE13507", "GSE32894", "GSE31684", "GSE48075")
OSCOH <- c(TRAIN, VALID, "IMvigor210", "GSE176307")    # 有生存数据的队列
N_RAND <- 1000

mdl    <- readRDS(file.path(DIR_PRO, "RLRS_model.rds"))
rlrs   <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
scores <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
load_cohort <- function(n) {
  o <- readRDS(file.path(DIR_PRO, paste0(n, ".rds")))
  if (n == TRAIN) o$clin <- o$clin[o$clin$group == "Tumor", ]
  o
}
coh <- lapply(stats::setNames(OSCOH, OSCOH), load_cohort)
sig <- mdl$genes; k <- length(sig)
log_msg("RLRS 基因数 k = ", k, "：", paste(sig, collapse = ", "))

# =============================================================================
# 1. 随机基因集基准
# =============================================================================
cidx <- function(lp, d) survival::concordance(survival::Surv(OS.time, OS) ~ lp,
                                              data = data.frame(d, lp = as.numeric(lp)), reverse = TRUE)$concordance
ridge_cindex <- function(genes, seed) {
  inp <- lapply(coh[c(TRAIN, VALID)], function(o) build_mime_input(o$expr, o$clin, genes))
  x <- as.matrix(inp[[TRAIN]][, -(1:3)]); y <- survival::Surv(inp[[TRAIN]]$OS.time, inp[[TRAIN]]$OS)
  set.seed(seed)
  cv <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 0, nfolds = 10)
  b <- as.matrix(stats::coef(cv, s = "lambda.min"))[, 1]
  vapply(VALID, function(n) cidx(as.matrix(inp[[n]][, names(b), drop = FALSE]) %*% b, inp[[n]]), numeric(1))
}

# 基因全集：5 个队列都有、且 TCGA 中有表达（平均 log2(TPM+1) ≥ 1）
univ <- Reduce(intersect, lapply(coh[c(TRAIN, VALID)], function(o) rownames(o$expr)))
univ <- univ[rowMeans(coh[[TRAIN]]$expr[univ, ]) >= 1]
log_msg("随机抽样的基因全集：", length(univ), " 个")

f_uc <- file.path(DIR_PRO, "TCGA_unicox_allgenes.rds")
if (!file.exists(f_uc)) {
  tr <- coh[[TRAIN]]; ok <- !is.na(tr$clin$OS.time) & tr$clin$OS.time > 0
  ex <- zscore_rows(tr$expr[univ, tr$clin$ID[ok]])
  sv <- survival::Surv(tr$clin$OS.time[ok], tr$clin$OS[ok])
  p_uc <- vapply(univ, function(g) summary(survival::coxph(sv ~ ex[g, ]))$coefficients[1, 5], numeric(1))
  saveRDS(p_uc, f_uc)
}
p_uc <- readRDS(f_uc)
pool_B <- names(p_uc)[p_uc < PARAM$unicox_p]
log_msg("零分布 B 的基因池（TCGA 单因素 Cox P < ", PARAM$unicox_p, "）：", length(pool_B), " 个")

f_null <- file.path(DIR_PRO, "random_signature_null.rds")
if (!file.exists(f_null)) {
  set.seed(PARAM$seed)
  setsA <- replicate(N_RAND, sample(univ, k), simplify = FALSE)
  setsB <- replicate(N_RAND, sample(pool_B, k), simplify = FALSE)
  run <- function(sets, tag) {
    out <- matrix(NA_real_, length(sets), length(VALID), dimnames = list(NULL, VALID))
    for (i in seq_along(sets)) {
      out[i, ] <- tryCatch(ridge_cindex(sets[[i]], PARAM$seed + i), error = function(e) NA_real_)
      if (i %% 100 == 0) log_msg(tag, "：", i, "/", length(sets))
    }
    out
  }
  obs <- ridge_cindex(sig, PARAM$seed)
  nullA <- run(setsA, "零分布 A"); nullB <- run(setsB, "零分布 B")
  saveRDS(list(setsA = setsA, setsB = setsB, nullA = nullA, nullB = nullB, obs = obs), f_null)
}
nl <- readRDS(f_null)

# 实际的 RLRS（Mime1 最优模型）在外部队列的 C-index
ci_best <- unlist(mdl$cindex[mdl$cindex$Model == mdl$model, VALID])
summ <- function(null, label) {
  m <- rowMeans(null); m <- m[!is.na(m)]
  data.frame(null = label, n_valid_draws = length(m), null_mean = mean(m), null_sd = stats::sd(m),
             null_q95 = stats::quantile(m, 0.95),
             RLRS_ridge_mean = mean(nl$obs), p_RLRS_ridge = (1 + sum(m >= mean(nl$obs))) / (length(m) + 1),
             RLRS_best_mean = mean(ci_best), p_RLRS_best = (1 + sum(m >= mean(ci_best))) / (length(m) + 1),
             row.names = NULL)
}
bench <- rbind(summ(nl$nullA, "A: random expressed genes"), summ(nl$nullB, "B: random TCGA-prognostic genes"))
print(bench, row.names = FALSE)
utils::write.csv(bench, file.path(OUT, "random_signature_benchmark.csv"), row.names = FALSE)
per_coh <- do.call(rbind, lapply(VALID, function(n) data.frame(
  cohort = n, RLRS_ridge = nl$obs[n], RLRS_best = ci_best[n],
  pctl_A = mean(nl$nullA[, n] < nl$obs[n], na.rm = TRUE), pctl_B = mean(nl$nullB[, n] < nl$obs[n], na.rm = TRUE))))
utils::write.csv(per_coh, file.path(OUT, "random_signature_benchmark_by_cohort.csv"), row.names = FALSE)

dd <- rbind(data.frame(null = "A: random expressed genes", C = rowMeans(nl$nullA)),
            data.frame(null = "B: random TCGA-prognostic genes", C = rowMeans(nl$nullB)))
p_null <- ggplot2::ggplot(dd[!is.na(dd$C), ], ggplot2::aes(C)) +
  ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = "grey40", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = mean(nl$obs), colour = "#D7301F", linewidth = 0.8) +
  ggplot2::geom_vline(xintercept = mean(ci_best), colour = "#2C7FB8", linewidth = 0.8, linetype = 2) +
  ggplot2::facet_wrap(~ null, ncol = 1, scales = "free_y") +
  ggplot2::labs(x = sprintf("Mean external C-index (%d cohorts), %d-gene ridge Cox", length(VALID), k), y = "Random gene sets",
                subtitle = "Red: RLRS genes refitted with ridge Cox; blue dashed: final RLRS model") +
  ggplot2::theme_classic(base_size = 10)
save_gg(p_null, file.path(OUT, "Fig_random_signature_benchmark.pdf"), 5.5, 5)

# =============================================================================
# 2. 增殖评分与分期校正
# =============================================================================
adv_stage <- function(n, c0) {     # 1 = 较晚期：TCGA 为 III–IV 期；混合队列为肌层浸润（T2–4）；MIBC 队列为 T3–4
  if (n == TRAIN) return(ifelse(c0$stage %in% c("III", "IV"), 1L, ifelse(c0$stage %in% c("I", "II"), 0L, NA)))
  if (!"T" %in% colnames(c0)) return(rep(NA_integer_, nrow(c0)))
  hi <- if (n %in% c("GSE13507", "GSE32894")) c("2", "3", "4") else c("3", "4")
  ifelse(c0$T %in% hi, 1L, ifelse(is.na(c0$T), NA, 0L))
}
rows <- list()
for (n in OSCOH) {
  c0 <- coh[[n]]$clin
  d <- merge(data.frame(ID = c0$ID, OS.time = c0$OS.time, OS = c0$OS, adv = adv_stage(n, c0)),
             rlrs[rlrs$cohort == n, c("ID", "RLRS")], by = "ID")
  d <- merge(d, scores[[n]][, c("ID", "Prolif", "RL_core")], by = "ID")
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
  d$RLRS_z <- as.numeric(scale(d$RLRS)); d$Prolif_z <- as.numeric(scale(d$Prolif))
  S <- survival::Surv(d$OS.time, d$OS)
  m0 <- survival::coxph(S ~ Prolif_z, data = d)
  m1 <- survival::coxph(S ~ Prolif_z + RLRS_z, data = d)
  mr <- survival::coxph(S ~ RLRS_z, data = d)
  r <- data.frame(cohort = n, n = nrow(d), events = sum(d$OS),
                  rho_RLRS_Prolif = stats::cor(d$RLRS, d$Prolif, method = "spearman"),
                  rho_RLRS_RLcore = stats::cor(d$RLRS, d$RL_core, method = "spearman"),
                  C_Prolif = summary(m0)$concordance[1], C_RLRS = summary(mr)$concordance[1],
                  C_Prolif_RLRS = summary(m1)$concordance[1])
  a1 <- summary(m1)$coefficients
  r$HR_RLRS_adjProlif <- exp(a1["RLRS_z", 1]); r$p_RLRS_adjProlif <- a1["RLRS_z", 5]
  r$LRT_p_add_RLRS <- stats::anova(m0, m1)[2, "Pr(>|Chi|)"]
  if (sum(!is.na(d$adv)) > 30 && length(unique(stats::na.omit(d$adv))) == 2) {
    m2 <- survival::coxph(S ~ RLRS_z + Prolif_z + adv, data = d)
    a2 <- summary(m2)$coefficients
    r$HR_RLRS_adjProlifStage <- exp(a2["RLRS_z", 1]); r$p_RLRS_adjProlifStage <- a2["RLRS_z", 5]
  } else { r$HR_RLRS_adjProlifStage <- NA; r$p_RLRS_adjProlifStage <- NA }
  rows[[n]] <- r
}
adj <- do.call(rbind, rows); print(adj, row.names = FALSE, digits = 3)
utils::write.csv(adj, file.path(OUT, "Table_S6_proliferation_adjusted.csv"), row.names = FALSE)

save_session(STEP)
log_msg("06b 完成")

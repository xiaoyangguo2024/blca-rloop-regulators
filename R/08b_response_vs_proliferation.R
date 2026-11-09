# =============================================================================
# 08b_response_vs_proliferation.R —— 治疗反应与增殖评分（对应 06b 对预后做的特异性检验）
#   免疫治疗（IMvigor210、GSE176307）：客观缓解（CR/PR vs SD/PD）
#   新辅助化疗（GSE169455，新辅助 125 例）：pCR、降期
#   1. 增殖评分的 AUC，方向固定（评分越高 → 越可能缓解），与 08 步对 RLRS 的做法相同
#   2. Logistic 回归：RLRS、增殖评分各自单独，以及 RLRS 校正增殖评分；似然比检验 RLRS 是否在增殖之外增加信息
#   3. IMvigor210：RLRS 校正 TMB（log2(TMB + 1)）
#   RLRS_z 用 RLRS_all_cohorts.rds 里队列内标准化的值（与 08 步相同），增殖评分同样在整个队列内标准化
# 输出：results/08_tme_therapy/Response_vs_proliferation.csv
# =============================================================================
source("R/00_config.R")
need(c("pROC"))
OUT <- res_dir("08_tme_therapy")

rlrs   <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
scores <- readRDS(file.path(DIR_PRO, "ssgsea_scores_all_cohorts.rds"))
log_msg("ssGSEA 评分队列：", paste(names(scores), collapse = ", "))

add_scores <- function(cl, n) {
  s <- scores[[n]][, c("ID", "Prolif")]
  s$Prolif_z <- as.numeric(scale(s$Prolif))
  e <- merge(cl, rlrs[rlrs$cohort == n, c("ID", "RLRS", "RLRS_z")], by = "ID")
  merge(e, s, by = "ID")
}
auc_fixed <- function(y, s) {
  roc <- pROC::roc(y, s, levels = c("No", "Yes"), direction = "<", quiet = TRUE)
  as.numeric(pROC::ci.auc(roc))
}
one <- function(cohort, endpoint, x, tmb = FALSE) {
  a_r <- auc_fixed(x$outcome, x$RLRS); a_p <- auc_fixed(x$outcome, x$Prolif)
  g0 <- stats::glm(outcome ~ Prolif_z, data = x, family = stats::binomial())
  g1 <- stats::glm(outcome ~ Prolif_z + RLRS_z, data = x, family = stats::binomial())
  gr <- stats::glm(outcome ~ RLRS_z, data = x, family = stats::binomial())
  cp <- summary(g0)$coefficients; c1 <- summary(g1)$coefficients; cr <- summary(gr)$coefficients
  r <- data.frame(cohort = cohort, endpoint = endpoint, n = nrow(x), n_yes = sum(x$outcome == "Yes"),
                  rho_RLRS_Prolif = stats::cor(x$RLRS, x$Prolif, method = "spearman"),
                  AUC_RLRS = a_r[2], AUC_RLRS_lower = a_r[1], AUC_RLRS_upper = a_r[3],
                  AUC_Prolif = a_p[2], AUC_Prolif_lower = a_p[1], AUC_Prolif_upper = a_p[3],
                  p_wilcox_Prolif = stats::wilcox.test(Prolif ~ outcome, x)$p.value,
                  OR_RLRS = exp(cr["RLRS_z", 1]), p_RLRS = cr["RLRS_z", 4],
                  OR_Prolif = exp(cp["Prolif_z", 1]), p_Prolif = cp["Prolif_z", 4],
                  OR_RLRS_adjProlif = exp(c1["RLRS_z", 1]), p_RLRS_adjProlif = c1["RLRS_z", 4],
                  OR_Prolif_adjRLRS = exp(c1["Prolif_z", 1]), p_Prolif_adjRLRS = c1["Prolif_z", 4],
                  LRT_p_add_RLRS = stats::anova(g0, g1, test = "LRT")[2, "Pr(>Chi)"],
                  n_TMB = NA, OR_RLRS_adjTMB = NA, p_RLRS_adjTMB = NA)
  if (tmb) {
    t <- x[!is.na(x$TMB), ]
    t$lTMB <- log2(t$TMB + 1)
    ct <- summary(stats::glm(outcome ~ RLRS_z + lTMB, data = t, family = stats::binomial()))$coefficients
    r$n_TMB <- nrow(t); r$OR_RLRS_adjTMB <- exp(ct["RLRS_z", 1]); r$p_RLRS_adjTMB <- ct["RLRS_z", 4]
  }
  r
}

rows <- list()
for (n in c("IMvigor210", "GSE176307")) {
  e <- add_scores(readRDS(file.path(DIR_PRO, paste0(n, ".rds")))$clin, n)
  x <- e[e$response %in% c("CR/PR", "SD/PD"), ]
  x$outcome <- factor(ifelse(x$response == "CR/PR", "Yes", "No"), levels = c("No", "Yes"))
  rows[[n]] <- one(n, "objective response", x, tmb = "TMB" %in% colnames(x) && sum(!is.na(x$TMB)) > 40)
}
e <- add_scores(readRDS(file.path(DIR_PRO, "GSE169455.rds"))$clin, "GSE169455")
e <- e[e$setting == "neoadjuvant", ]
for (ep in c("pCR", "downstage")) {
  x <- e[!is.na(e[[ep]]), ]
  x$outcome <- factor(ifelse(x[[ep]] == 1, "Yes", "No"), levels = c("No", "Yes"))
  rows[[ep]] <- one("GSE169455", ep, x)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab, digits = 3)
utils::write.csv(tab, file.path(OUT, "Response_vs_proliferation.csv"), row.names = FALSE)

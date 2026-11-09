# =============================================================================
# 07b_nomogram_validation.R —— 列线图 C-index 的 bootstrap 乐观校正（rms::validate，1000 次）
#   数据和模型与 07 步相同（RLRS + 年龄 + 分期，TCGA-BLCA 完整病例）；07 步只报告了表观 C-index。
# 输出：results/07_clinical_nomogram/Nomogram_optimism_corrected.csv
# =============================================================================
source("R/00_config.R")
need(c("survival", "rms"))
OUT <- res_dir("07_clinical_nomogram")

tum  <- readRDS(file.path(DIR_PRO, "TCGA_tumor_annot.rds"))
rlrs <- readRDS(file.path(DIR_PRO, "RLRS_all_cohorts.rds"))
d <- merge(tum, rlrs[rlrs$cohort == "TCGA_BLCA", c("ID", "RLRS")], by = "ID")
d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS), ]
d$time_y    <- d$OS.time / 365.25
d$stage_bin <- factor(d$stage_bin, levels = c("I-II", "III-IV"))
dn <- d[stats::complete.cases(d[, c("time_y", "OS", "RLRS", "age", "stage_bin")]), ]
dd <- rms::datadist(dn); options(datadist = "dd")
fit <- rms::cph(Surv(time_y, OS) ~ RLRS + age + stage_bin, data = dn, x = TRUE, y = TRUE, surv = TRUE)

set.seed(PARAM$seed)
v <- rms::validate(fit, method = "boot", B = 1000)
print(v)
dxy <- v["Dxy", ]
# rms 对 cph 报告的 Dxy 符号取决于版本（线性预测值与生存时间负相关）；C = 0.5 + |Dxy| / 2
c_app <- 0.5 + abs(dxy[["index.orig"]]) / 2
c_cor <- 0.5 + abs(dxy[["index.corrected"]]) / 2
res <- data.frame(n = nrow(dn), events = sum(dn$OS), B = 1000, C_apparent = c_app,
                  C_optimism = c_app - c_cor, C_corrected = c_cor)
print(res)
utils::write.csv(res, file.path(OUT, "Nomogram_optimism_corrected.csv"), row.names = FALSE)

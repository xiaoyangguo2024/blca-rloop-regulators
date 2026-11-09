# =============================================================================
# 03c_copykat_normal_control.R —— copyKAT 阴性对照：癌旁/正常样本的上皮细胞会不会也被判为非整倍体？
#   与 03 步完全相同的设置（免疫细胞作二倍体参照，同一数据集补足到 200 个），对象换成正常样本的上皮细胞。
#   若正常上皮大多被判为二倍体，说明 03 步的非整倍体判定不是"上皮 vs 免疫细胞"表达差异造成的假象。
#   结果只用于评估 copyKAT 判定的特异性，不改变 03 步的恶性细胞定义。
# 输入：data/processed/scRNA_integrated.rds（03 步）
# 输出：results/03_scRNA/copykat_normal_control.csv；copyKAT 预测文件放在 results/03_scRNA/copykat_control/
# 内存：与 03 步相同，先取出计数子矩阵再释放大对象，一次只跑一个样本
# =============================================================================
source("R/00_config.R")
need(c("Seurat", "Matrix"))
OUT <- res_dir("03_scRNA")
REF_TYPES <- c("T_cell", "NK", "B_cell", "Plasma", "Myeloid")
N_REF <- 200
ck_dir <- file.path(OUT, "copykat_control"); dir.create(ck_dir, showWarnings = FALSE)
ck_file <- function(smp) file.path(ck_dir, paste0(smp, "_copykat_prediction.txt"))

sc <- readRDS(file.path(DIR_PRO, "scRNA_integrated.rds"))
normals <- unique(sc$sample[sc$tissue == "Normal"])
plan <- list()
for (smp in normals) {
  epi <- colnames(sc)[sc$sample == smp & sc$celltype == "Epithelial"]
  if (length(epi) < 50) { log_msg(smp, "：上皮细胞 ", length(epi), " 个，不足 50，跳过"); next }
  own <- colnames(sc)[sc$sample == smp & sc$celltype %in% REF_TYPES]
  set.seed(PARAM$seed)
  ref <- if (length(own) >= N_REF) sample(own, N_REF) else {
    ds <- sc$dataset[match(smp, sc$sample)]
    pool <- setdiff(colnames(sc)[sc$dataset == ds & sc$celltype %in% REF_TYPES], own)
    c(own, sample(pool, min(length(pool), N_REF - length(own))))
  }
  plan[[smp]] <- list(epi = epi, ref = ref, n_own = length(own))
}
todo <- names(plan)[!file.exists(vapply(names(plan), ck_file, ""))]
if (length(todo)) {
  cnt <- LayerData(sc, assay = "RNA", layer = "counts")
  f_in <- stats::setNames(file.path(ck_dir, paste0(todo, "_input_counts.rds")), todo)
  for (smp in todo) {
    m <- cnt[, c(plan[[smp]]$epi, plan[[smp]]$ref)]
    saveRDS(m[Matrix::rowSums(m) > 0, ], f_in[[smp]])
  }
  rm(cnt, m)
}
rm(sc); invisible(gc())

if (length(todo)) {
  need("copykat")
  for (smp in todo) {
    raw <- as.matrix(readRDS(f_in[[smp]]))
    colnames(raw) <- make.names(colnames(raw))
    ref <- plan[[smp]]$ref
    log_msg("copyKAT 阴性对照：", smp, "（", length(plan[[smp]]$epi), " 个正常上皮 + ", length(ref), " 个参照细胞）")
    old <- setwd(ck_dir)
    ok <- tryCatch({
      invisible(copykat::copykat(rawmat = raw, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1, sam.name = smp,
                                 distance = "euclidean", norm.cell.names = make.names(ref),
                                 output.seg = "FALSE", plot.genes = "FALSE", genome = "hg20",
                                 n.cores = PARAM$copykat_cores))
      TRUE
    }, error = function(e) { log_msg("copyKAT 失败：", smp, "：", conditionMessage(e)); FALSE })
    setwd(old)
    rm(raw); invisible(gc())
    if (ok && file.exists(ck_file(smp))) {
      unlink(f_in[[smp]])
      unlink(list.files(ck_dir, pattern = paste0("^", smp, "_copykat_.*(gene_by_cell|bin_by_cell)\\.txt$"), full.names = TRUE))
    }
  }
}

res <- do.call(rbind, lapply(names(plan), function(smp) {
  if (!file.exists(ck_file(smp))) return(NULL)
  pr <- utils::read.delim(ck_file(smp)); epi <- make.names(plan[[smp]]$epi)
  p <- pr$copykat.pred[pr$cell.names %in% epi]
  data.frame(sample = smp, n_epithelial = length(epi), n_reference = length(plan[[smp]]$ref),
             n_aneuploid = sum(p == "aneuploid"), n_diploid = sum(p == "diploid"), n_undefined = sum(!p %in% c("aneuploid", "diploid")),
             aneuploid_frac = round(mean(p == "aneuploid"), 3))
}))
print(res, row.names = FALSE)
utils::write.csv(res, file.path(OUT, "copykat_normal_control.csv"), row.names = FALSE)
save_session("03c_copykat_normal_control")
log_msg("03c 完成")

# =============================================================================
# 03_scRNA.R —— 单细胞：质控 → Harmony 整合 → 注释 → copyKAT 识别恶性细胞 → UCell R-loop 评分
#               → 样本内中位数分组 → 伪 bulk 差异分析（edgeR）→ Hallmark → CellChat
# 数据（2026-11-07 核对，均为原始 UMI 计数的稠密表）：
#   GSE135337：BC1–BC7 肿瘤（Ta/T1/T1/T2/T3/Ta/Ta）+ BCN 癌旁
#   GSE222315：P1–P9 肿瘤（均为 T1）+ P2–P5 配对癌旁
#   注意：单细胞样本以非肌层浸润性（NMIBC）为主，MIBC 的空间分布由 09 步的 GSE319536（22 张 MIBC 切片）补充
# 输出：data/processed/scRNA_integrated.rds、scRNA_malignant.rds；results/03_scRNA/
# 耗时：copyKAT 每个样本 0.5–3 小时；结果会缓存，重跑时自动跳过已完成的样本
# =============================================================================
source("R/00_config.R")
need(c("Seurat", "harmony", "UCell", "edgeR", "data.table", "Matrix", "ggplot2", "patchwork", "ggpubr"))
STEP <- "03_scRNA"; OUT <- res_dir(STEP)
GS <- readRDS(file.path(DIR_GS, "Rloop_genesets.rds"))

MALIGNANT_RULE <- "all_tumor_epi"  # "copykat"：肿瘤样本中非整倍体的上皮细胞；"all_tumor_epi"：肿瘤样本中全部上皮细胞
# 2026-11-08 改为 all_tumor_epi：03c 阴性对照中正常样本的上皮细胞也有 77%–98% 被 copyKAT 判为非整倍体
# （BCN 0.983、P4-NAT 0.879、P5-NAT 0.771），而 P2、P9 两个肿瘤样本的上皮细胞几乎全判为二倍体，
# 说明以免疫细胞为参照时 copyKAT 的判定主要反映上皮与免疫细胞的表达差异，不能区分恶性细胞。
# copyKAT 结果仍计算并写入 copykat_by_sample.csv（补充表），按 copyKAT 定义的旧结果另存在 sensitivity_copykat_rule/
USE_GSE190888  <- FALSE          # 可选验证集（3 例肿瘤；LH0826 为腺性膀胱炎，不纳入）

# ---- 样本表（样本名用 "-" 不用 "_"，因为 Seurat 伪 bulk 会把 "_" 当分隔符）------
SC_SAMPLES <- rbind(
  data.frame(dataset = "GSE135337",
             gsm = c("GSM4006644", "GSM4006645", "GSM4006646", "GSM4006647", "GSM4006648", "GSM4751267", "GSM4751268", "GSM5329919"),
             sample = c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7", "BCN"),
             tissue = c(rep("Tumor", 7), "Normal"),
             T_stage = c("Ta", "T1", "T1", "T2", "T3", "Ta", "Ta", NA)),
  data.frame(dataset = "GSE222315",
             gsm = c("GSM6919778", "GSM6919779", "GSM6919780", "GSM6919781", "GSM6919782", "GSM6919783", "GSM6919784",
                     "GSM6919785", "GSM6919786", "GSM7898117", "GSM7898118", "GSM7898119", "GSM7898120"),
             sample = c("P1-BCa", "P2-NAT", "P2-BCa", "P3-NAT", "P3-BCa", "P4-NAT", "P4-BCa", "P5-NAT", "P5-BCa",
                        "P6-BCa", "P7-BCa", "P8-BCa", "P9-BCa"),
             tissue = ifelse(grepl("NAT", c("P1-BCa", "P2-NAT", "P2-BCa", "P3-NAT", "P3-BCa", "P4-NAT", "P4-BCa", "P5-NAT",
                                            "P5-BCa", "P6-BCa", "P7-BCa", "P8-BCa", "P9-BCa")), "Normal", "Tumor"),
             T_stage = "T1")
)
SC_SAMPLES$file <- with(SC_SAMPLES, ifelse(
  dataset == "GSE135337",
  paste0(gsm, "_", sample, "_gene_cell_exprs_table.", ifelse(sample == "BCN", "xls", "txt"), ".gz"),
  paste0(gsm, "_", sub("-", "_", sub("^P", "p", sample)), "_expression.",
         ifelse(gsm %in% c("GSM7898117", "GSM7898118", "GSM7898119", "GSM7898120"), "xls", "txt"), ".gz")))
if (USE_GSE190888) SC_SAMPLES <- rbind(SC_SAMPLES, data.frame(
  dataset = "GSE190888", gsm = c("GSM5733425", "GSM5733426", "GSM5733427"),
  sample = c("LHBC0803", "LHBC0911", "LHBC0917"), tissue = "Tumor", T_stage = NA,
  file = c("GSM5733425_LHBC0803-expression_matrix.csv.gz", "GSM5733426_LHBC0911-expression_matrix.csv.gz",
           "GSM5733427_LHBC0917-expression_matrix.csv.gz")))
utils::write.csv(SC_SAMPLES, file.path(OUT, "scRNA_samples.csv"), row.names = FALSE)

# ---- 读稠密计数表：自动识别基因符号列，同名基因计数相加 --------------------------
read_dense <- function(f, max_cells = Inf, seed = PARAM$seed) {
  # 先读前 50 行判断哪些列是细胞（数值列），超过上限时只读随机抽中的列，避免整表读入占满内存
  d0 <- data.table::fread(f, nrows = 50, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  num <- which(!vapply(d0, is.character, logical(1)))
  if (length(num) > max_cells) { set.seed(seed); num <- sort(sample(num, max_cells)) }
  keep <- unname(sort(c(which(vapply(d0, is.character, logical(1))), num))); rm(d0)   # 带名字的整数向量会被 fread 当成列类型
  d <- data.table::fread(f, select = keep, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  is_chr <- vapply(d, is.character, logical(1))
  chr <- which(is_chr)
  sym_j <- chr[vapply(chr, function(j) mean(grepl("^ENSG", utils::head(d[[j]], 200))) < 0.5, logical(1))][1]
  if (is.na(sym_j)) sym_j <- chr[1]
  m <- as.matrix(d[, !is_chr, drop = FALSE]); rownames(m) <- NULL
  m <- rowsum(m, d[[sym_j]], reorder = FALSE)
  chk <- m[, seq_len(min(200, ncol(m))), drop = FALSE]
  if (any(chk %% 1 != 0)) log_msg("【核对】", basename(f), " 含非整数，可能不是原始计数")
  m
}

# =============================================================================
# 1. 读入、质控、整合、聚类
# =============================================================================
f_base <- file.path(DIR_PRO, "scRNA_base.rds")
if (!file.exists(f_base)) {
  objs <- lapply(seq_len(nrow(SC_SAMPLES)), function(i) {
    s <- SC_SAMPLES[i, ]
    m <- read_dense(gsm_suppl(s$gsm, s$file, subdir = s$dataset), max_cells = PARAM$sc_max_cells_per_sample)
    if (ncol(m) > PARAM$sc_max_cells_per_sample) {
      set.seed(PARAM$seed); m <- m[, sort(sample(ncol(m), PARAM$sc_max_cells_per_sample))]
    }
    colnames(m) <- paste0(s$sample, "_", make.unique(sub("[.-][0-9]+$", "", colnames(m))))
    o <- CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE), project = s$sample,
                            min.cells = 3, min.features = PARAM$sc_min_features)
    o$sample <- s$sample; o$dataset <- s$dataset; o$tissue <- s$tissue; o$T_stage <- s$T_stage
    log_msg(s$sample, "：", ncol(o), " 个细胞")
    o
  })
  sc <- merge(objs[[1]], objs[-1]); rm(objs); invisible(gc())
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  save_gg(VlnPlot(sc, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", pt.size = 0, ncol = 1),
          file.path(OUT, "QC_before_filter.pdf"), 12, 10)
  sc <- subset(sc, subset = nFeature_RNA > PARAM$sc_min_features & nFeature_RNA < PARAM$sc_max_features &
                 percent.mt < PARAM$sc_max_mt)
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc, nfeatures = 2000)
  sc <- ScaleData(sc)
  sc <- RunPCA(sc, npcs = PARAM$sc_npcs, verbose = FALSE)
  sc <- harmony::RunHarmony(sc, group.by.vars = "sample")
  sc <- RunUMAP(sc, reduction = "harmony", dims = 1:PARAM$sc_npcs)
  sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:PARAM$sc_npcs)
  sc <- FindClusters(sc, resolution = PARAM$sc_resolution)
  sc <- JoinLayers(sc)
  saveRDS(sc, f_base)
}
sc <- readRDS(f_base)
# 之后的步骤都不用 scale.data（约 0.5 GB），删掉以节省内存
if ("scale.data" %in% SeuratObject::Layers(sc[["RNA"]])) sc[["RNA"]]$scale.data <- NULL
print(table(sc$sample))

# =============================================================================
# 2. 细胞类型注释：marker 打分自动给标签 → DotPlot 人工核对
#    人工修正方法：把 cluster_annotation_auto.csv 另存为 cluster_annotation_manual.csv，改 celltype 列后重跑
# =============================================================================
mk <- lapply(MARKERS, function(g) intersect(g, rownames(sc)))
mk <- mk[lengths(mk) >= 2]
sc <- AddModuleScore(sc, features = mk, name = "mk")
mk_cols <- paste0("mk", seq_along(mk))
cl_mean <- stats::aggregate(sc@meta.data[, mk_cols], by = list(cluster = sc$seurat_clusters), FUN = mean)
ann <- data.frame(cluster = cl_mean$cluster, celltype_auto = names(mk)[apply(as.matrix(cl_mean[, -1]), 1, which.max)])
f_manual <- file.path(OUT, "cluster_annotation_manual.csv")
if (file.exists(f_manual)) {
  man <- utils::read.csv(f_manual)
  ann$celltype <- man$celltype[match(ann$cluster, man$cluster)]
  log_msg("使用人工注释：", f_manual)
} else {
  ann$celltype <- ann$celltype_auto
  utils::write.csv(ann, file.path(OUT, "cluster_annotation_auto.csv"), row.names = FALSE)
  log_msg("【核对】自动注释已写入 cluster_annotation_auto.csv，请对照 DotPlot 检查")
}
sc$celltype <- ann$celltype[match(as.character(sc$seurat_clusters), as.character(ann$cluster))]
for (cc in mk_cols) sc[[cc]] <- NULL

save_gg(DotPlot(sc, features = unique(unlist(mk)), group.by = "seurat_clusters") + RotatedAxis(),
        file.path(OUT, "DotPlot_markers_by_cluster.pdf"), 18, 10)
save_gg(DotPlot(sc, features = unique(unlist(mk)), group.by = "celltype") + RotatedAxis(),
        file.path(OUT, "Fig2_DotPlot_markers_by_celltype.pdf"), 16, 6)
save_gg(DimPlot(sc, group.by = "celltype", label = TRUE, repel = TRUE) +
          DimPlot(sc, group.by = "sample") + DimPlot(sc, group.by = "tissue"),
        file.path(OUT, "Fig2_UMAP_celltype_sample_tissue.pdf"), 21, 6)

# =============================================================================
# 3. copyKAT：每个肿瘤样本单独运行，以免疫细胞作二倍体参照
#   内存（8 GB 机器）：不要从持有整个 Seurat 对象（约 2.5 GB）的进程 fork 多个 copyKAT——
#   每个子进程都会逐渐占用大量内存，曾导致整机卡死。现在的做法：
#   先把每个样本需要的计数子矩阵存成小文件，把 sc 暂存到磁盘并释放，再一次只跑一个样本，最后重新读回 sc。
#   去掉全零的基因行不改变 copyKAT 结果（它的第一步只保留在 > 5% 细胞中检出的基因）
# =============================================================================
REF_TYPES <- c("T_cell", "NK", "B_cell", "Plasma", "Myeloid")
N_REF <- 200     # 每次 copyKAT 的二倍体参照细胞数
ck_dir <- file.path(OUT, "copykat"); dir.create(ck_dir, showWarnings = FALSE)
ck_file <- function(smp) file.path(ck_dir, paste0(smp, "_copykat_prediction.txt"))
ck_plan <- list()
for (smp in SC_SAMPLES$sample[SC_SAMPLES$tissue == "Tumor"]) {
  epi <- colnames(sc)[sc$sample == smp & sc$celltype == "Epithelial"]
  if (length(epi) < 50) { log_msg(smp, "：上皮细胞不足 50 个，跳过 copyKAT"); next }
  # 二倍体参照：本样本免疫细胞优先，不足 N_REF 个时用同一数据集其他样本的免疫细胞补足。
  # 不给足参照时 copyKAT 会进入逐细胞 GMM 基线估计（实测约 2–3 分钟/细胞，1500 个细胞跑不完）
  own <- colnames(sc)[sc$sample == smp & sc$celltype %in% REF_TYPES]
  set.seed(PARAM$seed)
  ref <- if (length(own) >= N_REF) sample(own, N_REF) else {
    ds <- sc$dataset[match(smp, sc$sample)]
    pool <- setdiff(colnames(sc)[sc$dataset == ds & sc$celltype %in% REF_TYPES], own)
    c(own, sample(pool, min(length(pool), N_REF - length(own))))
  }
  ck_plan[[smp]] <- list(epi = epi, ref = ref, n_own = length(own))
}
utils::write.csv(data.frame(sample = names(ck_plan), n_epithelial = vapply(ck_plan, function(p) length(p$epi), 1L),
                            n_reference = vapply(ck_plan, function(p) length(p$ref), 1L),
                            n_reference_same_sample = vapply(ck_plan, function(p) p$n_own, 1L)),
                 file.path(OUT, "copykat_reference_cells.csv"), row.names = FALSE)

todo <- names(ck_plan)[!file.exists(vapply(names(ck_plan), ck_file, ""))]
if (length(todo)) {
  need("copykat")
  cnt <- LayerData(sc, assay = "RNA", layer = "counts")
  f_in <- stats::setNames(file.path(ck_dir, paste0(todo, "_input_counts.rds")), todo)
  for (smp in todo) {
    m <- cnt[, c(ck_plan[[smp]]$epi, ck_plan[[smp]]$ref)]
    saveRDS(m[Matrix::rowSums(m) > 0, ], f_in[[smp]])
  }
  rm(cnt, m)
  f_tmp <- file.path(DIR_PRO, "scRNA_annotated_tmp.rds")
  saveRDS(sc, f_tmp, compress = FALSE); rm(sc); invisible(gc())
  for (smp in todo) {
    raw <- as.matrix(readRDS(f_in[[smp]]))
    colnames(raw) <- make.names(colnames(raw))
    ref <- ck_plan[[smp]]$ref
    log_msg("copyKAT：", smp, "（", length(ck_plan[[smp]]$epi), " 个上皮 + ", length(ref), " 个参照细胞，其中本样本 ",
            ck_plan[[smp]]$n_own, " 个）")
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
      big <- list.files(ck_dir, pattern = paste0("^", smp, "_copykat_.*(gene_by_cell|bin_by_cell)\\.txt$"), full.names = TRUE)
      unlink(big)                                   # 中间文件，每个 0.15–0.35 GB
    }
  }
  sc <- readRDS(f_tmp); unlink(f_tmp)
}
ck <- do.call(rbind, lapply(names(ck_plan), function(smp) {
  if (!file.exists(ck_file(smp))) return(NULL)
  pr <- utils::read.delim(ck_file(smp)); epi <- ck_plan[[smp]]$epi
  out <- data.frame(cell = epi[match(pr$cell.names, make.names(epi))], copykat = pr$copykat.pred)   # 只保留本样本的上皮细胞
  out[!is.na(out$cell), ]
}))
sc$copykat <- if (is.null(ck)) NA else ck$copykat[match(colnames(sc), ck$cell)]

is_te <- sc$celltype == "Epithelial" & sc$tissue == "Tumor"
sc$cell_class <- sc$celltype
sc$cell_class[sc$celltype == "Epithelial" & sc$tissue == "Normal"] <- "Normal_epi"
sc$cell_class[is_te] <- if (MALIGNANT_RULE == "all_tumor_epi") "Malignant" else
  ifelse(sc$copykat[is_te] %in% "aneuploid", "Malignant",
         ifelse(sc$copykat[is_te] %in% "diploid", "Tumor_diploid_epi", "Tumor_epi_undefined"))
ck_lab <- ifelse(sc$copykat[is_te] %in% "aneuploid", "aneuploid",
                 ifelse(sc$copykat[is_te] %in% "diploid", "diploid", "not_analysed_or_undefined"))
ck_tab <- as.data.frame.matrix(table(sc$sample[is_te],
                                     factor(ck_lab, levels = c("aneuploid", "diploid", "not_analysed_or_undefined"))))
ck_tab$aneuploid_frac <- round(ck_tab$aneuploid / rowSums(ck_tab), 3)
print(ck_tab)
utils::write.csv(ck_tab, file.path(OUT, "copykat_by_sample.csv"))
if (any(ck_tab$aneuploid_frac < 0.1))
  log_msg("【核对】部分样本非整倍体比例 < 10%（Ta 期肿瘤常见近二倍体）：",
          paste(rownames(ck_tab)[ck_tab$aneuploid_frac < 0.1], collapse = ", "),
          "；可考虑 MALIGNANT_RULE <- \"all_tumor_epi\" 做敏感性分析")
epi_obj <- subset(sc, subset = celltype == "Epithelial")
save_gg(DimPlot(epi_obj, group.by = "cell_class") + DimPlot(epi_obj, group.by = "sample"),
        file.path(OUT, "Fig2_UMAP_epithelial_copykat.pdf"), 14, 6)
rm(epi_obj)

# =============================================================================
# 4. R-loop 评分（UCell，基于基因排序，受测序深度和样本组成影响小）
# =============================================================================
sc <- UCell::AddModuleScore_UCell(sc, features = GS[c("RL_core", "RL_all", "RL_promote", "RL_suppress")],
                                  maxRank = 2500, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
sc$RHI <- as.numeric(scale(sc$RL_promote_UCell) - scale(sc$RL_suppress_UCell))

save_gg(FeaturePlot(sc, c("RL_core_UCell", "RHI"), order = TRUE),
        file.path(OUT, "Fig2_UMAP_Rloop_score.pdf"), 12, 5)
save_gg(VlnPlot(sc, "RL_core_UCell", group.by = "cell_class", pt.size = 0) + NoLegend(),
        file.path(OUT, "Fig2_Vln_Rloop_by_cellclass.pdf"), 10, 5)

# 恶性 vs 正常上皮：以样本为单位比较（避免把细胞当独立样本）
md <- sc@meta.data
agg <- stats::aggregate(cbind(RL_core_UCell, RHI) ~ sample + cell_class,
                        data = md[md$cell_class %in% c("Malignant", "Tumor_diploid_epi", "Normal_epi"), ], FUN = mean)
cmp <- Filter(function(x) all(x %in% agg$cell_class),
              list(c("Malignant", "Normal_epi"), c("Malignant", "Tumor_diploid_epi")))
p <- ggpubr::ggboxplot(agg, x = "cell_class", y = "RL_core_UCell", add = "jitter", color = "cell_class") +
  ggpubr::stat_compare_means(comparisons = cmp, method = "wilcox.test") +
  ggplot2::labs(x = NULL, y = "R-loop score (per-sample mean)")
save_gg(p, file.path(OUT, "Fig2_Malignant_vs_normal_epi_per_sample.pdf"), 5, 5)
utils::write.csv(agg, file.path(OUT, "Rloop_score_per_sample_class.csv"), row.names = FALSE)

# =============================================================================
# 5. 恶性细胞：每个样本内部按 RL_core 中位数分高/低两组
# =============================================================================
mal <- subset(sc, subset = cell_class == "Malignant")
n_mal <- table(mal$sample)
keep_s <- names(n_mal)[n_mal >= 2 * PARAM$sc_min_cells_pb]
log_msg("纳入高/低分组的样本：", paste(keep_s, collapse = ", "))
mal <- subset(mal, subset = sample %in% keep_s)
grp <- character(ncol(mal))
for (s in keep_s) {
  i <- mal$sample == s
  grp[i] <- ifelse(mal$RL_core_UCell[i] >= stats::median(mal$RL_core_UCell[i]), "High", "Low")
}
mal$RL_group <- factor(grp, levels = c("Low", "High"))

# ---- 伪 bulk + edgeR（~ sample + group，配对设计）------------------------------
pb <- AggregateExpression(mal, assays = "RNA", group.by = c("sample", "RL_group"), return.seurat = FALSE)$RNA
pb <- as.matrix(pb)
pb_smp <- sub("_[^_]+$", "", colnames(pb)); pb_grp <- factor(sub("^.*_", "", colnames(pb)), levels = c("Low", "High"))
y <- edgeR::DGEList(pb, group = pb_grp)
y <- y[edgeR::filterByExpr(y, group = pb_grp), , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y)
design <- stats::model.matrix(~ factor(pb_smp) + pb_grp)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
deg_pb <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = "pb_grpHigh"), n = Inf)$table
deg_pb$gene <- rownames(deg_pb)
utils::write.csv(deg_pb, file.path(OUT, "pseudobulk_DE_High_vs_Low.csv"), row.names = FALSE)
sc_deg <- deg_pb$gene[deg_pb$FDR < PARAM$sc_deg_fdr & abs(deg_pb$logFC) > PARAM$sc_deg_lfc]
writeLines(sc_deg, file.path(OUT, "sc_DEG_genes.txt"))
log_msg("伪 bulk 差异基因：", length(sc_deg), " 个（FDR < ", PARAM$sc_deg_fdr, "，|log2FC| > ", PARAM$sc_deg_lfc, "）")

deg_pb$sig <- ifelse(deg_pb$FDR < PARAM$sc_deg_fdr & deg_pb$logFC > PARAM$sc_deg_lfc, "Up",
                     ifelse(deg_pb$FDR < PARAM$sc_deg_fdr & deg_pb$logFC < -PARAM$sc_deg_lfc, "Down", "NS"))
top_lab <- utils::head(deg_pb[deg_pb$sig != "NS", ][order(deg_pb$FDR[deg_pb$sig != "NS"]), ], 20)
p <- ggplot2::ggplot(deg_pb, ggplot2::aes(logFC, -log10(FDR), color = sig)) +
  ggplot2::geom_point(size = 0.6) +
  ggplot2::scale_color_manual(values = c(Up = "#D7301F", Down = "#2C7FB8", NS = "grey75")) +
  ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = gene), size = 3, max.overlaps = 30) +
  ggplot2::theme_classic() + ggplot2::labs(title = "Malignant cells: R-loop high vs low (pseudobulk)")
save_gg(p, file.path(OUT, "Fig3_volcano_pseudobulk.pdf"), 6, 5)

# 补充：细胞层面 Wilcoxon（只作参考，结论以伪 bulk 为准）
Idents(mal) <- "RL_group"
wil <- FindMarkers(mal, ident.1 = "High", ident.2 = "Low", logfc.threshold = 0.25, min.pct = 0.1)
wil$gene <- rownames(wil)
utils::write.csv(wil, file.path(OUT, "cell_level_wilcox_High_vs_Low.csv"), row.names = FALSE)

# =============================================================================
# 6. Hallmark 通路：UCell 打分 → 样本内 High−Low 均值差 → 配对 Wilcoxon
# =============================================================================
H <- get_msig("H")
mal <- UCell::AddModuleScore_UCell(mal, features = H, name = "_UCell", BPPARAM = BiocParallel::SerialParam())
hcols <- paste0(names(H), "_UCell")
hm <- stats::aggregate(mal@meta.data[, hcols], by = list(sample = mal$sample, group = mal$RL_group), FUN = mean)
hall <- do.call(rbind, lapply(hcols, function(h) {
  w <- stats::reshape(hm[, c("sample", "group", h)], idvar = "sample", timevar = "group", direction = "wide")
  d <- w[[paste0(h, ".High")]] - w[[paste0(h, ".Low")]]
  data.frame(pathway = sub("^HALLMARK_", "", sub("_UCell$", "", h)), mean_diff = mean(d),
             p = tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA))
}))
hall$FDR <- stats::p.adjust(hall$p, "BH")
hall <- hall[order(hall$mean_diff), ]
utils::write.csv(hall, file.path(OUT, "Hallmark_High_vs_Low_paired.csv"), row.names = FALSE)
show <- rbind(utils::head(hall, 10), utils::tail(hall, 10))
show$pathway <- factor(show$pathway, levels = show$pathway)
p <- ggplot2::ggplot(show, ggplot2::aes(mean_diff, pathway, fill = mean_diff > 0)) + ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(FDR < 0.05, "*", "")), hjust = -0.3) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#2C7FB8"), guide = "none") +
  ggplot2::theme_classic() + ggplot2::labs(x = "UCell score difference (High - Low, per-sample mean)", y = NULL)
save_gg(p, file.path(OUT, "Fig3_Hallmark_High_vs_Low.pdf"), 7, 6)

# =============================================================================
# 7. CellChat：R-loop 高/低恶性细胞与免疫细胞的通讯差异（只用肿瘤样本）
# =============================================================================
if (requireNamespace("CellChat", quietly = TRUE)) {
  others <- colnames(sc)[sc$tissue == "Tumor" &
                           sc$cell_class %in% c("T_cell", "NK", "B_cell", "Plasma", "Myeloid", "DC", "Mast",
                                                "Fibroblast", "Myofib_SMC", "Endothelial")]
  cc_cells <- c(colnames(mal), others)
  meta <- data.frame(labels = c(paste0("Malig_", mal$RL_group), sc$cell_class[match(others, colnames(sc))]),
                     samples = factor(sc$sample[match(cc_cells, colnames(sc))]), row.names = cc_cells)
  cellchat <- CellChat::createCellChat(object = LayerData(sc, assay = "RNA", layer = "data")[, cc_cells],
                                       meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean")
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(DIR_PRO, "scRNA_cellchat.rds"))

  comm <- CellChat::subsetCommunication(cellchat)
  utils::write.csv(comm[comm$source %in% c("Malig_High", "Malig_Low"), ],
                   file.path(OUT, "CellChat_from_malignant.csv"), row.names = FALSE)
  net <- data.frame(source = rep(rownames(cellchat@net$count), ncol(cellchat@net$count)),
                    target = rep(colnames(cellchat@net$count), each = nrow(cellchat@net$count)),
                    count = as.vector(cellchat@net$count), weight = as.vector(cellchat@net$weight))
  utils::write.csv(net[net$source %in% c("Malig_High", "Malig_Low"), ],
                   file.path(OUT, "CellChat_count_weight_from_malignant.csv"), row.names = FALSE)
  targets <- intersect(c("T_cell", "NK", "Myeloid", "DC", "B_cell"), levels(cellchat@idents))
  p <- CellChat::netVisual_bubble(cellchat, sources.use = c("Malig_High", "Malig_Low"),
                                  targets.use = targets, remove.isolate = TRUE)
  save_gg(p, file.path(OUT, "Fig3_CellChat_bubble.pdf"), 10, 12)
} else log_msg("未安装 CellChat，跳过通讯分析")

# =============================================================================
# 保存
# =============================================================================
saveRDS(sc, file.path(DIR_PRO, "scRNA_integrated.rds"))
saveRDS(mal, file.path(DIR_PRO, "scRNA_malignant.rds"))
meta_out <- sc@meta.data[, intersect(c("sample", "dataset", "tissue", "T_stage", "seurat_clusters", "celltype", "copykat",
                                       "cell_class", "RL_core_UCell", "RL_all_UCell", "RHI"), colnames(sc@meta.data))]
meta_out$RL_group <- mal$RL_group[match(rownames(meta_out), colnames(mal))]
utils::write.csv(meta_out, file.path(OUT, "cell_metadata.csv"))
save_session(STEP)
log_msg("03 完成")

# =============================================================================
# 09b_spatial_check.R —— 空间转录组的补充检查（在 09 步之后运行）
#   1. spot 测序深度：各评分与 log10(UMI) 的相关；校正 log10(UMI)（及增殖评分）后的秩偏相关。
#      09 步发现几乎所有程序在 spot 水平都两两正相关（包括增殖与 IFN-γ），需要排除深度/细胞数造成的假相关。
#   2. 核心基因与 RL_core 的 spot 水平相关：校正 log10(UMI)、以及校正 log10(UMI) + 增殖评分后的偏相关。
#   3. 导出精简的 spot 表（坐标、评分、核心基因表达），供 11 步画 Figure 9，作图时不必再读入 Seurat 对象。
# 输入：data/processed/spatial/BCxx.rds（09 步生成）、results/09_hub_gene_spatial/SHAP_gene_ranking.csv
# 输出：results/09_hub_gene_spatial/Spatial_depth_check_by_slice.csv、Spatial_hub_partial_by_slice.csv
#       data/processed/spatial_spots_slim.rds
# 内存：一次只读入一张切片（每张约 45 MB 的 rds），用完即释放。
# =============================================================================
source("R/00_config.R")
need(c("Seurat"))
STEP <- "09_hub_gene_spatial"; OUT <- res_dir(STEP)
SP <- sprintf("BC%02d", 1:22)
hubs <- utils::head(utils::read.csv(file.path(OUT, "SHAP_gene_ranking.csv"))$gene, 3)

# 秩偏相关：x、y 的秩分别对 z（一列或多列）的秩回归，取残差的 Pearson 相关
prho <- function(x, y, z) {
  Z <- as.data.frame(lapply(as.data.frame(z), rank))
  rx <- stats::resid(stats::lm(rank(x) ~ ., data = Z)); ry <- stats::resid(stats::lm(rank(y) ~ ., data = Z))
  stats::cor(rx, ry)
}
sp_cor <- function(x, y) stats::cor(x, y, method = "spearman")

dep <- list(); hub <- list(); slim <- list()
for (smp in SP) {
  f <- file.path(DIR_PRO, "spatial", paste0(smp, ".rds"))
  if (!file.exists(f)) { log_msg("缺少 ", f, "，跳过"); next }
  ob <- readRDS(f)
  m <- ob@meta.data
  lu <- log10(m$nCount_Spatial)
  dep[[smp]] <- data.frame(sample = smp, n_spots = nrow(m), median_UMI = stats::median(m$nCount_Spatial),
    rho_UMI_RLcore = sp_cor(lu, m$RL_core_UCell), rho_UMI_Prolif = sp_cor(lu, m$Prolif),
    rho_UMI_IFNg = sp_cor(lu, m$IFN_gamma_UCell), rho_UMI_GEP18 = sp_cor(lu, m$GEP18_UCell),
    rho_RLcore_Prolif = sp_cor(m$RL_core_UCell, m$Prolif),
    prho_RLcore_Prolif_given_UMI = prho(m$RL_core_UCell, m$Prolif, lu),
    rho_RLcore_IFNg = sp_cor(m$RL_core_UCell, m$IFN_gamma_UCell),
    prho_RLcore_IFNg_given_UMI = prho(m$RL_core_UCell, m$IFN_gamma_UCell, lu),
    prho_RLcore_IFNg_given_Prolif_UMI = prho(m$RL_core_UCell, m$IFN_gamma_UCell, cbind(m$Prolif, lu)),
    rho_Prolif_IFNg = sp_cor(m$Prolif, m$IFN_gamma_UCell),
    prho_Prolif_IFNg_given_UMI = prho(m$Prolif, m$IFN_gamma_UCell, lu),
    rho_RLcore_GEP18 = sp_cor(m$RL_core_UCell, m$GEP18_UCell),
    prho_RLcore_GEP18_given_Prolif_UMI = prho(m$RL_core_UCell, m$GEP18_UCell, cbind(m$Prolif, lu)))
  hs <- intersect(hubs, rownames(ob))
  ex <- as.matrix(SeuratObject::LayerData(ob, assay = "Spatial", layer = "data")[hs, rownames(m), drop = FALSE])
  for (g in hs) hub[[paste(smp, g)]] <- data.frame(sample = smp, gene = g, pct_detected = mean(ex[g, ] > 0),
    rho = sp_cor(m$RL_core_UCell, ex[g, ]),
    prho_given_UMI = prho(m$RL_core_UCell, ex[g, ], lu),
    prho_given_Prolif_UMI = prho(m$RL_core_UCell, ex[g, ], cbind(m$Prolif, lu)))
  # 坐标：Seurat v5（VisiumV2）返回 x、y、cell；旧版（VisiumV1）返回 imagerow、imagecol
  xy <- as.data.frame(SeuratObject::GetTissueCoordinates(ob))
  if ("cell" %in% names(xy)) rownames(xy) <- xy$cell
  if (all(c("imagerow", "imagecol") %in% names(xy))) { xy$x <- xy$imagerow; xy$y <- xy$imagecol }
  slim[[smp]] <- data.frame(sample = smp, spot = rownames(m), x = xy[rownames(m), "x"], y = xy[rownames(m), "y"],
                            UMI = m$nCount_Spatial, RL_core = m$RL_core_UCell, Prolif = m$Prolif,
                            IFN_gamma = m$IFN_gamma_UCell, t(ex), check.names = FALSE, row.names = NULL)
  rm(ob, m, ex); invisible(gc())
}
dep <- do.call(rbind, dep); hub <- do.call(rbind, hub)
utils::write.csv(dep, file.path(OUT, "Spatial_depth_check_by_slice.csv"), row.names = FALSE)
utils::write.csv(hub, file.path(OUT, "Spatial_hub_partial_by_slice.csv"), row.names = FALSE)
saveRDS(do.call(rbind, slim), file.path(DIR_PRO, "spatial_spots_slim.rds"))

md <- function(v) sprintf("%.2f (%.2f to %.2f)", stats::median(v), min(v), max(v))
log_msg("切片数：", nrow(dep), "；spot 中位 UMI 的范围：", round(min(dep$median_UMI)), "–", round(max(dep$median_UMI)))
for (v in setdiff(names(dep), c("sample", "n_spots", "median_UMI"))) log_msg(v, "：中位数（范围）", md(dep[[v]]))
for (g in hubs) {
  h <- hub[hub$gene == g, ]
  log_msg(g, "：检出 spot 比例 ", md(h$pct_detected), "；rho ", md(h$rho), "；校正 UMI ", md(h$prho_given_UMI),
          "；校正 UMI + 增殖 ", md(h$prho_given_Prolif_UMI), "；偏相关 > 0 的切片 ", sum(h$prho_given_Prolif_UMI > 0), "/", nrow(h))
}
writeLines(capture.output(utils::sessionInfo()), file.path(OUT, "sessionInfo_09b.txt"))
log_msg("09b 完成")

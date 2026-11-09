# =============================================================================
# compat_rfsrc_var_select.R —— 给 Mime1 补回 randomForestSRC::var.select()
#   Mime1 的 RSF 组合模型调用 var.select(object = fit, conservative = "high") 选变量；
#   randomForestSRC 3.4.0 删除了 var.select()（NEWS.md），本环境是 3.10.0。
#   R/vendor/var.select.rfsrc_3.3.1.R 是 randomForestSRC 3.3.1（CRAN 存档，GPL ≥ 3，作者 Ishwaran & Kogalur）
#   中 R/var.select.rfsrc.R 的原样副本，未做任何修改；函数在 randomForestSRC 命名空间里执行，
#   使用该包现有的内部函数（max.subtree、get.varselect.* 等），算法与 Mime1 开发时的版本一致。
# =============================================================================
local({
  e <- new.env(parent = asNamespace("randomForestSRC"))
  sys.source(file.path("R", "vendor", "var.select.rfsrc_3.3.1.R"), envir = e)
  assign("var.select", e$var.select.rfsrc, envir = globalenv())
  assign("var.select.rfsrc", e$var.select.rfsrc, envir = globalenv())
})

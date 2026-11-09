# Vendored file

`var.select.rfsrc_3.3.1.R` is an unmodified copy of `R/var.select.rfsrc.R` from randomForestSRC 3.3.1 by Hemant Ishwaran and Udaya B. Kogalur (CRAN archive: https://cran.r-project.org/src/contrib/Archive/randomForestSRC/). randomForestSRC is distributed under the GNU General Public License, version 3 or later. This file keeps that licence (full text in `GPL-3`) and is not covered by the MIT licence of the rest of the repository.

`R/compat_rfsrc_var_select.R` loads the file inside the randomForestSRC namespace, because Mime1 calls `var.select()`, which was removed in randomForestSRC 3.4.0.

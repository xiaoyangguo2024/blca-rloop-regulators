#!/bin/bash
# 依次运行一个或多个步骤，日志写到 logs/<步骤>.log；某一步失败就停下
# 用法（在 BLCA_Rloop/ 目录下）：
#   bash run_step.sh 02_download_bulk
#   bash run_step.sh 04_bulk_score 05_candidates 06_ml_signature
# 默认使用项目自带的 conda 环境 ./conda_env（见 README）；没有这个环境时用系统里的 Rscript
cd "$(dirname "$0")" || exit 1
mkdir -p logs
export PATH="$HOME/.local/bin:$PATH"
if [ -x ./conda_env/bin/Rscript ] && command -v micromamba >/dev/null 2>&1; then
  RUN=(micromamba run -p ./conda_env Rscript)
elif [ -x ./conda_env/bin/Rscript ]; then
  RUN=(./conda_env/bin/Rscript)
else
  RUN=(Rscript)
fi
export LANG=${LANG:-en_US.UTF-8}
# 缓存放在项目内（msigdbr 的基因集缓存已在 env/cache 里，不必重新下载）
export XDG_CACHE_HOME="$PWD/env/cache/xdg_cache" XDG_DATA_HOME="$PWD/env/cache/xdg_data" XDG_CONFIG_HOME="$PWD/env/cache/xdg_config" CONDA_PKGS_DIRS="$PWD/env/cache/conda_pkgs"
for s in "$@"; do
  s=${s%.R}
  echo "[$(date '+%F %T')] 开始 $s"
  if ! "${RUN[@]}" "R/$s.R" > "logs/$s.log" 2>&1; then
    echo "[$(date '+%F %T')] 失败：$s（见 logs/$s.log 末尾）"
    tail -20 "logs/$s.log"
    exit 1
  fi
  echo "[$(date '+%F %T')] 完成 $s"
done

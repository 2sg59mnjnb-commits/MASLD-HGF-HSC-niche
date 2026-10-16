#!/usr/bin/env bash
# =============================================================================
# MR 数据准备：把 iCloud 上被驱逐的文件下载到本地，并校验完整性
#
# 背景：本机开启了 iCloud「优化存储」，项目里绝大多数大文件只有占位符
#       （stat 显示 blocks=0），读取会挂起。此外 eQTLGen 与 FinnGen 两个文件
#       本身是残缺/损坏的下载。
#
# 用法：  cd ~/Desktop/sjfx3 && bash MR_数据准备.sh
# =============================================================================
set -u
cd "$(dirname "$0")"

hr() { printf '%s\n' "------------------------------------------------------------"; }

# ---- 工具：判断文件是否已在本地 -------------------------------------------
is_local() { [ "$(stat -f %b "$1" 2>/dev/null || echo 0)" -gt 0 ]; }

check() {   # check <路径> <期望字节数下限>
  local f="$1" min="${2:-1}" sz
  [ -e "$f" ] || { printf '  ✗ 不存在      %s\n' "$f"; return 1; }
  sz=$(stat -f %z "$f")
  if ! is_local "$f"; then printf '  ✗ iCloud 占位  %s\n' "$f"; return 1; fi
  if [ "$sz" -lt "$min" ]; then printf '  ✗ 尺寸不足(%s) %s\n' "$sz" "$f"; return 1; fi
  printf '  ✓ 本地就绪    %s (%s bytes)\n' "$f" "$sz"
}

hr; echo "[1/4] 触发 iCloud 下载"; hr
FILES=(
  "data/MR/GTEx_Analysis_v8_eQTL(1).tar"
  "data/MR/GCST90091033_buildGRCh37.tsv.gz"
  "data/MR/1kg.v3.tgz"
  "data/MR/cis-eQTL-SMR_20191212.tar.gz"
  "g1000_eur.zip"
)
for f in "${FILES[@]}"; do
  [ -e "$f" ] && { echo "  → $f"; brctl download "$f" 2>/dev/null; }
done

echo "  等待 iCloud 下载完成（每 15 秒检查一次，最多 40 分钟）…"
for i in $(seq 1 160); do
  sleep 15; pending=0
  for f in "${FILES[@]}"; do [ -e "$f" ] && ! is_local "$f" && pending=$((pending+1)); done
  [ "$pending" -eq 0 ] && { echo "  全部就绪"; break; }
  [ $((i % 8)) -eq 0 ] && echo "    [$((i*15))s] 仍有 $pending 个文件在下载"
done

hr; echo "[2/4] 重新下载 eQTLGen 完整版（你本地那份只有 1,204 个基因，是残缺的）"; hr
EQTL_URL="https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
EQTL_OUT="data/MR/eQTLGen_cis_full.txt.gz"
echo "  支持断点续传；若中断可重复运行本脚本"
curl -L --fail --retry 5 --retry-delay 10 -C - -o "$EQTL_OUT" "$EQTL_URL" \
     -w "  HTTP %{http_code} | 本次下载 %{size_download} bytes | 平均 %{speed_download} B/s\n" || true
if gzip -t "$EQTL_OUT" 2>/dev/null; then echo "  ✓ gzip 校验通过"
else echo "  ✗ gzip 校验失败，请重跑本脚本继续下载"; fi

hr; echo "[3/4] 重新下载 FinnGen R12 肝纤维化/肝硬化 GWAS"; hr
FG_URL="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_K11_FIBROCHIRLIV.gz"
FG_OUT="data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"
echo "  注意：你现有的 'finngen_R12_K11_FIBROCHIRLIV(1).gz' 字节流已损坏"
echo "        （染色体列出现乱码，可读 25 万行后 CRC 失败），必须重下"
curl -L --fail --retry 5 --retry-delay 10 -C - -o "$FG_OUT" "$FG_URL" \
     -w "  HTTP %{http_code} | 本次下载 %{size_download} bytes | 平均 %{speed_download} B/s\n" || true
if gzip -t "$FG_OUT" 2>/dev/null; then echo "  ✓ gzip 校验通过"
else echo "  ✗ gzip 校验失败，请重跑本脚本"; fi

hr; echo "[4/4] 最终校验"; hr
check "data/MR/eQTLGen_cis_full.txt.gz"                400000000
check "data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"     100000000
check "data/MR/GTEx_Analysis_v8_eQTL(1).tar"            700000000
check "g1000_eur.zip"                                   400000000

cat <<'EOF'

全部就绪后，运行：
    Rscript MR_pipeline_corrected.R

说明：eQTLGen 完整版约 500 MB，实测约 64 KB/s → 需 1-2 小时。
      建议挂着跑，或用有线网络 / 换时段重试。
EOF

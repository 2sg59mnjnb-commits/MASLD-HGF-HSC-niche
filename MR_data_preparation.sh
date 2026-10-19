#!/usr/bin/env bash
# =============================================================================
# MR data preparation: download iCloud-evicted files locally and verify integrity
#
# Background: this machine has iCloud 'Optimise Storage' enabled, so most large project
#       files are placeholders only (stat reports blocks=0) and reads hang. In addition,
#       the eQTLGen and FinnGen files themselves are incomplete/corrupt downloads.
#
# Usage:  cd ~/Desktop/sjfx3 && bash MR_data_preparation.sh
# =============================================================================
set -u
cd "$(dirname "$0")"

hr() { printf '%s\n' "------------------------------------------------------------"; }

# ---- Helper: is the file already local? ------------------------------------
is_local() { [ "$(stat -f %b "$1" 2>/dev/null || echo 0)" -gt 0 ]; }

check() {   # check <path> <minimum expected bytes>
  local f="$1" min="${2:-1}" sz
  [ -e "$f" ] || { printf '  x missing      %s\n' "$f"; return 1; }
  sz=$(stat -f %z "$f")
  if ! is_local "$f"; then printf '  x iCloud stub  %s\n' "$f"; return 1; fi
  if [ "$sz" -lt "$min" ]; then printf '  x too small(%s) %s\n' "$sz" "$f"; return 1; fi
  printf '  v ready        %s (%s bytes)\n' "$f" "$sz"
}

hr; echo "[1/4] Triggering iCloud download"; hr
FILES=(
  "data/MR/GTEx_Analysis_v8_eQTL(1).tar"
  "data/MR/GCST90091033_buildGRCh37.tsv.gz"
  "data/MR/1kg.v3.tgz"
  "data/MR/cis-eQTL-SMR_20191212.tar.gz"
  "g1000_eur.zip"
)
for f in "${FILES[@]}"; do
  [ -e "$f" ] && { echo "  -> $f"; brctl download "$f" 2>/dev/null; }
done

echo "  Waiting for iCloud downloads (checked every 15 s, up to 40 min)..."
for i in $(seq 1 160); do
  sleep 15; pending=0
  for f in "${FILES[@]}"; do [ -e "$f" ] && ! is_local "$f" && pending=$((pending+1)); done
  [ "$pending" -eq 0 ] && { echo "  all ready"; break; }
  [ $((i % 8)) -eq 0 ] && echo "    [$((i*15))s] $pending file(s) still downloading"
done

hr; echo "[2/4] Re-downloading the full eQTLGen release (the local copy has only 1,204 genes)"; hr
EQTL_URL="https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
EQTL_OUT="data/MR/eQTLGen_cis_full.txt.gz"
echo "  Resumable; re-run this script if interrupted"
curl -L --fail --retry 5 --retry-delay 10 -C - -o "$EQTL_OUT" "$EQTL_URL" \
     -w "  HTTP %{http_code} | downloaded %{size_download} bytes | mean %{speed_download} B/s\n" || true
if gzip -t "$EQTL_OUT" 2>/dev/null; then echo "  v gzip check passed"
else echo "  x gzip check failed -- re-run this script to continue"; fi

hr; echo "[3/4] Re-downloading FinnGen R12 liver fibrosis/cirrhosis GWAS"; hr
FG_URL="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_K11_FIBROCHIRLIV.gz"
FG_OUT="data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"
echo "  Note: the existing 'finngen_R12_K11_FIBROCHIRLIV(1).gz' byte stream is corrupt"
echo "        (garbled chromosome values, CRC failure after ~250k lines); re-download required"
curl -L --fail --retry 5 --retry-delay 10 -C - -o "$FG_OUT" "$FG_URL" \
     -w "  HTTP %{http_code} | downloaded %{size_download} bytes | mean %{speed_download} B/s\n" || true
if gzip -t "$FG_OUT" 2>/dev/null; then echo "  v gzip check passed"
else echo "  x gzip check failed -- re-run this script"; fi

hr; echo "[4/4] Final check"; hr
check "data/MR/eQTLGen_cis_full.txt.gz"                400000000
check "data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"     100000000
check "data/MR/GTEx_Analysis_v8_eQTL(1).tar"            700000000
check "g1000_eur.zip"                                   400000000

cat <<'EOF'

Once everything is ready, run:
    Rscript MR_pipeline_corrected.R

Note: the full eQTLGen release is ~500 MB per gene block; at ~64 KB/s this takes 1-2 h.
      It is best to leave it running, or use a wired connection / retry at another time.
EOF

#!/bin/bash
# CRESSENT Benchmarking Script - Simple Version (fixed)
set -Eeuo pipefail

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

BENCHMARK_DIR="benchmark_genomoviridae_results_rep_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BENCHMARK_DIR"

LOG_FILE="$BENCHMARK_DIR/benchmark.log"

print_status() { echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"; }
print_header() {
  echo -e "${YELLOW}========================================${NC}" | tee -a "$LOG_FILE"
  echo -e "${YELLOW}$1${NC}" | tee -a "$LOG_FILE"
  echo -e "${YELLOW}========================================${NC}" | tee -a "$LOG_FILE"
}

# Run a command with timing; $3 is a *base path* without extension
run_with_timing() {
  local description="$1"; local command="$2"; local base="$3"
  mkdir -p "$(dirname "$base")"
  print_status "Running: $description"
  /usr/bin/time -f "\n[METRICS]\nReal time: %E\nUser time: %U seconds\nSystem time: %S seconds\nMax memory: %M KB\nAverage memory: %K KB" \
    -o "${base}.time" \
    bash -c "$command" > "${base}.log" 2>&1 || {
      print_status "${RED}FAILED:${NC} $description (see ${base}.log)"
      return 1
    }
  # Show metrics
  [[ -f "${base}.time" ]] && cat "${base}.time" | tee -a "$LOG_FILE"
}

# Datasets
GENOME="./test/genomoviridae/genomo_reann_1a_virus_sequences.fna"
PROTEINS="./test/genomoviridae/genomo_reann_1a_virus_AA.faa"
CAPS="./test/genomoviridae/genomo_reann_1a_virus_AA.caps.faa"
REPS="./test/genomoviridae/genomo_reann_1a_virus_AA.reps.faa"

# Create subdirs we'll use
mkdir -p "$BENCHMARK_DIR/rep/cluster" "$BENCHMARK_DIR/rep/align" "$BENCHMARK_DIR/rep/tree" "$BENCHMARK_DIR/rep/motif"

# Sizes
print_status "Dataset sizes:"
{
  echo "dataset (Naryaviridae):"
  echo "  - Genome sequences: $(grep -c '>' "$GENOME" 2>/dev/null || echo 0)"
  echo "  - Protein sequences: $(grep -c '>' "$PROTEINS" 2>/dev/null || echo 0)"
  echo "  - Caps sequences: $(grep -c '>' "$CAPS" 2>/dev/null || echo 0)"
  echo "  - Reps sequences: $(grep -c '>' "$REPS" 2>/dev/null || echo 0)"
} | tee -a "$LOG_FILE"


############# REPS ###################3
print_header "BENCHMARKING: REP PROTEINS - Cluster > Align > Tree > Motif Discovery"
# 1) CLUSTERING 
print_header "BENCHMARKING: CLUSTERING MODULE"
run_with_timing "Clustering - Reps" \
  "cressent cluster -i '$REPS' -o '$BENCHMARK_DIR/rep/cluster' -t 40 --min_ani 95.0" \
  "$BENCHMARK_DIR/rep/cluster/cluster"

# 2) ALIGNMENT (Rep or Cap)
print_header "BENCHMARKING: ALIGNMENT MODULE"
run_with_timing "Alignment - Rep dataset" \
  "cressent align -i '$REPS' -o '$BENCHMARK_DIR/rep/align' -t 40 --db_family 'Genomoviridae' --protein_type reps --db_path /fs/project/PAS1117/ricardo/cressent/DB" \
  "$BENCHMARK_DIR/rep/align/align"

# 3) TREE BUILDING (use aligned file from align step)
print_header "BENCHMARKING: TREE BUILDING MODULE"
REP_ALIGNED=$(ls -1 "$BENCHMARK_DIR"/rep/align/*aligned_trimmed*.fasta 2>/dev/null | head -1 || true)
if [[ -z "${REP_ALIGNED:-}" ]]; then
  REP_ALIGNED=$(ls -1 "$BENCHMARK_DIR"/rep/align/*aligned*.fasta 2>/dev/null | head -1 || true)
fi
if [[ -n "${REP_ALIGNED:-}" ]]; then
  run_with_timing "Tree Building - Rep dataset" \
    "cressent build_tree -i '$REP_ALIGNED' -o '$BENCHMARK_DIR/rep/tree'" \
    "$BENCHMARK_DIR/rep/tree/tree"
else
  print_status "${YELLOW}Skipping tree: no aligned file found in $BENCHMARK_DIR/rep/align${NC}"
fi

# 4) MOTIF DISCOVERY (Rep proteins)
print_header "BENCHMARKING: MOTIF DISCOVERY MODULE"
run_with_timing "Motif Discovery - Rep dataset" \
  "cressent motif_disc -i '$REPS' -o '$BENCHMARK_DIR/rep/motif' --scanprosite" \
  "$BENCHMARK_DIR/rep/motif/motif"

# SUMMARY (Markdown)
print_header "GENERATING SUMMARY REPORT"
REPORT_FILE="$BENCHMARK_DIR/benchmark_summary.md"
{
  echo "# CRESSENT Benchmark Report"
  echo
  echo "**Date**: $(date)"
  echo "**System**: $(uname -a)"
  echo "**CPU**: $(nproc) cores"
  echo "**Memory**: $(free -h | awk 'NR==2{print $2}')"
  echo
  echo "## Dataset Information"
  echo
  echo "| Dataset | Type | Sequences |"
  echo "|---------|------|-----------|"
  echo "| Genomoviridae | genomes: $(grep -c '>' "$GENOME" 2>/dev/null || echo 0), proteins: $(grep -c '>' "$PROTEINS" 2>/dev/null || echo 0) |"
  echo
  echo "## Performance Metrics"
  for f in $(find "$BENCHMARK_DIR" -type f -name '*.time' | sort); do
    ds=$(basename "$(dirname "$f")")        # cap
    mod=$(basename "${f%.time}")            # cluster|align|tree|motif
    echo
    echo "### Module: ${mod^} — Dataset: $ds"
    echo '```'
    cat "$f"
    echo '```'
  done
} > "$REPORT_FILE"

# CSV summary
CSV_FILE="$BENCHMARK_DIR/benchmark_metrics.csv"
echo "Module,Dataset,Real_Time,User_Time,System_Time,Max_Memory_KB" > "$CSV_FILE"
while IFS= read -r -d '' f; do
  ds=$(basename "$(dirname "$f")")
  mod=$(basename "${f%.time}")
  real_time=$(grep "Real time:" "$f" | awk '{print $3}')
  user_time=$(grep "User time:" "$f" | awk '{print $3}')
  sys_time=$(grep "System time:" "$f" | awk '{print $3}')
  max_mem=$(grep "Max memory:" "$f" | awk '{print $3}')
  echo "$mod,$ds,$real_time,$user_time,$sys_time,$max_mem" >> "$CSV_FILE"
done < <(find "$BENCHMARK_DIR" -type f -name '*.time' -print0)

print_header "BENCHMARKING COMPLETE"
print_status "Results saved to: $BENCHMARK_DIR"
print_status "Summary report: $REPORT_FILE"
print_status "CSV metrics saved to: $CSV_FILE"

#!/usr/bin/env bash
# run.sh (repo root)
#
# EXP-2 (thesis §6.6): Synthetic scalability + phase breakdown (Dim=2 locked)
# -------------------------------------------------------------------------
# This script runs the full d=2 grid described in your §6.6 (A/B/C):
#   (A) vary alpha_out
#   (B) vary N (|R|=|S|=N)
#   (C) vary t
#
# It executes 7 model variants:
#   ours/{enum_sampling,sampling,adaptive}
#   range_tree/{enum_sampling,sampling,adaptive}
#   kd_tree/sampling
#
# Features:
#   - Parallel execution with a conservative, memory-aware default.
#   - Robust failure handling:
#       * a failed run does NOT stop the whole experiment
#       * failures are automatically retried (MAX_RETRIES)
#       * rerunning this script will resume and only rerun unfinished/failed tasks
#   - Per-task logs + per-task CSV; final merged CSVs per sweep.
#
# Notes:
#   - This repo build is Dim=2 only.
#   - Synthetic data generation uses the local Alacarte source only.
#
# Usage:
#   chmod +x run.sh
#   ./run.sh
#
# Common overrides (env vars):
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   CLEAN_BUILD=0|1
#   JOBS=8                     # build parallelism
#   THREADS=1                  # passed to --threads (internal baseline threads)
#   REPEATS=3                  # per-task repeats inside sjs_run
#   SEED=1                     # base seed for sampling (rep adds +rep)
#   GEN_SEED=1                 # dataset generation seed
#   GEN=alacarte_rectgen       # dataset generator name (aliases supported)
#   RECTGEN_SCRIPT=tools/alacarte_rectgen_generate.py
#   AUDIT_PAIRS=2000000
#   AUDIT_SEED=1
#
#   # Framework knobs (Chapter 6.4.3):
#   BUDGET_B=10000000          # B (maps to --j_star)
#   W_SMALL=1024               # full-cache threshold (cfg.run.extra["w_small"])
#
#   # Enum safety cap (0 disables cap; beware OOM for dense joins)
#   ENUM_CAP=0
#
#   # Concurrency & robustness
#   MAX_PARALLEL=4             # hard cap on concurrent tasks
#   PAR_ALPHA=8                # optional per-sweep override
#   PAR_N=2                    # optional per-sweep override
#   PAR_T_SMALL=4              # optional per-sweep override
#   PAR_T_LARGE=1              # optional per-sweep override
#   MAX_PARALLEL_CAP=12        # default cap used when MAX_PARALLEL is unset
#   MEM_PER_TASK_GB=12         # heuristic memory per task when MAX_PARALLEL is unset
#   MEM_RESERVE_GB=8           # reserved memory for OS / page cache
#   MAX_RETRIES=2
#   RETRY_BACKOFF_SEC=2        # linear backoff between retry rounds
#   TIMEOUT_SEC=0              # 0 disables timeout; if >0 uses `timeout` when available
#
#   # Parameter grids
#   ALPHAS_A="0.1 0.3 1 3 10 30 100 300 1000"
#   N_FIXED=1000000
#   T_FIXED=100000
#   ALPHA_FIXED=200
#   NS_B="100000 200000 500000 1000000 2000000 5000000"
#   TS_C="100000 300000 1000000 3000000 10000000 30000000 100000000 300000000"

set -Eeuo pipefail
IFS=$' \t\n'

trap 'echo -e "[run.sh][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

log()  { echo -e "[run.sh] $*"; }
warn() { echo -e "[run.sh][WARN] $*" >&2; }

task_count() {
  # Count runnable tasks in a TSV (exclude blank lines and comment/header lines starting with '#').
  local tf="$1"
  [[ -f "$tf" ]] || { echo 0; return 0; }
  awk 'BEGIN{c=0} !/^[[:space:]]*($|#)/ {c++} END{print c}' "$tf"
}


nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  else
    echo 4
  fi
}

mem_gb_safe() {
  # Best-effort physical memory size (GB).
  if [[ -r /proc/meminfo ]]; then
    awk '/MemTotal:/ {printf("%.0f", $2/1024/1024)}' /proc/meminfo
  else
    echo 16
  fi
}

build_subdir_for_type() {
  case "$1" in
    Release|release) echo "release" ;;
    Debug|debug) echo "debug" ;;
    RelWithDebInfo|relwithdebinfo) echo "relwithdebinfo" ;;
    MinSizeRel|minsizerel) echo "minsizerel" ;;
    *)
      echo "$(echo "$1" | tr '[:upper:]' '[:lower:]')" ;;
  esac
}

find_exe() {
  local build_dir="$1"
  local name="$2"
  local c
  for c in \
    "${build_dir}/${name}" \
    "${build_dir}/apps/${name}" \
    "${build_dir}/src/apps/${name}" \
    ; do
    if [[ -x "${c}" ]]; then
      echo "${c}"
      return 0
    fi
  done
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" && -x "${p}" ]] || return 1
  echo "${p}"
}

cmake_cache_internal_value() {
  local cache="$1"
  local key="$2"
  [[ -f "$cache" ]] || return 1
  local line
  line="$(grep -E "^${key}:INTERNAL=" "$cache" | head -n 1 || true)"
  [[ -n "$line" ]] || return 1
  echo "${line#*=}"
}

ensure_cmake_configured() {
  local cache="${BUILD_DIR}/CMakeCache.txt"

  if [[ ! -f "$cache" ]]; then
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  local cache_home=""
  local cache_dir=""
  cache_home="$(cmake_cache_internal_value "$cache" "CMAKE_HOME_DIRECTORY" || true)"
  cache_dir="$(cmake_cache_internal_value "$cache" "CMAKE_CACHEFILE_DIR" || true)"

  # If repository path changed (e.g., directory renamed/moved), old absolute
  # paths in CMakeCache can break cmake --build. Recreate this build dir once.
  if [[ -n "$cache_home" && "$cache_home" != "$ROOT" ]]; then
    warn "Detected stale CMake cache source path:"
    warn "  cache: $cache_home"
    warn "  root : $ROOT"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  if [[ -n "$cache_dir" && "$cache_dir" != "$BUILD_DIR" ]]; then
    warn "Detected stale CMake cache build path:"
    warn "  cache: $cache_dir"
    warn "  build: $BUILD_DIR"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi
}

sanitize_token() {
  # Make a token safe for filenames.
  # - replace '.' with 'p' (e.g., 0.1 -> 0p1)
  # - replace '-' with 'm'
  # - replace other non [A-Za-z0-9_] with '_'
  echo "$1" \
    | sed -e 's/-/m/g' -e 's/\./p/g' -e 's/[^A-Za-z0-9_]/_/g'
}

supports_wait_n() {
  # Avoid probing by executing `wait -n` because that can accidentally reap jobs.
  local major="${BASH_VERSINFO[0]:-0}"
  local minor="${BASH_VERSINFO[1]:-0}"
  (( major > 4 || (major == 4 && minor >= 3) ))
}

wait_for_slot() {
  local max_parallel="$1"
  while true; do
    local running
    running="$(jobs -pr | wc -l | tr -d '[:space:]')"
    [[ -z "${running}" ]] && running=0
    if (( running < max_parallel )); then
      return 0
    fi
    # Wait for at least one job to finish.
    if supports_wait_n; then
      wait -n || true
    else
      sleep 0.2
    fi
  done
}

csv_col_idx() {
  local file="$1"
  local col="$2"
  awk -F',' -v c="$col" 'NR==1{for(i=1;i<=NF;i++){if($i==c){print i; exit}}}' "$file"
}

csv_all_ok() {
  local file="$1"
  local ok_idx
  ok_idx="$(csv_col_idx "$file" "ok")"
  [[ -n "${ok_idx}" ]] || return 2
  # all ok==1, and at least one data row
  awk -F',' -v k="${ok_idx}" 'NR==1{next} {n++; if($k != 1) bad=1} END{ if(n==0) exit 3; exit(bad?1:0) }' "$file"
}

csv_row_count() {
  local file="$1"
  awk 'END{print NR-1}' "$file" 2>/dev/null || echo 0
}

merge_csv_dir() {
  local dir="$1"
  local out_csv="$2"
  local first=""
  first="$(find "$dir" -type f -name run.csv | sort | head -n 1 || true)"
  if [[ -z "$first" ]]; then
    warn "No run.csv files found under: $dir"
    return 0
  fi
  mkdir -p "$(dirname "$out_csv")"
  head -n 1 "$first" > "$out_csv"
  # Append all data rows from all run.csv files.
  while IFS= read -r f; do
    tail -n +2 "$f" >> "$out_csv" || true
  done < <(find "$dir" -type f -name run.csv | sort)
}

build_run_signature() {
  {
    echo "build_type=$BUILD_TYPE"
    echo "threads=$THREADS"
    echo "repeats=$REPEATS"
    echo "seed=$SEED"
    echo "gen=$GEN"
    echo "gen_seed=$GEN_SEED"
    echo "rectgen_script=$RECTGEN_SCRIPT"
    echo "audit_pairs=$AUDIT_PAIRS"
    echo "audit_seed=$AUDIT_SEED"
    echo "budget_B=$BUDGET_B"
    echo "w_small=$W_SMALL"
    echo "enum_cap=$ENUM_CAP"
    echo "alphas_A=$ALPHAS_A"
    echo "N_fixed=$N_FIXED"
    echo "t_fixed=$T_FIXED"
    echo "alpha_fixed=$ALPHA_FIXED"
    echo "Ns_B=$NS_B"
    echo "ts_C=$TS_C"
  } | cksum | awk '{print $1 "-" $2}'
}

# ----------------------------
# Resolve repo root
# ----------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------
# Parameters
# ----------------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
CLEAN_TEMP="${CLEAN_TEMP:-0}"
# Build parallelism: default is conservative to avoid OOM on large machines
# (override with JOBS=... if you want maximum build speed).
if [[ -n "${JOBS:-}" ]]; then
  JOBS="$JOBS"
else
  # Derive a safe default from memory.
  # Rule of thumb: ~1 compile job per GB is usually OK for C++ builds.
  local_nproc="$(nproc_safe)"
  local_mem_gb="$(mem_gb_safe)"
  JOBS="$local_mem_gb"
  [[ -z "$JOBS" || "$JOBS" -lt 1 ]] && JOBS=1
  (( JOBS > local_nproc )) && JOBS="$local_nproc"
  # Hard cap to avoid accidental oversubscription on huge machines.
  (( JOBS > 16 )) && JOBS=16
fi

THREADS="${THREADS:-1}"
REPEATS="${REPEATS:-3}"
SEED="${SEED:-1}"
GEN_SEED="${GEN_SEED:-1}"
GEN="${GEN:-alacarte_rectgen}"
RECTGEN_SCRIPT="${RECTGEN_SCRIPT:-$ROOT/tools/alacarte_rectgen_generate.py}"
AUDIT_PAIRS="${AUDIT_PAIRS:-2000000}"
AUDIT_SEED="${AUDIT_SEED:-$GEN_SEED}"

BUDGET_B="${BUDGET_B:-10000000}"
W_SMALL="${W_SMALL:-1024}"
ENUM_CAP="${ENUM_CAP:-0}"

MAX_RETRIES="${MAX_RETRIES:-2}"
RETRY_BACKOFF_SEC="${RETRY_BACKOFF_SEC:-2}"
TIMEOUT_SEC="${TIMEOUT_SEC:-0}"

MAX_PARALLEL_CAP="${MAX_PARALLEL_CAP:-12}"
MEM_PER_TASK_GB="${MEM_PER_TASK_GB:-12}"
MEM_RESERVE_GB="${MEM_RESERVE_GB:-8}"

ALPHAS_A="${ALPHAS_A:-0.1 0.3 1 3 10 30 100 300 1000}"
N_FIXED="${N_FIXED:-1000000}"
T_FIXED="${T_FIXED:-100000}"
ALPHA_FIXED="${ALPHA_FIXED:-200}"
NS_B="${NS_B:-100000 200000 500000 1000000 2000000 5000000}"
TS_C="${TS_C:-100000 300000 1000000 3000000 10000000 30000000 100000000 300000000}"

# 7 model variants
MODELS=(
  "ours enum_sampling"
  "ours sampling"
  "ours adaptive"
  "range_tree enum_sampling"
  "range_tree sampling"
  "range_tree adaptive"
  "kd_tree sampling"
)

# ----------------------------
# Directories
# ----------------------------
BUILD_SUBDIR="$(build_subdir_for_type "$BUILD_TYPE")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"

TEMP_ROOT="${ROOT}/run/temp/exp2_synth_d2"
OUT_ROOT="${TEMP_ROOT}/out"
LOG_ROOT="${TEMP_ROOT}/logs"
STATUS_ROOT="${TEMP_ROOT}/status"
MANIFEST_DIR="${TEMP_ROOT}/manifest"

RESULT_ROOT="${ROOT}/results/raw/exp2_synth_d2"

if [[ "$CLEAN_TEMP" == "1" ]]; then
  warn "CLEAN_TEMP=1: removing temp dir: $TEMP_ROOT"
  rm -rf "$TEMP_ROOT"
fi

mkdir -p "$TEMP_ROOT" "$OUT_ROOT" "$LOG_ROOT" "$STATUS_ROOT" "$MANIFEST_DIR"

# ----------------------------
# Build (if needed)
# ----------------------------
log "Repo root  : $ROOT"
log "Build type : $BUILD_TYPE"
log "Build dir  : $BUILD_DIR"
log "Temp dir   : $TEMP_ROOT"
log "Results dir: $RESULT_ROOT"

if [[ "$CLEAN_BUILD" == "1" ]]; then
  log "Cleaning build dir: $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

ensure_cmake_configured

log "Building... (JOBS=$JOBS)"
cmake --build "$BUILD_DIR" -j "$JOBS"

SJS_RUN="$(find_exe "$BUILD_DIR" sjs_run || true)"
SJS_GEN="$(find_exe "$BUILD_DIR" sjs_gen_dataset || true)"
if [[ -z "$SJS_RUN" ]]; then
  echo "[run.sh][FATAL] Could not locate sjs_run under $BUILD_DIR" >&2
  exit 2
fi
log "Using sjs_run: $SJS_RUN"
if [[ -n "$SJS_GEN" ]]; then
  log "Found sjs_gen_dataset: $SJS_GEN"
fi

# ----------------------------
# Parallelism defaults
# ----------------------------
NPROC="$(nproc_safe)"
MEM_GB="$(mem_gb_safe)"

if [[ -n "${MAX_PARALLEL:-}" ]]; then
  MAX_PARALLEL_USER="$MAX_PARALLEL"
else
  # Auto parallelism uses both CPU and memory limits, then caps hard to avoid
  # over-committing on very large hosts.
  cpu_per_task="$THREADS"
  (( cpu_per_task < 1 )) && cpu_per_task=1
  cpu_cap=$(( NPROC / cpu_per_task ))
  (( cpu_cap < 1 )) && cpu_cap=1
  if (( cpu_cap > 1 )); then
    cpu_cap=$(( cpu_cap - 1 ))  # leave one core for system + I/O
  fi

  mem_budget=$(( MEM_GB - MEM_RESERVE_GB ))
  if (( mem_budget < MEM_PER_TASK_GB )); then
    mem_budget="$MEM_PER_TASK_GB"
  fi
  mem_cap=$(( mem_budget / MEM_PER_TASK_GB ))
  (( mem_cap < 1 )) && mem_cap=1

  MAX_PARALLEL_USER="$cpu_cap"
  (( MAX_PARALLEL_USER > mem_cap )) && MAX_PARALLEL_USER="$mem_cap"
  (( MAX_PARALLEL_USER > MAX_PARALLEL_CAP )) && MAX_PARALLEL_USER="$MAX_PARALLEL_CAP"
fi
MAX_PARALLEL="$MAX_PARALLEL_USER"

PAR_ALPHA="${PAR_ALPHA:-$MAX_PARALLEL}"
PAR_N="${PAR_N:-$(( MAX_PARALLEL > 2 ? 2 : MAX_PARALLEL ))}"
PAR_T_SMALL="${PAR_T_SMALL:-$(( MAX_PARALLEL > 4 ? 4 : MAX_PARALLEL ))}"
PAR_T_LARGE="${PAR_T_LARGE:-1}"

(( PAR_ALPHA < 1 )) && PAR_ALPHA=1
(( PAR_N < 1 )) && PAR_N=1
(( PAR_T_SMALL < 1 )) && PAR_T_SMALL=1
(( PAR_T_LARGE < 1 )) && PAR_T_LARGE=1

RUN_SIGNATURE="$(build_run_signature)"

log "Host: nproc=$NPROC, mem~${MEM_GB}GB"
log "Parallelism: max=$MAX_PARALLEL, alpha=$PAR_ALPHA, N=$PAR_N, t_small=$PAR_T_SMALL, t_large=$PAR_T_LARGE"
log "Auto-parallel knobs: cap=$MAX_PARALLEL_CAP, mem_per_task_gb=$MEM_PER_TASK_GB, mem_reserve_gb=$MEM_RESERVE_GB"
log "Run signature: $RUN_SIGNATURE"

# ----------------------------
# Core task runner
# ----------------------------
run_one_task() {
  local task_id="$1"
  local sweep="$2"
  local N="$3"
  local alpha="$4"
  local t="$5"
  local method="$6"
  local variant="$7"
  local force="$8"
  local attempt="${9:-0}"

  local ok_flag="$STATUS_ROOT/${task_id}.ok"
  local fail_flag="$STATUS_ROOT/${task_id}.fail"
  local exit_flag="$STATUS_ROOT/${task_id}.exit"
  local meta_flag="$STATUS_ROOT/${task_id}.meta"

  local out_dir="$OUT_ROOT/${sweep}/${task_id}"
  local log_file="$LOG_ROOT/${task_id}.log"
  local csv_file="$out_dir/run.csv"

  if [[ "$force" != "1" && -f "$ok_flag" && -f "$csv_file" && -f "$meta_flag" ]]; then
    # Only resume if both output CSV and run-signature match.
    if grep -Fxq "run_signature=$RUN_SIGNATURE" "$meta_flag" \
      && csv_all_ok "$csv_file" >/dev/null 2>&1; then
      echo 0 > "$exit_flag"
      return 0
    fi
  fi

  rm -f "$ok_flag" "$fail_flag" "$exit_flag" "$meta_flag"
  rm -rf "$out_dir"
  mkdir -p "$out_dir"

  local ds_name
  ds_name="d2_${GEN}_n${N}_a$(sanitize_token "$alpha")_gs${GEN_SEED}"

  local -a cmd=(
    "$SJS_RUN"
    "--dataset_source=synthetic"
    "--dataset=$ds_name"
    "--dim=2"
    "--gen=$GEN"
    "--n_r=$N" "--n_s=$N"
    "--alpha=$alpha"
    "--gen_seed=$GEN_SEED"
    "--method=$method"
    "--variant=$variant"
    "--t=$t"
    "--seed=$SEED"
    "--repeats=$REPEATS"
    "--threads=$THREADS"
    "--j_star=$BUDGET_B"
    "--enum_cap=$ENUM_CAP"
    "--write_samples=0"
    "--verify=0"
    "--out_dir=$out_dir"
    "--results_file=$csv_file"
    "--log_level=info"
    "--log_timestamp=1"
    "--rectgen_script=$RECTGEN_SCRIPT"
    "--audit_pairs=$AUDIT_PAIRS"
    "--audit_seed=$AUDIT_SEED"
    "--w_small=$W_SMALL"
  )

  {
    echo "run_signature=$RUN_SIGNATURE"
    echo "attempt=$attempt"
    echo "task_id=$task_id"
    echo "sweep=$sweep"
    echo "N=$N"
    echo "alpha=$alpha"
    echo "t=$t"
    echo "method=$method"
    echo "variant=$variant"
    echo "cmd=${cmd[*]}"
  } > "$meta_flag"

  local rc=0
  if (( TIMEOUT_SEC > 0 )) && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SEC" "${cmd[@]}" >"$log_file" 2>&1 || rc=$?
  else
    "${cmd[@]}" >"$log_file" 2>&1 || rc=$?
  fi
  echo "$rc" > "$exit_flag"

  # Validate result CSV.
  if [[ "$rc" -ne 0 ]]; then
    echo "nonzero_exit=$rc" >> "$meta_flag"
    touch "$fail_flag"
    return 0
  fi
  if [[ ! -f "$csv_file" ]]; then
    echo "missing_csv=1" >> "$meta_flag"
    touch "$fail_flag"
    return 0
  fi
  local rows
  rows="$(csv_row_count "$csv_file")"
  if [[ "$rows" -lt "$REPEATS" ]]; then
    echo "short_csv_rows=$rows" >> "$meta_flag"
    touch "$fail_flag"
    return 0
  fi
  if ! csv_all_ok "$csv_file" >/dev/null 2>&1; then
    echo "not_all_ok=1" >> "$meta_flag"
    touch "$fail_flag"
    return 0
  fi

  touch "$ok_flag"
  return 0
}

run_task_file_with_parallel() {
  local task_file="$1"
  local parallel="$2"
  local attempt="$3"
  local force="$4"

  log "Running tasks: $(basename "$task_file") (attempt=$attempt, parallel=$parallel, force=$force)"

  while IFS=$'\t' read -r task_id sweep N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    # Skip comments
    [[ "$task_id" =~ ^# ]] && continue
    wait_for_slot "$parallel"
    (run_one_task "$task_id" "$sweep" "$N" "$alpha" "$t" "$method" "$variant" "$force" "$attempt") &
  done < "$task_file"

  wait || true
}

collect_failures_from_task_file() {
  local task_file="$1"
  local out_fail_list="$2"
  : > "$out_fail_list"

  while IFS=$'\t' read -r task_id sweep N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    if [[ ! -f "$STATUS_ROOT/${task_id}.ok" ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "$sweep" "$N" "$alpha" "$t" "$method" "$variant" >> "$out_fail_list"
    fi
  done < "$task_file"
}

write_manifest_csv() {
  local task_file="$1"
  local out_csv="$2"
  mkdir -p "$(dirname "$out_csv")"
  echo "task_id,sweep,N,alpha,t,method,variant,out_dir,log_file" > "$out_csv"
  while IFS=$'\t' read -r task_id sweep N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    local out_dir="$OUT_ROOT/${sweep}/${task_id}"
    local log_file="$LOG_ROOT/${task_id}.log"
    echo "${task_id},${sweep},${N},${alpha},${t},${method},${variant},${out_dir},${log_file}" >> "$out_csv"
  done < "$task_file"
}

# ----------------------------
# Generate task files
# ----------------------------
TASK_A="$MANIFEST_DIR/tasks_A_alpha.tsv"
TASK_B="$MANIFEST_DIR/tasks_B_N.tsv"
TASK_C_SMALL="$MANIFEST_DIR/tasks_C_t_small.tsv"
TASK_C_LARGE="$MANIFEST_DIR/tasks_C_t_large.tsv"

gen_task_id() {
  local sweep="$1"; local N="$2"; local alpha="$3"; local t="$4"; local method="$5"; local variant="$6"
  echo "${sweep}__n${N}__a$(sanitize_token "$alpha")__t${t}__${method}__${variant}"
}

log "Generating task manifests..."

{
  echo "# task_id\tsweep\tN\talpha\tt\tmethod\tvariant"
  for alpha in $ALPHAS_A; do
    for mv in "${MODELS[@]}"; do
      read -r method variant <<< "$mv"
      task_id="$(gen_task_id "A_alpha" "$N_FIXED" "$alpha" "$T_FIXED" "$method" "$variant")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "A_alpha" "$N_FIXED" "$alpha" "$T_FIXED" "$method" "$variant"
    done
  done
} > "$TASK_A"

{
  echo "# task_id\tsweep\tN\talpha\tt\tmethod\tvariant"
  for N in $NS_B; do
    for mv in "${MODELS[@]}"; do
      read -r method variant <<< "$mv"
      task_id="$(gen_task_id "B_N" "$N" "$ALPHA_FIXED" "$T_FIXED" "$method" "$variant")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "B_N" "$N" "$ALPHA_FIXED" "$T_FIXED" "$method" "$variant"
    done
  done
} > "$TASK_B"

# t sweep: split small vs large to avoid memory blowups.
{
  echo "# task_id\tsweep\tN\talpha\tt\tmethod\tvariant"
  for t in $TS_C; do
    # Heuristic threshold: >3M samples tends to become memory-heavy.
    if [[ "$t" -le 3000000 ]]; then
      for mv in "${MODELS[@]}"; do
        read -r method variant <<< "$mv"
        task_id="$(gen_task_id "C_t" "$N_FIXED" "$ALPHA_FIXED" "$t" "$method" "$variant")"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "C_t" "$N_FIXED" "$ALPHA_FIXED" "$t" "$method" "$variant"
      done
    fi
  done
} > "$TASK_C_SMALL"

{
  echo "# task_id\tsweep\tN\talpha\tt\tmethod\tvariant"
  for t in $TS_C; do
    if [[ "$t" -gt 3000000 ]]; then
      for mv in "${MODELS[@]}"; do
        read -r method variant <<< "$mv"
        task_id="$(gen_task_id "C_t" "$N_FIXED" "$ALPHA_FIXED" "$t" "$method" "$variant")"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "C_t" "$N_FIXED" "$ALPHA_FIXED" "$t" "$method" "$variant"
      done
    fi
  done
} > "$TASK_C_LARGE"

write_manifest_csv "$TASK_A" "$MANIFEST_DIR/manifest_A_alpha.csv"
write_manifest_csv "$TASK_B" "$MANIFEST_DIR/manifest_B_N.csv"
write_manifest_csv "$TASK_C_SMALL" "$MANIFEST_DIR/manifest_C_t_small.csv"
write_manifest_csv "$TASK_C_LARGE" "$MANIFEST_DIR/manifest_C_t_large.csv"

# Task counts (sanity check)
CNT_A="$(task_count "$TASK_A")"
CNT_B="$(task_count "$TASK_B")"
CNT_CS="$(task_count "$TASK_C_SMALL")"
CNT_CL="$(task_count "$TASK_C_LARGE")"
TOTAL_TASKS=$(( CNT_A + CNT_B + CNT_CS + CNT_CL ))
log "Task counts: A_alpha=$CNT_A, B_N=$CNT_B, C_t_small=$CNT_CS, C_t_large=$CNT_CL (total=$TOTAL_TASKS)"
if (( TOTAL_TASKS == 0 )); then
  warn "No tasks generated. Please check your grids (ALPHAS_A/NS_B/TS_C) and generator settings."
  exit 2
fi


# ----------------------------
# Execute sweeps (with retries)
# ----------------------------

run_with_retries() {
  local task_file="$1"
  local parallel="$2"
  local name="$3"

  # Attempt 0: normal run (resume supported)
  run_task_file_with_parallel "$task_file" "$parallel" 0 0

  local tmp_fail="$MANIFEST_DIR/${name}_FAILURES.tsv"
  collect_failures_from_task_file "$task_file" "$tmp_fail"
  local failed_now
  failed_now="$(task_count "$tmp_fail")"

  local attempt=1
  while (( failed_now > 0 && attempt <= MAX_RETRIES )); do
    warn "$name: retrying failures (attempt $attempt / $MAX_RETRIES, failed=$failed_now)"
    if (( RETRY_BACKOFF_SEC > 0 )); then
      sleep $(( RETRY_BACKOFF_SEC * attempt ))
    fi
    # Force rerun: delete previous outputs for these failed tasks.
    run_task_file_with_parallel "$tmp_fail" "$parallel" "$attempt" 1
    collect_failures_from_task_file "$task_file" "$tmp_fail"
    failed_now="$(task_count "$tmp_fail")"
    attempt=$((attempt + 1))
  done

  if (( failed_now > 0 )); then
    warn "$name: still has failures after retries: $failed_now tasks"
  else
    log "$name: all tasks OK"
  fi
}

if (( CNT_A > 0 )); then
  run_with_retries "$TASK_A" "$PAR_ALPHA" "A_alpha"
else
  warn "A_alpha: 0 tasks (skipped)"
fi

# For large-N points, reduce parallelism further to be safe.
if (( CNT_B > 0 )); then
  run_with_retries "$TASK_B" "$PAR_N" "B_N"
else
  warn "B_N: 0 tasks (skipped)"
fi

# t sweep: small t points can be parallel; large t points are serial by default.
if (( CNT_CS > 0 )); then
  run_with_retries "$TASK_C_SMALL" "$PAR_T_SMALL" "C_t_small"
else
  warn "C_t_small: 0 tasks (skipped)"
fi
if (( CNT_CL > 0 )); then
  run_with_retries "$TASK_C_LARGE" "$PAR_T_LARGE" "C_t_large"
else
  warn "C_t_large: 0 tasks (skipped)"
fi

# ----------------------------
# Post: merge CSVs + write failure report
# ----------------------------

mkdir -p "$MANIFEST_DIR/merged"
merge_csv_dir "$OUT_ROOT/A_alpha" "$MANIFEST_DIR/merged/A_alpha_merged.csv"
merge_csv_dir "$OUT_ROOT/B_N"     "$MANIFEST_DIR/merged/B_N_merged.csv"
merge_csv_dir "$OUT_ROOT/C_t"     "$MANIFEST_DIR/merged/C_t_merged.csv"

# Global failures (any task without .ok)
GLOBAL_FAIL="$MANIFEST_DIR/FAILURES.tsv"
{
  echo "# task_id\tsweep\tN\talpha\tt\tmethod\tvariant"
  for tf in "$TASK_A" "$TASK_B" "$TASK_C_SMALL" "$TASK_C_LARGE"; do
    while IFS=$'\t' read -r task_id sweep N alpha t method variant; do
      [[ -z "$task_id" ]] && continue
      [[ "$task_id" =~ ^# ]] && continue
      if [[ ! -f "$STATUS_ROOT/${task_id}.ok" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$task_id" "$sweep" "$N" "$alpha" "$t" "$method" "$variant"
      fi
    done < "$tf"
  done
} > "$GLOBAL_FAIL"

FAIL_N=0
if [[ -s "$GLOBAL_FAIL" ]]; then
  # subtract header
  FAIL_N=$(( $(wc -l < "$GLOBAL_FAIL" | tr -d '[:space:]') - 1 ))
  (( FAIL_N < 0 )) && FAIL_N=0
fi

{
  echo "# EXP-2 synthetic (Dim=2) runner output"
  echo "timestamp=$(date -Is)"
  echo "build_type=$BUILD_TYPE"
  echo "clean_temp=$CLEAN_TEMP"
  echo "threads=$THREADS"
  echo "repeats=$REPEATS"
  echo "seed=$SEED"
  echo "gen=$GEN"
  echo "gen_seed=$GEN_SEED"
  echo "rectgen_script=$RECTGEN_SCRIPT"
  echo "audit_pairs=$AUDIT_PAIRS"
  echo "audit_seed=$AUDIT_SEED"
  echo "budget_B(j_star)=$BUDGET_B"
  echo "w_small=$W_SMALL"
  echo "enum_cap=$ENUM_CAP"
  echo "run_signature=$RUN_SIGNATURE"
  echo "alphas_A=$ALPHAS_A"
  echo "N_fixed=$N_FIXED"
  echo "t_fixed=$T_FIXED"
  echo "alpha_fixed=$ALPHA_FIXED"
  echo "Ns_B=$NS_B"
  echo "ts_C=$TS_C"
  echo "max_parallel=$MAX_PARALLEL"
  echo "par_alpha=$PAR_ALPHA"
  echo "par_N=$PAR_N"
  echo "par_t_small=$PAR_T_SMALL"
  echo "par_t_large=$PAR_T_LARGE"
  echo "max_parallel_cap=$MAX_PARALLEL_CAP"
  echo "mem_per_task_gb=$MEM_PER_TASK_GB"
  echo "mem_reserve_gb=$MEM_RESERVE_GB"
  echo "max_retries=$MAX_RETRIES"
  echo "retry_backoff_sec=$RETRY_BACKOFF_SEC"
  echo "timeout_sec=$TIMEOUT_SEC"
  echo "failures=$FAIL_N"
  echo
  echo "Outputs:"
  echo "  - Per-task logs : $LOG_ROOT/<task_id>.log"
  echo "  - Per-task CSV  : $OUT_ROOT/<sweep>/<task_id>/run.csv"
  echo "  - Merged CSVs   : $MANIFEST_DIR/merged/*_merged.csv"
  echo "  - Failures list : $GLOBAL_FAIL"
} > "$MANIFEST_DIR/RUN_SUMMARY.txt"

# Sync temp -> results (always; even with failures, to preserve partial outputs)
rm -rf "$RESULT_ROOT"
mkdir -p "$(dirname "$RESULT_ROOT")"
cp -a "$TEMP_ROOT" "$RESULT_ROOT"

log "Synced results to: $RESULT_ROOT"
log "Merged CSVs:"
log "  $RESULT_ROOT/manifest/merged/A_alpha_merged.csv"
log "  $RESULT_ROOT/manifest/merged/B_N_merged.csv"
log "  $RESULT_ROOT/manifest/merged/C_t_merged.csv"

if (( FAIL_N > 0 )); then
  warn "There are still $FAIL_N failed tasks after retries. See: $RESULT_ROOT/manifest/FAILURES.tsv"
  exit 1
fi

log "All tasks completed successfully."

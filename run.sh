#!/usr/bin/env bash
# run.sh (repo root)
#
# EXP-2* (thesis §6.6 revised): Synthetic scalability in an SJS-favorable regime
# -----------------------------------------------------------------------------
# This script implements the *new* experiment plan in 《新实验计划.md》:
#   (A*) vary alpha_out   @ N=1e6,  t=1e7,   families F0+F1
#   (B*) vary N           @ alpha=200, t=1e7, family F0 (F1 optional)
#   (C*) vary t           @ N=1e6, alpha=200, family F0 (F1 optional)
#   (D*) family sensitivity @ N=1e6, alpha=200, t=1e7, families F0+F1
#   (E)  budget diagnosis @ N=1e6, alpha=200, t in {1e5,3e5}, sweep (B,w_small)
#   (F)  high-density alpha sweep (optional)
#
# Key changes vs the original run.sh:
#   - Throughput mode defaults: t=10,000,000; TH config: B=5e7, w_small=512
#   - Two data families (F0/F1) are supported by generating *preset* rectgen scripts
#     without modifying repo sources.
#   - Datasets are generated ONCE (via sjs_gen_dataset) and reused by all model runs
#     (dataset_source=binary). This saves massive wall-clock time.
#   - To keep the same CSV schema as synthetic runs, we patch each task's run.csv
#     to inject the dataset's gen_report_json.
#
# Usage:
#   chmod +x run.sh
#   ./run.sh
#
# Common overrides (env vars):
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   CLEAN_BUILD=0|1
#   CLEAN_TEMP=0|1
#   JOBS=8
#   THREADS=1
#   REPEATS=3
#   SEED=1
#   GEN_SEED=1
#   AUDIT_PAIRS=2000000
#   AUDIT_SEED=1
#
#   # Run subsets
#   RUN_A=1 RUN_B=1 RUN_C=1 RUN_D=1 RUN_E=1 RUN_F=1
#   INCLUDE_B_F1=0
#   INCLUDE_C_F1=0
#
#   # Framework knobs
#   B_TH=50000000
#   W_SMALL_TH=512
#   BUDGETS_E="10000000 20000000 50000000 100000000"
#   W_SMALLS_E="256 512 1024"
#
#   # Enum safety cap (0 disables; beware OOM for very dense joins)
#   ENUM_CAP=0
#
#   # Parallelism & robustness
#   MAX_PARALLEL=2
#   MAX_RETRIES=2
#   TIMEOUT_SEC=0
#
#   # Parameter grids
#   ALPHAS_A="0.1 0.3 1 3 10 30 100 300 1000"
#   NS_B="100000 200000 500000 1000000 2000000 5000000"
#   TS_C="1000000 3000000 10000000 30000000 100000000 300000000 1000000000"
#   ALPHAS_F="1000 3000 10000"
#
#   # Logging + retention
#   RUN_LOG_TO_STDOUT=0
#   RESULT_COPY_DATASETS=0
#   PRUNE_TEMP_DATASETS=1
#
# Notes:
#   - This repository build is Dim=2 only.
#   - F0/F1 datasets differ ONLY by Alacarte generator parameters:
#       F0: volume_dist=fixed,     volume_cv=0.25, shape_sigma=0.0
#       F1: volume_dist=lognormal, volume_cv=1.0,  shape_sigma=1.0
#   - If you are tight on RAM, consider lowering MAX_PARALLEL and/or setting ENUM_CAP.

set -Eeuo pipefail
IFS=$' \t\n'

timestamp_now() { date -Is; }

log()  { printf '[%s][run.sh][INFO] %s\n' "$(timestamp_now)" "$*"; }
warn() { printf '[%s][run.sh][WARN] %s\n' "$(timestamp_now)" "$*" >&2; }
fatal_handler() {
  local line="$1"
  local cmd="$2"
  printf '[%s][run.sh][FATAL] Failed at line %s: %s\n' "$(timestamp_now)" "$line" "$cmd" >&2
}
trap 'fatal_handler "${LINENO}" "${BASH_COMMAND}"' ERR

# ----------------------------
# Small helpers (adapted from the original run.sh)
# ----------------------------

task_count() {
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
    *) echo "$(echo "$1" | tr '[:upper:]' '[:lower:]')" ;;
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

  if [[ -n "$cache_home" && "$cache_home" != "$ROOT" ]]; then
    warn "Detected stale CMake cache source path: $cache_home (root is $ROOT)"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  if [[ -n "$cache_dir" && "$cache_dir" != "$BUILD_DIR" ]]; then
    warn "Detected stale CMake cache build path: $cache_dir (build is $BUILD_DIR)"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi
}

sanitize_token() {
  echo "$1" \
    | sed -e 's/-/m/g' -e 's/\./p/g' -e 's/[^A-Za-z0-9_]/_/g'
}

supports_wait_n() {
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
  while IFS= read -r f; do
    tail -n +2 "$f" >> "$out_csv" || true
  done < <(find "$dir" -type f -name run.csv | sort)
}

copy_dir_if_exists() {
  local src="$1"
  local dst="$2"
  [[ -d "$src" ]] || return 0
  cp -a "$src" "$dst"
}

dir_size_human() {
  local path="$1"
  [[ -e "$path" ]] || { echo "0"; return 0; }
  du -sh "$path" 2>/dev/null | awk '{print $1}'
}

# ----------------------------
# Resolve repo root
# ----------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------
# Parameters (defaults follow the new plan)
# ----------------------------
EXP_TAG="${EXP_TAG:-exp2_sjs_adv_d2}"

BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
CLEAN_TEMP="${CLEAN_TEMP:-0}"

THREADS="${THREADS:-1}"
REPEATS="${REPEATS:-3}"
SEED="${SEED:-1}"
GEN_SEED="${GEN_SEED:-1}"

# Generator
GEN="${GEN:-alacarte_rectgen}"
RECTGEN_BASE_SCRIPT="${RECTGEN_BASE_SCRIPT:-$ROOT/tools/alacarte_rectgen_generate.py}"
ALACARTE_MODULE="${ALACARTE_MODULE:-$ROOT/Alacarte/alacarte_rectgen.py}"
AUDIT_PAIRS="${AUDIT_PAIRS:-2000000}"
AUDIT_SEED="${AUDIT_SEED:-$GEN_SEED}"

# Families (Alacarte params)
F0_VOLUME_DIST="${F0_VOLUME_DIST:-fixed}"
F0_VOLUME_CV="${F0_VOLUME_CV:-0.25}"
F0_SHAPE_SIGMA="${F0_SHAPE_SIGMA:-0.0}"

F1_VOLUME_DIST="${F1_VOLUME_DIST:-lognormal}"
F1_VOLUME_CV="${F1_VOLUME_CV:-1.0}"
F1_SHAPE_SIGMA="${F1_SHAPE_SIGMA:-1.0}"

# Framework knobs
B_TH="${B_TH:-50000000}"          # TH config: B=5e7
W_SMALL_TH="${W_SMALL_TH:-512}"   # TH config: w_small=512
ENUM_CAP="${ENUM_CAP:-0}"

# Robustness
MAX_RETRIES="${MAX_RETRIES:-2}"
RETRY_BACKOFF_SEC="${RETRY_BACKOFF_SEC:-2}"
TIMEOUT_SEC="${TIMEOUT_SEC:-0}"

# Auto-parallelism (conservative defaults; adjust if you know your RAM headroom)
MAX_PARALLEL_CAP="${MAX_PARALLEL_CAP:-8}"
MEM_PER_TASK_GB="${MEM_PER_TASK_GB:-24}"
MEM_RESERVE_GB="${MEM_RESERVE_GB:-8}"

# Which sweeps to run
RUN_A="${RUN_A:-1}"
RUN_B="${RUN_B:-1}"
RUN_C="${RUN_C:-1}"
RUN_D="${RUN_D:-1}"
RUN_E="${RUN_E:-1}"
RUN_F="${RUN_F:-1}"

INCLUDE_B_F1="${INCLUDE_B_F1:-0}"
INCLUDE_C_F1="${INCLUDE_C_F1:-0}"

# Grids
ALPHAS_A="${ALPHAS_A:-0.1 0.3 1 3 10 30 100 300 1000}"
N_A="${N_A:-1000000}"
T_A="${T_A:-10000000}"

NS_B="${NS_B:-100000 200000 500000 1000000 2000000 5000000}"
ALPHA_B="${ALPHA_B:-200}"
T_B="${T_B:-10000000}"

TS_C="${TS_C:-1000000 3000000 10000000 30000000 100000000 300000000 1000000000}"
N_C="${N_C:-1000000}"
ALPHA_C="${ALPHA_C:-200}"

N_D="${N_D:-1000000}"
ALPHA_D="${ALPHA_D:-200}"
T_D="${T_D:-10000000}"

TS_E="${TS_E:-100000 300000}"
BUDGETS_E="${BUDGETS_E:-10000000 20000000 50000000 100000000}"
W_SMALLS_E="${W_SMALLS_E:-256 512 1024}"

N_F="${N_F:-1000000}"
T_F="${T_F:-1000000}"
ALPHAS_F="${ALPHAS_F:-1000 3000 10000}"

# Split thresholds (to keep OOM risk low)
N_LARGE_THRESHOLD="${N_LARGE_THRESHOLD:-2000000}"
T_SMALL_MAX="${T_SMALL_MAX:-10000000}"

# Logging + result retention
RUN_LOG_TO_STDOUT="${RUN_LOG_TO_STDOUT:-0}"
RESULT_COPY_DATASETS="${RESULT_COPY_DATASETS:-0}"
PRUNE_TEMP_DATASETS="${PRUNE_TEMP_DATASETS:-1}"

# Build parallelism
if [[ -n "${JOBS:-}" ]]; then
  JOBS="$JOBS"
else
  local_nproc="$(nproc_safe)"
  local_mem_gb="$(mem_gb_safe)"
  JOBS="$local_mem_gb"
  [[ -z "$JOBS" || "$JOBS" -lt 1 ]] && JOBS=1
  (( JOBS > local_nproc )) && JOBS="$local_nproc"
  (( JOBS > 16 )) && JOBS=16
fi

# ----------------------------
# Model sets
# ----------------------------
MODELS_MAIN=(
  "ours enum_sampling"
  "ours sampling"
  "ours adaptive"
  "range_tree enum_sampling"
  "range_tree sampling"
  "range_tree adaptive"
  "kd_tree sampling"
)

MODELS_E_VAR=(
  "ours sampling"
  "ours adaptive"
  "range_tree sampling"
  "range_tree adaptive"
)

MODELS_F=(
  "ours sampling"
  "ours adaptive"
  "range_tree sampling"
  "range_tree adaptive"
  "kd_tree sampling"
)

# ----------------------------
# Directories
# ----------------------------
BUILD_SUBDIR="$(build_subdir_for_type "$BUILD_TYPE")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"

TEMP_ROOT="${ROOT}/run/temp/${EXP_TAG}"
OUT_ROOT="${TEMP_ROOT}/out"
LOG_ROOT="${TEMP_ROOT}/logs"
STATUS_ROOT="${TEMP_ROOT}/status"
MANIFEST_DIR="${TEMP_ROOT}/manifest"

DATA_ROOT="${TEMP_ROOT}/datasets"
DATA_STATUS_ROOT="${TEMP_ROOT}/dataset_status"
DATA_LOG_ROOT="${TEMP_ROOT}/dataset_logs"
PRESET_DIR="${TEMP_ROOT}/rectgen_presets"

RESULT_ROOT="${ROOT}/results/raw/${EXP_TAG}"
RUN_LOG_DIR="${ROOT}/log"
RUN_LOG_STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG_FILE="${RUN_LOG_DIR}/run_${EXP_TAG}_${RUN_LOG_STAMP}_pid$$.log"

setup_run_logger() {
  mkdir -p "$RUN_LOG_DIR"
  if [[ "$RUN_LOG_TO_STDOUT" == "1" ]]; then
    exec > >(tee -a "$RUN_LOG_FILE") 2>&1
  else
    exec >>"$RUN_LOG_FILE" 2>&1
  fi
}

setup_run_logger
log "Run log initialized: $RUN_LOG_FILE"

if [[ "$CLEAN_TEMP" == "1" ]]; then
  warn "CLEAN_TEMP=1: removing temp dir: $TEMP_ROOT"
  rm -rf "$TEMP_ROOT"
fi

mkdir -p "$TEMP_ROOT" "$OUT_ROOT" "$LOG_ROOT" "$STATUS_ROOT" "$MANIFEST_DIR" \
         "$DATA_ROOT" "$DATA_STATUS_ROOT" "$DATA_LOG_ROOT" "$PRESET_DIR"

# ----------------------------
# Build
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
if [[ -z "$SJS_GEN" ]]; then
  echo "[run.sh][FATAL] Could not locate sjs_gen_dataset under $BUILD_DIR (needed for dataset caching)" >&2
  exit 2
fi

log "Using sjs_run        : $SJS_RUN"
log "Using sjs_gen_dataset: $SJS_GEN"

# ----------------------------
# Parallelism defaults
# ----------------------------
NPROC="$(nproc_safe)"
MEM_GB="$(mem_gb_safe)"

if [[ -n "${MAX_PARALLEL:-}" ]]; then
  MAX_PARALLEL_USER="$MAX_PARALLEL"
else
  cpu_per_task="$THREADS"
  (( cpu_per_task < 1 )) && cpu_per_task=1
  cpu_cap=$(( NPROC / cpu_per_task ))
  (( cpu_cap < 1 )) && cpu_cap=1
  (( cpu_cap > 1 )) && cpu_cap=$(( cpu_cap - 1 ))

  mem_budget=$(( MEM_GB - MEM_RESERVE_GB ))
  (( mem_budget < MEM_PER_TASK_GB )) && mem_budget="$MEM_PER_TASK_GB"
  mem_cap=$(( mem_budget / MEM_PER_TASK_GB ))
  (( mem_cap < 1 )) && mem_cap=1

  MAX_PARALLEL_USER="$cpu_cap"
  (( MAX_PARALLEL_USER > mem_cap )) && MAX_PARALLEL_USER="$mem_cap"
  (( MAX_PARALLEL_USER > MAX_PARALLEL_CAP )) && MAX_PARALLEL_USER="$MAX_PARALLEL_CAP"
fi
MAX_PARALLEL="$MAX_PARALLEL_USER"
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1

# Per-sweep parallelism (override with env if needed)
PAR_A="${PAR_A:-$MAX_PARALLEL}"
PAR_B_SMALL="${PAR_B_SMALL:-$(( MAX_PARALLEL > 2 ? 2 : MAX_PARALLEL ))}"
PAR_B_LARGE="${PAR_B_LARGE:-1}"
PAR_C_SMALL="${PAR_C_SMALL:-$(( MAX_PARALLEL > 3 ? 3 : MAX_PARALLEL ))}"
PAR_C_LARGE="${PAR_C_LARGE:-1}"
PAR_D="${PAR_D:-$MAX_PARALLEL}"
PAR_E="${PAR_E:-$MAX_PARALLEL}"
PAR_F="${PAR_F:-1}"

log "Host: nproc=$NPROC, mem~${MEM_GB}GB"
log "Parallelism: max=$MAX_PARALLEL (A=$PAR_A, B_small=$PAR_B_SMALL, B_large=$PAR_B_LARGE, C_small=$PAR_C_SMALL, C_large=$PAR_C_LARGE, D=$PAR_D, E=$PAR_E, F=$PAR_F)"
log "Auto-parallel knobs: cap=$MAX_PARALLEL_CAP, mem_per_task_gb=$MEM_PER_TASK_GB, mem_reserve_gb=$MEM_RESERVE_GB"

# ----------------------------
# Build run signature (for resumability)
# ----------------------------
build_run_signature() {
  {
    echo "exp_tag=$EXP_TAG"
    echo "build_type=$BUILD_TYPE"
    echo "threads=$THREADS"
    echo "repeats=$REPEATS"
    echo "seed=$SEED"
    echo "gen=$GEN"
    echo "gen_seed=$GEN_SEED"
    echo "audit_pairs=$AUDIT_PAIRS"
    echo "audit_seed=$AUDIT_SEED"
    echo "rectgen_base_script=$RECTGEN_BASE_SCRIPT"
    echo "alacarte_module=$ALACARTE_MODULE"

    echo "F0=$F0_VOLUME_DIST,$F0_VOLUME_CV,$F0_SHAPE_SIGMA"
    echo "F1=$F1_VOLUME_DIST,$F1_VOLUME_CV,$F1_SHAPE_SIGMA"

    echo "B_TH=$B_TH"
    echo "W_SMALL_TH=$W_SMALL_TH"
    echo "enum_cap=$ENUM_CAP"

    echo "RUN_A=$RUN_A RUN_B=$RUN_B RUN_C=$RUN_C RUN_D=$RUN_D RUN_E=$RUN_E RUN_F=$RUN_F"
    echo "INCLUDE_B_F1=$INCLUDE_B_F1 INCLUDE_C_F1=$INCLUDE_C_F1"

    echo "ALPHAS_A=$ALPHAS_A N_A=$N_A T_A=$T_A"
    echo "NS_B=$NS_B ALPHA_B=$ALPHA_B T_B=$T_B"
    echo "TS_C=$TS_C N_C=$N_C ALPHA_C=$ALPHA_C"
    echo "N_D=$N_D ALPHA_D=$ALPHA_D T_D=$T_D"
    echo "TS_E=$TS_E BUDGETS_E=$BUDGETS_E W_SMALLS_E=$W_SMALLS_E"
    echo "ALPHAS_F=$ALPHAS_F N_F=$N_F T_F=$T_F"

    echo "N_LARGE_THRESHOLD=$N_LARGE_THRESHOLD"
    echo "T_SMALL_MAX=$T_SMALL_MAX"
  } | cksum | awk '{print $1 "-" $2}'
}

RUN_SIGNATURE="$(build_run_signature)"
log "Run signature: $RUN_SIGNATURE"

# ----------------------------
# Prepare rectgen preset scripts for F0/F1
# ----------------------------

if [[ ! -f "$RECTGEN_BASE_SCRIPT" ]]; then
  echo "[run.sh][FATAL] RECTGEN_BASE_SCRIPT not found: $RECTGEN_BASE_SCRIPT" >&2
  exit 2
fi
if [[ ! -f "$ALACARTE_MODULE" ]]; then
  echo "[run.sh][FATAL] ALACARTE_MODULE not found: $ALACARTE_MODULE" >&2
  exit 2
fi

make_rectgen_preset() {
  local src="$1"; local dst="$2"; local vol_dist="$3"; local vol_cv="$4"; local shape_sigma="$5"; local tag="$6"

  python3 - "$src" "$dst" "$vol_dist" "$vol_cv" "$shape_sigma" "$tag" <<'PY'
import re
import sys
from pathlib import Path

src, dst, vol_dist, vol_cv, shape_sigma, tag = sys.argv[1:7]
text = Path(src).read_text(encoding='utf-8')

def sub_one(pattern, repl, name):
    global text
    new, n = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit(f"[preset][FATAL] Expected to patch 1 occurrence of {name}, got {n}")
    text = new

sub_one(r'^\s*volume_dist\s*=\s*"[^"]+"\s*$', f'    volume_dist = "{vol_dist}"  # preset:{tag}', 'volume_dist')
sub_one(r'^\s*volume_cv\s*=\s*[-+0-9.eE]+\s*$', f'    volume_cv = {vol_cv}  # preset:{tag}', 'volume_cv')
sub_one(r'^\s*shape_sigma\s*=\s*[-+0-9.eE]+\s*$', f'    shape_sigma = {shape_sigma}  # preset:{tag}', 'shape_sigma')

Path(dst).write_text(text, encoding='utf-8')
PY

  chmod +x "$dst"
}

RECTGEN_F0="$PRESET_DIR/alacarte_rectgen_F0.py"
RECTGEN_F1="$PRESET_DIR/alacarte_rectgen_F1.py"

make_rectgen_preset "$RECTGEN_BASE_SCRIPT" "$RECTGEN_F0" "$F0_VOLUME_DIST" "$F0_VOLUME_CV" "$F0_SHAPE_SIGMA" "F0"
make_rectgen_preset "$RECTGEN_BASE_SCRIPT" "$RECTGEN_F1" "$F1_VOLUME_DIST" "$F1_VOLUME_CV" "$F1_SHAPE_SIGMA" "F1"

log "Rectgen preset scripts:"
log "  F0: $RECTGEN_F0"
log "  F1: $RECTGEN_F1"

rectgen_for_family() {
  local fam="$1"
  case "$fam" in
    F0) echo "$RECTGEN_F0" ;;
    F1) echo "$RECTGEN_F1" ;;
    *) echo "$RECTGEN_F0" ;;
  esac
}

# ----------------------------
# Dataset cache helpers
# ----------------------------

ds_name_for() {
  local fam="$1"; local N="$2"; local alpha="$3"
  echo "d2_${fam}_n${N}_a$(sanitize_token "$alpha")_gs${GEN_SEED}"
}

ds_dir_for() {
  local fam="$1"; echo "$DATA_ROOT/$fam"
}

ds_paths_for() {
  local fam="$1"; local N="$2"; local alpha="$3"
  local name; name="$(ds_name_for "$fam" "$N" "$alpha")"
  local dir; dir="$(ds_dir_for "$fam")"
  echo "$dir/${name}_R.bin"$'\t'"$dir/${name}_S.bin"$'\t'"$dir/${name}_gen_report.json"$'\t'"$name"
}

patch_run_csv_gen_report() {
  local csv_file="$1"
  local gen_json_path="$2"
  if [[ ! -f "$csv_file" || ! -f "$gen_json_path" ]]; then
    return 0
  fi
  python3 - "$csv_file" "$gen_json_path" <<'PY'
import csv
import sys
from pathlib import Path

csv_path, gen_path = sys.argv[1], sys.argv[2]

gen_json = Path(gen_path).read_text(encoding='utf-8').strip()
if not gen_json:
    gen_json = "{}"

with open(csv_path, newline='') as f:
    reader = csv.reader(f)
    rows = list(reader)
if not rows:
    sys.exit(0)
header = rows[0]
try:
    idx = header.index('gen_report_json')
except ValueError:
    sys.exit(0)

out_rows = [header]
for row in rows[1:]:
    if len(row) < len(header):
        row = row + [''] * (len(header) - len(row))
    row[idx] = gen_json
    out_rows.append(row)

with open(csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(out_rows)
PY
}

ensure_dataset() {
  local fam="$1"; local N="$2"; local alpha="$3"

  local paths
  paths="$(ds_paths_for "$fam" "$N" "$alpha")"
  local path_r path_s rep name
  IFS=$'\t' read -r path_r path_s rep name <<< "$paths"

  mkdir -p "$(ds_dir_for "$fam")"

  local ok_flag="$DATA_STATUS_ROOT/${name}.ok"
  local meta_flag="$DATA_STATUS_ROOT/${name}.meta"
  local log_file="$DATA_LOG_ROOT/${name}.log"

  if [[ -f "$ok_flag" && -f "$path_r" && -f "$path_s" && -f "$rep" ]]; then
    if grep -Fxq "run_signature=$RUN_SIGNATURE" "$meta_flag" 2>/dev/null; then
      return 0
    fi
  fi

  local lock_dir="$DATA_STATUS_ROOT/${name}.lock"
  while ! mkdir "$lock_dir" 2>/dev/null; do
    if [[ -f "$ok_flag" && -f "$path_r" && -f "$path_s" && -f "$rep" ]]; then
      return 0
    fi
    sleep 0.5
  done

  rm -f "$ok_flag" "$meta_flag"

  {
    echo "run_signature=$RUN_SIGNATURE"
    echo "family=$fam"
    echo "N=$N"
    echo "alpha=$alpha"
    echo "gen_seed=$GEN_SEED"
    echo "gen=$GEN"
    echo "rectgen_script=$(rectgen_for_family "$fam")"
    echo "alacarte_module=$ALACARTE_MODULE"
    echo "audit_pairs=$AUDIT_PAIRS"
    echo "audit_seed=$AUDIT_SEED"
  } > "$meta_flag"

  log "[dataset] Generating $name (fam=$fam N=$N alpha=$alpha)"

  local -a gen_cmd=(
    "$SJS_GEN"
    "--dataset_source=synthetic"
    "--gen=$GEN"
    "--dataset=$name"
    "--dim=2"
    "--n_r=$N" "--n_s=$N"
    "--alpha=$alpha"
    "--gen_seed=$GEN_SEED"
    "--out_dir=$(ds_dir_for "$fam")"
    "--write_csv=0"
    "--rectgen_script=$(rectgen_for_family "$fam")"
    "--alacarte_module=$ALACARTE_MODULE"
    "--audit_pairs=$AUDIT_PAIRS"
    "--audit_seed=$AUDIT_SEED"
  )

  {
    echo "[$(timestamp_now)][run.sh][dataset] start name=$name fam=$fam N=$N alpha=$alpha"
    echo "[$(timestamp_now)][run.sh][dataset] cmd=${gen_cmd[*]}"
  } > "$log_file"

  local rc=0
  if (( TIMEOUT_SEC > 0 )) && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SEC" "${gen_cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  else
    "${gen_cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  fi
  echo "[$(timestamp_now)][run.sh][dataset] end rc=$rc name=$name" >> "$log_file"

  if [[ "$rc" -ne 0 ]]; then
    warn "[dataset] FAILED rc=$rc name=$name (see $log_file)"
    rmdir "$lock_dir" || true
    return 1
  fi

  if [[ ! -f "$path_r" || ! -f "$path_s" || ! -f "$rep" ]]; then
    warn "[dataset] FAILED: missing outputs for $name"
    warn "  R=$path_r"
    warn "  S=$path_s"
    warn "  rep=$rep"
    rmdir "$lock_dir" || true
    return 1
  fi

  touch "$ok_flag"
  rmdir "$lock_dir" || true
  return 0
}

# ----------------------------
# Task registry (for pre-generating datasets)
# ----------------------------

declare -A DATASET_SPEC
DATASET_ORDER=()

register_dataset() {
  local fam="$1"; local N="$2"; local alpha="$3"
  local name
  name="$(ds_name_for "$fam" "$N" "$alpha")"
  if [[ -z "${DATASET_SPEC[$name]+x}" ]]; then
    DATASET_SPEC[$name]="${fam}"$'\t'"${N}"$'\t'"${alpha}"
    DATASET_ORDER+=("$name")
  fi
}

# ----------------------------
# Task runner
# ----------------------------

run_one_task() {
  local task_id="$1"
  local sweep="$2"
  local fam="$3"
  local N="$4"
  local alpha="$5"
  local t="$6"
  local method="$7"
  local variant="$8"
  local j_star="$9"
  local w_small="${10}"
  local force="${11}"
  local attempt="${12:-0}"

  local ok_flag="$STATUS_ROOT/${task_id}.ok"
  local fail_flag="$STATUS_ROOT/${task_id}.fail"
  local exit_flag="$STATUS_ROOT/${task_id}.exit"
  local meta_flag="$STATUS_ROOT/${task_id}.meta"

  local out_dir="$OUT_ROOT/${sweep}/${task_id}"
  local log_file="$LOG_ROOT/${task_id}.log"
  local csv_file="$out_dir/run.csv"

  if [[ "$force" != "1" && -f "$ok_flag" && -f "$csv_file" && -f "$meta_flag" ]]; then
    if grep -Fxq "run_signature=$RUN_SIGNATURE" "$meta_flag" \
      && csv_all_ok "$csv_file" >/dev/null 2>&1; then
      echo 0 > "$exit_flag"
      return 0
    fi
  fi

  rm -f "$ok_flag" "$fail_flag" "$exit_flag" "$meta_flag"
  rm -rf "$out_dir"
  mkdir -p "$out_dir"

  if ! ensure_dataset "$fam" "$N" "$alpha"; then
    echo "dataset_generation_failed=1" > "$meta_flag"
    touch "$fail_flag"
    echo 1 > "$exit_flag"
    return 0
  fi

  local paths
  paths="$(ds_paths_for "$fam" "$N" "$alpha")"
  local path_r path_s rep_path ds_name
  IFS=$'\t' read -r path_r path_s rep_path ds_name <<< "$paths"

  local -a cmd=(
    "$SJS_RUN"
    "--dataset_source=binary"
    "--dataset=$ds_name"
    "--dim=2"
    "--path_r=$path_r" "--path_s=$path_s"
    "--method=$method"
    "--variant=$variant"
    "--t=$t"
    "--seed=$SEED"
    "--repeats=$REPEATS"
    "--threads=$THREADS"
    "--j_star=$j_star"
    "--enum_cap=$ENUM_CAP"
    "--write_samples=0"
    "--verify=0"
    "--out_dir=$out_dir"
    "--results_file=$csv_file"
    "--log_level=info"
    "--log_timestamp=1"
    "--w_small=$w_small"
  )

  {
    echo "run_signature=$RUN_SIGNATURE"
    echo "attempt=$attempt"
    echo "task_id=$task_id"
    echo "sweep=$sweep"
    echo "family=$fam"
    echo "dataset=$ds_name"
    echo "N=$N"
    echo "alpha=$alpha"
    echo "t=$t"
    echo "method=$method"
    echo "variant=$variant"
    echo "j_star=$j_star"
    echo "w_small=$w_small"
    echo "path_r=$path_r"
    echo "path_s=$path_s"
    echo "gen_report=$rep_path"
    echo "cmd=${cmd[*]}"
  } > "$meta_flag"

  {
    echo "[$(timestamp_now)][run.sh][task] start task_id=$task_id sweep=$sweep family=$fam N=$N alpha=$alpha t=$t method=$method variant=$variant attempt=$attempt"
    echo "[$(timestamp_now)][run.sh][task] cmd=${cmd[*]}"
  } > "$log_file"

  local rc=0
  if (( TIMEOUT_SEC > 0 )) && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SEC" "${cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  else
    "${cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  fi
  echo "[$(timestamp_now)][run.sh][task] end task_id=$task_id rc=$rc" >> "$log_file"
  echo "$rc" > "$exit_flag"

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

  patch_run_csv_gen_report "$csv_file" "$rep_path" || true

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

  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant j_star w_small; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    wait_for_slot "$parallel"
    (run_one_task "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant" "$j_star" "$w_small" "$force" "$attempt") &
  done < "$task_file"

  wait || true
}

collect_failures_from_task_file() {
  local task_file="$1"
  local out_fail_list="$2"
  : > "$out_fail_list"

  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant j_star w_small; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    if [[ ! -f "$STATUS_ROOT/${task_id}.ok" ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant" "$j_star" "$w_small" \
        >> "$out_fail_list"
    fi
  done < "$task_file"
}

write_manifest_csv() {
  local task_file="$1"
  local out_csv="$2"
  mkdir -p "$(dirname "$out_csv")"
  echo "task_id,sweep,family,N,alpha,t,method,variant,j_star,w_small,out_dir,log_file" > "$out_csv"
  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant j_star w_small; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    local out_dir="$OUT_ROOT/${sweep}/${task_id}"
    local log_file="$LOG_ROOT/${task_id}.log"
    echo "${task_id},${sweep},${fam},${N},${alpha},${t},${method},${variant},${j_star},${w_small},${out_dir},${log_file}" >> "$out_csv"
  done < "$task_file"
}

run_with_retries() {
  local task_file="$1"
  local parallel="$2"
  local name="$3"

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

# ----------------------------
# Generate task manifests (and register datasets)
# ----------------------------

gen_task_id() {
  local sweep="$1"; local fam="$2"; local N="$3"; local alpha="$4"; local t="$5"; local method="$6"; local variant="$7"; local j_star="$8"; local w_small="$9"
  echo "${sweep}__${fam}__n${N}__a$(sanitize_token "$alpha")__t${t}__${method}__${variant}__B${j_star}__w${w_small}"
}

log "Generating task manifests..."

TASK_A="$MANIFEST_DIR/tasks_A_alpha.tsv"
TASK_B_SMALL="$MANIFEST_DIR/tasks_B_N_small.tsv"
TASK_B_LARGE="$MANIFEST_DIR/tasks_B_N_large.tsv"
TASK_C_SMALL="$MANIFEST_DIR/tasks_C_t_small.tsv"
TASK_C_LARGE="$MANIFEST_DIR/tasks_C_t_large.tsv"
TASK_D="$MANIFEST_DIR/tasks_D_family.tsv"
TASK_E="$MANIFEST_DIR/tasks_E_budget.tsv"
TASK_F="$MANIFEST_DIR/tasks_F_high_alpha.tsv"

TASK_HEADER="# task_id\tsweep\tfamily\tN\talpha\tt\tmethod\tvariant\tj_star\tw_small"

: > "$TASK_A"; echo -e "$TASK_HEADER" > "$TASK_A"
: > "$TASK_B_SMALL"; echo -e "$TASK_HEADER" > "$TASK_B_SMALL"
: > "$TASK_B_LARGE"; echo -e "$TASK_HEADER" > "$TASK_B_LARGE"
: > "$TASK_C_SMALL"; echo -e "$TASK_HEADER" > "$TASK_C_SMALL"
: > "$TASK_C_LARGE"; echo -e "$TASK_HEADER" > "$TASK_C_LARGE"
: > "$TASK_D"; echo -e "$TASK_HEADER" > "$TASK_D"
: > "$TASK_E"; echo -e "$TASK_HEADER" > "$TASK_E"
: > "$TASK_F"; echo -e "$TASK_HEADER" > "$TASK_F"

# (A*) alpha sweep, F0+F1
if [[ "$RUN_A" == "1" ]]; then
  for alpha in $ALPHAS_A; do
    for fam in F0 F1; do
      register_dataset "$fam" "$N_A" "$alpha"
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="A_alpha_${fam}"
        task_id="$(gen_task_id "$sweep" "$fam" "$N_A" "$alpha" "$T_A" "$method" "$variant" "$B_TH" "$W_SMALL_TH")"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$task_id" "$sweep" "$fam" "$N_A" "$alpha" "$T_A" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
          >> "$TASK_A"
      done
    done
  done
fi

# (B*) N sweep, default F0; F1 optional
if [[ "$RUN_B" == "1" ]]; then
  families_b=("F0")
  if [[ "$INCLUDE_B_F1" == "1" ]]; then
    families_b+=("F1")
  fi

  for fam in "${families_b[@]}"; do
    for N in $NS_B; do
      register_dataset "$fam" "$N" "$ALPHA_B"
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="B_N_${fam}"
        task_id="$(gen_task_id "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant" "$B_TH" "$W_SMALL_TH")"
        if (( N >= N_LARGE_THRESHOLD )); then
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
            >> "$TASK_B_LARGE"
        else
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
            >> "$TASK_B_SMALL"
        fi
      done
    done
  done
fi

# (C*) t sweep, default F0; F1 optional
if [[ "$RUN_C" == "1" ]]; then
  families_c=("F0")
  if [[ "$INCLUDE_C_F1" == "1" ]]; then
    families_c+=("F1")
  fi

  for fam in "${families_c[@]}"; do
    register_dataset "$fam" "$N_C" "$ALPHA_C"
    for t in $TS_C; do
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="C_t_${fam}"
        task_id="$(gen_task_id "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant" "$B_TH" "$W_SMALL_TH")"
        if (( t <= T_SMALL_MAX )); then
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
            >> "$TASK_C_SMALL"
        else
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
            >> "$TASK_C_LARGE"
        fi
      done
    done
  done
fi

# (D*) family sensitivity (fixed point), F0+F1
if [[ "$RUN_D" == "1" ]]; then
  for fam in F0 F1; do
    register_dataset "$fam" "$N_D" "$ALPHA_D"
    for mv in "${MODELS_MAIN[@]}"; do
      read -r method variant <<< "$mv"
      sweep="D_family_${fam}"
      task_id="$(gen_task_id "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant" "$B_TH" "$W_SMALL_TH")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
        >> "$TASK_D"
    done
  done
fi

# (E) budget diagnosis (B, w_small sweep) @ small t; F0 only
if [[ "$RUN_E" == "1" ]]; then
  fam="F0"
  register_dataset "$fam" "$N_D" "$ALPHA_D"  # ensure base dataset exists

  for t in $TS_E; do
    for B in $BUDGETS_E; do
      for ws in $W_SMALLS_E; do
        for mv in "${MODELS_E_VAR[@]}"; do
          read -r method variant <<< "$mv"
          sweep="E_budget_${fam}"
          task_id="$(gen_task_id "$sweep" "$fam" "$N_D" "$ALPHA_D" "$t" "$method" "$variant" "$B" "$ws")"
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N_D" "$ALPHA_D" "$t" "$method" "$variant" "$B" "$ws" \
            >> "$TASK_E"
        done
      done
    done

    sweep="E_budget_${fam}"
    task_id="$(gen_task_id "$sweep" "$fam" "$N_D" "$ALPHA_D" "$t" "kd_tree" "sampling" "$B_TH" "$W_SMALL_TH")"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$task_id" "$sweep" "$fam" "$N_D" "$ALPHA_D" "$t" "kd_tree" "sampling" "$B_TH" "$W_SMALL_TH" \
      >> "$TASK_E"
  done
fi

# (F) high-density alpha sweep (optional), F0 only; limited model set
if [[ "$RUN_F" == "1" ]]; then
  fam="F0"
  for alpha in $ALPHAS_F; do
    register_dataset "$fam" "$N_F" "$alpha"
    for mv in "${MODELS_F[@]}"; do
      read -r method variant <<< "$mv"
      sweep="F_high_alpha_${fam}"
      task_id="$(gen_task_id "$sweep" "$fam" "$N_F" "$alpha" "$T_F" "$method" "$variant" "$B_TH" "$W_SMALL_TH")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N_F" "$alpha" "$T_F" "$method" "$variant" "$B_TH" "$W_SMALL_TH" \
        >> "$TASK_F"
    done
  done
fi

write_manifest_csv "$TASK_A" "$MANIFEST_DIR/manifest_A_alpha.csv"
write_manifest_csv "$TASK_B_SMALL" "$MANIFEST_DIR/manifest_B_N_small.csv"
write_manifest_csv "$TASK_B_LARGE" "$MANIFEST_DIR/manifest_B_N_large.csv"
write_manifest_csv "$TASK_C_SMALL" "$MANIFEST_DIR/manifest_C_t_small.csv"
write_manifest_csv "$TASK_C_LARGE" "$MANIFEST_DIR/manifest_C_t_large.csv"
write_manifest_csv "$TASK_D" "$MANIFEST_DIR/manifest_D_family.csv"
write_manifest_csv "$TASK_E" "$MANIFEST_DIR/manifest_E_budget.csv"
write_manifest_csv "$TASK_F" "$MANIFEST_DIR/manifest_F_high_alpha.csv"

CNT_A="$(task_count "$TASK_A")"
CNT_BS="$(task_count "$TASK_B_SMALL")"
CNT_BL="$(task_count "$TASK_B_LARGE")"
CNT_CS="$(task_count "$TASK_C_SMALL")"
CNT_CL="$(task_count "$TASK_C_LARGE")"
CNT_D="$(task_count "$TASK_D")"
CNT_E="$(task_count "$TASK_E")"
CNT_F="$(task_count "$TASK_F")"
TOTAL_TASKS=$(( CNT_A + CNT_BS + CNT_BL + CNT_CS + CNT_CL + CNT_D + CNT_E + CNT_F ))

log "Task counts: A=$CNT_A, B_small=$CNT_BS, B_large=$CNT_BL, C_small=$CNT_CS, C_large=$CNT_CL, D=$CNT_D, E=$CNT_E, F=$CNT_F (total=$TOTAL_TASKS)"
if (( TOTAL_TASKS == 0 )); then
  warn "No tasks generated. Check RUN_* flags and grids."
  exit 2
fi

# ----------------------------
# Pre-generate all required datasets (once)
# ----------------------------

log "Pre-generating datasets (unique count=${#DATASET_ORDER[@]})..."
for name in "${DATASET_ORDER[@]}"; do
  IFS=$'\t' read -r fam N alpha <<< "${DATASET_SPEC[$name]}"
  if ! ensure_dataset "$fam" "$N" "$alpha"; then
    warn "Dataset generation failed for $name. Inspect: $DATA_LOG_ROOT/${name}.log"
    exit 3
  fi
done
log "All datasets are ready."

# ----------------------------
# Execute sweeps
# ----------------------------

if (( CNT_A > 0 )); then
  run_with_retries "$TASK_A" "$PAR_A" "A_alpha"
fi
if (( CNT_BS > 0 )); then
  run_with_retries "$TASK_B_SMALL" "$PAR_B_SMALL" "B_N_small"
fi
if (( CNT_BL > 0 )); then
  run_with_retries "$TASK_B_LARGE" "$PAR_B_LARGE" "B_N_large"
fi
if (( CNT_CS > 0 )); then
  run_with_retries "$TASK_C_SMALL" "$PAR_C_SMALL" "C_t_small"
fi
if (( CNT_CL > 0 )); then
  run_with_retries "$TASK_C_LARGE" "$PAR_C_LARGE" "C_t_large"
fi
if (( CNT_D > 0 )); then
  run_with_retries "$TASK_D" "$PAR_D" "D_family"
fi
if (( CNT_E > 0 )); then
  run_with_retries "$TASK_E" "$PAR_E" "E_budget"
fi
if (( CNT_F > 0 )); then
  run_with_retries "$TASK_F" "$PAR_F" "F_high_alpha"
fi

# ----------------------------
# Post: merge CSVs + write failure report
# ----------------------------

MERGED_DIR="$MANIFEST_DIR/merged"
mkdir -p "$MERGED_DIR"

while IFS= read -r d; do
  sweep="$(basename "$d")"
  merge_csv_dir "$d" "$MERGED_DIR/${sweep}_merged.csv"
done < <(find "$OUT_ROOT" -mindepth 1 -maxdepth 1 -type d | sort)

GLOBAL_FAIL="$MANIFEST_DIR/FAILURES.tsv"
{
  echo -e "$TASK_HEADER"
  for tf in "$TASK_A" "$TASK_B_SMALL" "$TASK_B_LARGE" "$TASK_C_SMALL" "$TASK_C_LARGE" "$TASK_D" "$TASK_E" "$TASK_F"; do
    while IFS=$'\t' read -r task_id sweep fam N alpha t method variant j_star w_small; do
      [[ -z "$task_id" ]] && continue
      [[ "$task_id" =~ ^# ]] && continue
      if [[ ! -f "$STATUS_ROOT/${task_id}.ok" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant" "$j_star" "$w_small"
      fi
    done < "$tf"
  done
} > "$GLOBAL_FAIL"

FAIL_N=0
if [[ -s "$GLOBAL_FAIL" ]]; then
  FAIL_N=$(( $(wc -l < "$GLOBAL_FAIL" | tr -d '[:space:]') - 1 ))
  (( FAIL_N < 0 )) && FAIL_N=0
fi

{
  echo "# EXP-2* synthetic runner output (SJS-favorable plan)"
  echo "timestamp=$(date -Is)"
  echo "exp_tag=$EXP_TAG"
  echo "build_type=$BUILD_TYPE"
  echo "threads=$THREADS"
  echo "repeats=$REPEATS"
  echo "seed=$SEED"
  echo "gen=$GEN"
  echo "gen_seed=$GEN_SEED"
  echo "audit_pairs=$AUDIT_PAIRS"
  echo "audit_seed=$AUDIT_SEED"
  echo "rectgen_base_script=$RECTGEN_BASE_SCRIPT"
  echo "alacarte_module=$ALACARTE_MODULE"
  echo "rectgen_f0=$RECTGEN_F0"
  echo "rectgen_f1=$RECTGEN_F1"
  echo "F0=$F0_VOLUME_DIST,$F0_VOLUME_CV,$F0_SHAPE_SIGMA"
  echo "F1=$F1_VOLUME_DIST,$F1_VOLUME_CV,$F1_SHAPE_SIGMA"
  echo "B_TH=$B_TH"
  echo "W_SMALL_TH=$W_SMALL_TH"
  echo "enum_cap=$ENUM_CAP"
  echo "run_signature=$RUN_SIGNATURE"
  echo "run_log_file=$RUN_LOG_FILE"
  echo "result_copy_datasets=$RESULT_COPY_DATASETS"
  echo "prune_temp_datasets=$PRUNE_TEMP_DATASETS"
  echo
  echo "RUN_A=$RUN_A (ALPHAS_A=$ALPHAS_A N_A=$N_A T_A=$T_A)"
  echo "RUN_B=$RUN_B (NS_B=$NS_B ALPHA_B=$ALPHA_B T_B=$T_B INCLUDE_B_F1=$INCLUDE_B_F1)"
  echo "RUN_C=$RUN_C (TS_C=$TS_C N_C=$N_C ALPHA_C=$ALPHA_C INCLUDE_C_F1=$INCLUDE_C_F1)"
  echo "RUN_D=$RUN_D (N_D=$N_D ALPHA_D=$ALPHA_D T_D=$T_D)"
  echo "RUN_E=$RUN_E (TS_E=$TS_E BUDGETS_E=$BUDGETS_E W_SMALLS_E=$W_SMALLS_E)"
  echo "RUN_F=$RUN_F (ALPHAS_F=$ALPHAS_F N_F=$N_F T_F=$T_F)"
  echo
  echo "Parallelism: max=$MAX_PARALLEL (A=$PAR_A, B_small=$PAR_B_SMALL, B_large=$PAR_B_LARGE, C_small=$PAR_C_SMALL, C_large=$PAR_C_LARGE, D=$PAR_D, E=$PAR_E, F=$PAR_F)"
  echo "Failures=$FAIL_N"
  echo
  echo "Outputs (temp workspace):"
  echo "  - Root run log    : $RUN_LOG_FILE"
  echo "  - Per-task logs   : $LOG_ROOT/<task_id>.log"
  echo "  - Per-task CSV    : $OUT_ROOT/<sweep>/<task_id>/run.csv"
  echo "  - Merged CSVs     : $MERGED_DIR/*_merged.csv"
  echo "  - Failures list   : $GLOBAL_FAIL"
  echo "  - Dataset logs    : $DATA_LOG_ROOT/*.log"
  if [[ "$RESULT_COPY_DATASETS" == "1" ]]; then
    echo "  - Dataset binaries: $DATA_ROOT/<F0|F1>/*.bin, *_gen_report.json"
  else
    echo "  - Dataset binaries: [not copied to results by default]"
  fi
} > "$MANIFEST_DIR/RUN_SUMMARY.txt"

TMP_DATASET_SIZE="$(dir_size_human "$DATA_ROOT")"
TMP_OUT_SIZE="$(dir_size_human "$OUT_ROOT")"
TMP_TASK_LOG_SIZE="$(dir_size_human "$LOG_ROOT")"
TMP_DATASET_LOG_SIZE="$(dir_size_human "$DATA_LOG_ROOT")"

rm -rf "$RESULT_ROOT"
mkdir -p "$RESULT_ROOT"

copy_dir_if_exists "$MANIFEST_DIR" "$RESULT_ROOT/manifest"
copy_dir_if_exists "$OUT_ROOT" "$RESULT_ROOT/out"
copy_dir_if_exists "$LOG_ROOT" "$RESULT_ROOT/logs"
copy_dir_if_exists "$DATA_LOG_ROOT" "$RESULT_ROOT/dataset_logs"
copy_dir_if_exists "$STATUS_ROOT" "$RESULT_ROOT/status"
copy_dir_if_exists "$DATA_STATUS_ROOT" "$RESULT_ROOT/dataset_status"
copy_dir_if_exists "$PRESET_DIR" "$RESULT_ROOT/rectgen_presets"
cp -a "$RUN_LOG_FILE" "$RESULT_ROOT/run.log"

if [[ "$RESULT_COPY_DATASETS" == "1" ]]; then
  copy_dir_if_exists "$DATA_ROOT" "$RESULT_ROOT/datasets"
else
  log "Skipping dataset binaries in result snapshot (RESULT_COPY_DATASETS=0)."
fi

if [[ "$PRUNE_TEMP_DATASETS" == "1" ]]; then
  log "PRUNE_TEMP_DATASETS=1: removing temporary dataset binaries under $DATA_ROOT"
  rm -rf "$DATA_ROOT"
fi

RESULT_SIZE="$(dir_size_human "$RESULT_ROOT")"

log "Synced curated results to: $RESULT_ROOT"
log "Space usage: temp datasets=$TMP_DATASET_SIZE, temp out=$TMP_OUT_SIZE, temp task_logs=$TMP_TASK_LOG_SIZE, temp dataset_logs=$TMP_DATASET_LOG_SIZE, result_total=$RESULT_SIZE"
log "Merged CSVs under: $RESULT_ROOT/manifest/merged/"
log "Root run log: $RUN_LOG_FILE"

if (( FAIL_N > 0 )); then
  warn "There are still $FAIL_N failed tasks after retries. See: $RESULT_ROOT/manifest/FAILURES.tsv"
  exit 1
fi

log "All tasks completed successfully." 

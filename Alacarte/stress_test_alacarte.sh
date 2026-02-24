#!/usr/bin/env bash
# stress_test_alacarte.sh  (INDUSTRIAL / 最严苛摔打版)
#
# 目标：
#   1) 只保留“最严苛”一档压力测试（工业级门禁）
#   2) 将全部运行记录写入同目录 run.log（同时打印到终端）
#
# 用法：
#   bash stress_test_alacarte.sh
#
# 可选环境变量：
#   PYTHON=python3.10     # 指定 Python 可执行文件（默认自动探测 python3/python）
#
# 输出：
#   - run.log：完整运行日志（覆盖写；旧 run.log 会备份为 run.YYYYmmdd_HHMMSS.log）
#
# 退出码：
#   0  全部测试通过
#   !=0 任一断言失败或异常

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

: "${PYTHON:=}"

if [[ -z "${PYTHON}" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
  elif command -v python >/dev/null 2>&1; then
    PYTHON=python
  else
    echo "[FATAL] python3/python not found in PATH"
    exit 127
  fi
fi

if [[ ! -f "alacarte_rectgen.py" ]]; then
  echo "[FATAL] alacarte_rectgen.py not found in ${DIR}"
  exit 2
fi

# -----------------------------
# Logging (industrial gate)
# -----------------------------
LOG="${DIR}/run.log"
if [[ -f "${LOG}" ]]; then
  TS="$(date +%Y%m%d_%H%M%S)"
  cp "${LOG}" "${DIR}/run.${TS}.log" || true
fi
: > "${LOG}"

# Tee EVERYTHING (stdout+stderr) to run.log
exec > >(tee -a "${LOG}") 2>&1

START_ISO="$(date -Is)"
echo "============================================================"
echo "[INFO] alacarte_rectgen INDUSTRIAL stress test"
echo "[INFO] START=${START_ISO}"
echo "[INFO] DIR=${DIR}"
echo "[INFO] LOG=${LOG}"
echo "[INFO] PYTHON=${PYTHON}"
echo "[INFO] USER=$(whoami)  HOST=$(hostname)"
echo "============================================================"

# Helpful env snapshot for reproducibility
echo "[ENV] uname: $(uname -a)"
if command -v git >/dev/null 2>&1 && git -C "${DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "[ENV] git: $(git -C "${DIR}" rev-parse --short HEAD)  dirty=$(git -C "${DIR}" status --porcelain | wc -l | tr -d ' ') files"
fi
echo "[ENV] ulimit -a:"
ulimit -a || true
echo "============================================================"

# Stabilize thread usage (reproducibility + avoid BLAS oversubscription noise)
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

# Better error context in logs
trap 'echo "[ERROR] line ${LINENO}: command failed: ${BASH_COMMAND}"' ERR

"${PYTHON}" - <<'PY'
import os, sys, time, math, platform
import numpy as np
import alacarte_rectgen as ar

# -----------------------------
# Industrial configuration (single mode)
# -----------------------------
CFG = {
    # Coverage solver accuracy target (relative)
    "TUNE_TOL_REL": 0.005,          # 0.5% target for solver convergence
    "TUNE_MAX_ITER": 45,

    # Monte Carlo sizes
    "TUNE_SAMPLES_BASE": 1_200_000, # used for most cases
    "TUNE_SAMPLES_DEEP": 2_500_000, # used for hardest/critical cases
    "RECHECK_SAMPLES_BASE": 1_800_000,
    "RECHECK_SAMPLES_DEEP": 3_500_000,

    # Pair-sampling sanity checks
    "PAIR_MAX_PAIRS": 10_000_000,   # cap (avoid runaway)
    "PAIR_MIN_PAIRS": 500_000,      # minimum when we decide to do it
    "PAIR_TARGET_HITS": 2_000,      # aim for this many expected intersections (if feasible)

    # Large-scale stress
    "LARGE_N_TOTAL": 2_000_000,     # total N=|R|+|S|
    "LARGE_ALPHA_TARGETS": [0.1, 10.0, 1000.0, 250_000.0],  # includes near-saturation regime
}

np.seterr(divide="raise", invalid="raise", over="raise", under="ignore")

def log(msg: str) -> None:
    print(msg, flush=True)

def assert_true(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)

def rel_err(a: float, b: float) -> float:
    if b == 0:
        return abs(a - b)
    return abs(a - b) / abs(b)

def abs_err(a: float, b: float) -> float:
    return abs(a - b)

def tol_expected(alpha: float) -> float:
    """
    Industrial tolerance for expected-alpha checks.
    Mix relative + absolute to be meaningful at small alpha.
    """
    return max(alpha * 0.02, 0.002)  # 2% rel OR 0.002 abs (stricter than prior script)

def tol_recheck(alpha: float) -> float:
    """
    Slightly looser than primary check, but still tight.
    """
    return max(alpha * 0.03, 0.004)  # 3% rel OR 0.004 abs

def make_universe(d: int, kind: str) -> np.ndarray:
    """
    Produce a (d,2) universe for different stress regimes.
    """
    if kind == "unit":
        return np.array([[0.0, 1.0]] * d, dtype=np.float64)
    if kind == "anisotropic":
        # mix large/small spans
        spans = []
        for k in range(d):
            if k == 0:
                spans.append([0.0, 10.0])
            elif k == 1:
                spans.append([0.0, 1.0])
            else:
                spans.append([0.0, 0.1])
        return np.array(spans, dtype=np.float64)
    if kind == "shifted":
        spans = []
        for k in range(d):
            if k == 0:
                spans.append([-1000.0, 1000.0])
            elif k == 1:
                spans.append([-1.0, 3.0])
            else:
                spans.append([-7.0, 5.0])
        return np.array(spans, dtype=np.float64)
    if kind == "huge":
        return np.array([[0.0, 1e9]] * d, dtype=np.float64)
    if kind == "tiny":
        spans = []
        for k in range(d):
            if k == 0:
                spans.append([0.0, 1e-6])
            elif k == 1:
                spans.append([0.0, 2e-6])
            else:
                spans.append([0.0, 5e-7])
        return np.array(spans, dtype=np.float64)
    raise ValueError(f"unknown universe kind: {kind}")

def check_boxset_invariants(B: ar.BoxSet, name: str) -> None:
    assert_true(isinstance(B.lower, np.ndarray), f"{name}.lower not ndarray")
    assert_true(isinstance(B.upper, np.ndarray), f"{name}.upper not ndarray")
    assert_true(isinstance(B.universe, np.ndarray), f"{name}.universe not ndarray")
    n, d = B.lower.shape
    assert_true(B.upper.shape == (n, d), f"{name}.upper shape mismatch: {B.upper.shape} vs {(n,d)}")
    assert_true(B.universe.shape == (d,2), f"{name}.universe shape mismatch: {B.universe.shape} vs {(d,2)}")

    loU = B.universe[:,0].astype(np.float64)
    hiU = B.universe[:,1].astype(np.float64)
    spans = hiU - loU

    lower = B.lower.astype(np.float64)
    upper = B.upper.astype(np.float64)

    assert_true(np.all(np.isfinite(lower)), f"{name}.lower contains NaN/Inf")
    assert_true(np.all(np.isfinite(upper)), f"{name}.upper contains NaN/Inf")

    # inside universe (with tiny eps for float rounding)
    eps = 1e-9
    assert_true(np.all(lower >= loU - eps), f"{name}.lower out of universe (below)")
    assert_true(np.all(upper <= hiU + eps), f"{name}.upper out of universe (above)")

    lens = upper - lower
    assert_true(np.all(lens > 0), f"{name} has non-positive side length(s)")
    assert_true(np.all(lens < spans * (1.0 + 1e-9)), f"{name} has length exceeding universe span(s)")
    assert_true(np.all(lower < upper), f"{name} violates lower<upper")

def check_invalid_params() -> None:
    log("[TEST] invalid parameter handling (hard failures must raise)")
    nR, nS = 100, 200
    alpha_max = (nR*nS)/(nR+nS)
    try:
        ar.make_rectangles_R_S(nR=nR, nS=nS, alpha_out=alpha_max*1.01, d=2, seed=0, tune_samples=50_000)
        raise AssertionError("Expected ValueError for alpha_out > alpha_max, but none raised")
    except ValueError:
        pass

    try:
        ar.make_rectangles_R_S(nR=10, nS=10, alpha_out=1.0, d=2, universe=np.array([0,1,2]), seed=0, tune_samples=10_000)
        raise AssertionError("Expected ValueError for bad universe shape")
    except ValueError:
        pass

    try:
        ar.make_rectangles_R_S(nR=10, nS=10, alpha_out=1.0, d=2, universe=np.array([[0.0,0.0],[0.0,1.0]]), seed=0, tune_samples=10_000)
        raise AssertionError("Expected ValueError for bad universe bounds")
    except ValueError:
        pass

def check_interval_prob_formula() -> None:
    log("[TEST] interval_overlap_prob formula vs Monte Carlo (randomized 1D regression)")
    rng = np.random.default_rng(42)
    W = 1.0
    # multiple random (a,b) including edge regimes
    pairs = [(0.2,0.35),(0.9,0.05),(0.6,0.6),(1e-6,0.8),(0.999999,1e-6)]
    for _ in range(10):
        a = float(rng.uniform(1e-6, 0.999999))
        b = float(rng.uniform(1e-6, 0.999999))
        pairs.append((a,b))

    M = 1_500_000  # strong Monte Carlo
    rng2 = np.random.default_rng(123)
    for (a,b) in pairs:
        p_formula = float(ar.interval_overlap_prob(np.array([a]), np.array([b]), W)[0])
        x = rng2.random(M) * (W - a)
        y = rng2.random(M) * (W - b)
        hit = np.sum(np.maximum(x,y) < np.minimum(x+a, y+b))
        p_mc = hit / M
        err = abs(p_mc - p_formula)
        # Very tight tolerance thanks to large M
        assert_true(err < 0.003, f"interval_overlap_prob deviates too much: a={a} b={b} err={err}")
    log("    ok")

def check_reproducibility() -> None:
    log("[TEST] reproducibility (same seed => bitwise identical; different seed => different)")
    cfg = dict(nR=8000, nS=12000, alpha_out=10.0, d=2, volume_dist="fixed",
               shape_sigma=0.6, seed=123,
               tune_samples=CFG["TUNE_SAMPLES_BASE"], tune_tol_rel=CFG["TUNE_TOL_REL"], tune_max_iter=CFG["TUNE_MAX_ITER"])
    R1, S1, info1 = ar.make_rectangles_R_S(**cfg)
    R2, S2, info2 = ar.make_rectangles_R_S(**cfg)

    assert_true(info1["coverage"] == info2["coverage"], "coverage differs under same seed/config")
    assert_true(np.array_equal(R1.lower, R2.lower) and np.array_equal(R1.upper, R2.upper), "R differs under same seed")
    assert_true(np.array_equal(S1.lower, S2.lower) and np.array_equal(S1.upper, S2.upper), "S differs under same seed")

    cfg2 = dict(cfg)
    cfg2["seed"] = 124
    R3, S3, info3 = ar.make_rectangles_R_S(**cfg2)
    assert_true(not np.array_equal(R1.lower, R3.lower), "R unexpectedly identical under different seeds")
    assert_true(not np.array_equal(S1.lower, S3.lower), "S unexpectedly identical under different seeds")

def check_alpha_zero_behavior() -> None:
    log("[TEST] alpha_out=0 behavior (no crash; strict invariants hold; near-zero expected alpha)")
    R, S, info = ar.make_rectangles_R_S(
        nR=5000, nS=5000, alpha_out=0.0, d=2, seed=0,
        tune_samples=200_000, tune_tol_rel=0.10
    )
    check_boxset_invariants(R, "R_alpha0")
    check_boxset_invariants(S, "S_alpha0")
    assert_true(info["coverage"] == 0.0, "alpha=0 should return coverage=0.0")
    assert_true(info["alpha_expected_est"] < 1e-2, f"alpha=0 should yield tiny expected alpha, got {info['alpha_expected_est']}")

def estimate_pairs_budget(p_est: float) -> int:
    """
    Choose #pairs for pair-sampling sanity check.
    We aim for PAIR_TARGET_HITS expected hits but cap to PAIR_MAX_PAIRS.
    """
    if p_est <= 0:
        return 0
    target = int(math.ceil(CFG["PAIR_TARGET_HITS"] / p_est))
    target = max(target, CFG["PAIR_MIN_PAIRS"])
    target = min(target, CFG["PAIR_MAX_PAIRS"])
    return int(target)

def check_monotonicity_small_fixed() -> None:
    log("[TEST] monotonicity sanity (alpha_expected must increase with coverage; fixed lengths => deterministic)")
    # Choose deterministic setting: fixed volumes, shape_sigma=0 => lengths constant => no MC noise.
    nR, nS, d = 5000, 5000, 2
    U = make_universe(d, "unit")
    # Evaluate alpha for increasing coverages
    Cs = [0.01, 0.02, 0.05, 0.1, 0.2]
    alphas = []
    for C in Cs:
        a, p = ar.estimate_alpha_expected(
            nR=nR, nS=nS, alpha_out=0.0, d=d, universe=U,
            coverage=C, volume_dist="fixed", volume_cv=0.25, shape_sigma=0.0,
            num_samples=50_000, seed=777  # deterministic anyway
        )
        alphas.append(a)
    # Strictly increasing (allow tiny epsilon)
    for i in range(1, len(alphas)):
        assert_true(alphas[i] >= alphas[i-1] - 1e-9, f"non-monotone alpha_expected: C={Cs[i-1]}->{Cs[i]} {alphas[i-1]}->{alphas[i]}")
    log("    ok")

def solve_and_check(case_name: str, *, nR: int, nS: int, alpha: float, d: int,
                    universe_kind: str | None,
                    volume_dist: str, volume_cv: float, shape_sigma: float,
                    seed: int, deep: bool = False,
                    do_pair_check: bool = True) -> None:
    alpha_max = (nR*nS)/(nR+nS)
    assert_true(alpha <= alpha_max + 1e-12, f"{case_name}: alpha_target exceeds alpha_max")

    U = None if universe_kind is None else make_universe(d, universe_kind)

    tune_samples = CFG["TUNE_SAMPLES_DEEP"] if deep else CFG["TUNE_SAMPLES_BASE"]
    recheck_samples = CFG["RECHECK_SAMPLES_DEEP"] if deep else CFG["RECHECK_SAMPLES_BASE"]
    tune_tol = CFG["TUNE_TOL_REL"] if not deep else min(CFG["TUNE_TOL_REL"], 0.004)

    t0 = time.time()
    R, S, info = ar.make_rectangles_R_S(
        nR=nR, nS=nS, alpha_out=alpha,
        d=d, universe=U,
        volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        tune_samples=tune_samples, tune_tol_rel=tune_tol, tune_max_iter=CFG["TUNE_MAX_ITER"],
        seed=seed, dtype=np.float32
    )
    dt = time.time() - t0

    check_boxset_invariants(R, f"{case_name}.R")
    check_boxset_invariants(S, f"{case_name}.S")

    C = float(info["coverage"])

    # Primary expected-alpha check (recomputed with base recheck samples, not relying on info field)
    a1, p1 = ar.estimate_alpha_expected(
        nR=nR, nS=nS, alpha_out=alpha, d=d, universe=U,
        coverage=C, volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        num_samples=recheck_samples, seed=seed + 99991
    )

    err1 = abs_err(a1, alpha)
    tol1 = tol_expected(alpha)
    log(f"    [{case_name}] nR={nR} nS={nS} d={d} U={universe_kind or 'default'} dist={volume_dist} cv={volume_cv} shape_sigma={shape_sigma}")
    log(f"        alpha_target={alpha:.6g} alpha_max={alpha_max:.6g}")
    log(f"        solved coverage C={C:.6g}  expected_alpha(recheck)={a1:.6g}  abs_err={err1:.3g} tol={tol1:.3g}  p_est={p1:.6g}  time={dt:.2f}s")
    assert_true(math.isfinite(C) and C >= 0.0, f"{case_name}: coverage must be finite and >=0")
    assert_true(err1 <= tol1, f"{case_name}: expected alpha too far (abs_err={err1}, tol={tol1})")

    # Optional pair-sampling sanity check on realized sets (only if meaningful)
    if do_pair_check:
        # choose sample budget based on expected p1
        m_pairs = estimate_pairs_budget(p1)
        if m_pairs <= 0:
            log(f"        realized_check: skipped (p_est≈0)")
        else:
            # If expected hits is tiny, this check becomes noisy. We will *log* and only enforce if hits>=200.
            exp_hits = m_pairs * p1
            alpha_hat, p_hat = ar.estimate_alpha_by_pair_sampling(R, S, num_pairs=m_pairs, seed=seed + 12345)
            # predicted relative std for p_hat (small p approx): 1/sqrt(M*p)
            rel_std = 1.0 / math.sqrt(max(exp_hits, 1.0))
            # allow 6-sigma (very strict) + a tiny floor
            tol_rel = max(6.0 * rel_std, 0.02)
            # compare in alpha-space
            rerr = rel_err(alpha_hat, alpha) if alpha > 0 else abs(alpha_hat)
            log(f"        realized_check: pairs={m_pairs} exp_hits≈{exp_hits:.1f}  alpha_hat≈{alpha_hat:.6g}  rel_err≈{rerr:.3g}  tol_rel≈{tol_rel:.3g}")
            if exp_hits >= 200.0:
                assert_true(rerr <= tol_rel, f"{case_name}: realized alpha check failed (rel_err={rerr}, tol_rel={tol_rel}, exp_hits={exp_hits})")
            else:
                log("        realized_check: not enforced (expected hits < 200 => estimator too noisy for gate)")

def check_core_suite() -> None:
    log("[TEST] core industrial suite (diverse + harsh)")
    # Hand-picked hard cases (deterministic)
    cases = [
        # Small alpha / huge N regime (rare intersections) -> tests numeric stability
        dict(name="rare_p_unit_fixed",    nR=2000,  nS=2000,  d=2, U="unit",      dist="fixed",      cv=0.25, ss=0.0,  alpha=0.1,  deep=True,  pair=True),
        dict(name="rare_p_shift_logn",    nR=5000,  nS=3000,  d=2, U="shifted",   dist="lognormal",  cv=0.50, ss=0.6,  alpha=0.1,  deep=True,  pair=False),

        # Moderate alpha / skew sizes
        dict(name="skew_sizes_exp",       nR=800,   nS=12000, d=2, U="unit",      dist="exponential",cv=0.25, ss=0.6,  alpha=15.0, deep=False, pair=False),

        # Higher dimensions (stress d>2)
        dict(name="d3_aniso_normal",      nR=3000,  nS=3000,  d=3, U="anisotropic",dist="normal",   cv=0.25, ss=0.6,  alpha=20.0, deep=False, pair=False),
        dict(name="d5_unit_fixed",        nR=2500,  nS=2500,  d=5, U="unit",      dist="fixed",      cv=0.25, ss=0.0,  alpha=5.0,  deep=False, pair=False),
        dict(name="d8_unit_logn",         nR=2000,  nS=2000,  d=8, U="unit",      dist="lognormal",  cv=0.25, ss=1.8,  alpha=2.0,  deep=False, pair=False),

        # Extreme universes (huge/tiny) + aspect ratio extremes
        dict(name="huge_coords_fixed",    nR=4000,  nS=4000,  d=2, U="huge",      dist="fixed",      cv=0.25, ss=0.0,  alpha=80.0, deep=False, pair=True),
        dict(name="tiny_coords_normal",   nR=4000,  nS=4000,  d=2, U="tiny",      dist="normal",     cv=0.50, ss=1.8,  alpha=30.0, deep=False, pair=False),

        # Near saturation on small sets (tests high coverage, denominators near zero)
        dict(name="near_sat_fixed",       nR=1500,  nS=1500,  d=2, U="unit",      dist="fixed",      cv=0.25, ss=0.0,  alpha=0.98*((1500*1500)/(3000)), deep=True, pair=True),

        # High alpha but still within experiments range
        dict(name="alpha_1000_fixed",     nR=5000,  nS=5000,  d=2, U=None,        dist="fixed",      cv=0.25, ss=0.0,  alpha=1000.0, deep=True, pair=True),
    ]

    seed0 = 2025
    for i, c in enumerate(cases):
        solve_and_check(
            c["name"],
            nR=c["nR"], nS=c["nS"], alpha=float(c["alpha"]), d=c["d"],
            universe_kind=c["U"],
            volume_dist=c["dist"], volume_cv=float(c["cv"]), shape_sigma=float(c["ss"]),
            seed=seed0 + i*17,
            deep=bool(c["deep"]),
            do_pair_check=bool(c["pair"]),
        )

def check_fuzz_suite() -> None:
    log("[TEST] deterministic fuzz suite (random but reproducible configs)")
    rng = np.random.default_rng(2026)
    kinds = ["unit","anisotropic","shifted","huge","tiny"]
    d_choices = [1,2,3,5]  # keep moderate for runtime/memory
    dists = ["fixed","exponential","normal","lognormal"]

    for t in range(20):
        d = int(rng.choice(d_choices))
        nR = int(rng.integers(500, 6000))
        nS = int(rng.integers(500, 6000))
        alpha_max = (nR*nS)/(nR+nS)

        # sample alpha across wide range but keep feasible; bias to small/moderate
        r = float(rng.random())
        if r < 0.40:
            alpha = float(rng.uniform(0.0, min(1.0, alpha_max)))
        elif r < 0.85:
            alpha = float(rng.uniform(1.0, min(100.0, alpha_max)))
        else:
            alpha = float(rng.uniform(100.0, min(alpha_max, 800.0)))

        dist = str(rng.choice(dists))
        cv = float(rng.choice([0.10, 0.25, 0.50]))
        ss = float(rng.choice([0.0, 0.6, 1.8]))
        U = str(rng.choice(kinds))

        # For fuzz we still enforce the strong expected-alpha gate, but skip realized sampling to keep runtime controlled.
        solve_and_check(
            f"fuzz_{t:02d}",
            nR=nR, nS=nS, alpha=alpha, d=d,
            universe_kind=U,
            volume_dist=dist, volume_cv=cv, shape_sigma=ss,
            seed=7000 + t*31,
            deep=False,
            do_pair_check=False
        )

def check_large_scale() -> None:
    log("[TEST] large-scale generation (stress memory/time) + expected-alpha gate")
    N = int(CFG["LARGE_N_TOTAL"])
    nR = N // 2
    nS = N - nR
    d = 2
    U = None

    for alpha in CFG["LARGE_ALPHA_TARGETS"]:
        alpha_max = (nR*nS)/(nR+nS)
        if alpha > alpha_max:
            log(f"    skip alpha={alpha} (exceeds alpha_max={alpha_max})")
            continue

        # Large scale uses deep solver but moderate recheck to keep runtime bounded.
        solve_and_check(
            f"large_N{N}_alpha{alpha:g}",
            nR=nR, nS=nS, alpha=float(alpha), d=d,
            universe_kind=None,
            volume_dist="fixed", volume_cv=0.25, shape_sigma=0.0,
            seed=9000 + int(alpha) % 997,
            deep=True,
            do_pair_check=(alpha >= 1000.0)  # only meaningful when p is not tiny
        )

def main() -> None:
    log("------------------------------------------------------------")
    log(f"[ENV] python={sys.version.split()[0]} numpy={np.__version__} platform={platform.platform()}")
    log("------------------------------------------------------------")

    t0 = time.time()

    check_invalid_params()
    check_interval_prob_formula()
    check_reproducibility()
    check_alpha_zero_behavior()
    check_monotonicity_small_fixed()
    check_core_suite()
    check_fuzz_suite()
    check_large_scale()

    dt = time.time() - t0
    log("------------------------------------------------------------")
    log(f"[PASS] ALL INDUSTRIAL TESTS PASSED ✅  (total_time={dt:.2f}s)")
    log("------------------------------------------------------------")

if __name__ == "__main__":
    main()
PY

END_ISO="$(date -Is)"
echo "============================================================"
echo "[PASS] bash script completed successfully ✅"
echo "[INFO] END=${END_ISO}"
echo "============================================================"

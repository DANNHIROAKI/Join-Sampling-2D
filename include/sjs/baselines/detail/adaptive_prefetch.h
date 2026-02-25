#pragma once
// sjs/baselines/detail/adaptive_prefetch.h
//
// Helper utilities for Framework III (Adaptive) baselines (SJS v3 §3.5).
//
// This header is intentionally lightweight and header-only:
//  - Poisson tail score approximation for slot valuation
//  - A min-heap container for budgeted prefetch slots
//
// IMPORTANT: The scoring only affects performance (how many samples we prefetch
// in pass 1). It must be monotone decreasing in r to preserve the “prefix
// per event” invariant, but it does NOT affect distribution correctness.

#include "sjs/core/assert.h"
#include "sjs/core/types.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>

namespace sjs {
namespace baselines {
namespace detail {

// Estimate total join size W during a streaming sweep.
// SJS v3 suggests extrapolating from the prefix average.
inline long double EstimateTotalW(u32 total_events, u32 i, long double w_sofar) noexcept {
  if (total_events == 0 || i == 0) return w_sofar;
  // Simple extrapolation: W_hat = W_sofar * (M / i).
  const long double M = static_cast<long double>(total_events);
  const long double ii = static_cast<long double>(i);
  return (ii > 0.0L) ? (w_sofar * (M / ii)) : w_sofar;
}

// Poisson survival function approximation: Pr(Poisson(mu) >= r).
// For small mu we compute an exact CDF prefix; for large mu we use a normal
// approximation with continuity correction.
inline double PoissonSurvival(double mu, u32 r) noexcept {
  if (r <= 0) return 1.0;
  if (!(mu > 0.0) || !std::isfinite(mu)) return 0.0;

  // Small-mu exact prefix sum.
  if (mu <= 50.0 && r <= 200U) {
    const double p0 = std::exp(-mu);
    double p = p0;
    double cdf = p;
    for (u32 k = 1; k < r; ++k) {
      p *= mu / static_cast<double>(k);
      cdf += p;
    }
    double sf = 1.0 - cdf;
    if (sf < 0.0) sf = 0.0;
    if (sf > 1.0) sf = 1.0;
    return sf;
  }

  // Normal approximation: X ~ N(mu, mu). Use continuity correction.
  const double sigma = std::sqrt(mu);
  if (!(sigma > 0.0) || !std::isfinite(sigma)) return 0.0;
  const double z = (static_cast<double>(r) - 0.5 - mu) / sigma;
  // survival = 0.5 * erfc(z / sqrt(2))
  const double sf = 0.5 * std::erfc(z * 0.70710678118654752440);  // 1/sqrt(2)
  if (sf < 0.0) return 0.0;
  if (sf > 1.0) return 1.0;
  return sf;
}

// Slot value score for event i and slot index r (1-based) under the Poisson
// approximation described in SJS v3 §3.5.3.
inline double SlotScorePoisson(u64 w_i, u32 total_events, u32 i, long double w_sofar, u64 t, u32 r) noexcept {
  if (r == 0) return 1.0;
  if (w_i == 0 || t == 0) return 0.0;
  const long double W_hat = EstimateTotalW(total_events, i, w_sofar);
  if (!(W_hat > 0.0L)) return 0.0;
  const long double p_hat = static_cast<long double>(w_i) / W_hat;
  const long double mu_hat = static_cast<long double>(t) * p_hat;
  const double mu = static_cast<double>(mu_hat);
  return PoissonSurvival(mu, r);
}

struct PrefetchSlot {
  double score = 0.0;
  u32 sid = 0;
};

struct PrefetchSlotGreater {
  bool operator()(const PrefetchSlot& a, const PrefetchSlot& b) const noexcept {
    // min-heap by score
    if (a.score != b.score) return a.score > b.score;
    return a.sid > b.sid;
  }
};

// A min-heap of prefetch slots with capacity. Each pushed element corresponds
// to one unit of sample-cache budget.
class PrefetchHeap {
 public:
  void Clear() {
    while (!h_.empty()) h_.pop();
  }

  usize Size() const noexcept { return h_.size(); }
  bool Empty() const noexcept { return h_.empty(); }

  double MinScore() const noexcept {
    if (h_.empty()) return -std::numeric_limits<double>::infinity();
    return h_.top().score;
  }

  const PrefetchSlot& Top() const noexcept { return h_.top(); }

  void Push(PrefetchSlot s) { h_.push(s); }
  PrefetchSlot PopMin() {
    SJS_DASSERT(!h_.empty());
    const PrefetchSlot t = h_.top();
    h_.pop();
    return t;
  }

 private:
  std::priority_queue<PrefetchSlot, std::vector<PrefetchSlot>, PrefetchSlotGreater> h_;
};

}  // namespace detail
}  // namespace baselines
}  // namespace sjs

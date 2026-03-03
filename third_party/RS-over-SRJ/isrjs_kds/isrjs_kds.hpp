#pragma once
// isrjs_kds/isrjs_kds.hpp (modified)
//
// This file adapts the RS-over-SRJ KDS range-sampling baseline to the
// Join-Sampling-2D repository:
//   - point range join  ->  rectangle intersection join
//   - 2D KD-tree         ->  4D KD-tree on the standard embedding
//   - global datasets     ->  sjs::Dataset / sjs::IBaseline API
//   - per-sample query    ->  per-r batching (slot grouping) for performance
//
// Geometry / semantics alignment
// ------------------------------
// Dataset rectangles are half-open: [lo, hi)
// Intersection predicate is strict overlap (touching edges is NOT intersection):
//   r ∩ s != ∅  <=>  L(r) < R(s) && L(s) < R(r)  (per axis)
//
// Embedding for S rectangles:
//   p(s) = (s.lo.x, s.lo.y, s.hi.x, s.hi.y) ∈ R^4
// For a fixed r, s intersects r iff:
//   s.lo.x < r.hi.x
//   s.lo.y < r.hi.y
//   s.hi.x > r.lo.x
//   s.hi.y > r.lo.y
//
// Our KDTree implementation uses CLOSED ranges (<=, >=). We therefore convert
// strict inequalities using nextafter():
//   x < a  <=>  x <= nextDown(a)
//   x > a  <=>  x >= nextUp(a)
// (for floating-point T).

#include "../utils/kdtree.hpp"
#include "../utils/weighted_sampling.hpp"  // alias wrapper (uses sjs::AliasTable + sjs::Rng)

#include "sjs/baselines/baseline_api.h"
#include "sjs/join/join_enumerator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace rs_over_srj {

namespace detail {

// Wrapper enumerator: deterministic plane sweep on axis 0.
// This is only used for optional verification / enumerator API.

template <int Dim, class T>
class PlaneSweepEnumerator final : public IJoinEnumerator {
 public:
  PlaneSweepEnumerator(const Relation<Dim, T>* R, const Relation<Dim, T>* S, int axis = 0)
      : r_(R), s_(S) {
    join::PlaneSweepOptions opt;
    opt.axis = axis;
    opt.side_order = join::SideTieBreak::RBeforeS;
    stream_ = std::make_unique<join::PlaneSweepJoinStream<Dim, T>>(*r_, *s_, opt);
  }

  void Reset() override { stream_->Reset(); }
  bool Next(PairId* out) override { return stream_->Next(out); }
  const join::JoinStats& Stats() const noexcept override { return stream_->Stats(); }

 private:
  const Relation<Dim, T>* r_;
  const Relation<Dim, T>* s_;
  std::unique_ptr<join::PlaneSweepJoinStream<Dim, T>> stream_;
};

// nextafter helpers (for strict -> closed conversion)

template <class T>
inline T NextUp(T x) {
  if (!std::isfinite(x)) return x;
  return std::nextafter(x, std::numeric_limits<T>::infinity());
}

template <class T>
inline T NextDown(T x) {
  if (!std::isfinite(x)) return x;
  return std::nextafter(x, -std::numeric_limits<T>::infinity());
}

}  // namespace detail

// --------------------------
// RSOverSRJKDTreeSamplingBaseline
// --------------------------

template <int Dim, class T = Scalar>
class RSOverSRJKDTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RSOverSRJKDTreeSamplingBaseline is currently specialized for Dim==2 (embedding is 4D)");

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int D = 2 * Dim;  // embedding dimension (4)

  Method method() const noexcept override { return Method::RSOverSRJ; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "kd_tree_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    w_r_.clear();
    alias_r_.Clear();
    tree_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)err;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    // Build KD-tree on embedded points of S.
    {
      auto _ = phases ? phases->Scoped("build_points") : PhaseRecorder::ScopedPhase(nullptr, "");

      std::vector<KDPoint<D, T>> pts;
      pts.reserve(ds.S.Size());

      for (usize i = 0; i < ds.S.boxes.size(); ++i) {
        const auto& b = ds.S.boxes[i];
        KDPoint<D, T> p;
        p.coord[0] = b.lo.v[0];
        p.coord[1] = b.lo.v[1];
        p.coord[2] = b.hi.v[0];
        p.coord[3] = b.hi.v[1];
        p.idx = static_cast<std::uint32_t>(i);  // store S index as payload
        pts.push_back(p);
      }

      tree_.Build(std::move(pts));
    }

    w_r_.assign(ds.R.Size(), 0ULL);
    built_ = true;
    weights_valid_ = false;
    W_ = 0;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic

    if (!built_ || !ds_) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (ds_->R.Size() == 0 || ds_->S.Size() == 0) {
      std::fill(w_r_.begin(), w_r_.end(), 0ULL);
      W_ = 0;
      weights_valid_ = true;
      alias_r_.Clear();
      *out = MakeExactCount(0);
      return true;
    }

    const T inf = std::numeric_limits<T>::infinity();

    u64 W = 0;
    for (usize i = 0; i < ds_->R.boxes.size(); ++i) {
      const auto& r = ds_->R.boxes[i];

      // Closed 4D query (strict overlap converted via nextafter)
      std::array<T, D> qlo;
      std::array<T, D> qhi;
      qlo[0] = -inf;
      qlo[1] = -inf;
      qhi[0] = detail::NextDown<T>(r.hi.v[0]);
      qhi[1] = detail::NextDown<T>(r.hi.v[1]);

      qlo[2] = detail::NextUp<T>(r.lo.v[0]);
      qlo[3] = detail::NextUp<T>(r.lo.v[1]);
      qhi[2] = inf;
      qhi[3] = inf;

      const u64 w = static_cast<u64>(tree_.Count(qlo, qhi));
      w_r_[i] = w;

      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }
    }

    W_ = W;
    weights_valid_ = true;

    if (!alias_r_.Build(w_r_, err)) {
      if (err && err->empty()) *err = "RSOverSRJKDTreeSamplingBaseline::Count: failed to build alias";
      return false;
    }

    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: null rng/out";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: t too large for u32";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    // Ensure weights are ready.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    // --------------------------
    // Phase 2: plan slots (draw ridx for each output slot)
    // --------------------------
    struct SlotAssign {
      u32 ridx;
      u32 slot;
    };

    std::vector<SlotAssign> asg;
    asg.reserve(static_cast<usize>(t));

    {
      auto _ = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u32 slot = 0; slot < t; ++slot) {
        const u32 ridx = static_cast<u32>(alias_r_.Sample(rng));
        asg.push_back(SlotAssign{ridx, slot});
      }

      std::sort(asg.begin(), asg.end(), [](const SlotAssign& a, const SlotAssign& b) {
        if (a.ridx < b.ridx) return true;
        if (b.ridx < a.ridx) return false;
        return a.slot < b.slot;
      });
    }

    out->pairs.assign(static_cast<usize>(t), PairId{});

    // --------------------------
    // Phase 3: fill slots by batching per ridx
    // --------------------------
    {
      auto _ = phases ? phases->Scoped("phase3_fill") : PhaseRecorder::ScopedPhase(nullptr, "");

      const T inf = std::numeric_limits<T>::infinity();

      std::vector<typename KDTree<D, T>::SearchItem> items;
      items.reserve(256);

      std::vector<u64> weights;
      weights.reserve(256);

      usize ptr = 0;
      while (ptr < asg.size()) {
        const u32 ridx = asg[ptr].ridx;
        usize end = ptr;
        while (end < asg.size() && asg[end].ridx == ridx) ++end;
        const u32 k = static_cast<u32>(end - ptr);

        // Build query for this r.
        const auto& r = ds_->R.boxes[static_cast<usize>(ridx)];
        std::array<T, D> qlo;
        std::array<T, D> qhi;
        qlo[0] = -inf;
        qlo[1] = -inf;
        qhi[0] = detail::NextDown<T>(r.hi.v[0]);
        qhi[1] = detail::NextDown<T>(r.hi.v[1]);
        qlo[2] = detail::NextUp<T>(r.lo.v[0]);
        qlo[3] = detail::NextUp<T>(r.lo.v[1]);
        qhi[2] = inf;
        qhi[3] = inf;

        // Decompose query result set into blocks + single points.
        tree_.Search(qlo, qhi, &items);

        if (items.empty()) {
          // This should not happen if Count() and Search() are consistent.
          if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: empty query result for a sampled r";
          return false;
        }

        // Build alias over blocks.
        weights.clear();
        weights.reserve(items.size());
        for (const auto& it : items) {
          if (it.fully_contained) {
            const auto& nd = tree_.GetNode(it.node);
            weights.push_back(static_cast<u64>(nd.Size()));
          } else {
            weights.push_back(1ULL);
          }
        }

        alias alias_q;
        if (!alias_q.Build(weights, err)) {
          if (err && err->empty()) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: failed to build per-query alias";
          return false;
        }

        const Id rid = ds_->R.GetId(static_cast<usize>(ridx));

        // Draw k conditional samples from S(r).
        for (u32 j = 0; j < k; ++j) {
          const u32 slot = asg[ptr + j].slot;

          const usize ii = alias_q.Sample(rng);
          if (ii >= items.size()) {
            if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Sample: alias sampled out of range";
            return false;
          }

          const auto& it = items[ii];
          std::uint32_t s_idx = 0;

          if (it.fully_contained) {
            const auto& nd = tree_.GetNode(it.node);
            const u64 off = rng->UniformU64(static_cast<u64>(nd.Size()));
            const std::uint32_t pos = nd.l + static_cast<std::uint32_t>(off);
            s_idx = tree_.Points()[static_cast<usize>(pos)].idx;
          } else {
            s_idx = tree_.NodePointPayload(it.node);
          }

          out->pairs[static_cast<usize>(slot)] = PairId{rid, ds_->S.GetId(static_cast<usize>(s_idx))};
        }

        ptr = end;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;

    if (!built_ || !ds_) {
      if (err) *err = "RSOverSRJKDTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::PlaneSweepEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0);
  }

 private:
    const DatasetT* ds_ = nullptr;
    bool built_ = false;

    KDTree<D, T> tree_;

    std::vector<u64> w_r_;
    alias alias_r_;
    u64 W_ = 0;
    bool weights_valid_ = false;
};

}  // namespace rs_over_srj
}  // namespace baselines
}  // namespace sjs

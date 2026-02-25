#pragma once
// sjs/baselines/range_tree/enum_sampling.h
//
// RangeTree baseline (Variant::EnumSampling) — RT-Enumerate (SJS v3 Framework I).
//
// Framework I:
//   1) Enumerate *all* join pairs once and materialize them into memory.
//   2) Draw t i.i.d. uniform samples with replacement by indexing that array.
//
// Notes
// -----
// * This is a correctness / constant-factor reference. It is intentionally
//   Theta(|J|) in time and space and may be infeasible when |J| is huge.
// * Enumeration is deterministic and reuses the same sweep + dynamic range-tree
//   machinery as the Sampling/Adaptive variants.

#include "sjs/baselines/range_tree/sampling.h"

#include "sjs/baselines/detail/vector_join_enumerator.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace range_tree {

template <int Dim, class T = Scalar>
class RangeTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeEnumSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "range_tree_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    pairs_cached_ = false;
    pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)err;
    Reset();
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    ds_ = &ds;
    built_ = true;
    return true;
  }

  // Exact |J| for Framework I.
  //
  // For strict SJS v3 alignment (single materialization), we materialize all
  // join pairs here and cache them for Sample().
  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic

    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    if (!pairs_cached_) {
      auto scoped = phases ? phases->Scoped("phase1_enumerate_materialize")
                           : PhaseRecorder::ScopedPhase(nullptr, "");

      pairs_.clear();
      const u64 cap = cfg.run.enum_cap;

      detail::RangeJoinEnumerator<Dim, T> stream(&ds_->R, &ds_->S, /*axis=*/0, join::SideTieBreak::RBeforeS);
      PairId p;
      while (stream.Next(&p)) {
        pairs_.push_back(p);
        if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
          pairs_.clear();
          pairs_cached_ = false;
          if (err) *err = "RangeTreeEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
          return false;
        }
      }

      pairs_cached_ = true;
    }

    if (out) *out = MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  // Draw samples by uniform indexing into the cached join vector.
  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    if (!pairs_cached_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    const u64 W = static_cast<u64>(pairs_.size());
    if (W == 0) return true;

    {
      auto scoped = phases ? phases->Scoped("phase2_resample") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u32 i = 0; i < t; ++i) {
        const u64 j = rng->UniformU64(W);
        out->pairs[static_cast<usize>(i)] = pairs_[static_cast<usize>(j)];
      }
    }

    return true;
  }

  // Provide a deterministic enumerator.
  //
  // If we already materialized the join (Count ran), return an enumerator over
  // the cached vector to avoid re-running a sweep.
  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    if (pairs_cached_) {
      return std::make_unique<baselines::detail::VectorJoinEnumerator>(&pairs_);
    }
    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  bool pairs_cached_ = false;
  std::vector<PairId> pairs_;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs

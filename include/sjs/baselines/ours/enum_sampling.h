#pragma once
// sjs/baselines/ours/enum_sampling.h
//
// "Our Method" — Enumerate+Sampling variant (SJS v3 Framework I).
//
// Framework I:
//   1) Enumerate *all* join pairs once and materialize them into memory.
//   2) Draw t i.i.d. uniform samples with replacement by indexing that array.
//
// Memory: Theta(|J|).

#include "sjs/baselines/ours/sampling.h"  // reuse 2D preprocessing + enumerator

#include "sjs/baselines/detail/vector_join_enumerator.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {

template <int Dim, class T = Scalar>
class OursEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "ours_enum+sampling"; }

  void Reset() override {
    ctx_.Reset();
    built_ = false;
    pairs_cached_ = false;
    pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if (!ctx_.Build(ds, phases, err)) return false;

    built_ = true;
    return true;
  }

  // Exact join cardinality for Framework I.
  //
  // For strict SJS v3 alignment (single materialization), we *materialize* all
  // join pairs here and cache them for Sample().
  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    if (!pairs_cached_) {
      auto scoped = phases ? phases->Scoped("phase1_enumerate_materialize")
                           : PhaseRecorder::ScopedPhase(nullptr, "");

      pairs_.clear();
      const u64 cap = cfg.run.enum_cap;

      // Deterministic enumerator based on the same sweep+REPORT primitives.
      detail::OursReportJoinEnumerator2D<Dim, T> it(&ctx_);
      PairId p;
      while (it.Next(&p)) {
        pairs_.push_back(p);
        if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
          pairs_.clear();
          pairs_cached_ = false;
          if (err) *err = "OursEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
          return false;
        }
      }

      pairs_cached_ = true;
    }

    if (out) *out = MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  // Draw samples by uniform indexing into the cached join vector.
  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: run.t too large (must fit in u32)";
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
      out->pairs.resize(t);
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[i] = pairs_[static_cast<usize>(idx)];
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

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    if (pairs_cached_) {
      return std::make_unique<baselines::detail::VectorJoinEnumerator>(&pairs_);
    }
    return std::make_unique<detail::OursReportJoinEnumerator2D<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};

  bool pairs_cached_{false};
  std::vector<PairId> pairs_;

  detail::Ours2DContext<Dim, T> ctx_;
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs

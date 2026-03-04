#pragma once
// baselines/runners/enum_sampling_runner.h
//
// Runner for Variant::EnumSampling (Enumerate + Sampling).
//
// SJS v3 alignment
// ---------------
// In SJS v3, the EnumSampling variant corresponds to "Framework I":
//   - enumerate/materialize the full join result J
//   - draw t i.i.d. uniform samples with replacement by uniform indexing
//
// To keep responsibilities clear, the runner enforces a simple, consistent
// protocol and delegates the concrete enumeration/materialization strategy to
// the baseline implementation:
//   Reset -> Build -> Count -> Sample
//
// Notes:
//  - Baselines are expected to return exact counts (or fail) for this variant.
//  - cfg.run.enum_cap is treated as a safety cap: if a baseline cannot
//    materialize the full join under the cap, it should return an error rather
//    than producing biased samples.

#include "baselines/baseline_api.h"
#include "core/assert.h"
#include "core/logging.h"

#include <string>

namespace sjs {
namespace baselines {

template <int Dim, class T = Scalar>
bool RunEnumSamplingOnce(IBaseline<Dim, T>* baseline,
                         const Dataset<Dim, T>& dataset,
                         const Config& cfg,
                         u64 seed,
                         RunReport* out,
                         std::string* err = nullptr,
                         RunnerReusePolicy policy = RunnerReusePolicy{}) {
  if (!baseline) {
    if (err) *err = "RunEnumSamplingOnce: baseline is null";
    return false;
  }
  if (!out) {
    if (err) *err = "RunEnumSamplingOnce: out is null";
    return false;
  }
  if (!policy.build_before_run && policy.reset_before_run) {
    if (err) *err = "RunEnumSamplingOnce: invalid RunnerReusePolicy (reset=true with build=false)";
    return false;
  }

  out->ok = false;
  out->error.clear();
  out->method = baseline->method();
  out->variant = Variant::EnumSampling;
  out->baseline_name = std::string(baseline->Name());
  out->dataset_name = dataset.name;
  out->seed = seed;
  out->t = cfg.run.t;
  out->count = CountResult{};
  out->samples.Clear();
  out->used_enumeration = true;
  out->enumeration_truncated = false;
  out->enumeration_cap = cfg.run.enum_cap;
  out->enumeration_pairs_pass1 = 0;
  out->enumeration_pairs_pass2 = 0;
  out->enum_stats_pass1.Reset();
  out->enum_stats_pass2.Reset();
  out->adaptive_branch.clear();
  out->adaptive_pilot_pairs = 0;
  out->note.clear();
  out->phases.Clear();

  std::string local_err;
  // SJS v3: separate RNG streams across phases.
  Rng rng_count(DeriveSeed(seed, 1));
  Rng rng_sample(DeriveSeed(seed, 2));

  if (policy.reset_before_run) baseline->Reset();

  if (policy.build_before_run) {
    auto _ = out->phases.Scoped("run_build");
    if (!baseline->Build(dataset, cfg, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_count");
    if (!baseline->Count(cfg, &rng_count, &out->count, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }
  if (!out->count.exact) {
    out->error = "RunEnumSamplingOnce: EnumSampling requires exact count/materialization";
    if (err) *err = out->error;
    return false;
  }

  {
    auto _ = out->phases.Scoped("run_sample");
    if (!baseline->Sample(cfg, &rng_sample, &out->samples, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // Defensive validation.
  {
    std::string v_err;
    if (!out->samples.Validate(&v_err)) {
      out->error = "SampleSet validation failed: " + v_err;
      if (err) *err = out->error;
      return false;
    }
  }

  const u64 rounded_count = out->count.RoundedU64();
  if (rounded_count == 0ULL && !out->samples.Empty()) {
    out->error = "RunEnumSamplingOnce: baseline returned non-empty samples on empty join";
    if (err) *err = out->error;
    return false;
  }
  if (rounded_count > 0ULL &&
      out->samples.Size() != static_cast<usize>(cfg.run.t)) {
    out->error = "RunEnumSamplingOnce: baseline returned sample size != t on non-empty join";
    if (err) *err = out->error;
    return false;
  }

  // Populate a few convenience fields.
  if (out->count.exact) {
    out->enumeration_pairs_pass1 = out->count.RoundedU64();
  }

  out->ok = true;
  return true;
}

}  // namespace baselines
}  // namespace sjs

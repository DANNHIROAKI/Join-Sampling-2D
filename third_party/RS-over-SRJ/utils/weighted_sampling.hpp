#pragma once
// utils/weighted_sampling.hpp (modified)
//
// This file originally contained an independent alias-table implementation with
// its own RNG (pcg32). In the Join-Sampling-2D repository, we already have a
// project-wide alias table (sjs::sampling::AliasTable) and RNG (sjs::Rng).
//
// To ensure reproducibility across baselines (same seed => same sampling stream)
// and to avoid mixing RNG sources, we wrap the project implementation here.

#include "core/rng.h"
#include "core/types.h"
#include "sampling/alias_table.h"

#include <string>
#include <vector>

// Keep the original class name `alias` to minimize downstream edits.
class alias {
 public:
  alias() = default;

  explicit alias(const std::vector<sjs::u64>& weights) {
    Build(weights);
  }

  void Clear() { table_.Clear(); }

  sjs::usize Size() const noexcept { return table_.Size(); }
  bool Empty() const noexcept { return table_.Empty(); }

  // Build from u64 weights.
  bool Build(const std::vector<sjs::u64>& weights, std::string* err = nullptr) {
    return table_.BuildFromU64(sjs::Span<const sjs::u64>(weights), err);
  }

  // Sample one index using the provided RNG.
  sjs::usize Sample(sjs::Rng* rng) const noexcept {
    return table_.Sample(rng);
  }

 private:
  sjs::sampling::AliasTable table_;
};
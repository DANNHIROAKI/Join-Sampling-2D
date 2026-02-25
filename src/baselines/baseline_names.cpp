// src/baselines/baseline_names.cpp
//
// Baseline registry + help strings.
//
// Keeping this in a .cpp (not a header) reduces compile-time and keeps the
// registry centralized for CLI/help/error messaging.

#include "baselines/baseline_factory_2d.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace sjs {
namespace baselines {

Span<const BaselineSpec2D> BaselineRegistry2D() noexcept {
  // IMPORTANT: Keep this list consistent with CreateBaseline2D(...) in
  // baseline_factory_2d.cpp.
  //
  // Canonical key convention: "<method>/<variant>".
  static const BaselineSpec2D kReg[] = {
      // Ours
      {Method::Ours, Variant::Sampling, "ours/sampling",
       "Our method (Framework II): 2-pass sweep + event-block primitives, exact i.i.d. uniform"},
      {Method::Ours, Variant::EnumSampling, "ours/enum_sampling",
       "Our method (Framework I): materialize full join then uniform index sampling"},
      {Method::Ours, Variant::Adaptive, "ours/adaptive",
       "Our method (Framework III): budgeted full-cache + prefetch sample-cache"},

      // Range tree baseline
      {Method::RangeTree, Variant::Sampling, "range_tree/sampling",
       "Comparison method 1 (Framework II): plane sweep + orthogonal range-tree event-block primitives"},
      {Method::RangeTree, Variant::EnumSampling, "range_tree/enum_sampling",
       "Range-tree baseline (Framework I): materialize full join then uniform index sampling"},
      {Method::RangeTree, Variant::Adaptive, "range_tree/adaptive",
       "Range-tree baseline (Framework III): budgeted full-cache + prefetch sample-cache"},

      // KD-tree baseline (comparison method 2, SJS v3 Ch.5) — sampling-only
      {Method::KDTree, Variant::Sampling, "kd_tree/sampling",
       "KD-tree on 2d-dimensional embedding: exact COUNT + range SAMPLE join (sampling only)"},
  };

  return Span<const BaselineSpec2D>(kReg, sizeof(kReg) / sizeof(kReg[0]));
}

bool IsBaselineSupported2D(Method method, Variant variant) noexcept {
  const auto reg = BaselineRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    if (reg[i].method == method && reg[i].variant == variant) return true;
  }
  return false;
}

std::string BaselineHelp2D() {
  std::ostringstream oss;

  oss << "  Methods (canonical):\n";
  // Gather unique method names in registry order.
  {
    std::vector<Method> seen;
    const auto reg = BaselineRegistry2D();
    seen.reserve(reg.size());
    for (usize i = 0; i < reg.size(); ++i) {
      const Method m = reg[i].method;
      if (std::find(seen.begin(), seen.end(), m) == seen.end()) seen.push_back(m);
    }
    for (Method m : seen) {
      oss << "    - " << ToString(m) << "\n";
    }
  }

  oss << "\n  Variants (canonical):\n";
  oss << "    - sampling\n";
  oss << "    - enum_sampling\n";
  oss << "    - adaptive\n";

  oss << "\n  Supported (Dim=2) method/variant combinations:\n";
  const auto reg = BaselineRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    const auto& r = reg[i];
    oss << "    - " << r.key << "  (" << ToString(r.method) << ", " << ToString(r.variant) << ")";
    if (!r.desc.empty()) {
      oss << ": " << r.desc;
    }
    oss << "\n";
  }

  oss << "\n  Example:\n";
  oss << "    --method=ours --variant=sampling --t=100000 --seed=1\n";

  return oss.str();
}

}  // namespace baselines
}  // namespace sjs

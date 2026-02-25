// src/baselines/baseline_factory_2d.cpp
//
// Implementation of the Dim=2 baseline factory.
//
// This translation unit intentionally includes baseline headers.
// Keeping this centralized means apps only include baseline_factory_2d.h.

#include "baselines/baseline_factory_2d.h"

#include "sjs/core/types.h"  // ParseMethod/ParseVariant

// --------------------------
// Baseline implementations (each provides 3 variants)
// --------------------------
#include "sjs/baselines/ours/adaptive.h"
#include "sjs/baselines/ours/enum_sampling.h"
#include "sjs/baselines/ours/sampling.h"

#include "sjs/baselines/range_tree/adaptive.h"
#include "sjs/baselines/range_tree/enum_sampling.h"
#include "sjs/baselines/range_tree/sampling.h"

// Comparison method 2 (SJS v3 Ch.5): kd-tree on 2d-dimensional embedding.
#include "sjs/baselines/kd_tree/sampling.h"

#include <sstream>

namespace sjs {
namespace baselines {

namespace {

inline std::unique_ptr<IBaseline<2, Scalar>> MakeUnsupported(std::string_view method,
                                                            std::string_view variant,
                                                            std::string* err) {
  if (err) {
    std::ostringstream oss;
    oss << "Unsupported baseline for Dim=2: method='" << method << "' variant='" << variant << "'.\n";
    oss << "Available baselines:\n" << BaselineHelp2D();
    *err = oss.str();
  }
  return nullptr;
}

}  // namespace

std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(Method method,
                                                       Variant variant,
                                                       std::string* err) {
  // Keep consistent with BaselineRegistry2D() in baseline_names.cpp.

  switch (method) {
    case Method::Ours:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<ours::OursSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<ours::OursEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<ours::OursAdaptiveBaseline<2, Scalar>>();
        default:
          break;
      }
      break;

    case Method::RangeTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<range_tree::RangeTreeSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<range_tree::RangeTreeEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<range_tree::RangeTreeAdaptiveBaseline<2, Scalar>>();
        default:
          break;
      }
      break;

    case Method::KDTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<kd_tree::KDTree2DRangeSamplingBaseline<2, Scalar>>();
        default:
          break;
      }
      break;

    case Method::Unknown:
    default:
      break;
  }

  return MakeUnsupported(ToString(method), ToString(variant), err);
}

std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(std::string_view method,
                                                       std::string_view variant,
                                                       std::string* err) {
  Method m = Method::Unknown;
  Variant v = Variant::Sampling;

  if (!ParseMethod(method, &m) || m == Method::Unknown) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown method: '" << method << "'.\n";
      oss << "Available baselines:\n" << BaselineHelp2D();
      *err = oss.str();
    }
    return nullptr;
  }
  if (!ParseVariant(variant, &v)) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown variant: '" << variant << "'.\n";
      oss << "Allowed variants: sampling | enum_sampling | adaptive.\n";
      oss << "Available baselines:\n" << BaselineHelp2D();
      *err = oss.str();
    }
    return nullptr;
  }

  return CreateBaseline2D(m, v, err);
}

}  // namespace baselines
}  // namespace sjs

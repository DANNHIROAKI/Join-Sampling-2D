#pragma once
// sjs/baselines/kd_tree/sampling.h
//
// KD-tree baseline (Comparison method 2 in SJS v3, Chapter 5):
//   kd-tree on the 2d-dimensional embedding, supporting exact range COUNT and
//   range SAMPLE, used to produce i.i.d. uniform (with replacement) samples on
//   the full join J.
//
// This baseline has a single variant: Variant::Sampling.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/join/join_enumerator.h"
#include "sjs/geometry/box.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace kd_tree {

namespace detail {

inline bool ParseU64(std::string_view s, u64* out) {
  if (!out) return false;
  if (s.empty()) return false;
  try {
    std::size_t idx = 0;
    const std::string tmp(s);
    const unsigned long long v = std::stoull(tmp, &idx, 10);
    if (idx != tmp.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline u64 ExtraU64Or(const Config& cfg, std::string_view key, u64 def) {
  auto it = cfg.run.extra.find(std::string(key));
  if (it == cfg.run.extra.end()) return def;
  u64 v = def;
  if (!ParseU64(it->second, &v)) return def;
  return v;
}

// --------------------------
// Static kd-tree for open-bound orthogonal ranges in D dimensions.
// --------------------------

template <int D, class T>
struct KDPoint {
  std::array<T, D> x{};
  u32 payload = 0;  // e.g., index into the right relation
};

template <int D, class T>
struct OpenQuery {
  std::array<bool, D> has_lo{};
  std::array<bool, D> has_hi{};
  std::array<T, D> lo{};
  std::array<T, D> hi{};

  bool Contains(const KDPoint<D, T>& p) const noexcept {
    for (int j = 0; j < D; ++j) {
      const T v = p.x[static_cast<usize>(j)];
      if (has_hi[static_cast<usize>(j)] && !(v < hi[static_cast<usize>(j)])) return false;
      if (has_lo[static_cast<usize>(j)] && !(v > lo[static_cast<usize>(j)])) return false;
    }
    return true;
  }
};

template <int D, class T>
class StaticKDTree {
 public:
  struct BuildOptions {
    u32 leaf_size = 32;
  };

  void Clear() {
    pts_.clear();
    nodes_.clear();
    root_ = kNull;
    opt_ = BuildOptions{};
  }

  bool Empty() const noexcept { return root_ == kNull; }
  usize Size() const noexcept { return pts_.size(); }

  void Build(std::vector<KDPoint<D, T>>&& points, BuildOptions opt = BuildOptions{}) {
    Clear();
    pts_ = std::move(points);
    opt_ = opt;
    if (pts_.empty()) return;
    nodes_.reserve(pts_.size() * 2);
    root_ = BuildRec(0U, static_cast<u32>(pts_.size()), 0);
  }

  u64 Count(const OpenQuery<D, T>& q) const {
    if (root_ == kNull) return 0;
    return CountRec(root_, q);
  }

  // Sample k points i.i.d. uniformly (with replacement) from pts ∩ q.
  // Returns false iff k>0 and the query set is empty.
  bool Sample(const OpenQuery<D, T>& q,
              u32 k,
              Rng* rng,
              std::vector<u32>* out_payload,
              std::string* err = nullptr) const {
    if (!out_payload) {
      if (err) *err = "StaticKDTree::Sample: out_payload is null";
      return false;
    }
    out_payload->clear();
    if (k == 0) return true;
    if (!rng) {
      if (err) *err = "StaticKDTree::Sample: rng is null";
      return false;
    }
    if (root_ == kNull) {
      if (err) *err = "StaticKDTree::Sample: empty tree";
      return false;
    }

    std::vector<SampleBlock> blocks;
    blocks.reserve(64);
    std::vector<u32> hit_pool;
    hit_pool.reserve(256);

    CollectBlocks(root_, q, &blocks, &hit_pool);

    if (blocks.empty()) {
      if (err) *err = "StaticKDTree::Sample: query set is empty";
      return false;
    }

    std::vector<u64> weights;
    weights.reserve(blocks.size());
    for (const auto& b : blocks) weights.push_back(b.weight);

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(weights), err)) {
      if (err && err->empty()) *err = "StaticKDTree::Sample: failed to build alias";
      return false;
    }

    out_payload->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const u32 bi = static_cast<u32>(alias.Sample(rng));
      if (bi >= static_cast<u32>(blocks.size())) {
        if (err) *err = "StaticKDTree::Sample: alias sampled out of range";
        return false;
      }
      const SampleBlock& b = blocks[static_cast<usize>(bi)];
      if (b.contained) {
        const Node& nd = nodes_[b.node];
        const u64 off = rng->UniformU64(static_cast<u64>(nd.r - nd.l));
        const u32 pos = nd.l + static_cast<u32>(off);
        out_payload->push_back(pts_[static_cast<usize>(pos)].payload);
      } else {
        const u64 off = rng->UniformU64(static_cast<u64>(b.hit_len));
        const u32 pos = hit_pool[static_cast<usize>(b.hit_off + static_cast<u32>(off))];
        out_payload->push_back(pts_[static_cast<usize>(pos)].payload);
      }
    }
    return true;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    u32 left = kNull;
    u32 right = kNull;
    u32 l = 0;
    u32 r = 0;
    u8 split_dim = 0;
    std::array<T, D> mn{};
    std::array<T, D> mx{};

    u32 Size() const noexcept { return r - l; }
    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  std::vector<KDPoint<D, T>> pts_;
  std::vector<Node> nodes_;
  u32 root_ = kNull;
  BuildOptions opt_{};

  // A sampling block per SJS v3 §5.5:
  //  - contained node -> whole subtree slice is inside the open query
  //  - boundary leaf  -> only some points in leaf match (tracked by hit_pool)
  struct SampleBlock {
    bool contained = false;
    u32 node = kNull;
    u32 hit_off = 0;
    u32 hit_len = 0;
    u64 weight = 0;
  };

  static bool DisjointDim(const Node& nd, const OpenQuery<D, T>& q, int j) noexcept {
    const usize u = static_cast<usize>(j);
    const bool has_lo = q.has_lo[u];
    const bool has_hi = q.has_hi[u];
    if (has_hi && !has_lo) {
      return !(nd.mn[u] < q.hi[u]);  // mn >= hi
    }
    if (has_lo && !has_hi) {
      return !(q.lo[u] < nd.mx[u]);  // mx <= lo
    }
    if (has_lo && has_hi) {
      return (!(nd.mn[u] < q.hi[u])) || (!(q.lo[u] < nd.mx[u]));
    }
    return false;
  }

  static bool ContainedDim(const Node& nd, const OpenQuery<D, T>& q, int j) noexcept {
    const usize u = static_cast<usize>(j);
    const bool has_lo = q.has_lo[u];
    const bool has_hi = q.has_hi[u];
    if (has_hi && !has_lo) {
      return nd.mx[u] < q.hi[u];
    }
    if (has_lo && !has_hi) {
      return q.lo[u] < nd.mn[u];
    }
    if (has_lo && has_hi) {
      return (q.lo[u] < nd.mn[u]) && (nd.mx[u] < q.hi[u]);
    }
    return true;  // unbounded
  }

  static bool IsDisjoint(const Node& nd, const OpenQuery<D, T>& q) noexcept {
    for (int j = 0; j < D; ++j) {
      if (DisjointDim(nd, q, j)) return true;
    }
    return false;
  }

  static bool IsContained(const Node& nd, const OpenQuery<D, T>& q) noexcept {
    for (int j = 0; j < D; ++j) {
      if (!ContainedDim(nd, q, j)) return false;
    }
    return true;
  }

  u32 BuildRec(u32 l, u32 r, int depth) {
    const u32 count = r - l;
    SJS_DASSERT(count > 0);

    const u32 me = static_cast<u32>(nodes_.size());
    nodes_.push_back(Node{});
    Node& nd = nodes_.back();
    nd.l = l;
    nd.r = r;
    nd.split_dim = static_cast<u8>(depth % D);

    if (count <= opt_.leaf_size) {
      ComputeBounds(nd);
      return me;
    }

    const u32 mid = l + count / 2;
    const int axis = depth % D;
    std::nth_element(pts_.begin() + l, pts_.begin() + mid, pts_.begin() + r,
                     [&](const KDPoint<D, T>& a, const KDPoint<D, T>& b) {
                       const T av = a.x[static_cast<usize>(axis)];
                       const T bv = b.x[static_cast<usize>(axis)];
                       if (av < bv) return true;
                       if (bv < av) return false;
                       return a.payload < b.payload;
                     });

    // Recurse.
    const u32 left_n = BuildRec(l, mid, depth + 1);
    const u32 right_n = BuildRec(mid, r, depth + 1);
    nodes_[me].left = left_n;
    nodes_[me].right = right_n;

    // Combine bounds.
    const Node& L = nodes_[left_n];
    const Node& R = nodes_[right_n];
    for (int j = 0; j < D; ++j) {
      const usize u = static_cast<usize>(j);
      nodes_[me].mn[u] = std::min(L.mn[u], R.mn[u]);
      nodes_[me].mx[u] = std::max(L.mx[u], R.mx[u]);
    }
    return me;
  }

  void ComputeBounds(Node& nd) {
    SJS_DASSERT(nd.r > nd.l);
    // Initialize to first point.
    const KDPoint<D, T>& p0 = pts_[static_cast<usize>(nd.l)];
    for (int j = 0; j < D; ++j) {
      const usize u = static_cast<usize>(j);
      nd.mn[u] = p0.x[u];
      nd.mx[u] = p0.x[u];
    }
    for (u32 i = nd.l + 1; i < nd.r; ++i) {
      const KDPoint<D, T>& p = pts_[static_cast<usize>(i)];
      for (int j = 0; j < D; ++j) {
        const usize u = static_cast<usize>(j);
        nd.mn[u] = std::min(nd.mn[u], p.x[u]);
        nd.mx[u] = std::max(nd.mx[u], p.x[u]);
      }
    }
  }

  u64 CountRec(u32 node, const OpenQuery<D, T>& q) const {
    const Node& nd = nodes_[node];
    if (IsDisjoint(nd, q)) return 0;
    if (IsContained(nd, q)) return static_cast<u64>(nd.Size());
    if (nd.IsLeaf()) {
      u64 c = 0;
      for (u32 i = nd.l; i < nd.r; ++i) {
        if (q.Contains(pts_[static_cast<usize>(i)])) ++c;
      }
      return c;
    }
    return CountRec(nd.left, q) + CountRec(nd.right, q);
  }

  void CollectBlocks(u32 root,
                     const OpenQuery<D, T>& q,
                     std::vector<SampleBlock>* blocks,
                     std::vector<u32>* hit_pool) const {
    blocks->clear();
    hit_pool->clear();

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root);

    while (!stack.empty()) {
      const u32 ni = stack.back();
      stack.pop_back();
      const Node& nd = nodes_[ni];

      if (IsDisjoint(nd, q)) continue;
      if (IsContained(nd, q)) {
        blocks->push_back(SampleBlock{true, ni, 0, 0, static_cast<u64>(nd.Size())});
        continue;
      }

      if (nd.IsLeaf()) {
        const u32 off = static_cast<u32>(hit_pool->size());
        for (u32 i = nd.l; i < nd.r; ++i) {
          if (q.Contains(pts_[static_cast<usize>(i)])) hit_pool->push_back(i);
        }
        const u32 len = static_cast<u32>(hit_pool->size()) - off;
        if (len > 0) {
          blocks->push_back(SampleBlock{false, kNull, off, len, static_cast<u64>(len)});
        }
      } else {
        stack.push_back(nd.left);
        stack.push_back(nd.right);
      }
    }
  }
};

// Adapter enumerator for tests / tooling.
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

}  // namespace detail

// --------------------------
// KDTree2DRangeSamplingBaseline (Dim==2)
// --------------------------

template <int Dim, class T = Scalar>
class KDTree2DRangeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "KDTree2DRangeSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int D = 2 * Dim;  // embedding dimension

  Method method() const noexcept override { return Method::KDTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "kd_tree_2drc_sampling"; }

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

    // Build point set P from S.
    {
      auto _ = phases ? phases->Scoped("build_points") : PhaseRecorder::ScopedPhase(nullptr, "");
      std::vector<detail::KDPoint<D, T>> pts;
      pts.reserve(ds.S.Size());
      for (usize i = 0; i < ds.S.boxes.size(); ++i) {
        const auto& b = ds.S.boxes[i];
        detail::KDPoint<D, T> p;
        p.x[0] = b.lo.v[0];
        p.x[1] = b.lo.v[1];
        p.x[2] = b.hi.v[0];
        p.x[3] = b.hi.v[1];
        p.payload = static_cast<u32>(i);
        pts.push_back(p);
      }

      typename detail::StaticKDTree<D, T>::BuildOptions opt;
      u64 leaf = detail::ExtraU64Or(cfg, "kd_leaf_size", 32);
      if (leaf < 1) leaf = 1;
      if (leaf > static_cast<u64>(std::numeric_limits<u32>::max())) {
        leaf = static_cast<u64>(std::numeric_limits<u32>::max());
      }
      opt.leaf_size = static_cast<u32>(leaf);
      tree_.Build(std::move(pts), opt);
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
      if (err) *err = "KDTree2DRangeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTree2DRangeSamplingBaseline::Count: out is null";
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

    u64 W = 0;
    for (usize i = 0; i < ds_->R.boxes.size(); ++i) {
      const auto& r = ds_->R.boxes[i];
      detail::OpenQuery<D, T> q;
      // (-inf, R_i(r)) for L-coordinates
      q.has_hi[0] = true; q.hi[0] = r.hi.v[0];
      q.has_hi[1] = true; q.hi[1] = r.hi.v[1];
      // (L_i(r), +inf) for R-coordinates
      q.has_lo[2] = true; q.lo[2] = r.lo.v[0];
      q.has_lo[3] = true; q.lo[3] = r.lo.v[1];

      const u64 w = tree_.Count(q);
      w_r_[i] = w;
      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "KDTree2DRangeSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }
    }

    W_ = W;
    weights_valid_ = true;
    if (!alias_r_.BuildFromU64(Span<const u64>(w_r_), err)) return false;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "KDTree2DRangeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "KDTree2DRangeSamplingBaseline::Sample: null rng/out";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "KDTree2DRangeSamplingBaseline::Sample: t too large for u32";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    // Assign output slots to left objects r by w(r)/|J|.
    struct SlotAssign { u32 ridx; u32 slot; };
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

    // For each distinct r, sample within Q(r) in batch.
    {
      auto _ = phases ? phases->Scoped("phase3_fill") : PhaseRecorder::ScopedPhase(nullptr, "");

      usize ptr = 0;
      while (ptr < asg.size()) {
        const u32 ridx = asg[ptr].ridx;
        usize end = ptr;
        while (end < asg.size() && asg[end].ridx == ridx) ++end;
        const u32 k = static_cast<u32>(end - ptr);

        const auto& r = ds_->R.boxes[static_cast<usize>(ridx)];
        detail::OpenQuery<D, T> q;
        q.has_hi[0] = true; q.hi[0] = r.hi.v[0];
        q.has_hi[1] = true; q.hi[1] = r.hi.v[1];
        q.has_lo[2] = true; q.lo[2] = r.lo.v[0];
        q.has_lo[3] = true; q.lo[3] = r.lo.v[1];

        std::vector<u32> picked;
        picked.reserve(static_cast<usize>(k));

        // Conditional range sampling on the right relation (SJS v3 §5.5).
        if (!tree_.Sample(q, k, rng, &picked, err)) return false;
        if (picked.size() != k) {
          if (err) *err = "KDTree2DRangeSamplingBaseline::Sample: SAMPLE returned wrong count";
          return false;
        }

        const Id rid = ds_->R.GetId(static_cast<usize>(ridx));
        for (u32 j = 0; j < k; ++j) {
          const u32 slot = asg[ptr + j].slot;
          const u32 s_idx = picked[static_cast<usize>(j)];
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
      if (err) *err = "KDTree2DRangeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::PlaneSweepEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0);
  }

 private:

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  detail::StaticKDTree<D, T> tree_;

  std::vector<u64> w_r_;
  sampling::AliasTable alias_r_;
  u64 W_ = 0;
  bool weights_valid_ = false;
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs

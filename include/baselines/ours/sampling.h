#pragma once
// baselines/ours/sampling.h
//
// "Our Method" — Sampling variant (3-phase uniform sampler).
//
// This implements our main method as specified in SJS v3 (Ch. 2–3):
//
//  - Common preprocessing (all variants share):
//      * Build plane-sweep events on axis 0 with half-open semantics:
//          - sort by coordinate
//          - END before START at the same coordinate
//          - START tie-break by SideTieBreak (R before S by default)
//      * For each side (R/S) build per-pattern indices (skeleton) and keep the
//        active set empty.
//
//  - Sampling variant (Ours):
//      Phase 1: single sweep, for each START event e with box q compute exact
//               w_e^A and w_e^B and w_e = w_e^A + w_e^B. Let W = sum_e w_e = |J|.
//      Phase 2: build event-level alias on (w_e). For each output slot j=1..t:
//               sample event E_j ~ w_e/W; sample pattern G_j in {A,B} with
//               Pr(G_j=A|E_j=e)=w_e^A/w_e; then append slot j into S_e^G.
//      Phase 3: second sweep. When visiting START event e for q:
//               for each pattern g, let k = |S_e^g|. If k>0, call
//                 Sample_g(q*, k)
//               to get k i.i.d. uniform "other-side" rectangles, then fill
//               Ans at those slots, keeping pair order as (R,S).
//
// Dim notes:
//  - The repository currently instantiates Dim=2. This implementation is
//    specialized for Dim=2 but keeps the template parameter so it integrates
//    cleanly with the baseline factory.

#include "baselines/baseline_api.h"

#include "core/assert.h"
#include "join/join_enumerator.h"  // optional (kept for debugging / comparison)
#include "join/sweep_events.h"
#include "sampling/alias_table.h"

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
namespace ours {

namespace detail {

// -----------------------------------------------------------------------------
// Pattern definitions for 2D (y-axis)
// -----------------------------------------------------------------------------
// For half-open intervals [lo,hi):
//   [a0,a1) intersects [b0,b1)  <=>  a0 < b1 AND b0 < a1.
// When sweep axis (x) overlap is guaranteed by active-set membership, a pair
// intersects iff y-intervals intersect.
//
// Partition by comparing the lower endpoints (ties go to Pattern A):
//   Pattern A: y_r_lo <= y_q_lo < y_r_hi
//   Pattern B: y_q_lo < y_r_lo < y_q_hi
// These sets are disjoint and cover all y-intersections.

// -----------------------------------------------------------------------------
// RangePointSegTree (Pattern B): dynamic point range count/sample/report
// -----------------------------------------------------------------------------
// We maintain active points keyed by y_r_lo (rank on the compressed ylo-domain).
//
// Key design goal (per SJS v3):
//   SampleRange(l,r,k) should cost O(log m + k) for fixed query (l,r), by
//   building a small alias table once and drawing k samples from it.
//
// Technique:
//   - Segment tree over ranks [0, P), P = next power-of-two >= m.
//   - Insert a point into *all* nodes on its leaf-to-root path.
//   - For a query range [l,r), compute the canonical segment-tree cover nodes.
//     These cover nodes are disjoint and each leaf in the range belongs to
//     exactly one cover node; therefore each point in [l,r) appears in exactly
//     one cover node bucket (even though it is stored on multiple ancestors).
//   - To sample uniformly from all points in [l,r): choose a cover node with
//     probability proportional to its bucket size, then choose uniformly from
//     that bucket.
//
// Deletion is O(log m) using swap-delete + backrefs, same as the stabbing tree.


class RangePointSegTree {
 public:
  void Init(u32 num_handles, u32 num_ranks) {
    n_handles_ = num_handles;
    m_ = num_ranks;

    // Next power-of-two.
    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    nodes_.assign(static_cast<usize>(2 * p_), Node{});

    // log2(p_) (p_ is power-of-two)
    leaf_depth_ = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++leaf_depth_;

    // Path length leaf->root inclusive.
    max_refs_ = leaf_depth_ + 1;

    // Store only the position within each node bucket.
    // backref is derived from node depth: backref = leaf_depth - depth(node).
    pos_in_node_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), 0U);

    rank_of_handle_.assign(static_cast<usize>(n_handles_), kInvalidRank);
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    leaf_depth_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    pos_in_node_.clear();
    rank_of_handle_.clear();
  }

  // Keep the skeleton but drop all active points.
  void ResetActive() {
    for (auto& n : nodes_) n.items.clear();
    std::fill(rank_of_handle_.begin(), rank_of_handle_.end(), kInvalidRank);
  }

  void Insert(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);
    SJS_DASSERT(rank_of_handle_[handle] == kInvalidRank);

    rank_of_handle_[handle] = rank;

    u32 idx = rank + p_;
    while (idx > 0) {
      AddToNode(handle, idx);
      idx >>= 1;
    }
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 rank = rank_of_handle_[handle];
    if (rank == kInvalidRank) return;

    SJS_DASSERT(rank < m_);
    u32 idx = rank + p_;
    while (idx > 0) {
      const u32 backref = BackrefOfNode(idx);
      const u32 pos = pos_in_node_[PlacementIndex(handle, backref)];
      RemoveFromNode(idx, pos);
      idx >>= 1;
    }

    rank_of_handle_[handle] = kInvalidRank;
  }

  // Allocation-free exact count of active handles in ranks [l, r).
  u64 CountRange(u32 l, u32 r) const {
    if (r <= l || m_ == 0) return 0ULL;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return 0ULL;

    u32 L = l + p_;
    u32 R = r + p_;

    u64 total = 0;
    while (L < R) {
      if (L & 1U) total += static_cast<u64>(nodes_[L++].items.size());
      if (R & 1U) total += static_cast<u64>(nodes_[--R].items.size());
      L >>= 1;
      R >>= 1;
    }
    return total;
  }

  // Sample k handles uniformly from ranks in [l, r).
  // Returns false iff the range is empty.
  bool SampleRange(u32 l, u32 r, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (k == 0) return true;

    if (r <= l || m_ == 0) return false;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return false;

    // Canonical cover nodes in left-to-right order (size <= 2*log2(p_)+2).
    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(l, r, &cover);

    // Filter non-empty buckets (weights > 0).
    std::array<u32, kMaxCoverNodes> nodes;
    std::array<u64, kMaxCoverNodes> weights;
    u32 nz = 0;
    u64 total = 0;

    for (u32 i = 0; i < cover_sz; ++i) {
      const u32 node = cover[i];
      const u64 w = static_cast<u64>(nodes_[node].items.size());
      if (w == 0) continue;
      SJS_DASSERT(nz < kMaxCoverNodes);
      nodes[nz] = node;
      weights[nz] = w;
      total += w;
      ++nz;
    }

    if (total == 0 || nz == 0) return false;

    // For heavy k, use a stack-based alias table to avoid O(k * nz).
    // For light k, a prefix scan is usually faster.
    static constexpr u32 kAliasMinK = 64;
    if (nz > 1 && k >= kAliasMinK) {
      sampling::AliasSmall<kMaxCoverNodes> alias;
      (void)alias.BuildFromU64(Span<const u64>(weights.data(), nz), /*err=*/nullptr);

      for (u32 i = 0; i < k; ++i) {
        const u32 bi = static_cast<u32>(alias.Sample(rng));
        const u32 node = nodes[bi];
        const auto& bucket = nodes_[node].items;
        SJS_DASSERT(!bucket.empty());
        const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
        out->push_back(bucket[pos]);
      }
      return true;
    }

    // Light-k fallback: prefix scan selection.
    for (u32 i = 0; i < k; ++i) {
      const u64 x = rng->UniformU64(total);  // in [0,total)
      u64 cum = 0;
      u32 bi = 0;
      for (; bi < nz; ++bi) {
        cum += weights[bi];
        if (x < cum) break;
      }
      if (bi >= nz) bi = nz - 1;  // defensive

      const u32 node = nodes[bi];
      const auto& bucket = nodes_[node].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos]);
    }
    return true;
  }

  // Report all handles in [l, r) in a deterministic left-to-right cover order.
  void ReportRange(u32 l, u32 r, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l || m_ == 0) return;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return;

    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(l, r, &cover);

    for (u32 i = 0; i < cover_sz; ++i) {
      const u32 node = cover[i];
      const auto& bucket = nodes_[node].items;
      for (const u32 h : bucket) out->push_back(h);
    }
  }

 private:
  struct Node {
    std::vector<u32> items;  // handle-only
  };

  static constexpr u32 kInvalidRank = std::numeric_limits<u32>::max();
  static constexpr usize kMaxCoverNodes = 128;

  static u32 FloorLog2(u32 x) noexcept {
    SJS_DASSERT(x > 0);
#if defined(__GNUC__) || defined(__clang__)
    return 31U - static_cast<u32>(__builtin_clz(x));
#else
    u32 d = 0;
    while (x > 1) {
      x >>= 1;
      ++d;
    }
    return d;
#endif
  }

  u32 BackrefOfNode(u32 node) const noexcept {
    // backref = distance from leaf to node along the path.
    // depth(root)=0, depth(leaf)=leaf_depth_, so backref = leaf_depth - depth(node).
    const u32 d = FloorLog2(node);
    SJS_DASSERT(d <= leaf_depth_);
    return leaf_depth_ - d;
  }

  usize PlacementIndex(u32 handle, u32 backref) const {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node) {
    SJS_DASSERT(node < nodes_.size());

    const u32 backref = BackrefOfNode(node);
    SJS_DASSERT(backref < max_refs_);

    auto& bucket = nodes_[node].items;
    const u32 pos = static_cast<u32>(bucket.size());

    pos_in_node_[PlacementIndex(handle, backref)] = pos;
    bucket.push_back(handle);
  }

  void RemoveFromNode(u32 node, u32 pos) {
    SJS_DASSERT(node < nodes_.size());
    auto& bucket = nodes_[node].items;
    SJS_DASSERT(!bucket.empty());
    SJS_DASSERT(pos < bucket.size());

    const u32 backref = BackrefOfNode(node);

    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const u32 swapped_handle = bucket[last_pos];
      bucket[pos] = swapped_handle;
      pos_in_node_[PlacementIndex(swapped_handle, backref)] = pos;
    }
    bucket.pop_back();
  }

  // Canonical cover decomposition of [l,r) (0<=l<=r<=m_).
  // The returned nodes are disjoint and ordered left-to-right.
  u32 DecomposeOrdered(u32 l, u32 r, std::array<u32, kMaxCoverNodes>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l) return 0;

    u32 L = l + p_;
    u32 R = r + p_;

    std::array<u32, kMaxCoverNodes> left;
    std::array<u32, kMaxCoverNodes> right;
    u32 lsz = 0;
    u32 rsz = 0;

    while (L < R) {
      if (L & 1U) {
        SJS_DASSERT(lsz < kMaxCoverNodes);
        left[lsz++] = L++;
      }
      if (R & 1U) {
        SJS_DASSERT(rsz < kMaxCoverNodes);
        right[rsz++] = --R;
      }
      L >>= 1;
      R >>= 1;
    }

    const u32 out_sz = lsz + rsz;
    SJS_DASSERT(out_sz <= kMaxCoverNodes);

    u32 out_i = 0;
    for (u32 i = 0; i < lsz; ++i) (*out)[out_i++] = left[i];
    for (u32 i = 0; i < rsz; ++i) (*out)[out_i++] = right[rsz - 1 - i];

    SJS_DASSERT(out_i == out_sz);
    return out_sz;
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 leaf_depth_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;               // size 2*p_
  std::vector<u32> pos_in_node_;          // size n_handles_*max_refs_ (position only)
  std::vector<u32> rank_of_handle_;       // per-handle rank (needed to recompute nodes on erase)
};



class StabbingSegTree {
 public:
  struct Item {
    u32 handle = 0;
    u32 backref = 0;  // index into handle's placement list
  };

  void Init(u32 num_handles, u32 num_points) {
    n_handles_ = num_handles;
    m_ = num_points;

    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    nodes_.assign(static_cast<usize>(2 * p_), Node{});

    // log2(p_)
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;

    // Interval range decomposition touches <= 2*log2(p_) nodes.
    max_refs_ = static_cast<u32>(2 * logp + 4);

    // Store only the position within each node bucket.
    pos_in_node_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), 0U);

    lo_rank_.assign(static_cast<usize>(n_handles_), kInvalidRank);
    hi_rank_.assign(static_cast<usize>(n_handles_), kInvalidRank);
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    pos_in_node_.clear();
    lo_rank_.clear();
    hi_rank_.clear();
  }

  // Keep the skeleton but drop all active intervals.
  void ResetActive() {
    for (auto& n : nodes_) n.items.clear();
    std::fill(lo_rank_.begin(), lo_rank_.end(), kInvalidRank);
    std::fill(hi_rank_.begin(), hi_rank_.end(), kInvalidRank);
  }

  // Insert interval [L, R) on ranks.
  void Insert(u32 handle, u32 L, u32 R) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(lo_rank_[handle] == kInvalidRank && hi_rank_[handle] == kInvalidRank);

    if (R <= L || m_ == 0) return;
    L = std::min(L, m_);
    R = std::min(R, m_);
    if (R <= L) return;

    lo_rank_[handle] = L;
    hi_rank_[handle] = R;

    u32 l = L + p_;
    u32 r = R + p_;

    u32 backref = 0;
    while (l < r) {
      if (l & 1U) {
        AddToNode(handle, l, backref);
        ++backref;
        ++l;
      }
      if (r & 1U) {
        --r;
        AddToNode(handle, r, backref);
        ++backref;
      }
      l >>= 1;
      r >>= 1;
    }
    SJS_DASSERT(backref <= max_refs_);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    if (lo_rank_[handle] == kInvalidRank) return;

    const u32 L = lo_rank_[handle];
    const u32 R = hi_rank_[handle];
    SJS_DASSERT(L != kInvalidRank && R != kInvalidRank);

    u32 l = L + p_;
    u32 r = R + p_;

    u32 backref = 0;
    while (l < r) {
      if (l & 1U) {
        const u32 pos = pos_in_node_[PlacementIndex(handle, backref)];
        RemoveFromNode(l, pos);
        ++backref;
        ++l;
      }
      if (r & 1U) {
        --r;
        const u32 pos = pos_in_node_[PlacementIndex(handle, backref)];
        RemoveFromNode(r, pos);
        ++backref;
      }
      l >>= 1;
      r >>= 1;
    }

    SJS_DASSERT(backref <= max_refs_);
    lo_rank_[handle] = kInvalidRank;
    hi_rank_[handle] = kInvalidRank;
  }

  u64 Count(u32 q) const {
    if (m_ == 0 || q >= m_) return 0ULL;
    u64 total = 0;
    u32 idx = q + p_;
    while (idx > 0) {
      total += static_cast<u64>(nodes_[idx].items.size());
      idx >>= 1;
    }
    return total;
  }

  bool Sample(u32 q, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (k == 0) return true;
    if (m_ == 0 || q >= m_) return false;

    // Collect non-empty buckets on root-to-leaf path (allocation-free).
    std::array<u32, kMaxPathNodes> path_nodes;
    std::array<u64, kMaxPathNodes> weights;
    u32 n = 0;

    u32 idx = q + p_;
    u64 total = 0;
    while (idx > 0) {
      const u64 w = static_cast<u64>(nodes_[idx].items.size());
      if (w > 0) {
        SJS_DASSERT(n < kMaxPathNodes);
        path_nodes[n] = idx;
        weights[n] = w;
        total += w;
        ++n;
      }
      idx >>= 1;
    }

    if (total == 0 || n == 0) return false;

    // Heavy-k: stack alias.
    static constexpr u32 kAliasMinK = 64;
    if (n > 1 && k >= kAliasMinK) {
      sampling::AliasSmall<kMaxPathNodes> alias;
      (void)alias.BuildFromU64(Span<const u64>(weights.data(), n), /*err=*/nullptr);

      for (u32 i = 0; i < k; ++i) {
        const u32 bi = static_cast<u32>(alias.Sample(rng));
        const u32 node = path_nodes[bi];
        const auto& bucket = nodes_[node].items;
        SJS_DASSERT(!bucket.empty());
        const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
        out->push_back(bucket[pos].handle);
      }
      return true;
    }

    // Light-k: prefix scan.
    for (u32 i = 0; i < k; ++i) {
      const u64 x = rng->UniformU64(total);
      u64 cum = 0;
      u32 bi = 0;
      for (; bi < n; ++bi) {
        cum += weights[bi];
        if (x < cum) break;
      }
      if (bi >= n) bi = n - 1;  // defensive

      const u32 node = path_nodes[bi];
      const auto& bucket = nodes_[node].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos].handle);
    }
    return true;
  }

  void Report(u32 q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (m_ == 0 || q >= m_) return;

    u32 idx = q + p_;
    while (idx > 0) {
      const auto& bucket = nodes_[idx].items;
      for (const auto& it : bucket) out->push_back(it.handle);
      idx >>= 1;
    }
  }

 private:
  struct Node {
    std::vector<Item> items;
  };

  static constexpr u32 kInvalidRank = std::numeric_limits<u32>::max();
  static constexpr usize kMaxPathNodes = 64;

  usize PlacementIndex(u32 handle, u32 backref) const {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node, u32 backref) {
    SJS_DASSERT(node < nodes_.size());
    SJS_DASSERT(backref < max_refs_);

    auto& bucket = nodes_[node].items;
    const u32 pos = static_cast<u32>(bucket.size());

    pos_in_node_[PlacementIndex(handle, backref)] = pos;
    bucket.push_back(Item{handle, backref});
  }

  void RemoveFromNode(u32 node, u32 pos) {
    SJS_DASSERT(node < nodes_.size());
    auto& bucket = nodes_[node].items;
    SJS_DASSERT(!bucket.empty());
    SJS_DASSERT(pos < bucket.size());

    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const Item swapped = bucket[last_pos];
      bucket[pos] = swapped;
      pos_in_node_[PlacementIndex(swapped.handle, swapped.backref)] = pos;
    }
    bucket.pop_back();
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<u32> pos_in_node_;      // size n_handles_*max_refs_ (position only)
  std::vector<u32> lo_rank_;          // stored interval endpoints for erase (kInvalidRank => inactive)
  std::vector<u32> hi_rank_;
};

// -----------------------------------------------------------------------------
// ActiveIndex2D: per-side wrapper for both patterns
// -----------------------------------------------------------------------------

class ActiveIndex2D {
 public:
  void Init(u32 num_handles, u32 num_ylo_ranks) {
    n_ = num_handles;
    m_ = num_ylo_ranks;
    stab_.Init(num_handles, num_ylo_ranks);
    pts_.Init(num_handles, num_ylo_ranks);
  }

  void Clear() {
    n_ = 0;
    m_ = 0;
    stab_.Clear();
    pts_.Clear();
  }

  void ResetActive() {
    stab_.ResetActive();
    pts_.ResetActive();
  }

  void Insert(u32 handle, u32 ylo_rank, u32 yhi_lb_rank) {
    SJS_DASSERT(handle < n_);
    SJS_DASSERT(ylo_rank < m_);
    // yhi_lb_rank is in [0, m_].
    stab_.Insert(handle, ylo_rank, yhi_lb_rank);
    pts_.Insert(handle, ylo_rank);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_);
    stab_.Erase(handle);
    pts_.Erase(handle);
  }

  // Pattern A: y_r contains y_q_lo.
  u64 CountA(u32 y_q_lo_rank) const { return stab_.Count(y_q_lo_rank); }

  // Pattern B: y_r_lo in (y_q_lo, y_q_hi).
  u64 CountB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank) const {
    if (m_ == 0) return 0ULL;
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.CountRange(l, r);
  }

  bool SampleA(u32 y_q_lo_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    return stab_.Sample(y_q_lo_rank, k, rng, out);
  }

  bool SampleB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.SampleRange(l, r, k, rng, out);
  }

  void ReportA(u32 y_q_lo_rank, std::vector<u32>* out) const { stab_.Report(y_q_lo_rank, out); }

  void ReportB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    pts_.ReportRange(l, r, out);
  }

 private:
  u32 n_{0};
  u32 m_{0};
  StabbingSegTree stab_;
  RangePointSegTree pts_;
};

// -----------------------------------------------------------------------------
// CountIndex2D: count-only active-set index (two Fenwick trees)
//
// Used only in Phase 1 (Count) to compute exact wa/wb without maintaining the
// heavier sampling/report structures.
//
//   - Pattern A count: stabbing count via difference Fenwick
//       Insert interval [ylo, yhi_lb):  diff[ylo]+=1, diff[yhi_lb]-=1
//       CountA(q_ylo) = prefix_sum(diff, q_ylo)
//
//   - Pattern B count: point range count via point Fenwick
//       Insert point at ylo: pts[ylo]+=1
//       CountB(q_ylo, q_yhi_lb) = sum(pts, (q_ylo+1 .. q_yhi_lb))
// -----------------------------------------------------------------------------

class Fenwick {
 public:
  void Init(u32 n) {
    n_ = n;
    bit_.assign(static_cast<usize>(n_ + 1), 0);
  }

  void Clear() {
    n_ = 0;
    bit_.clear();
  }

  void Reset() { std::fill(bit_.begin(), bit_.end(), 0); }

  // Add delta at index i (0-based). Requires i < n_.
  void Add(u32 i, i64 delta) {
    SJS_DASSERT(i < n_);
    // BIT uses 1-based internally.
    for (u32 x = i + 1; x <= n_; x += x & (~x + 1)) {
      bit_[x] += delta;
    }
  }

  // Prefix sum over [0..i] inclusive. If i >= n_, returns SumAll().
  i64 PrefixSumI64(u32 i) const {
    if (n_ == 0) return 0;
    if (i >= n_) i = n_ - 1;
    i64 res = 0;
    for (u32 x = i + 1; x > 0; x &= (x - 1)) {
      res += bit_[x];
    }
    return res;
  }

  i64 RangeSumI64(u32 l, u32 r) const {
    if (r <= l || n_ == 0) return 0;
    if (l >= n_) return 0;
    if (r > n_) r = n_;
    if (r <= l) return 0;
    const i64 a = PrefixSumI64(r - 1);
    const i64 b = (l == 0) ? 0 : PrefixSumI64(l - 1);
    return a - b;
  }

 private:
  u32 n_{0};
  std::vector<i64> bit_;  // index 0 unused
};

class CountIndex2D {
 public:
  void Init(u32 num_ylo_ranks) {
    m_ = num_ylo_ranks;
    // diff has length m_+1 (we update at yhi_lb which can be == m_).
    diff_.Init(m_ + 1);
    pts_.Init(m_);
    ResetActive();
  }

  void Clear() {
    m_ = 0;
    diff_.Clear();
    pts_.Clear();
  }

  void ResetActive() {
    diff_.Reset();
    pts_.Reset();
  }

  // Insert rectangle with (ylo_rank, yhi_lb_rank) into active set.
  void Insert(u32 ylo_rank, u32 yhi_lb_rank) {
    if (m_ == 0) return;
    SJS_DASSERT(ylo_rank < m_);
    SJS_DASSERT(yhi_lb_rank <= m_);

    pts_.Add(ylo_rank, +1);
    diff_.Add(ylo_rank, +1);
    diff_.Add(yhi_lb_rank, -1);
  }

  void Erase(u32 ylo_rank, u32 yhi_lb_rank) {
    if (m_ == 0) return;
    SJS_DASSERT(ylo_rank < m_);
    SJS_DASSERT(yhi_lb_rank <= m_);

    pts_.Add(ylo_rank, -1);
    diff_.Add(ylo_rank, -1);
    diff_.Add(yhi_lb_rank, +1);
  }

  // Pattern A: y_r contains y_q_lo.
  u64 CountA(u32 y_q_lo_rank) const {
    if (m_ == 0 || y_q_lo_rank >= m_) return 0ULL;
    const i64 v = diff_.PrefixSumI64(y_q_lo_rank);
    return (v <= 0) ? 0ULL : static_cast<u64>(v);
  }

  // Pattern B: y_r_lo in (y_q_lo, y_q_hi).
  u64 CountB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank) const {
    if (m_ == 0) return 0ULL;
    u32 l = y_q_lo_rank + 1;
    if (l > m_) l = m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    const i64 v = pts_.RangeSumI64(l, r);
    return (v <= 0) ? 0ULL : static_cast<u64>(v);
  }

 private:
  u32 m_{0};
  Fenwick diff_;  // size m_+1
  Fenwick pts_;   // size m_
};

// -----------------------------------------------------------------------------
// Shared 2D preprocessing context (events + ranks + active indices)
// -----------------------------------------------------------------------------

template <int Dim, class T>
class Ours2DContext {
 public:
  using DatasetT = Dataset<Dim, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;

    events_.clear();
    start_id_of_event_.clear();
    start_event_pos_.clear();

    y_coords_.clear();
    ylo_rank_r_.clear();
    yhi_lb_rank_r_.clear();
    ylo_rank_s_.clear();
    yhi_lb_rank_s_.clear();

    active_r_.Clear();
    active_s_.Clear();
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    if constexpr (Dim != 2) {
      if (err) *err = "Ours2DContext: only Dim=2 is implemented";
      return false;
    }

    ds_ = &ds;

    const usize nR = ds.R.Size();
    const usize nS = ds.S.Size();
    if (nR > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "Ours2DContext: relation size exceeds u32";
      return false;
    }

    // 1) Events
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // 2) START id mapping (dense 0..|E|-1)
    start_id_of_event_.assign(events_.size(), -1);
    start_event_pos_.clear();
    start_event_pos_.reserve(nR + nS);
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_event_pos_.size());
        start_event_pos_.push_back(i);
      }
    }

    // 3) y-domain from all y-lower endpoints
    {
      auto scoped = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_coords_.clear();
      y_coords_.reserve(nR + nS);
      for (const auto& b : ds.R.boxes) y_coords_.push_back(b.lo.v[1]);
      for (const auto& b : ds.S.boxes) y_coords_.push_back(b.lo.v[1]);
      std::sort(y_coords_.begin(), y_coords_.end());
      y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());
    }

    if (y_coords_.empty()) {
      if (err) *err = "Ours2DContext: empty y-domain";
      return false;
    }

    const u32 m = static_cast<u32>(y_coords_.size());

    auto lb_rank = [&](T v) -> u32 {
      const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), v);
      return static_cast<u32>(std::distance(y_coords_.begin(), it));
    };

    // 4) Precompute ranks per box
    {
      auto scoped = phases ? phases->Scoped("build_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");

      ylo_rank_r_.resize(nR);
      yhi_lb_rank_r_.resize(nR);
      for (usize i = 0; i < nR; ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_r_[i] = lo;
        yhi_lb_rank_r_[i] = lb_rank(b.hi.v[1]);
      }

      ylo_rank_s_.resize(nS);
      yhi_lb_rank_s_.resize(nS);
      for (usize i = 0; i < nS; ++i) {
        const auto& b = ds.S.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_s_[i] = lo;
        yhi_lb_rank_s_[i] = lb_rank(b.hi.v[1]);
      }
    }

    // 5) Build active indices skeleton
    {
      auto scoped = phases ? phases->Scoped("build_active_indices") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.Init(static_cast<u32>(nR), m);
      active_s_.Init(static_cast<u32>(nS), m);
      active_r_.ResetActive();
      active_s_.ResetActive();
    }

    built_ = true;
    return true;
  }

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  const std::vector<join::Event>& events() const noexcept { return events_; }
  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  usize num_start_events() const noexcept { return start_event_pos_.size(); }

  u32 num_ylo_ranks() const noexcept { return static_cast<u32>(y_coords_.size()); }


  const std::vector<u32>& ylo_rank_r() const noexcept { return ylo_rank_r_; }
  const std::vector<u32>& yhi_lb_rank_r() const noexcept { return yhi_lb_rank_r_; }
  const std::vector<u32>& ylo_rank_s() const noexcept { return ylo_rank_s_; }
  const std::vector<u32>& yhi_lb_rank_s() const noexcept { return yhi_lb_rank_s_; }

  ActiveIndex2D& active_r() noexcept { return active_r_; }
  ActiveIndex2D& active_s() noexcept { return active_s_; }
  const ActiveIndex2D& active_r() const noexcept { return active_r_; }
  const ActiveIndex2D& active_s() const noexcept { return active_s_; }

  void ResetActive() {
    active_r_.ResetActive();
    active_s_.ResetActive();
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;   // per event position (size events_.size())
  std::vector<usize> start_event_pos_;   // positions of START events in events_

  std::vector<T> y_coords_;              // unique y-lower values
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  ActiveIndex2D active_r_;
  ActiveIndex2D active_s_;
};

// -----------------------------------------------------------------------------
// Phase2 slot plan (shared by Sampling + Adaptive large branch)
// -----------------------------------------------------------------------------

struct SlotPlan2D {
  // Offsets are size E+1; slots arrays are size t.
  std::vector<u32> offset_a;
  std::vector<u32> offset_b;
  std::vector<u32> slots_a;
  std::vector<u32> slots_b;

  void Clear() {
    offset_a.clear();
    offset_b.clear();
    slots_a.clear();
    slots_b.clear();
  }
};

inline bool BuildSlotPlan2D(u32 t,
                           Rng* rng,
                           const std::vector<u64>& w_total,
                           const std::vector<u64>& w_a,
                           const std::vector<u64>& w_b,
                           SlotPlan2D* plan,
                           std::string* err) {
  SJS_DASSERT(rng != nullptr);
  SJS_DASSERT(plan != nullptr);
  plan->Clear();

  const usize E = w_total.size();
  if (E == 0) {
    if (err) *err = "BuildSlotPlan2D: empty event set";
    return false;
  }

  // Pack per-slot (eid,pattern) into one u32: info = (eid<<1) | pat.
  // This requires eid to fit in 31 bits.
  if (E > (static_cast<usize>(std::numeric_limits<u32>::max()) >> 1)) {
    if (err) *err = "BuildSlotPlan2D: too many events to pack slot_info";
    return false;
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(Span<const u64>(w_total), err)) {
    if (err && err->empty()) *err = "BuildSlotPlan2D: failed to build alias";
    return false;
  }

  // First pass: per-slot assignment + per-event counts.
  std::vector<u32> slot_info;
  slot_info.resize(t);

  std::vector<u32> cnt_a(E, 0U);
  std::vector<u32> cnt_b(E, 0U);

  for (u32 j = 0; j < t; ++j) {
    const u32 eid = static_cast<u32>(alias.Sample(rng));
    const u64 wa = w_a[eid];
    const u64 wb = w_b[eid];
    const u64 w = wa + wb;
    SJS_DASSERT(w == w_total[eid]);

    // Choose pattern conditional on the event.
    u32 pat = 0;  // 0=A, 1=B
    if (wa == 0) {
      pat = 1;
    } else if (wb == 0) {
      pat = 0;
    } else {
      const u64 r = rng->UniformU64(w);
      pat = (r < wa) ? 0U : 1U;
    }

    slot_info[j] = (eid << 1) | pat;

    if (pat == 0) {
      ++cnt_a[eid];
    } else {
      ++cnt_b[eid];
    }
  }

  // Prefix sums -> offsets.
  plan->offset_a.resize(E + 1);
  plan->offset_b.resize(E + 1);
  plan->offset_a[0] = 0;
  plan->offset_b[0] = 0;
  for (usize e = 0; e < E; ++e) {
    plan->offset_a[e + 1] = plan->offset_a[e] + cnt_a[e];
    plan->offset_b[e + 1] = plan->offset_b[e] + cnt_b[e];
  }

  const u32 total_a = plan->offset_a[E];
  const u32 total_b = plan->offset_b[E];
  SJS_DASSERT(total_a + total_b == t);

  plan->slots_a.resize(total_a);
  plan->slots_b.resize(total_b);

  // Second pass: stable-fill the slot indices into the flat arrays.
  std::vector<u32> cur_a = plan->offset_a;
  std::vector<u32> cur_b = plan->offset_b;

  for (u32 j = 0; j < t; ++j) {
    const u32 info = slot_info[j];
    const u32 eid = info >> 1;
    const u32 pat = info & 1U;
    if (pat == 0) {
      plan->slots_a[cur_a[eid]++] = j;
    } else {
      plan->slots_b[cur_b[eid]++] = j;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Wrapper enumerator around join::PlaneSweepJoinStream so it satisfies
// baselines::IJoinEnumerator (used by the experimental framework).
// -----------------------------------------------------------------------------

template <int Dim, class T>
class PlaneSweepEnumeratorWrapper final : public baselines::IJoinEnumerator {
 public:
  PlaneSweepEnumeratorWrapper(const Relation<Dim, T>& R,
                              const Relation<Dim, T>& S,
                              join::PlaneSweepOptions opt)
      : stream_(R, S, opt) {}

  void Reset() override { stream_.Reset(); }
  bool Next(PairId* out) override { return stream_.Next(out); }
  const join::JoinStats& Stats() const noexcept override { return stream_.Stats(); }

 private:
  join::PlaneSweepJoinStream<Dim, T> stream_;
};

// -----------------------------------------------------------------------------
// OursReportJoinEnumerator2D
// -----------------------------------------------------------------------------
// A deterministic join enumerator that follows the same sweep + pattern
// decomposition used by "Our Method".
//
// Why this exists:
//   The project-wide EnumSampling and Adaptive runners rely on baseline->Enumerate()
//   and may perform one or two passes over the join stream. Returning a generic
//   PlaneSweep enumerator here can be catastrophically slow on adversarial cases
//   (e.g., large active sets), even if the baseline's Count/Sample logic is fast.
//
// This enumerator uses the same ActiveIndex2D (stabbing + range-point segtree)
// to *directly* report true intersections, avoiding O(|active|) scans.

template <int Dim, class T>
class OursReportJoinEnumerator2D final : public baselines::IJoinEnumerator {
 public:
  using Ctx = Ours2DContext<Dim, T>;

  explicit OursReportJoinEnumerator2D(Ctx* ctx) : ctx_(ctx) {
    static_assert(Dim == 2, "OursReportJoinEnumerator2D is only implemented for Dim=2");
    SJS_DASSERT(ctx_ != nullptr);
    Reset();
  }

  void Reset() override {
    SJS_DASSERT(ctx_ != nullptr);
    ctx_->ResetActive();

    pos_ = 0;
    stage_ = Stage::Scan;
    tmp_.clear();
    tmp_i_ = 0;

    active_r_ = 0;
    active_s_ = 0;

    stats_.Reset();
    stats_.num_events = static_cast<u64>(ctx_->events().size());
  }

  bool Next(PairId* out) override {
    SJS_DASSERT(out != nullptr);
    SJS_DASSERT(ctx_ != nullptr);

    const auto* ds = ctx_->dataset();
    SJS_DASSERT(ds != nullptr);

    auto& ar = ctx_->active_r();
    auto& as = ctx_->active_s();

    const auto& events = ctx_->events();
    const auto& ylo_r = ctx_->ylo_rank_r();
    const auto& yhi_r = ctx_->yhi_lb_rank_r();
    const auto& ylo_s = ctx_->ylo_rank_s();
    const auto& yhi_s = ctx_->yhi_lb_rank_s();

    while (true) {
      // If we're in the middle of emitting pairs for the current START event.
      if (stage_ == Stage::Emit) {
        if (tmp_i_ < tmp_.size()) {
          const u32 oh = tmp_[tmp_i_++];
          if (q_is_r_) {
            *out = PairId{ds->R.GetId(q_idx_), ds->S.GetId(oh)};
          } else {
            *out = PairId{ds->R.GetId(oh), ds->S.GetId(q_idx_)};
          }
          ++stats_.output_pairs;
          ++stats_.candidate_checks;
          return true;
        }

        // Finished current buffer; advance to the next pattern or finalize START.
        if (pat_ == Pattern::A) {
          // Switch to Pattern B.
          pat_ = Pattern::B;
          tmp_.clear();
          const detail::ActiveIndex2D& other = q_is_r_ ? as : ar;
          other.ReportB(q_ylo_, q_yhi_, &tmp_);
          tmp_i_ = 0;
          continue;
        }

        // Finished Pattern B: insert q and move on.
        if (q_is_r_) {
          ar.Insert(q_handle_, q_ylo_, q_yhi_);
          ++active_r_;
          if (active_r_ > stats_.active_max_r) stats_.active_max_r = active_r_;
        } else {
          as.Insert(q_handle_, q_ylo_, q_yhi_);
          ++active_s_;
          if (active_s_ > stats_.active_max_s) stats_.active_max_s = active_s_;
        }

        stage_ = Stage::Scan;
        ++pos_;
        continue;
      }

      // Scan events until we find the next START (or reach end).
      if (pos_ >= events.size()) {
        // Leave ctx_ in a clean state for safety.
        ctx_->ResetActive();
        return false;
      }

      const auto& ev = events[pos_];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) {
          ar.Erase(handle);
          SJS_DASSERT(active_r_ > 0);
          --active_r_;
        } else {
          as.Erase(handle);
          SJS_DASSERT(active_s_ > 0);
          --active_s_;
        }
        ++pos_;
        continue;
      }

      // START event: prepare Pattern A buffer, but do NOT insert q yet.
      q_is_r_ = (ev.side == join::Side::R);
      q_idx_ = static_cast<u32>(ev.index);
      q_handle_ = handle;
      q_ylo_ = q_is_r_ ? ylo_r[ev.index] : ylo_s[ev.index];
      q_yhi_ = q_is_r_ ? yhi_r[ev.index] : yhi_s[ev.index];

      pat_ = Pattern::A;
      tmp_.clear();
      const detail::ActiveIndex2D& other = q_is_r_ ? as : ar;
      other.ReportA(q_ylo_, &tmp_);
      tmp_i_ = 0;
      stage_ = Stage::Emit;
      // Loop back to emit from tmp_.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  enum class Stage : u8 {
    Scan = 0,
    Emit = 1,
  };

  enum class Pattern : u8 {
    A = 0,
    B = 1,
  };

  Ctx* ctx_{nullptr};

  // Iteration state.
  usize pos_{0};
  Stage stage_{Stage::Scan};
  Pattern pat_{Pattern::A};

  // Current START event (q).
  bool q_is_r_{true};
  u32 q_idx_{0};
  u32 q_handle_{0};
  u32 q_ylo_{0};
  u32 q_yhi_{0};

  // Buffer of opposite-side handles for the current pattern.
  std::vector<u32> tmp_;
  usize tmp_i_{0};

  // Lightweight stats.
  join::JoinStats stats_;
  u64 active_r_{0};
  u64 active_s_{0};
};

}  // namespace detail

// -----------------------------------------------------------------------------
// OursSamplingBaseline (Sampling variant, Dim=2)
// -----------------------------------------------------------------------------

template <int Dim, class T = Scalar>
class OursSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "ours_sampling"; }

  void Reset() override {
    ctx_.Reset();
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
    count_r_.Clear();
    count_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if (!ctx_.Build(ds, phases, err)) {
      return false;
    }

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    // Count-only BIT indices for Phase 1.
    const u32 m = ctx_.num_ylo_ranks();
    count_r_.Init(m);
    count_s_.Init(m);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    // Count-only active sets (two BITs per side).
    count_r_.ResetActive();
    count_s_.ResetActive();

    u64 W = 0;

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) {
          count_r_.Erase(ylo_r[ev.index], yhi_r[ev.index]);
        } else {
          count_s_.Erase(ylo_s[ev.index], yhi_s[ev.index]);
        }
        continue;
      }

      // START event.
      const i32 sid_i32 = sid_of_pos[pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 sid = static_cast<u32>(sid_i32);

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
      const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

      const detail::CountIndex2D& other = q_is_r ? count_s_ : count_r_;

      const u64 wa = other.CountA(q_ylo);
      const u64 wb = other.CountB(q_ylo, q_yhi);
      const u64 w = wa + wb;

      w_a_[sid] = wa;
      w_b_[sid] = wb;
      w_total_[sid] = w;

      W += w;

      // Insert q into its side.
      if (q_is_r) {
        count_r_.Insert(q_ylo, q_yhi);
      } else {
        count_s_.Insert(q_ylo, q_yhi);
      }
    }

    W_ = W;
    weights_valid_ = true;

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    // Ensure Phase1 weights exist.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // Phase2: slot plan.
    detail::SlotPlan2D plan;
    {
      auto scoped = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!detail::BuildSlotPlan2D(t, rng, w_total_, w_a_, w_b_, &plan, err)) {
        if (err && err->empty()) *err = "OursSamplingBaseline::Sample: failed to build slot plan";
        return false;
      }
    }

    out->pairs.resize(t);

    // Phase3: second sweep.
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      ctx_.ResetActive();

      auto& ar = ctx_.active_r();
      auto& as = ctx_.active_s();

      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();

      const auto& ylo_r = ctx_.ylo_rank_r();
      const auto& yhi_r = ctx_.yhi_lb_rank_r();
      const auto& ylo_s = ctx_.ylo_rank_s();
      const auto& yhi_s = ctx_.yhi_lb_rank_s();

      std::vector<u32> sampled;

      const auto* ds = ctx_.dataset();
      SJS_DASSERT(ds != nullptr);

      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          if (ev.side == join::Side::R) {
            ar.Erase(handle);
          } else {
            as.Erase(handle);
          }
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        const bool q_is_r = (ev.side == join::Side::R);
        const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
        const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

        const detail::ActiveIndex2D& other = q_is_r ? as : ar;

        // Pattern A slots.
        {
          const u32 begin = plan.offset_a[sid];
          const u32 end = plan.offset_a[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleA(q_ylo, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursSamplingBaseline::Sample: SampleA failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_a[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) {
                out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              } else {
                out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
              }
            }
          }
        }

        // Pattern B slots.
        {
          const u32 begin = plan.offset_b[sid];
          const u32 end = plan.offset_b[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleB(q_ylo, q_yhi, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursSamplingBaseline::Sample: SampleB failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_b[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) {
                out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              } else {
                out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
              }
            }
          }
        }

        // Insert q into its side.
        if (q_is_r) {
          ar.Insert(handle, q_ylo, q_yhi);
        } else {
          as.Insert(handle, q_ylo, q_yhi);
        }
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    // IMPORTANT: return the report-based enumerator (not the generic plane sweep).
    // EnumSampling / Adaptive runners depend on Enumerate() and would otherwise
    // measure an unrelated and potentially adversarially slow algorithm.
    return std::make_unique<detail::OursReportJoinEnumerator2D<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};

  detail::Ours2DContext<Dim, T> ctx_;

  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  detail::CountIndex2D count_r_;
  detail::CountIndex2D count_s_;


  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs

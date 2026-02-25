#pragma once
// sjs/baselines/range_tree/adaptive.h
//
// RangeTree baseline — Adaptive variant (SJS v3 Framework III).
//
// Implements Algorithm 3.4 (Framework III) in SJS v3 while using Chapter 4's
// orthogonal range-tree event-block primitives:
//   - Full-cache small blocks (REPORT once) under budget
//   - Prefetch sample cache via a global heap of valued "slots" (SAMPLE once)
//   - Optional second sweep only for residual positions
//
// Knobs:
//   - Budget B: cfg.run.j_star (alias extra["budget"]).
//   - w_small: cfg.run.extra["w_small"] (default 0 -> no caching).

#include "sjs/baselines/range_tree/sampling.h"

#include "sjs/baselines/detail/adaptive_prefetch.h"

#include "sjs/core/assert.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace range_tree {

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

}  // namespace detail

template <int Dim, class T = Scalar>
class RangeTreeAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeAdaptiveBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "range_tree_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    cache_valid_ = false;
    W_ = 0;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    start_side_.clear();
    start_index_.clear();

    cached_.clear();
    cache_off_.clear();
    cache_len_.clear();
    cache_partners_.clear();

    prefetch_keep_.clear();
    prefetch_partners_.clear();

    budget_B_ = 0;
    budget_used_ = 0;
    w_small_ = 0;

    embed_.Reset();
    rt_r_.Clear();
    rt_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    // Events.
    {
      auto _ = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // START ids.
    {
      auto _ = phases ? phases->Scoped("build_start_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      start_id_of_event_.assign(events_.size(), kInvalidStartId);
      u32 sid = 0;
      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == join::EventKind::Start) start_id_of_event_[i] = sid++;
      }
      w_total_.assign(static_cast<usize>(sid), 0ULL);
      start_side_.assign(static_cast<usize>(sid), join::Side::R);
      start_index_.assign(static_cast<usize>(sid), 0U);
      cached_.assign(static_cast<usize>(sid), 0U);
      cache_off_.assign(static_cast<usize>(sid), 0ULL);
      cache_len_.assign(static_cast<usize>(sid), 0ULL);

      prefetch_keep_.assign(static_cast<usize>(sid), 0U);
      prefetch_partners_.assign(static_cast<usize>(sid), std::vector<u32>{});

      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind != join::EventKind::Start) continue;
        const u32 s = start_id_of_event_[i];
        if (s == kInvalidStartId) continue;
        start_side_[static_cast<usize>(s)] = events_[i].side;
        start_index_[static_cast<usize>(s)] = static_cast<u32>(events_[i].index);
      }
    }

    // Embedding + range trees.
    {
      std::vector<detail::RankPoint<K>> pts_r;
      std::vector<detail::RankPoint<K>> pts_s;
      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) return false;
      auto _ = phases ? phases->Scoped("build_range_trees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!rt_r_.Build(pts_r, err)) return false;
      if (!rt_s_.Build(pts_s, err)) return false;
    }

    built_ = true;
    weights_valid_ = false;
    cache_valid_ = false;
    W_ = 0;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_and_cache") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Budget and threshold.
    budget_B_ = detail::ExtraU64Or(cfg, "budget", cfg.run.j_star);
    w_small_ = detail::ExtraU64Or(cfg, "w_small", 0ULL);

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(cached_.begin(), cached_.end(), 0U);
    std::fill(cache_off_.begin(), cache_off_.end(), 0ULL);
    std::fill(cache_len_.begin(), cache_len_.end(), 0ULL);
    cache_partners_.clear();

    std::fill(prefetch_keep_.begin(), prefetch_keep_.end(), 0U);
    for (auto& v : prefetch_partners_) v.clear();

    // Prefetch is optional. If rng is null, we skip it.
    const bool enable_prefetch = (rng != nullptr) && (budget_B_ > 0) && (cfg.run.t > 0);
    const u64 base_prefetch_seed = enable_prefetch ? rng->NextU64() : 0ULL;

    baselines::detail::PrefetchHeap heap;
    u64 mem_full = 0;

    rt_r_.ResetToEmpty();
    rt_s_.ResetToEmpty();

    u64 W = 0;
    long double W_sofar = 0.0L;

    std::vector<u32> hits;
    hits.reserve(1024);

    const usize E = w_total_.size();

    for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
      const join::Event& e = events_[ev_i];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
        else rt_s_.Deactivate(static_cast<u32>(e.index));
        continue;
      }

      // START
      const u32 sid = start_id_of_event_[ev_i];
      SJS_DASSERT(sid != kInvalidStartId);
      const usize sid_u = static_cast<usize>(sid);
      const u32 i = sid + 1;

      const BoxT& q = (e.side == join::Side::R)
                          ? ds_->R.boxes[static_cast<usize>(e.index)]
                          : ds_->S.boxes[static_cast<usize>(e.index)];

      detail::RankBox<K> range;
      u64 w = 0;
      if (embed_.MakeQueryRange(q, &range)) {
        if (e.side == join::Side::R) w = rt_s_.Count(range);
        else w = rt_r_.Count(range);
      }
      w_total_[sid_u] = w;

      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "RangeTreeAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }
      W_sofar += static_cast<long double>(w);

      // Full-cache small blocks.
      const bool can_full_cache = (w_small_ > 0 && w > 0 && w <= w_small_ && (mem_full + w) <= budget_B_);
      if (can_full_cache) {
        if (!embed_.MakeQueryRange(q, &range)) {
          if (err) *err = "RangeTreeAdaptiveBaseline::Count: empty query range for cached block";
          return false;
        }
        cached_[sid_u] = 1U;
        cache_off_[sid_u] = static_cast<u64>(cache_partners_.size());
        hits.clear();
        if (e.side == join::Side::R) rt_s_.Report(range, &hits, /*st=*/nullptr);
        else rt_r_.Report(range, &hits, /*st=*/nullptr);
        cache_partners_.insert(cache_partners_.end(), hits.begin(), hits.end());
        const u64 off = cache_off_[sid_u];
        const u64 len = static_cast<u64>(cache_partners_.size()) - off;
        cache_len_[sid_u] = len;
        if (len != w) {
          if (err) *err = "RangeTreeAdaptiveBaseline::Count: cache REPORT size mismatch (bug)";
          return false;
        }
        mem_full += w;

        const u64 B_samp = (budget_B_ > mem_full) ? (budget_B_ - mem_full) : 0ULL;
        while (heap.Size() > static_cast<usize>(B_samp)) {
          const auto popped = heap.PopMin();
          const usize j = static_cast<usize>(popped.sid);
          if (j < prefetch_keep_.size() && prefetch_keep_[j] > 0) {
            --prefetch_keep_[j];
            if (!prefetch_partners_[j].empty()) prefetch_partners_[j].pop_back();
          }
        }
      } else if (enable_prefetch && w > 0) {
        const u64 B_samp = (budget_B_ > mem_full) ? (budget_B_ - mem_full) : 0ULL;
        if (B_samp > 0) {
          while (true) {
            const u32 r = prefetch_keep_[sid_u] + 1U;
            const double score = baselines::detail::SlotScorePoisson(w, static_cast<u32>(E), i, W_sofar, cfg.run.t, r);
            if (heap.Size() < static_cast<usize>(B_samp)) {
              heap.Push(baselines::detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[sid_u];
              continue;
            }
            if (score > heap.MinScore()) {
              heap.Push(baselines::detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[sid_u];
              const auto popped = heap.PopMin();
              const usize j = static_cast<usize>(popped.sid);
              if (j < prefetch_keep_.size() && prefetch_keep_[j] > 0) {
                --prefetch_keep_[j];
                if (!prefetch_partners_[j].empty()) prefetch_partners_[j].pop_back();
              }
              continue;
            }
            break;
          }

          const u32 s_keep = prefetch_keep_[sid_u];
          if (s_keep > 0 && embed_.MakeQueryRange(q, &range)) {
            const u64 ev_seed = DeriveSeed(base_prefetch_seed, 0xA11C3EULL, static_cast<u64>(sid));
            Rng rng_samp(DeriveSeed(ev_seed, 1));
            auto& dst = prefetch_partners_[sid_u];
            std::string local_err;
            const bool ok = (e.side == join::Side::R)
                                ? rt_s_.Sample(range, s_keep, &rng_samp, &dst, &local_err)
                                : rt_r_.Sample(range, s_keep, &rng_samp, &dst, &local_err);
            if (!ok || dst.size() != s_keep) {
              if (err) *err = local_err.empty() ? "RangeTreeAdaptiveBaseline::Count: prefetch Sample failed" : local_err;
              return false;
            }
          }
        }
      }

      // Activate current.
      if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
      else rt_s_.Activate(static_cast<u32>(e.index));
    }

    W_ = W;
    weights_valid_ = true;
    cache_valid_ = true;
    budget_used_ = mem_full + static_cast<u64>(heap.Size());
    if (budget_used_ > budget_B_) budget_used_ = budget_B_;

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: null rng/out";
      return false;
    }
    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: t too large for u32";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();
    if (t == 0) return true;

    if (!weights_valid_ || !cache_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/rng, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    // Phase 2: assign slots to START events by w_e/|J|.
    struct SlotAssign { u32 sid; u32 slot; };
    auto by_sid_slot = [](const SlotAssign& a, const SlotAssign& b) {
      if (a.sid < b.sid) return true;
      if (b.sid < a.sid) return false;
      return a.slot < b.slot;
    };

    sampling::AliasTable alias;
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) return false;
    }

    // Split RNG streams within Sample for determinism.
    const u64 seed_assign = rng->NextU64();
    const u64 seed_cache = rng->NextU64();
    const u64 seed_sweep = rng->NextU64();

    Rng rng_assign(seed_assign);

    std::vector<SlotAssign> asg;
    asg.reserve(static_cast<usize>(t));
    for (u32 slot = 0; slot < t; ++slot) {
      const u32 sid = static_cast<u32>(alias.Sample(&rng_assign));
      asg.push_back(SlotAssign{sid, slot});
    }
    std::sort(asg.begin(), asg.end(), by_sid_slot);

    out->pairs.assign(static_cast<usize>(t), PairId{});

    // Residual (event,slot) assignments for pass 2.
    std::vector<SlotAssign> asg_need;
    asg_need.reserve(asg.size());

    // Phase 2 fill: full cache -> prefetched samples -> residual.
    {
      auto _ = phases ? phases->Scoped("phase2_fill") : PhaseRecorder::ScopedPhase(nullptr, "");
      usize ptr = 0;
      while (ptr < asg.size()) {
        const u32 sid = asg[ptr].sid;
        usize end = ptr + 1;
        while (end < asg.size() && asg[end].sid == sid) ++end;

        const usize sid_u = static_cast<usize>(sid);
        const u32 k = static_cast<u32>(end - ptr);

        const join::Side side = start_side_[sid_u];
        const u32 q_idx = start_index_[sid_u];
        const Id q_id = (side == join::Side::R) ? ds_->R.GetId(static_cast<usize>(q_idx))
                                                : ds_->S.GetId(static_cast<usize>(q_idx));

        if (cached_[sid_u]) {
          const u64 len = cache_len_[sid_u];
          const u64 off = cache_off_[sid_u];
          if (len == 0) {
            if (err) *err = "RangeTreeAdaptiveBaseline::Sample: cached sid has empty cache";
            return false;
          }
          Rng crng(DeriveSeed(seed_cache, 0xCA11CEULL, static_cast<u64>(sid)));
          for (u32 j = 0; j < k; ++j) {
            const u32 slot = asg[ptr + j].slot;
            const u64 pick = crng.UniformU64(len);
            const u32 oh = cache_partners_[static_cast<usize>(off + pick)];
            if (side == join::Side::R) out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds_->S.GetId(static_cast<usize>(oh))};
            else out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(oh)), q_id};
          }
          ptr = end;
          continue;
        }

        const auto& pref = prefetch_partners_[sid_u];
        const u32 s = static_cast<u32>(pref.size());
        const u32 a = (s < k) ? s : k;
        for (u32 j = 0; j < a; ++j) {
          const u32 slot = asg[ptr + j].slot;
          const u32 oh = pref[static_cast<usize>(j)];
          if (side == join::Side::R) out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds_->S.GetId(static_cast<usize>(oh))};
          else out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(oh)), q_id};
        }
        for (u32 j = a; j < k; ++j) {
          asg_need.push_back(SlotAssign{sid, asg[ptr + j].slot});
        }
        ptr = end;
      }
    }

    if (asg_need.empty()) return true;

    // Build residual offsets.
    const usize E = w_total_.size();
    std::sort(asg_need.begin(), asg_need.end(), by_sid_slot);
    std::vector<u32> off;
    off.assign(E + 1, 0U);
    for (const auto& x : asg_need) off[static_cast<usize>(x.sid) + 1]++;
    for (usize i = 0; i < E; ++i) off[i + 1] += off[i];
    std::vector<u32> slots;
    slots.resize(asg_need.size());
    {
      std::vector<u32> cur = off;
      for (const auto& x : asg_need) {
        const u32 pos = cur[static_cast<usize>(x.sid)]++;
        slots[static_cast<usize>(pos)] = x.slot;
      }
    }

    // Pass 2: second sweep for residual positions.
    {
      auto _ = phases ? phases->Scoped("phase3_fill_residual") : PhaseRecorder::ScopedPhase(nullptr, "");
      rt_r_.ResetToEmpty();
      rt_s_.ResetToEmpty();

      for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
        const join::Event& e = events_[ev_i];
        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
          else rt_s_.Deactivate(static_cast<u32>(e.index));
          continue;
        }

        const u32 sid = start_id_of_event_[ev_i];
        SJS_DASSERT(sid != kInvalidStartId);
        const usize sid_u = static_cast<usize>(sid);
        const u32 begin = off[sid_u];
        const u32 end = off[sid_u + 1];
        const u32 k = end - begin;

        if (k > 0) {
          const BoxT& q = (e.side == join::Side::R)
                              ? ds_->R.boxes[static_cast<usize>(e.index)]
                              : ds_->S.boxes[static_cast<usize>(e.index)];
          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "RangeTreeAdaptiveBaseline::Sample: empty query range for START with residual slots";
            return false;
          }

          Rng ev_rng(DeriveSeed(seed_sweep, 0x5EEDULL, static_cast<u64>(sid)));
          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(k));
          if (e.side == join::Side::R) {
            if (!rt_s_.Sample(range, k, &ev_rng, &picked, err)) return false;
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = slots[static_cast<usize>(begin + j)];
              const u32 s_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                          ds_->S.GetId(static_cast<usize>(s_idx))};
            }
          } else {
            if (!rt_r_.Sample(range, k, &ev_rng, &picked, err)) return false;
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = slots[static_cast<usize>(begin + j)];
              const u32 r_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(r_idx)),
                                                          ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }
        }

        if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
        else rt_s_.Activate(static_cast<u32>(e.index));
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr u32 kInvalidStartId = std::numeric_limits<u32>::max();

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  std::vector<join::Event> events_;
  std::vector<u32> start_id_of_event_;

  std::vector<u64> w_total_;
  u64 W_ = 0;
  bool weights_valid_ = false;

  // sid -> (side,index) for direct cached fill.
  std::vector<join::Side> start_side_;
  std::vector<u32> start_index_;

  // Cache: partner indices in the opposite relation.
  std::vector<u8> cached_;
  std::vector<u64> cache_off_;
  std::vector<u64> cache_len_;
  std::vector<u32> cache_partners_;

  // Prefetch sample cache S_i (Framework III): for each START event i,
  // store a prefix of i.i.d. partner samples (with replacement).
  std::vector<u32> prefetch_keep_;
  std::vector<std::vector<u32>> prefetch_partners_;
  bool cache_valid_ = false;

  u64 budget_B_ = 0;
  u64 budget_used_ = 0;
  u64 w_small_ = 0;

  detail::GlobalRankEmbedding<Dim, T> embed_;
  detail::ActiveRangeTree<K> rt_r_;
  detail::ActiveRangeTree<K> rt_s_;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs

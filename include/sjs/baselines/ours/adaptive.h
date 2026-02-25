#pragma once
// sjs/baselines/ours/adaptive.h
//
// "Our Method" — Adaptive variant (SJS v3 Framework III).
//
// Implements Algorithm 3.4 (Adaptive Budgeted Join Sampling) in SJS v3:
//   - Pass 1 (sweep): compute exact per-event weights w_e (and pattern split
//     w_e^A/w_e^B), and cache under budget B:
//       * full cache C_i for small blocks (REPORT once, store all partners)
//       * prefetch sample cache S_i via a global heap of valued "slots"
//         (SAMPLE once per event, keep a prefix of i.i.d. samples)
//   - Planning: assign t output positions to event indices with Pr(I=i)=w_i/W.
//   - Fill: prefer C_i, then use prefetched samples S_i, then record residual.
//   - Pass 2 (optional): second sweep to generate residual samples only.
//
// Knobs (read from Config, aligned with SJS v3 §3.5):
//   - Budget B (max cached partner records): cfg.run.j_star (alias "budget").
//   - Small-block threshold w_small: cfg.run.extra["w_small"] (default 0).
//
// Setting w_small=0 disables caching and makes this variant behave like the
// Sampling variant (Framework II).  This does not change the output
// distribution; it only affects performance.

#include "sjs/baselines/ours/sampling.h"  // Ours2DContext + ActiveIndex2D + enumerator

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
namespace ours {

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
class OursAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "ours_adaptive"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();

    start_side_.clear();
    start_index_.clear();

    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
    W_ = 0;
    weights_valid_ = false;

    cached_.clear();
    cache_off_.clear();
    cache_len_.clear();
    cache_partners_.clear();

    prefetch_keep_.clear();
    prefetch_partners_.clear();

    cache_valid_ = false;
    budget_B_ = 0;
    budget_used_ = 0;
    w_small_ = 0;
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    Reset();
    if (!ctx_.Build(ds, phases, err)) return false;

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    start_side_.assign(E, join::Side::R);
    start_index_.assign(E, 0U);

    // Build sid -> (side,index) mapping.
    {
      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();
      for (usize pos = 0; pos < events.size(); ++pos) {
        if (events[pos].kind != join::EventKind::Start) continue;
        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);
        if (sid < static_cast<u32>(E)) {
          start_side_[static_cast<usize>(sid)] = events[pos].side;
          start_index_[static_cast<usize>(sid)] = static_cast<u32>(events[pos].index);
        }
      }
    }

    cached_.assign(E, 0U);
    cache_off_.assign(E, 0ULL);
    cache_len_.assign(E, 0ULL);
    cache_partners_.clear();

    prefetch_keep_.assign(E, 0U);
    prefetch_partners_.assign(E, std::vector<u32>{});

    prefetch_keep_.assign(E, 0U);
    prefetch_partners_.assign(E, {});

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_and_cache") : PhaseRecorder::ScopedPhase(nullptr, "");

    const usize E = w_total_.size();
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);
    std::fill(cached_.begin(), cached_.end(), 0U);
    std::fill(cache_off_.begin(), cache_off_.end(), 0ULL);
    std::fill(cache_len_.begin(), cache_len_.end(), 0ULL);
    cache_partners_.clear();

    std::fill(prefetch_keep_.begin(), prefetch_keep_.end(), 0U);
    for (auto& v : prefetch_partners_) v.clear();

    // Budget and threshold.
    budget_B_ = detail::ExtraU64Or(cfg, "budget", cfg.run.j_star);
    w_small_ = detail::ExtraU64Or(cfg, "w_small", 0ULL);

    // Prefetch is optional (performance only). If rng is null, we skip it.
    const bool enable_prefetch = (rng != nullptr) && (budget_B_ > 0) && (cfg.run.t > 0);

    // Split a base seed for per-event prefetch RNGs to avoid inter-event coupling.
    const u64 base_prefetch_seed = enable_prefetch ? rng->NextU64() : 0ULL;

    baselines::detail::PrefetchHeap heap;
    u64 mem_full = 0;

    ctx_.ResetActive();

    u64 W = 0;
    long double W_sofar = 0.0L;

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();
    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    std::vector<u32> tmp;
    tmp.reserve(1024);

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) ar.Erase(handle);
        else as.Erase(handle);
        continue;
      }

      // START
      const i32 sid_i32 = sid_of_pos[pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 sid = static_cast<u32>(sid_i32);
      if (sid >= static_cast<u32>(E)) {
        if (err) *err = "OursAdaptiveBaseline::Count: sid out of range";
        return false;
      }
      const u32 i = sid + 1;  // 1-based START index in event order

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
      const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];
      const detail::ActiveIndex2D& other = q_is_r ? as : ar;

      const u64 wa = other.CountA(q_ylo);
      const u64 wb = other.CountB(q_ylo, q_yhi);
      const u64 w = wa + wb;

      w_a_[static_cast<usize>(sid)] = wa;
      w_b_[static_cast<usize>(sid)] = wb;
      w_total_[static_cast<usize>(sid)] = w;

      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "OursAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }
      W_sofar += static_cast<long double>(w);

      // -----------------
      // Optional: full cache for small blocks (REPORT once).
      // -----------------
      const bool can_full_cache = (w_small_ > 0 && w > 0 && w <= w_small_ && (mem_full + w) <= budget_B_);
      if (can_full_cache) {
        cached_[static_cast<usize>(sid)] = 1U;
        cache_off_[static_cast<usize>(sid)] = static_cast<u64>(cache_partners_.size());

        tmp.clear();
        other.ReportA(q_ylo, &tmp);
        cache_partners_.insert(cache_partners_.end(), tmp.begin(), tmp.end());
        tmp.clear();
        other.ReportB(q_ylo, q_yhi, &tmp);
        cache_partners_.insert(cache_partners_.end(), tmp.begin(), tmp.end());

        const u64 off = cache_off_[static_cast<usize>(sid)];
        const u64 len = static_cast<u64>(cache_partners_.size()) - off;
        cache_len_[static_cast<usize>(sid)] = len;
        if (len != w) {
          if (err) *err = "OursAdaptiveBaseline::Count: cache REPORT size mismatch (bug)";
          return false;
        }

        mem_full += w;

        // Shrink prefetch heap capacity if needed.
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
        // -----------------
        // Prefetch sample cache under remaining capacity.
        // -----------------
        const u64 B_samp = (budget_B_ > mem_full) ? (budget_B_ - mem_full) : 0ULL;
        if (B_samp > 0) {
          while (true) {
            const u32 r = prefetch_keep_[static_cast<usize>(sid)] + 1U;
            const double score = baselines::detail::SlotScorePoisson(
                w,
                static_cast<u32>(E),
                i,
                W_sofar,
                cfg.run.t,
                r);

            if (heap.Size() < static_cast<usize>(B_samp)) {
              heap.Push(baselines::detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[static_cast<usize>(sid)];
              continue;
            }

            // Heap is full.
            if (score > heap.MinScore()) {
              heap.Push(baselines::detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[static_cast<usize>(sid)];
              const auto popped = heap.PopMin();
              const usize j = static_cast<usize>(popped.sid);
              if (j < prefetch_keep_.size() && prefetch_keep_[j] > 0) {
                --prefetch_keep_[j];
                if (!prefetch_partners_[j].empty()) prefetch_partners_[j].pop_back();
              }
              continue;
            }

            // Scores are non-increasing in r, so we can stop.
            break;
          }

          const u32 s_keep = prefetch_keep_[static_cast<usize>(sid)];
          if (s_keep > 0) {
            // Generate i.i.d. samples on K_e (with replacement) using the
            // same primitives as Framework II.
            const u64 ev_seed = DeriveSeed(base_prefetch_seed, 0xA11C3EULL, static_cast<u64>(sid));
            Rng rng_pat(DeriveSeed(ev_seed, 1));
            Rng rng_samp(DeriveSeed(ev_seed, 2));

            std::vector<u8> pat;
            pat.resize(s_keep);
            u32 kA = 0, kB = 0;
            for (u32 j = 0; j < s_keep; ++j) {
              u8 p = 0;
              if (wa == 0) p = 1;
              else if (wb == 0) p = 0;
              else p = (rng_pat.UniformU64(w) < wa) ? 0 : 1;
              pat[static_cast<usize>(j)] = p;
              if (p == 0) ++kA; else ++kB;
            }

            std::vector<u32> sampA;
            std::vector<u32> sampB;
            if (kA > 0) {
              if (!other.SampleA(q_ylo, kA, &rng_samp, &sampA) || sampA.size() != kA) {
                if (err) *err = "OursAdaptiveBaseline::Count: prefetch SampleA failed";
                return false;
              }
            }
            if (kB > 0) {
              if (!other.SampleB(q_ylo, q_yhi, kB, &rng_samp, &sampB) || sampB.size() != kB) {
                if (err) *err = "OursAdaptiveBaseline::Count: prefetch SampleB failed";
                return false;
              }
            }

            auto& dst = prefetch_partners_[static_cast<usize>(sid)];
            dst.clear();
            dst.reserve(s_keep);
            u32 ia = 0, ib = 0;
            for (u32 j = 0; j < s_keep; ++j) {
              if (pat[static_cast<usize>(j)] == 0) dst.push_back(sampA[static_cast<usize>(ia++)]);
              else dst.push_back(sampB[static_cast<usize>(ib++)]);
            }
          }
        }
      }

      // Insert q.
      if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
      else as.Insert(handle, q_ylo, q_yhi);
    }

    ctx_.ResetActive();

    W_ = W;
    weights_valid_ = true;
    cache_valid_ = true;
    budget_used_ = mem_full + static_cast<u64>(heap.Size());
    if (budget_used_ > budget_B_) budget_used_ = budget_B_;  // defensive clamp

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursAdaptiveBaseline::Sample: null rng/out";
      return false;
    }
    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursAdaptiveBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    if (t == 0) return true;

    // Ensure pass-1 state.
    if (!weights_valid_ || !cache_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/rng, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    // -----------------
    // Phase 2: assign t output positions to START events by exact weights.
    // -----------------
    struct SlotAssign {
      u32 sid;
      u32 slot;
    };
    auto by_sid_slot = [](const SlotAssign& a, const SlotAssign& b) {
      if (a.sid < b.sid) return true;
      if (b.sid < a.sid) return false;
      return a.slot < b.slot;
    };

    const usize E = w_total_.size();
    sampling::AliasTable alias;
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) {
        if (err && err->empty()) *err = "OursAdaptiveBaseline::Sample: failed to build alias";
        return false;
      }
    }

    const u64 seed_assign = rng->NextU64();
    const u64 seed_cache = rng->NextU64();
    const u64 seed_pat = rng->NextU64();
    const u64 seed_sweep = rng->NextU64();

    Rng rng_assign(seed_assign);

    std::vector<SlotAssign> asg;
    asg.reserve(static_cast<usize>(t));
    for (u32 j = 0; j < t; ++j) {
      const u32 sid = static_cast<u32>(alias.Sample(&rng_assign));
      asg.push_back(SlotAssign{sid, j});
    }
    std::sort(asg.begin(), asg.end(), by_sid_slot);

    out->pairs.assign(static_cast<usize>(t), PairId{});

    const auto* ds = ctx_.dataset();
    SJS_DASSERT(ds != nullptr);

    // Collect residual (event,slot) assignments split by pattern for pass 2.
    std::vector<SlotAssign> asg_a;
    std::vector<SlotAssign> asg_b;
    asg_a.reserve(asg.size());
    asg_b.reserve(asg.size());

    // Fill from full caches / prefetched samples; record residual.
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
        const Id q_id = (side == join::Side::R) ? ds->R.GetId(static_cast<usize>(q_idx))
                                                : ds->S.GetId(static_cast<usize>(q_idx));

        if (cached_[sid_u]) {
          const u64 len = cache_len_[sid_u];
          const u64 off = cache_off_[sid_u];
          if (len == 0) {
            if (err) *err = "OursAdaptiveBaseline::Sample: cached sid has empty cache";
            return false;
          }
          Rng crng(DeriveSeed(seed_cache, 0xCA11CEULL, static_cast<u64>(sid)));
          for (u32 j = 0; j < k; ++j) {
            const u32 slot = asg[ptr + j].slot;
            const u64 pick = crng.UniformU64(len);
            const u32 oh = cache_partners_[static_cast<usize>(off + pick)];
            if (side == join::Side::R) {
              out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds->S.GetId(static_cast<usize>(oh))};
            } else {
              out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(static_cast<usize>(oh)), q_id};
            }
          }
          ptr = end;
          continue;
        }

        // Not fully cached: use prefetched samples if any.
        const auto& pref = prefetch_partners_[sid_u];
        const u32 s = static_cast<u32>(pref.size());
        const u32 a = (s < k) ? s : k;
        for (u32 j = 0; j < a; ++j) {
          const u32 slot = asg[ptr + j].slot;
          const u32 oh = pref[static_cast<usize>(j)];
          if (side == join::Side::R) {
            out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds->S.GetId(static_cast<usize>(oh))};
          } else {
            out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(static_cast<usize>(oh)), q_id};
          }
        }

        if (a < k) {
          // Residual: decide pattern per slot in proportion to (w^A, w^B).
          const u64 wa = w_a_[sid_u];
          const u64 wb = w_b_[sid_u];
          const u64 w = wa + wb;
          if (w == 0) {
            if (err) *err = "OursAdaptiveBaseline::Sample: sampled an event with w_e==0 (unexpected)";
            return false;
          }
          Rng prng(DeriveSeed(seed_pat, 0xA7717ULL, static_cast<u64>(sid)));
          for (u32 j = a; j < k; ++j) {
            const u32 slot = asg[ptr + j].slot;
            u8 pat = 0;
            if (wa == 0) pat = 1;
            else if (wb == 0) pat = 0;
            else pat = (prng.UniformU64(w) < wa) ? 0 : 1;
            if (pat == 0) asg_a.push_back(SlotAssign{sid, slot});
            else asg_b.push_back(SlotAssign{sid, slot});
          }
        }

        ptr = end;
      }
    }

    if (asg_a.empty() && asg_b.empty()) {
      return true;  // one-pass completion
    }

    // -----------------
    // Build residual slot plan (A/B) for the second sweep.
    // -----------------
    detail::SlotPlan2D plan;
    plan.offset_a.assign(E + 1, 0U);
    plan.offset_b.assign(E + 1, 0U);

    std::sort(asg_a.begin(), asg_a.end(), by_sid_slot);
    std::sort(asg_b.begin(), asg_b.end(), by_sid_slot);
    for (const auto& x : asg_a) plan.offset_a[static_cast<usize>(x.sid) + 1]++;
    for (const auto& x : asg_b) plan.offset_b[static_cast<usize>(x.sid) + 1]++;
    for (usize i = 0; i < E; ++i) {
      plan.offset_a[i + 1] += plan.offset_a[i];
      plan.offset_b[i + 1] += plan.offset_b[i];
    }
    plan.slots_a.resize(asg_a.size());
    plan.slots_b.resize(asg_b.size());
    {
      std::vector<u32> cursor_a = plan.offset_a;
      for (const auto& x : asg_a) {
        const u32 pos = cursor_a[static_cast<usize>(x.sid)]++;
        plan.slots_a[static_cast<usize>(pos)] = x.slot;
      }
      std::vector<u32> cursor_b = plan.offset_b;
      for (const auto& x : asg_b) {
        const u32 pos = cursor_b[static_cast<usize>(x.sid)]++;
        plan.slots_b[static_cast<usize>(pos)] = x.slot;
      }
    }

    // -----------------
    // Pass 2: second sweep, generate residual samples only.
    // -----------------
    {
      auto _ = phases ? phases->Scoped("phase3_fill_residual") : PhaseRecorder::ScopedPhase(nullptr, "");
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
      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          if (ev.side == join::Side::R) ar.Erase(handle);
          else as.Erase(handle);
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);
        const usize sid_u = static_cast<usize>(sid);

        const bool q_is_r = (ev.side == join::Side::R);
        const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
        const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];
        const detail::ActiveIndex2D& other = q_is_r ? as : ar;

        const u64 ev_seed = DeriveSeed(seed_sweep, 0x5EEDULL, static_cast<u64>(sid));
        Rng rngA(DeriveSeed(ev_seed, 1));
        Rng rngB(DeriveSeed(ev_seed, 2));

        // Pattern A residual.
        {
          const u32 begin = plan.offset_a[sid_u];
          const u32 end = plan.offset_a[sid_u + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            if (!other.SampleA(q_ylo, k, &rngA, &sampled) || sampled.size() != k) {
              if (err) *err = "OursAdaptiveBaseline::Sample: residual SampleA failed";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_a[static_cast<usize>(begin + i)];
              const u32 oh = sampled[static_cast<usize>(i)];
              if (q_is_r) out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              else out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
            }
          }
        }

        // Pattern B residual.
        {
          const u32 begin = plan.offset_b[sid_u];
          const u32 end = plan.offset_b[sid_u + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            if (!other.SampleB(q_ylo, q_yhi, k, &rngB, &sampled) || sampled.size() != k) {
              if (err) *err = "OursAdaptiveBaseline::Sample: residual SampleB failed";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_b[static_cast<usize>(begin + i)];
              const u32 oh = sampled[static_cast<usize>(i)];
              if (q_is_r) out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              else out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
            }
          }
        }

        if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
        else as.Insert(handle, q_ylo, q_yhi);
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::OursReportJoinEnumerator2D<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};
  detail::Ours2DContext<Dim, T> ctx_;

  std::vector<join::Side> start_side_;
  std::vector<u32> start_index_;

  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;
  u64 W_{0};
  bool weights_valid_{false};

  std::vector<u8> cached_;
  std::vector<u64> cache_off_;
  std::vector<u64> cache_len_;
  std::vector<u32> cache_partners_;

  // Prefetch sample cache S_i (Framework III): for each START event i,
  // store a prefix of i.i.d. partner samples (with replacement).
  std::vector<u32> prefetch_keep_;
  std::vector<std::vector<u32>> prefetch_partners_;
  bool cache_valid_{false};

  u64 budget_B_{0};
  u64 budget_used_{0};
  u64 w_small_{0};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs

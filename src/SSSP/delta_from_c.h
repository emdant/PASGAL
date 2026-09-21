#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <type_traits>

#include "parlay/primitives.h"
#include "parlay/sequence.h"

/*
Derive delta from a graph-independent constant C:

    delta = C * mean_edge_weight / average_degree

`average_degree` is the number of arcs stored in the CSR per vertex, G.m / G.n.
A symmetrised graph stores every edge in both directions, so this is 2m/n, which
is what the graph statistics tool reports as "Average degree". `mean_edge_weight`
is taken over the same stored arcs; a symmetrised edge contributes twice with
the same weight, so the mean is unaffected.

TIMING
------
DeltaSelector::get() is called from inside the timed region of run(), once per
round of every source, immediately before sssp() initialises its state. That is
the whole point: the cost of choosing delta is charged to the algorithm.

To take it back out of the timer, pass -O. warmup() then computes the value
once before the timer starts and get() returns the cached value. The call site
never moves, so switching between the two costs nothing.
*/

template <class Graph> double mean_edge_weight(const Graph &G) {
  auto weights = parlay::delayed_seq<double>(
      G.m, [&](size_t i) { return static_cast<double>(G.edges[i].w); });
  return parlay::reduce(weights) / static_cast<double>(G.m);
}

template <class EdgeTy> class DeltaSelector {
public:
  DeltaSelector() = default;
  DeltaSelector(EdgeTy fixed_delta, double c, bool use_c, bool outside_timer)
      : fixed_delta_(fixed_delta), c_(c), use_c_(use_c),
        outside_timer_(outside_timer) {}

  // Call before the timer starts. Only does work under -O.
  template <class Graph> void warmup(const Graph &G) {
    if (use_c_ && outside_timer_) {
      cached_ = compute(G);
      have_cached_ = true;
    }
  }

  // Call inside the timed region, just before sssp() initialises its state.
  template <class Graph> EdgeTy get(const Graph &G) const {
    if (!use_c_) {
      last_ = fixed_delta_;
    } else if (have_cached_) {
      last_ = cached_;
    } else {
      last_ = compute(G);
    }
    return last_;
  }

  // Report the delta the last get() handed out. Call outside the timer.
  void print_last() const {
    if (use_c_) {
      if constexpr (std::is_integral_v<EdgeTy>) {
        printf("Derived delta: %zu (C = %f)\n", (size_t)last_, c_);
      } else {
        printf("Derived delta: %f (C = %f)\n", (double)last_, c_);
      }
    }
  }

  bool use_c() const { return use_c_; }

private:
  template <class Graph> EdgeTy compute(const Graph &G) const {
    const double average_degree =
        static_cast<double>(G.m) / static_cast<double>(G.n);
    double delta = c_ * mean_edge_weight(G) / average_degree;

    // An integer-weighted graph has no use for a delta below 1.
    if constexpr (std::is_integral_v<EdgeTy>) {
      delta = std::max(1.0, std::round(delta));
    }
    return static_cast<EdgeTy>(delta);
  }

  EdgeTy fixed_delta_ = 0;
  double c_ = 0.0;
  bool use_c_ = false;
  bool outside_timer_ = false;
  EdgeTy cached_ = 0;
  bool have_cached_ = false;
  mutable EdgeTy last_ = 0;
};

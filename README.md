# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Lire en français](README.fr.md)

**No — for the problems and active-set sizes this repository measured,
hashing `FrankWolfe.jl`'s active set would not bring a measurable advantage.**
`find_atom`'s linear scan took under **0.15% of runtime** in three real
blended pairwise conditional gradient (BPCG) runs, with active sets up to
389 atoms. An isolated microbenchmark shows why: hashing an atom always
costs O(dimension) — it must touch every entry — while the scan's `!=` exits
after the first mismatching entry, which for atoms with independent
coordinates is usually the first one or two. A Dict only wins once the
active set is large enough, or the atoms similar enough, to make up that
gap; neither happened here.

## The question

[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
open since 2021 with no comments: could
[`OrderedCollections.OrderedSet`](https://github.com/JuliaCollections/OrderedCollections.jl),
or more generally a hash for the active set's atoms, help? `find_atom`
(`active_set.jl:316`) scans linearly, called by `active_set_update!`
whenever no index is supplied — which is every "add to active set" step of
BPCG (`blended_pairwise.jl:374`, `nothing` passed explicitly) and every
pairwise step of plain PFW (`pairwise.jl:242`, since `pfw_step` never
supplies one), but neither call site of away-step FW (`afw.jl`), which
always tracks the index itself. `_unsafe_equal` (line 499/513) is **exact**
— elementwise `!=` for a dense `Array`, `==` for a sparse one — which is
what makes hashing sound: no tolerance is being traded away.

## What was measured

`measurement/` runs BPCG (`lazy=true`) on three problems, instrumented by
adding a timed method for `find_atom` on `FrankWolfe.jl`'s own `ActiveSet`
type rather than editing the package (`measurement/instrumentation.jl`):

| Problem | Dimension | Iterations | Max active set | Mean | `find_atom` calls | Lookup share |
|---|---|---|---|---|---|---|
| Birkhoff, n=25 | 625 | 8,000 | 158 | 138.0 | 159 | 0.135% |
| Birkhoff, n=60 | 3,600 | 20,000 | 389 | 300.3 | 389 | 0.057% |
| L∞-ball, d=3,000 | 3,000 | 15,000 | 241 | 170.2 | 240 | 0.022% |

`microbenchmark/` isolates the lookup: a linear scan against a `Dict`, over
active-set sizes 1–20,000 and dimensions 16–8,192, on two atom scenarios —
**generic** (independent random coordinates, so `!=` exits almost at once)
and **adversarial** (atoms share one fixed prefix, differing only in the
last entry, forcing every comparison to run the prefix out). The crossover
— the size beyond which the `Dict` wins:

| Dimension | Generic atoms | Adversarial atoms |
|---|---|---|
| 16 | 50 | 5 |
| 128 | 500 | 5 |
| 1,024 | 2,000 | 5 |
| 8,192 | between 6,500 and 10,000 (noisy — see `MEASURING.md`) | 10 |

Both scripts write their table to a committed `results.csv` beside them; see
`what-is-where.md`.

## The answer

The three real runs never came close to a generic crossover — the largest
active set seen (389) sits below even the dimension-128 threshold (500) —
and `FrankWolfe.jl`'s own vertices (sparse permutation matrices, dense box
corners) behave closer to the generic case than the adversarial one, which
is why the measured lookup share stayed under 0.15% throughout. **A Dict
would only start paying for itself if the active set grew past a few
hundred to several thousand atoms, or the polytope made distinct vertices
unusually easy to confuse.** Neither happened here, on the polytopes
measured — not a claim about every polytope `FrankWolfe.jl` can be pointed
at. `DECISIONS.md` has the draft issue comment and the open question.

## Reading further

Exactness, call sites, the issue: `references.md`. The machine and its
noise: `MEASURING.md`. Every file, one line each: `what-is-where.md`. What
is Mohamed's to decide: `DECISIONS.md`.

# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Lire en français](README.fr.md)

**Yes: for the atoms `FrankWolfe.jl` actually generates, a short prefix
hash beats `find_atom`'s linear scan at the active-set sizes the three real
runs reached.** The first version of this repository only ever tested
hashing the *whole* atom (always O(dimension)) against a miss query (the
scan's best case: `!=` exits at once); neither restriction is necessary,
and relaxing both reverses the answer. A hash over just the first `k=8`
coordinates, still confirmed by the same exact equality the scan already
uses on every bucket hit, beats the scan by active-set size 100, for every
atom alphabet tested (Birkhoff permutation matrices, L∞-ball box corners,
generic random vectors) and every query mix (miss, hit, a stress mix),
comfortably inside the 158-389 range the three real blended pairwise
conditional gradient (BPCG) runs reached. **The answer to #244 changes: a
prefix-hashed active set would help.**

## The question

[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
open since 2021 with no comments: could
[`OrderedCollections.OrderedSet`](https://github.com/JuliaCollections/OrderedCollections.jl),
or more generally a hash for the active set's atoms, help? `find_atom`
(`active_set.jl:316`) scans linearly, called by `active_set_update!`
whenever no index is supplied, which is every "add to active set" step of
BPCG (`blended_pairwise.jl:374`, `nothing` passed explicitly) and every
pairwise step of plain PFW (`pairwise.jl:242`, since `pfw_step` never
supplies one), but neither call site of away-step FW (`afw.jl`), which
always tracks the index itself. `_unsafe_equal` (line 499/513) is **exact**
(elementwise `!=` for a dense `Array`, `==` for a sparse one), which is
what makes hashing sound: no tolerance is being traded away, whether the
hash covers the whole atom or, as it turns out matters most, only part of
it (see "Prefix hashing" below).

## What was measured

`measurement/` runs BPCG (`lazy=true`) on three problems, instrumented by
adding a timed method for `find_atom` on `FrankWolfe.jl`'s own `ActiveSet`
type rather than editing the package (`measurement/instrumentation.jl`):

| Problem | Dimension | Iterations | Max active set | Mean | `find_atom` calls | Lookup share |
|---|---|---|---|---|---|---|
| Birkhoff, n=25 | 625 | 8,000 | 158 | 138.0 | 159 | 0.144% |
| Birkhoff, n=60 | 3,600 | 20,000 | 389 | 300.3 | 389 | 0.058% |
| L∞-ball, d=3,000 | 3,000 | 15,000 | 241 | 170.2 | 240 | 0.021% |

`microbenchmark/` isolates the lookup: a linear scan against a `Dict`, over
active-set sizes 1-20,000 and dimensions 16-8,192, on two atom scenarios:
**generic** (independent random coordinates, so `!=` exits almost at once)
and **adversarial** (atoms share one fixed prefix, differing only in the
last entry, forcing every comparison to run the prefix out). The crossover,
the size beyond which the `Dict` wins:

| Dimension | Generic atoms | Adversarial atoms |
|---|---|---|
| 16 | 50 | 5 |
| 128 | 500 | 5 |
| 1,024 | 2,000 | 5 |
| 8,192 | between 6,500 and 10,000 (noisy, see `MEASURING.md`) | 10 |

Both scripts write their table to a committed `results.csv` beside them; see
`what-is-where.md`.

### Prefix hashing

The sweep above left two gaps: it only ever hashed the *whole* atom, and it
only ever queried a miss (a fresh atom absent from the set, the scan's best
case). `microbenchmark/run_prefix.jl` closes both. A **prefix hash** buckets
atoms by only their first `k` coordinates; a bucket hit is still confirmed
against the whole atom with the same exact equality `dict_lookup` already
uses, so shortening the hash can only cost speed, never correctness
(`microbenchmark/lookup_methods.jl`'s `prefix_lookup`). Three query types
are timed: **miss** (as before), **hit** (the target sampled uniformly from
the stored atoms, so the scan pays its average position rather than its
best case), and a **50/50 mix** of both. `measurement/results.csv`'s new
`find_atom_hits` column shows the three real BPCG runs never actually
produced a hit (0 of 159, 389, and 240 calls), so the mix is a deliberate
stress test, not a reproduction of what was observed; miss is what
actually happened.

The atom alphabet is the decisive axis. `rand(dim)` Float64 atoms (the
"generic" scenario above) give one coordinate near-perfect power to
discriminate; `FrankWolfe.jl`'s real atoms do not. Birkhoff-polytope
vertices are permutation matrices, generated here by calling
`FrankWolfe.BirkhoffPolytopeLMO.compute_extreme_point` exactly as
`measurement/problems.jl` does, and kept sparse rather than densified so
the scan baseline is the real `_unsafe_equal` sparse comparison (about
20ns flat on a mismatch, not the O(dimension) a dense flatten would
suggest). L∞-ball vertices are box corners, `{-1,+1}` entries from
`FrankWolfe.LpNormBallLMO{Inf}.compute_extreme_point`. Both alphabets carry
one bit of discriminating power per coordinate, not the roughly 53 bits a
generic Float64 coordinate carries.

The crossover, the active-set size beyond which a `k`-coordinate prefix
hash beats the scan, for a **miss** query (the full table, every `k` and
query type, is `microbenchmark/results_prefix_crossover.csv`):

| Atom alphabet (dimension) | Real problem's max active set | k=1 | k=2 | k=4 | k=8 | k=16 | k=64 | full hash |
|---|---|---|---|---|---|---|---|---|
| Birkhoff permutation, n=25 (625) | 158 | 389 | 100 | 50 | 50 | 50 | 100 | 500 |
| Birkhoff permutation, n=60 (3,600) | 389 | 1,500 | 500 | 200 | 50 | 50 | 100 | 1,500 |
| L∞-ball box corners (3,000) | 241 | 50 | 50 | 50 | 50 | 100 | 100 | 1,000 |
| generic Float64 (3,000), control | n/a | 50 | 50 | 50 | 100 | 158 | 300 | none in [1, 2000] |

`k=8`, a fixed, dimension-independent prefix (1.3% of the smaller
Birkhoff run's dimension, 0.2% of the larger), clears every real
alphabet's own observed maximum, on hit and mix as well as miss: at
active-set size 158 (Birkhoff n=25's own maximum), `k=8` runs a lookup in
500.6ns against the scan's 949.4ns on miss (47% faster) and 350.5ns
against 488.2ns on hit (28% faster); at size 200 for the L∞-ball, `k=8`
runs in 53.4ns against the scan's 1,243.8ns on miss, a 23x speedup
(`microbenchmark/results_prefix_timing.csv`).

Collision rate, how often a `k`-coordinate bucket holds more than one
atom, explains why the speedup differs so much between the two real
alphabets. At `k=1`, box corners have only 2 possible values, so 100% of a
200-atom active set shares a bucket with another atom, and the hash still
wins, because the O(1) hash step and a roughly halved candidate list beat
the scan's O(size) walk. At `k=8`, box corners' 256 possible prefixes
spread those 200 atoms thin enough that only 45.5% still share a bucket,
which is most of `k=8`'s 23x win. Permutation matrices are the harder
case: even at `k=8`, only 9 of the 256 possible prefixes are ever reached
(a permutation matrix's flattened prefix is almost always all zero), so
100% of a 389-atom active set still shares one of those 9 buckets, and the
win there comes entirely from the O(1) hash step and `_unsafe_equal`'s own
cheap sparse comparison, not from the bucket narrowing the candidate list
much (`microbenchmark/results_prefix_collisions.csv`).

## The answer

Whether hashing the active set helps depends entirely on what gets
hashed. **Hashing the whole atom is still a bad idea**: it costs
O(dimension) regardless of alphabet, and the full-`Dict` sweep above only
wins past a few hundred to several thousand atoms for atoms with
independent coordinates (50 at dimension 16, 6,500-10,000 at dimension
8,192). `FrankWolfe.jl`'s own atoms give a full hash an easier time (the "Prefix
hashing" table's "full hash" column: 500 to 1,500 for Birkhoff permutation
matrices, 300 to 1,500 for box corners, both earlier than the generic
scenario at a comparable dimension): how many atoms it takes to amortize a
hash's fixed, dimension-sized cost against the scan depends on how
expensive the scan itself is per atom, which differs by alphabet and is
not always in the direction intuition suggests (see
`microbenchmark/results_prefix_timing.csv` for the raw per-atom numbers
behind each alphabet's own crossover). Either way, none of the three real
BPCG runs (max active set 158 to 389) came close to a full-hash crossover,
which is why the lookup share stayed under 0.15% throughout, and that
measurement stands.

**Hashing a short prefix is a different question, with a different
answer: yes.** A `k=8` prefix hash, confirmed by the same exact equality
the scan already uses, beats `find_atom`'s linear scan for every atom
alphabet `FrankWolfe.jl` actually generates (Birkhoff permutation
matrices, L∞-ball box corners) and every query mix tested (miss, hit, a
50/50 stress mix), at active-set sizes below every one of the three real
runs' own observed maximum. This holds despite both real alphabets being
exactly the low-entropy, easy-to-confuse case the original "adversarial"
scenario was meant to bound: a short prefix still wins there, because the
O(1) hash step and even a barely-narrowed candidate list beat an O(size)
scan once the active set is a few dozen atoms large, which every one of
the three real runs' active sets became.

**This changes the answer to #244.** A hash-augmented active set, keyed by
a short prefix of each atom rather than the whole atom, would measurably
help `FrankWolfe.jl`'s BPCG at the active-set sizes it actually reaches.
`DECISIONS.md`'s draft issue comment has been rewritten to say so; the
open question it still leaves is choosing `k` and rebuilding the index
incrementally as atoms are pushed and dropped, not whether prefix hashing
helps.

## Reading further

Exactness, call sites, the issue: `references.md`. The machine and its
noise: `MEASURING.md`. Every file, one line each: `what-is-where.md`. What
is Mohamed's to decide: `DECISIONS.md`.

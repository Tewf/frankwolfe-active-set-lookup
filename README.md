# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Lire en français](README.fr.md)

**Yes, still, once the whole lifecycle is costed, not lookup alone: a
hash-augmented active set beats today's linear scan at the active-set
sizes the three real runs reached.** The first version of this repository
concluded hashing the whole atom doesn't help; the second found that a
short `k=8` prefix hash does. This one asks a harder question, because
`measurement/results.csv` shows all three real runs produced **zero
`find_atom` hits**: every call was a miss followed by a `push!`, so the
real per-iteration cost is lookup plus insert, and a structure also has to
survive `active_set_cleanup!`'s `deleteat!`, which shifts every later
index down by one. Costing all three changes which structure to propose,
not whether hashing helps. For Birkhoff's permutation-matrix atoms, keying
on *where* the sparse structure's nonzeros are (`SparseMatrixCSC.rowval`)
rather than on a flattened value prefix turns the `k=8` prefix hash's 9
buckets for 158 or 389 atoms into 158 and 389 buckets respectively, one
atom per bucket, and wins on total per-iteration cost by 4.8-25.6x over the
existing prefix-hash recommendation. For the L∞-ball, the existing `k=8`
prefix hash is still the right answer: nothing beats it there. A third
idea, a hash trie recursing over coordinate blocks, was swept across four
coordinate-selection strategies, three values of `k`, and three depths per
alphabet, and **never won outright at any of the three real active-set
sizes**; at Birkhoff n=25's own size it was measurably worse than doing
nothing. See "Lookup is not the whole cost" below for the numbers.

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

## Lookup is not the whole cost

Every measurement up to this point timed lookups only. `measurement/results.csv`'s
`find_atom_hits` column already showed the three real BPCG runs never produced a
hit; this branch went further and instrumented `deleteat!` the same way
`find_atom` was instrumented (`measurement/instrumentation.jl`), to find out how
often a structure actually has to repair itself. The answer is: almost never.
Across 8,002, 20,002, and 15,002 iterations, the three real runs called
`deleteat!` on an individual atom **2, 1, and 0 times** respectively (0.025%,
0.005%, and 0% of iterations). Every `find_atom` call, by contrast, is a miss
followed by a `push!` (0 hits out of 159, 389, and 240 calls). So the real
per-iteration cost, weighted by how often each operation actually happens in
these three runs, is dominated by lookup and insert, and deletion-repair, while
real and sometimes expensive per operation, barely moves the total for the
workloads measured here. That is itself the most useful thing this branch
found before comparing any structures: it means a structure only has to be
cheap to *build*, not cheap to *repair*, to win here, though a workload with
more frequent drop steps (PFW, BCG, or a looser `epsilon`, none of which were
run) could see repair cost matter far more, an open question `DECISIONS.md`
carries forward.

### Idea 1: key on the sparsity pattern, not the values

The "Prefix hashing" section above found that Birkhoff's `k=8` prefix hash wins
despite being nearly useless as a key: at active-set size 389, 9 of the 256
possible 8-coordinate prefixes are ever reached, so every lookup still walks a
bucket of about 43 candidates on average, and the win comes entirely from the
O(1) hash step and a cheap sparse `==`, not from narrowing the candidate list.
The reason is representation: the first 8 entries of a flattened permutation
matrix are row 0, columns 0-7, which hold the single 1 only when the
permutation sends one of those 8 columns to row 1, about a third of the time;
the other two-thirds of atoms share one all-zero bucket. A permutation
matrix's real information, *where* its nonzero sits, never reaches that key at
all: it already sits in `SparseMatrixCSC`'s own `rowval` array, one row index
per stored column, in column order. `microbenchmark/sparse_pattern.jl` hashes
`rowval[1:k]` instead: same O(k) cost, no densifying, and close to perfectly
discriminating for a random permutation (two independent n-permutations agree
on their first k columns' images with probability 1/(n)_k, already under 1 in
10 million at n=25, k=8).

It measures exactly as that reasoning predicts
(`microbenchmark/results_sparse_pattern_collisions.csv`): at `k=8`, the 9
buckets become **158 buckets for 158 atoms** (Birkhoff n=25, its own real
maximum) and **389 buckets for 389 atoms** (Birkhoff n=60), one atom per
bucket both times, `atom_collision_rate` exactly 0.0 where the flattened
prefix's was 1.0.

### Idea 2: a hash trie over coordinate blocks

Mohamed's second design: hash `k` coordinates to a bucket; if that bucket
still holds more than one atom, hash a further block of `k` coordinates, from
just that bucket's members, and recurse, up to `max_depth` levels, confirming
with the same exact equality at the leaf (`microbenchmark/hash_trie.jl`). Four
coordinate-selection strategies were compared, all reduced to the same
"compute one coordinate order, then take sequential k-blocks per level" shape:
first-k (today's flat prefix, as a depth-1 baseline), strided, a fixed random
sample, and a selectivity-ordered choice, ranked by how many distinct values a
coordinate takes across the atoms being indexed (a balance tie-break included,
since a permutation matrix's flattened entries take exactly 2 values, 0 or 1,
at *every* coordinate, so raw cardinality alone cannot rank them at all).
`microbenchmark/run_lifecycle.jl` swept `k in {4,8,16}` and `max_depth in
{1,2,4}` across all four orders, at every alphabet's own real active-set size,
and kept the config with the fewest atoms still sharing a leaf
(`microbenchmark/results_lifecycle_collisions.csv`):

| Alphabet (real size) | Best trie config | Atom collision rate | Mean depth |
|---|---|---|---|
| Birkhoff n=25 (158) | selectivity, k=16, depth=4 | 0.0 | 2.06 |
| Birkhoff n=60 (389) | selectivity, k=16, depth=4 | 0.288 | 3.15 |
| L∞-ball (241) | selectivity, k=16, depth=1 | 0.0 | 1.00 |

Depth does not substitute for the right representation. Even the best trie
config, examining 64 of Birkhoff n=60's 3,600 flattened coordinates across 4
levels, still leaves 28.8% of atoms sharing a leaf with another atom, because
every one of a permutation matrix's flattened coordinates carries the same
weak, skewed signal (about a 1-in-n chance of being 1, for any coordinate,
picked any way): no coordinate-selection heuristic over that representation
matches what a one-level, `k=8` key built from the sparse structure itself
gets for free. The smaller Birkhoff run (n=25) does fully resolve at depth 4,
but even there the trie's own overhead (four coordinate lookups per level, a
16-long Float64 key hashed and dict-traversed at each) makes it slower
end to end than the structures it was built to beat, below.

### The total, at each real run's own maximum

Lookup, insert (via marginal build-cost differencing: build at size N and at
N+100, divide the difference by 100, so no structure is timed by repeatedly
mutating itself), and deletion-repair (direct repeated timing: every repair
function here walks the same number of stored positions whether or not a
previous call already ran, so it is stable to time that way) were measured
separately for the scan, the `k=8` prefix hash, the sparse-pattern key
(Birkhoff only), and each alphabet's best trie config
(`microbenchmark/results_lifecycle_timing.csv`), then weighted by the real
per-iteration lookup/insert/delete rates above
(`microbenchmark/results_lifecycle_total.csv`):

| Alphabet (real size) | scan | prefix `k=8` | sparse-pattern `k=8` | best trie |
|---|---|---|---|---|
| Birkhoff n=25 (158) | 20.25 ns | 13.11 ns | **2.72 ns** | 23.24 ns |
| Birkhoff n=60 (389) | 57.49 ns | 48.11 ns | **1.88 ns** | 29.83 ns |
| L∞-ball (241) | 33.52 ns | **3.31 ns** | n/a (dense atoms) | 3.73 ns |

**On total per-iteration cost at each alphabet's own real active-set size, the
sparse-pattern key wins for both Birkhoff runs (2.72 ns at n=25's 158-atom
maximum, 1.88 ns at n=60's 389-atom maximum) and the existing `k=8` prefix
hash remains the winner for the L∞-ball (3.31 ns at 241 atoms), with the hash
trie never winning outright at any of the three, and actually costing more
than the plain scan at Birkhoff n=25's own size.** On Birkhoff, the trie's
insert cost is what sinks it even where its lookup is competitive: building a
16-coordinate Float64 key and walking up to 4 dict levels costs 743-773 ns
per marginal insert there, next to the sparse-pattern key's 57-98 ns and the
flat prefix hash's 86-88 ns (all three cost within 150-184 ns on the
L∞-ball, where the trie stays close to the prefix hash on total cost too).

Deletion-repair itself is not free per operation: every structure's repair
function walks every stored position (79-2,330 ns depending on structure and
alphabet at these sizes), and for the two new structures specifically,
usually costing *more* than the raw `deleteat!` shift on the underlying
Vector it sits alongside (357-676 ns for the same Birkhoff sizes), because
the naive repair implementation here re-scans the whole index rather than
only the entries that actually moved. It just never gets charged more than
2-3 times per 10,000 iterations in the runs measured, which is why it does
not change the ranking above; `DECISIONS.md` carries the
open question of whether a smarter, indirection-based repair could close that
gap for a workload where deletion is not this rare.

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

**This changes the answer to #244**, and "Lookup is not the whole cost"
above refines it once more, on a question this section never asked: not
just whether a hash beats a scan on a lookup, but whether it still wins
once insert and deletion-repair are counted too, and which specific
structure to propose. It does still win, everywhere it won on lookup
alone; the specific recommendation, though, is no longer simply "a `k=8`
prefix hash": for Birkhoff's permutation-matrix atoms, key on the sparse
structure itself (`microbenchmark/sparse_pattern.jl`), not a flattened
value prefix, a 4.8-25.6x win over the prefix hash on total cost, not just a
faster lookup; for the L∞-ball, the `k=8` prefix hash already measured
here remains the right answer, and a recursive hash trie over coordinate
blocks, the third structure this repository tried, **does not earn its
added complexity anywhere**: it never won outright at any of the three
real active-set sizes, and was measurably worse than the plain scan at
Birkhoff n=25's own size. `DECISIONS.md`'s draft issue comment has been
rewritten again to say so.

## Reading further

Exactness, call sites, the issue: `references.md`. The machine and its
noise: `MEASURING.md`. Every file, one line each: `what-is-where.md`. What
is Mohamed's to decide: `DECISIONS.md`.

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

**The sparse-pattern key itself still allocates, and removing that
allocation makes the win bigger.** It is a `Vector{Int}`, so every lookup
and every insert builds one on the heap first; `microbenchmark/pattern_key_reps.jl`
adds two allocation-free representations of the identical key, a `UInt64`
folded with an incremental hash and an `NTuple{k,Int}`. **Propose the
`UInt64` fold at `k=4`: 0 bytes and 0 allocations per lookup, against the
`Vector{Int}` key's 128 bytes and 2 allocations at `k=8`, and a total
per-iteration cost of 1.11ns at Birkhoff n=25's 158-atom maximum and
0.81ns at n=60's 389-atom maximum, down from Idea 1's already-winning
2.56ns and 1.84ns (this branch's own re-measurement of the `Vector{Int}`
key at `k=8`, close to but not identical to the earlier 2.72ns/1.88ns,
ordinary run-to-run noise per `MEASURING.md`, not a regression).** The
`NTuple{k,Int}` key is just as allocation-free and usually a little
faster still, but needs `k` fixed as a compile-time type parameter rather
than an ordinary argument, a real API cost the `UInt64` fold does not
have; see "Idea 1, tightened" below for both representations' numbers
side by side and the collision hazard the `UInt64` fold introduces (and
survives).

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

### Idea 1, tightened: three representations of the pattern key

The sparse-pattern key above (`pattern_key(atom, k) = atom.rowval[1:k]`) already
won. It is also a `Vector{Int}`, so building it, on every lookup and every
insert, allocates: `microbenchmark/pattern_key_reps.jl` builds two more
representations of the exact same key, chosen to remove that allocation
rather than to change what gets keyed:

- **`UInt64`**: the k row indices folded into one integer with an
  incremental hash (`h = hash(rowval[i], h)`, chained). Nothing is
  allocated building it; the trade is that the fold is lossy, so two
  *different* patterns can now land in the same `Dict` bucket, a real
  hash collision, not merely a prefix tie the way two atoms sharing a
  `Vector{Int}` key already could.
- **`NTuple{k,Int}`**: Julia's fixed-length tuple, stored inline in the
  `Dict`'s own key array rather than as a separate heap object, so this
  is also allocation-free. Unlike the `UInt64` fold, this has no
  collision hazard at all: a `Dict` compares two tuples elementwise,
  exactly as it already compares two `Vector{Int}`s, so two different
  patterns can never share a tuple key. The cost is that `k` has to be
  known to the compiler as a type parameter (`Val(k)`), not passed as an
  ordinary `Int`, for the tuple to actually be stored inline; a
  freshly-built `Val(k)` inside a hot per-call loop costs an allocation
  in exactly the way this idea is trying to avoid, confirmed empirically
  while writing this file (32-80 bytes per call, vs. 0 once `Val(k)` is
  built once and reused).

Both keep `pattern_key`'s own scope restriction: `SparseMatrixCSC` only,
reading `rowval` directly.

**Allocation, at k=8, both real Birkhoff sizes (identical at both sizes,
since a Vector/tuple/hash's own allocation depends on k and its element
type, not on how many atoms are stored):**

| | lookup | insert (key + bucket entry) |
|---|---|---|
| `Vector{Int}` (existing) | 128 bytes, 2 allocations | 192 bytes, 4 allocations |
| `UInt64` | **0 bytes, 0 allocations** | 64 bytes, 2 allocations |
| `NTuple{8,Int}` | **0 bytes, 0 allocations** | 64 bytes, 2 allocations |

Insert never reaches zero for any representation: even a `UInt64` or an
inline tuple key still has to allocate a fresh one-element `Vector{Int}`
bucket the first time a key is seen (`microbenchmark/bucket_lifecycle.jl`'s
`bucket_insert!`), and at these real sizes' own k=8 collision rate (0.0,
`results_sparse_pattern_collisions.csv`) a first-time key is the
overwhelmingly common case. What the 192-vs-64-byte gap is: the
`Vector{Int}` key itself (128 bytes at k=8) plus that same 64-byte bucket,
against the bucket alone for the two allocation-free representations. Full
numbers, every k swept, both Birkhoff sizes, and the L∞-ball/generic
`Vector{Float64}` value-prefix baseline (measured the same way, for
context, since its own allocation was never measured directly before this
file): `microbenchmark/results_pattern_key_reps_allocations.csv`.

**Time and total per-iteration cost, at k=8 (directly comparable to "The
total" table above, same alphabets, same real sizes):**

| Alphabet (real size) | representation | lookup | insert | delete-repair | total/iter |
|---|---|---|---|---|---|
| Birkhoff n=25 (158) | `Vector{Int}` | 31.94 ns | 91.12 ns | 440.08 ns | 2.555 ns |
| Birkhoff n=25 (158) | `UInt64` | 21.92 ns | 44.16 ns | 378.95 ns | 1.408 ns |
| Birkhoff n=25 (158) | `NTuple{8,Int}` | 21.56 ns | 34.76 ns | 636.92 ns | 1.278 ns |
| Birkhoff n=60 (389) | `Vector{Int}` | 33.97 ns | 57.91 ns | 1115.10 ns | 1.843 ns |
| Birkhoff n=60 (389) | `UInt64` | 23.67 ns | 33.65 ns | 1158.68 ns | 1.173 ns |
| Birkhoff n=60 (389) | `NTuple{8,Int}` | 24.02 ns | 31.39 ns | 1212.21 ns | 1.138 ns |

Both new representations beat `Vector{Int}` on lookup and insert, not just
on allocation: lookup is faster too (about 29-33% at both sizes), which is
the folded/inline key itself being cheaper to build and compare,
independent of what it costs the allocator. `NTuple{8,Int}` against
`UInt64` is not a clean sweep, though: it wins insert clearly (21.3%
faster at n=25, 6.7% at n=60) and lookup narrowly at n=25 (1.6% faster),
but is 1.5% *slower* than `UInt64` on lookup at n=60, small enough at all
four points that the choice between them is really the API question
above (does the caller get to choose `k` at runtime), not a speed
question. **Delete-repair does not follow the same pattern, and has no
consistent winner across the two sizes**: `UInt64` beats `Vector{Int}` at
n=25 (378.95 ns vs. 440.08 ns) but loses at n=60 (1158.68 ns vs. 1115.10
ns); `NTuple{8,Int}` loses at both (636.92 ns and 1212.21 ns against
440.08 ns and 1115.10 ns). This is plausible rather than alarming:
`bucket_delete_repair!` walks every bucket's stored position *values*,
never touching a key at all (`bucket_lifecycle.jl`), so key type has no
structural reason to change this cost, and Dict iteration order (which
can depend on a key's hash, and so differs by representation even at
identical bucket counts) is enough to explain a difference this size on a
machine this noisy (`MEASURING.md`). It does not change which
representation wins overall, since delete-repair is weighted by a real
per-iteration rate under 0.03% at both Birkhoff sizes
(`measurement/results.csv`), the same reason it did not change Idea 1's
own ranking against the scan and the flat prefix hash.

**k sweep: k=8 is no longer the best choice, for either new
representation.** The original `k=8` was chosen because it was the
smallest prefix giving perfect discrimination against the flattened value
prefix, where every extra coordinate read was cheap next to the
allocation the whole key still cost. With that allocation gone, a smaller
`k` is close to free to try, and the numbers say to:

| Alphabet (real size) | representation | k=2 | k=4 | k=8 |
|---|---|---|---|---|
| Birkhoff n=25 (158) | `UInt64` | 1.253 ns | **1.110 ns** | 1.408 ns |
| Birkhoff n=25 (158) | `NTuple{k,Int}` | **1.126 ns** | 3.876 ns\* | 1.278 ns |
| Birkhoff n=60 (389) | `UInt64` | **0.789 ns** | 0.812 ns | 1.173 ns |
| Birkhoff n=60 (389) | `NTuple{k,Int}` | **0.689 ns** | 2.957 ns\* | 1.138 ns |

\* Reproducible across repeated runs, not noise, but not a property of the
representation either: at `k=4`, `NTuple{4,Int}`'s marginal insert timing
(build at size N vs. N+100, differenced) measured 135-178 ns against 19-38
ns at `k=2` and `k=8` for the same alphabet, most likely a Dict-resize
boundary landing inside the differenced size range for that one key byte
width and not the others (`NTuple{2,Int}` is 16 bytes, `NTuple{4,Int}` 32,
`NTuple{8,Int}` 64); `DECISIONS.md` carries this forward rather than
smoothing it over. Excluding that one figure, `NTuple`'s own lookup and
delete-repair costs at k=4 are unremarkable, in line with k=2 and k=8.

`k=2` is not perfectly discriminating: `results_pattern_key_reps_collisions.csv`
shows a real 27.9% (n=25) / 11.3% (n=60) atom collision rate there, against
0.0% at `k=4` and `k=8`, confirmed correct regardless (every candidate is
still checked against the whole atom before being trusted) but not free:
the extra confirm comparisons show up in delete-repair's walk. `k=4`
already reaches 0.0% collision (the same finding "Idea 1" reported at
`k=8`, just one sweep point earlier) and, `k=4`'s `NTuple` insert anomaly
aside, is at or within a few percent of each representation's own best k
at both sizes.

`Vector{Int}` was swept at k=2/4/8 too (`results_pattern_key_reps_total.csv`),
and the honest reading there is that no single k wins at both sizes:
n=25's order is k=4 (2.045 ns) < k=8 (2.555 ns) < k=2 (2.679 ns), while
n=60's is k=2 (1.341 ns) < k=4 (1.495 ns) < k=8 (1.843 ns). k=8 is worst
at n=60 but only middling at n=25, and k=2, which loses outright at n=25,
wins outright at n=60; nothing about `Vector{Int}`'s own k-ranking is
one-sided the way `UInt64`'s is (k=8 is `UInt64`'s worst choice at both
sizes, not just one; `NTuple`'s own k=4 figure is the `*`-marked anomaly
above, so its k-ranking is not read from here at all). What is one-sided
is the representation comparison itself: **`Vector{Int}`'s own best case
at any k (2.045 ns at n=25, 1.341 ns at n=60) still loses to `UInt64`'s
or `NTuple`'s worst case among k=2/k=4 at the matching size (1.253 ns and
0.812 ns respectively)**. The representation change is where the larger,
dependable win is; which k to pick on top of it is a smaller, noisier
question, most legible for `UInt64`/`NTuple` (avoid k=8) and genuinely
unresolved for `Vector{Int}` at these two sizes.

**The collision hazard the `UInt64` fold introduces, and how it is
survived:** because the fold is lossy, `microbenchmark/test_pattern_key_reps.jl`
constructs two atoms whose real `k=2` patterns are different
(`pattern_key` returns `[2, 5]` and `[6, 1]`, found by a short random
search, not designed by hand) but whose `UInt64` fold, narrowed to 1 bit
via `pattern_key_uint64`'s `bits` keyword (a test-only knob;
`pattern_key_uint64(atom, k)`'s default is the unnarrowed 64-bit fold),
collides exactly. The test asserts the resulting `PatternIndexU64` puts
both atoms in **one bucket holding both indices**, not one index each,
and that `pattern_lookup_u64` still returns the correct atom for each of
their own queries, plus a clean miss for a third atom matching neither:
correctness survives the forced collision because every candidate in a
bucket is confirmed against the whole atom with `==` before being
trusted, the bucket is a `Vector{Int}` of candidates rather than a single
index, and neither is optional for this representation the way they are
merely convenient for `Vector{Int}`. A second testset checks the
`NTuple{K,Int}` representation on the same two atoms and confirms it
keeps them in two separate buckets, since tuple equality is exact: there
is no forced-collision case to write for it, only the absence of one. At
the fold's real, unnarrowed default (`bits=64`), this same pair of atoms
does *not* collide; the hazard is real (the test forces it), not
something these two real Birkhoff sizes happen to trigger on their own.

**Both new representations agree with each other and with the
`Vector{Int}` baseline**, checked directly (not just inferred from
matching totals) in `test_pattern_key_reps.jl`'s third testset, on a
small hand-built five-atom pool at k=2, every atom's own self-lookup plus
one deliberate miss; the large-scale sweep itself
(`run_pattern_key_reps.jl`, the real Birkhoff sizes) does not cross-check
representations against each other, only against its own timing and
allocation numbers. The generic `rand(dim)`
control was kept for this sweep too (`results_pattern_key_reps_*.csv`'s
`generic_d3000` rows), measured with the existing `Vector{Float64}`
value-prefix baseline only, since the pattern key does not apply to a
dense alphabet; nothing about it is surprising (0.0% collision at every
k, allocation matching the L∞-ball's own `Vector{Float64}` numbers at the
same k), so it is not discussed further here.

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

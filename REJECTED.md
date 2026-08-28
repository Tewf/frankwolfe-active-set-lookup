# What was tried and refused, with the numbers that killed it

Not an appendix. If a future maintainer, or the author six months from now,
proposes any of these again, the numbers below are why they were not
shipped. Every claim here has a results CSV behind it; the file names are
given so a claim can be re-checked rather than taken on faith. Longer than the
80-line ceiling `CONTRIBUTING.md` sets, because each of four rejected methods needs its own number, not
just its name, to actually be useful to the next person who thinks of it.

The four methods below were refused against the folded sparse-pattern
key. That key is itself now second to the absence certificate wherever
the caller is a Frank-Wolfe step (`METHOD.md`), and stays the answer for a
caller that has only the atom; sections 5 and 6 are what was weighed
against the certificate.

## The first answer this repository gave was wrong, and here is why

Before any of the four methods below, this repository's very first draft
concluded hashing the active set does not help at all. That answer was
wrong, for three compounding reasons, all in the version of
`microbenchmark/run.jl`'s original sweep: it compared a full-atom hash's
worst case (`Dict` still has to hash every coordinate) against the linear
scan's *best* case (a fresh random query almost always fails on the first
coordinate, so `!=` exits immediately); it only ever measured a miss,
never a hit or a realistic mix; and it used dense, independent-coordinate
random vectors rather than `FrankWolfe.jl`'s own atom shapes, which are
far more confusable than a generic Float64 vector. Fixing any one of the
three would not have been enough; fixing all three is what turned "no"
into "yes, but not by hashing the whole atom" and eventually into the
method `README.md` and `METHOD.md` describe.

## 1. Hashing the whole atom

The first thing measured, and the idea the wrong answer above was drawn
from. A `Dict` keyed on the entire atom is sound (it is exactly as exact
as `_unsafe_equal`), but it pays `O(dimension)` to hash every query and
every stored atom, so it only wins once the scan's own comparable cost is
actually being paid. Against independent random coordinates, the
crossover ranges from 50 atoms at dimension 16 to somewhere between 6,500
and 10,000 at dimension 8,192 (`microbenchmark/results.csv`). Against
`FrankWolfe.jl`'s own atoms it is somewhat more competitive, since a full
hash still amortizes eventually, but the crossover is still 300 to 1,500
atoms depending on alphabet and k (`microbenchmark/results_prefix_crossover.csv`'s
"full hash" column). None of the three real BPCG runs measured here came
anywhere close: their active sets topped out at 158 to 389 atoms, an
order of magnitude below the closest crossover. Rejected: correct, but
never actually reached by a real run.

## 2. The value-prefix key applied to sparse atoms

`microbenchmark/lookup_methods.jl`'s `sparse_prefix` hashes a Birkhoff
atom's first k *flattened* coordinate values, the same recipe that works
for a dense atom. It still beats the scan (the O(1) hash step and even a
weak bucket narrowing win past a few dozen atoms), but it is a poor key
for this alphabet specifically: a flattened permutation matrix's first 8
entries are row 0, columns 0-7, which hold the single 1 only when the
permutation happens to send one of those columns to row 1, about a third
of the time. The other two-thirds of atoms share one all-zero bucket. At
k=8 and a 389-atom active set (Birkhoff n=60), this collapses into just 9
buckets, mean bucket size 43.2, `atom_collision_rate` 1.0
(`microbenchmark/results_prefix_collisions.csv`). Every lookup still walks
roughly 43 candidates on average; the win over the scan comes entirely
from the O(1) hash step and a cheap sparse `==`, not from the key
actually narrowing anything. Costed through the full lifecycle
(`microbenchmark/results_lifecycle_total.csv`), this key's total
per-iteration cost is 13.11ns at n=25 and 48.11ns at n=60, against the
sparse-pattern key's 2.72ns and 1.88ns, a 4.8-25.6x loss. Rejected for
sparse, permutation-like atoms specifically: it is still the right answer
for the L-infinity ball, where nothing beat it (3.31ns at 241 atoms,
same table).

## 3. A hash trie over coordinate blocks

Hash k coordinates to a bucket; if a bucket still holds more than one
atom, hash a further block of k coordinates drawn only from that bucket's
members, and recurse, up to `max_depth` levels
(`microbenchmark/hash_trie.jl`). Four coordinate-selection strategies were
compared (first-k, strided, a fixed random sample, and a
selectivity-ordered choice modeled on composite database indexes), swept
across `k in {4,8,16}` and `max_depth in {1,2,4}`, at each alphabet's own
real active-set size (`microbenchmark/results_lifecycle_collisions.csv`).
It looked right going in: depth lets a structure spend more work only
where atoms are genuinely confusable, instead of committing to one fixed
k for everyone. It never won outright at any of the three real sizes.
Even its best swept config for Birkhoff n=25 (selectivity order, k=16,
depth=4, the one config that fully resolves every collision there) was
measurably worse on total per-iteration cost than doing nothing: 23.24ns
against the plain scan's 20.25ns (`microbenchmark/results_lifecycle_total.csv`).
The reason is structural, not a tuning miss: every one of a permutation
matrix's flattened coordinates carries the same weak, skewed signal
(about a 1-in-n chance of being 1) regardless of which ones a level
picks, so no coordinate-*selection* strategy over that representation
substitutes for keying on the sparse structure itself, and the trie's own
overhead (building a multi-coordinate Float64 key and walking several
dict levels on every insert, 743-773ns per marginal insert at Birkhoff
sizes) outweighs whatever it saves on lookup. Rejected everywhere it was
tried; the underlying idea (recurse only where atoms are genuinely
confusable) is not ruled out for a polytope whose vertices are harder to
tell apart than either alphabet measured here, only unneeded for these
two.

## 4. `NTuple{k,Int}` keys

`microbenchmark/pattern_key_reps.jl`'s third representation of the exact
same pattern key (`rowval[1:k]`, stored as `NTuple{k,Int}` rather than a
`Vector{Int}` or a folded `UInt64`). It has a real advantage the shipped
`UInt64` fold does not: a `Dict` compares two tuples elementwise, so two
different patterns can never collide the way a folded hash can. It is
also usually a little faster: 21.3% faster on insert at Birkhoff n=25 and
6.7% at n=60, and 1.6% faster on lookup at n=25, though 1.5% *slower* on
lookup at n=60 (`microbenchmark/results_pattern_key_reps_total.csv`), so
not a clean sweep. The reason it was not shipped is an API cost, not a
speed or correctness one: staying allocation-free requires `k` to be
known to the compiler as a type parameter (`Val(k)`), and building a
fresh `Val(k)` from a runtime `Int` inside a hot per-call loop allocates
in exactly the way this representation exists to avoid (32-80 bytes per
call, confirmed empirically). `src/`'s own design requirement, that `k`
stay an ordinary keyword argument on every entry point rather than a
compile-time parameter a caller has to thread through their own API, is
incompatible with getting that speed for free. Rejected on API grounds,
not on a correctness or performance one: if a future integration is
comfortable fixing `k` once per index construction and carrying a `Val`
from there, this representation is a legitimate, slightly faster
alternative to the shipped `UInt64` fold, not a worse idea.

## 5. An identity supplied by the oracle

What the two implementations that solve this elsewhere do
(`references.md`): `copt` has each LMO return a hashable *vertex
representation* (`(sign, index)` for the L1 ball) and keys the active set
on it; `linearFW` builds a string from each 0-1 atom and maps it to its
position. The field's name for it is hash consing: one canonical identity
per distinct value, so equality becomes an identity test. It is exact and
O(1), and it was not measured here because there was nothing to gain by
measuring it: a dictionary probe on a small key costs what the folded key's
probe costs (27-42 ns per miss in `results_certificate_timing.csv`), and
the certificate decides a miss in one comparison with nothing to probe,
insert or repair. What decided against proposing it upstream is the size
of the change relative to that: every oracle in `FrankWolfe.jl` would have
to return a key beside its vertex, and the pairwise variant would only work
for oracles that do (which is exactly `copt`'s restriction), while the
certificate touches one line per call site and asks nothing of any oracle.
Where an oracle already supplies an identity, it composes: it is a
`fallback` for the certificate's tie, the same slot the folded key fills.

## 6. The fingerprint walk

The second method in `src/certificate.jl`: keep every `dot(g, a)` from the
argminmax loop and, on a lookup, compare only the atoms whose value equals
`dot(g, v)`. No argmin needed, no index, and one Float64 comparison per
atom instead of one `==`. Measured on a miss at 79.1, 160.5 and 104.5 ns at
the three real sizes (`results_certificate_timing.csv`): cheaper than the
scan by 3-13x and free of maintenance, but slower than the index's probe
and an order of magnitude slower than the certificate's one comparison,
because it walks the whole set. It is shipped, not headlined: the case
for it is a caller that has the values and no argmin, and no call site in
`FrankWolfe.jl` is in that position.

## Reading further

The methods that were kept, and why they are correct: `METHOD.md`. The
numbers: `README.md`. Every judgement call, including open questions
these rejections do not close: `DECISIONS.md`.

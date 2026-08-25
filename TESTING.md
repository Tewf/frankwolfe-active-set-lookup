# Testing the folded sparse-pattern key

`README.md`'s numbers are a measurement, not a proof: they say the folded
`UInt64` pattern key is fast, not that it is correct or that it stays
correct once atoms get inserted and removed. This file is the suite that
earns that trust, four files under `microbenchmark/`, each protecting one
invariant, plus what none of them cover yet.

Every test is seeded (`Random.Xoshiro`, an explicit `MASTER_SEED` per
file, distinct from every `run_*.jl` script's own `seed!`) and every
failure prints that seed plus enough context (alphabet, size, k, which
operation) to reproduce the exact case by rerunning the file. None of
these assert on timing: `MEASURING.md` already rules that out, and every
number here is a correctness count, not a nanosecond.

## What each test protects

**`test_equivalence.jl`: for any atom set and any query, the structure
returns exactly what the linear scan returns.** This is the one invariant
that actually matters; everything else in this repository is only worth
trusting once it holds. Checked as a property, not examples: random
seeds, active-set sizes from 0 to 500, k in {2,4,8,16}, all three
alphabets (Birkhoff permutation matrices routed to the folded pattern
key, L-infinity box corners and a generic dense control routed to the
value-prefix hash), queries that are present (first atom, last atom, a
few from the middle), and a query confirmed absent by the scan itself. A
second testset injects duplicate atoms and checks the structure agrees
with the scan on *which* index wins (the first occurrence): a bucketed
structure has no reason to preserve that for free, and this is the test
that would catch it silently returning the wrong duplicate.

**`test_lifecycle.jl`: after every single insert or delete, the structure
still agrees with a fresh scan over whatever the atoms vector currently
holds.** `active_set_cleanup!` calls `deleteat!`, which shifts every
later index down by one; nothing before this branch tested that a bucket
map's cached positions actually move with it. Randomised sequences of
~50/50 inserts and deletes (deliberately far more delete-heavy than the
real BPCG runs, which called `deleteat!` 0-2 times across 8,000-20,000
iterations, per `measurement/results.csv`: the point here is to exercise
the repair path hard, not to reproduce how rarely it fires), checked
after every operation against a fresh scan, with the specific operation
that broke it printed on failure. This is the test that stops someone
"optimising" `bucket_delete_repair!`'s current whole-index walk into an
indirection-based repair later and getting the shift direction or
boundary wrong.

**`test_fold_quality.jl`: the fold does not cluster.** At 389 real atoms
and the fold's real 64-bit width, no collision will ever be observed, so
this narrows the fold to 8/12/16/20/24 bits (via `pattern_key_uint64`'s
own `bits` keyword) and counts real collisions, at k=4 (the headline
recommendation) and n=60 (the larger real Birkhoff size), across five
independent trials of 8,000 atoms each, deduplicated on their real
pattern first so a masked collision is never a restatement of a genuine
tie in the data. Compared against two baselines: Julia's own `hash` over
the identical key, masked the same way, and the birthday approximation
`n^2 / 2^(b+1)`.

Both comparisons passed at every width tested. The birthday comparison is
the one actually doing the discriminating: fold and Julia's `hash`
turned out, while writing this test, to be *provably* related for a fixed
key length (`hash(pattern) - fold(pattern)` is an exact constant mod
2^64, confirmed across thousands of distinct real patterns), so masking
to any width preserves that relationship exactly and the two are
collision-isomorphic, not just empirically similar. Comparing against
`hash` is still a genuine sanity check (confirming Julia's own hash isn't
being reimplemented worse), it just is not the independent randomness
test it might look like; the birthday numbers (observed within roughly
0.96x-1.8x of predicted, worst case at the noisiest, lowest-count width,
24 bits) are the real evidence the fold doesn't cluster. A deliberately
broken fold (mixing only the first index instead of all k) was checked
against this suite while writing it and failed by orders of magnitude at
every width, confirming the assertion is not vacuous.

**`test_dispatch.jl`: sparse atoms route to the pattern key, dense atoms
route to the value-prefix hash, and both agree with the scan.** This is
the argument for the proposed API shape: `_unsafe_equal` itself already
dispatches this way (one method on `Array`, one on
`SparseArrays.AbstractSparseArray`), and `test_atom_generators.jl`'s
`route_build`/`route_lookup`/`route_scan` are that same routing, written
as ordinary Julia methods rather than an `if atom isa SparseMatrixCSC`
branch. Checked directly (`isa` on the built index, not just "the answer
came out right", since a routing bug that happened to still produce a
correct answer through the wrong structure would slip past a
correctness-only check) and under an interleaved mixed-alphabet batch, so
the dispatch is shown not to depend on call order.

## What was found

Nothing broke. All four files pass against the shipped implementation
without any change to `sparse_pattern.jl`, `pattern_key_reps.jl`, or
`bucket_lifecycle.jl`; the lifecycle test in particular is a regression
guard for changes not yet made; a deliberate mutation of
`bucket_delete_repair!` (disabling its shift, and separately reversing
its direction) was used to confirm it fails loudly rather than pass
vacuously. The one real finding is the fold/`hash` isomorphism described
above, which is a property of this Julia version's `hash` implementation,
not a bug.

## Running the suite

```
source ~/miniforge3/etc/profile.d/conda.sh && conda activate frankwolfe
julia --project=. microbenchmark/test_equivalence.jl
julia --project=. microbenchmark/test_lifecycle.jl
julia --project=. microbenchmark/test_dispatch.jl
julia --project=. microbenchmark/test_fold_quality.jl
```

All four run in `.github/workflows/ci.yml` on every push, alongside the
existing `test_soundness.jl` and `test_pattern_key_reps.jl`. Total added
runtime is about 20 seconds, `test_fold_quality.jl`'s ~13 seconds of real
atom generation being the only part that isn't sub-second.

## What is deliberately not tested yet

This is the part a maintainer should actually read before trusting the
suite above.

- **No real `FrankWolfe.jl` `ActiveSet`.** Every test here runs a
  standalone harness (`atoms::Vector`, a bucket-map index built and
  mutated by hand) that mirrors `find_atom`/`active_set_update!`/
  `active_set_cleanup!`'s shape, not the actual package types or a real
  BPCG run driving them. Nothing here would catch an integration bug in
  how a hashed index gets threaded through `ActiveSet`'s own fields.
- **PFW and BCG are untested.** `DECISIONS.md`'s own draft issue comment
  lists their call sites (`pairwise.jl:242`, `blended_cg.jl:358`,
  `corrective_frankwolfe.jl:240`, `block_coordinate_algorithms.jl:421`) as
  found by grep, never run; this suite inherits that gap unchanged. The
  lifecycle test's insert/delete mix is synthetic, not drawn from any of
  these algorithms' real call patterns.
- **Only three atom alphabets.** Birkhoff permutation matrices, L-infinity
  box corners, and a generic dense control, the same three
  `run_pattern_key_reps.jl` already swept. Any other polytope's LMO
  (simplex, general sparse structures beyond one-nonzero-per-row-and-column)
  is untested, and the pattern key's own scope note
  (`sparse_pattern.jl`) already flags that it needs exactly that
  structure to be discriminating.
- **No concurrency.** Every structure here is a plain `Dict`, mutated
  in-place, single-threaded. Nothing checks behaviour under concurrent
  reads and writes, and `Dict` itself gives no such guarantee.
- **Active-set sizes stop at 500.** The three real runs topped out at
  158-389 atoms; `test_equivalence.jl` goes up to 500 as a margin, not
  because anything here was run at the thousands-of-atoms scale a looser
  `epsilon` or a different algorithm could reach. Whether the folded key
  keeps its near-zero collision rate at, say, 50,000 atoms (where a
  64-bit fold's birthday bound is still comfortable but no longer
  astronomically so) is unmeasured.
- **The signed-zero hazard is not re-tested here.** It stays covered by
  the existing `test_soundness.jl`; this suite's property tests use only
  real atoms from `BirkhoffPolytopeLMO`/`LpNormBallLMO{Inf}`/`rand`, none
  of which ever produce a `-0.0` (confirmed in `test_soundness.jl`'s own
  comments), so the hazard never has a chance to fire here. The
  `DECISIONS.md`-recommended canonicalisation fix (`prefix .+ 0.0`) is
  still not applied to `lookup_methods.jl`'s shipped `PrefixIndex`; that
  gap is pre-existing and unrelated to this test suite, not something the
  property tests happen to paper over.

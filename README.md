# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Lire en français](README.fr.md)

Longer than 80 lines because a stranger's entry point needs the answer, a
usable code sample, and the reproduction commands together on one page,
not spread across links before anything has been said.

**Yes: a folded sparse-pattern key beats `FrankWolfe.jl`'s linear-scan
`find_atom` at the active-set sizes real runs actually reach, at 0.812ns
per iteration (lookup, insert, and occasional delete-repair, weighted by
how often each actually happens) against the scan's 57.49ns at Birkhoff
n=60 (k=4, the default; see "What was measured" below).** This answers
[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
open since 2021 with no replies. It took several reversals to get here,
including a first draft that concluded the opposite; `REJECTED.md` says
what was tried and refused, with the numbers, and is the most useful file
in this repository for anyone about to propose one of those again.

## The method, in three sentences

For a sparse atom (a Birkhoff-polytope permutation matrix, say), the
values are all 1.0 and carry no information, so the key hashes the first
`k` stored *positions* instead, folded into one `UInt64`. For a dense atom
(an L-infinity-ball box corner), there is no sparse structure to read, so
the key hashes the first `k` coordinate *values* instead. Either way the
key only ever picks a bucket: every candidate in it is confirmed against
the whole atom with exact equality before being trusted, so a collision
costs one comparison and never a wrong answer. `METHOD.md` has the full
argument, including a real soundness gap around the sign of zero and how
it is closed.

## Using it

```julia
include("src/ActiveSetLookup.jl")
using .ActiveSetLookup

# atoms is a Vector you already own (e.g. an ActiveSet's own atoms field).
index = build_index(atoms)              # k defaults to DEFAULT_K = 4
pos = lookup_atom(index, atoms, query)  # -1 if absent, mirrors find_atom
push_atom!(index, atoms, new_atom)      # keeps atoms and index in sync
delete_atom!(index, atoms, pos)         # repairs the index after deleteat!
```

`k` is an ordinary keyword everywhere (`build_index(atoms; k=8)`), never a
compile-time type parameter; `DEFAULT_K = 4` is measured, not guessed
(`src/keys.jl`). `build_index` dispatches on the atoms' element type to
pick the right key automatically: `SparseMatrixCSC` routes to the pattern
key, a dense `Array` to the value key. The whole module is three small
files, `src/keys.jl`, `src/confirm.jl`, `src/index.jl`, each readable on
its own.

## What was measured

`measurement/` runs blended pairwise conditional gradient (BPCG,
`lazy=true`) on three problems, instrumented without editing
`FrankWolfe.jl` itself:

| Problem | Dimension | Iterations | Max active set | `find_atom` calls | Hits |
|---|---|---|---|---|---|
| Birkhoff, n=25 | 625 | 8,000 | 158 | 159 | 0 |
| Birkhoff, n=60 | 3,600 | 20,000 | 389 | 389 | 0 |
| L-infinity ball, d=3,000 | 3,000 | 15,000 | 241 | 240 | 0 |

Every call was a miss followed by a `push!`, so the real per-iteration
cost is lookup plus insert, and `active_set_cleanup!`'s occasional
`deleteat!` (2, 1, and 0 times across the three runs) means a structure
also has to survive deletion. `microbenchmark/` costs each candidate
structure the same way, separately for lookup, insert, and delete-repair,
then weights them by these three runs' own call rates. The sparse-pattern
key wins both Birkhoff sizes on total per-iteration cost (2.72ns at
n=25's 158-atom maximum, 1.88ns at n=60's 389-atom maximum, before the
allocation-free `UInt64` fold made both faster still); the existing
value-prefix hash remains the right answer for the L-infinity ball
(3.31ns at 241 atoms), since it has no sparse structure to exploit.

## Reproducing the measurements

```
source ~/miniforge3/etc/profile.d/conda.sh && conda activate frankwolfe
julia --project=. measurement/run.jl                    # measurement/results.csv
julia --project=. microbenchmark/run.jl                 # the full-hash-vs-scan sweep
julia --project=. microbenchmark/run_prefix.jl          # the value-prefix-hash sweep
julia --project=. microbenchmark/run_lifecycle.jl       # sparse-pattern key + trie, full lifecycle
julia --project=. microbenchmark/run_pattern_key_reps.jl # UInt64/NTuple/Vector{Int} representations
julia --project=. test/test_public_api.jl               # this module's own correctness tests
```

Every script writes a committed `results*.csv` beside itself; `MEASURING.md`
states, once, how every timing here was taken and what it does not claim
(none are asserted by CI, only run end to end). `TESTING.md` covers the
six correctness test files under `microbenchmark/` and what each protects.

## Reading further

How the method works and why it is correct, including the signed-zero
subtlety: `METHOD.md`. What was tried and refused, with numbers:
`REJECTED.md`. What each test protects: `TESTING.md`. The machine and its
noise: `MEASURING.md`. Every judgement call, open questions, and the draft
issue comment: `DECISIONS.md`. Papers, call sites, the issue itself:
`references.md`. Every file, one line each: `what-is-where.md`.

# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Lire en français](README.fr.md)

Longer than 80 lines because a stranger's entry point needs the answer, a
usable code sample, and the reproduction commands together on one page,
not spread across links before anything has been said.

**The vertex the LMO just returned is not in the active set whenever
`<g,v> < <g,s>`, where `s` is the active atom the step has already found
best for the gradient `g`. That one comparison replaces `find_atom`'s
linear scan at every call site in `FrankWolfe.jl`, at about 10 ns against
the scan's 2,042 ns on Birkhoff n=60, with no index to build, insert into
or repair; a tie is settled by one comparison with `s`. On real runs it
decided every call: 788 of 788 in blended pairwise conditional gradients,
20,001 of 20,001 in pairwise Frank-Wolfe on Birkhoff n=60, where the
active set reaches 9,368 atoms and the scan is 6.4% of the run.** For a
caller that has only the atom, a folded structural hash with exact
confirmation beats the scan at every measured size and stays here as the
fall-back. This answers
[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
open since 2021 with no replies. It took several reversals to get here,
including a first draft that concluded hashing does not help and a second
that stopped at the hash; `REJECTED.md` says what was tried and refused,
with the numbers.

**Never opened FrankWolfe.jl?** Start with [`guide/`](guide/README.md):
both algorithms reimplemented in four short files with no dependency on
the library, the membership question answered three ways side by side, a
narrated run and a test of every claim.

## The method, in three sentences

Every Frank-Wolfe step that keeps an active set first minimises `<g, a>`
over the active atoms `a` (that is how it finds its local and away
vertices), then asks the LMO for the vertex `v` and computes `<g, v>` for
the dual gap; if `v` were already active, `<g, v>` would be one of the
values just minimised, so `<g, v> < <g, s>` proves it is not. When the
comparison fails, `v` ties the best atom `s`, and either `v == s`, which
one exact comparison settles (this is the pairwise Frank-Wolfe case, where
the LMO routinely returns an active vertex), or `v` ties a *different*
atom, which is the only case that searches at all and has probability zero
for a real-valued gradient. In blended pairwise conditional gradients the
FW step is taken only when `v` beats the whole active set, so there the
scan was never needed in exact arithmetic, and the certificate is the
floating-point-safe way to skip it; `METHOD.md` has both arguments and
what they rest on.

## Using it

```julia
include("src/ActiveSetLookup.jl")
using .ActiveSetLookup

# Inside a step. `atoms` is the active set's own Vector, `g` the gradient;
# the step has just minimised dot(g, a) over atoms (active_set_argminmax),
# so it holds `best` and `best_value`; the LMO returned `v`, and dot(g, v)
# was computed for the dual gap.
pos = certified_lookup(atoms, v, dot(g, v), best, best_value)  # -1 if absent, mirrors find_atom

# A caller that has only the atom keeps an index instead.
index = build_index(atoms)              # k defaults to DEFAULT_K = 4
pos = lookup_atom(index, atoms, query)  # -1 if absent
push_atom!(index, atoms, new_atom)      # keeps atoms and index in sync
delete_atom!(index, atoms, pos)         # repairs the index after deleteat!

# The index is also the natural fall-back for the certificate's rare tie.
pos = certified_lookup(atoms, v, dot(g, v), best, best_value;
                       fallback=(a, q) -> lookup_atom(index, a, q))
```

`certified_lookup` needs nothing of the atom's type: the same call serves
permutation matrices, box corners and anything else `dot` accepts. The
index dispatches on the atoms' element type (`SparseMatrixCSC` and
`SparseVector` to a key over stored *positions*, a dense `Array` to a key
over leading *values*), and every bucket hit is confirmed against the whole
atom before it is trusted, so a collision costs a comparison and never a
wrong answer. The module is four small files, `src/certificate.jl`,
`src/keys.jl`, `src/confirm.jl`, `src/index.jl`, each readable on its own.

## What was measured

`measurement/` runs four algorithms on three problems, instrumented
without editing `FrankWolfe.jl`: at every `find_atom` call it records what
the scan answered and, for the LMO's own vertex, what the certificate
would have decided, then checks the two agree.

| Algorithm | Problem | Max active set | `find_atom` calls | Hits | Scan's share of run | Certified absent | Tie, `v == s` | Tie, search needed |
|---|---|---|---|---|---|---|---|---|
| BPCG, lazy | Birkhoff n=25 | 158 | 159 | 0 | 0.07% | 159 | 0 | 0 |
| BPCG, lazy | Birkhoff n=60 | 389 | 389 | 0 | 0.06% | 389 | 0 | 0 |
| BPCG, lazy | L-inf ball d=3,000 | 241 | 240 | 0 | 0.02% | 240 | 0 | 0 |
| PFW | Birkhoff n=25 | 2,463 | 8,001 | 5,535 | 6.05% | 2,466 | 5,535 | 0 |
| PFW | Birkhoff n=60 | 9,368 | 20,001 | 10,628 | 6.41% | 9,373 | 10,628 | 0 |
| PFW | L-inf ball d=3,000 | 2,613 | 15,001 | 12,387 | 0.81% | 2,614 | 12,387 | 0 |
| PFW, lazy | Birkhoff n=25 | 181 | 204 | 0 | 0.17% | 204 | 0 | 0 |
| PFW, lazy | Birkhoff n=60 | 576 | 598 | 0 | 0.12% | 598 | 0 | 0 |
| PFW, lazy | L-inf ball d=3,000 | 295 | 323 | 0 | 0.03% | 323 | 0 | 0 |
| BCG | Birkhoff n=25 | 148 | 7,948 | 7,800 | 0.14% | 148 | 38 | 0 |
| BCG | Birkhoff n=60 | 623 | 19,988 | 19,365 | 0.14% | 623 | 126 | 0 |
| BCG | L-inf ball d=3,000 | 348 | 1,827 | 1,360 | 0.06% | 467 | 17 | 0 |

(`measurement/results.csv`, `measurement/results_algorithms.csv`; 8,000,
20,000 and 15,000 iterations, `epsilon=1e-9`.) The certificate never
contradicted the scan, and no call ever needed a search. Three things the
table says that the earlier BPCG-only measurement could not: the scan is
a real cost in non-lazy pairwise Frank-Wolfe, whose active set grows to
thousands of atoms (6.4% of a 99 s run at n=60); every one of its 53-83%
hits was the best atom itself, so one comparison found it; and BCG's
lookups are mostly for an atom it took out of its own active set a moment
earlier (`lp_separation_oracle` returns the atom without its position),
which is a separate, index-only waste.

Per call, on the same atoms in the same session
(`microbenchmark/results_certificate_timing.csv`, ns):

| Size | Case | Certificate | Folded key / prefix hash | Scan |
|---|---|---|---|---|
| Birkhoff n=25, 158 atoms | miss | 10.5 | 28.0 | 777.6 |
| | hit (best atom) | 87.4 | 106.9 | 200.3 |
| Birkhoff n=60, 389 atoms | miss | 10.3 | 27.4 | 2,041.8 |
| | hit (best atom) | 177.8 | 196.3 | 436.9 |
| L-inf ball d=3,000, 241 atoms | miss | 10.2 | 41.8 | 322.1 |
| | hit (best atom) | 1,392.7 | 1,448.1 | 1,619.9 |

A miss costs the certificate one Float64 comparison (the 10 ns is the
timer's floor for a call, the same for every column); a hit costs every
method the one exact comparison, which on a 3,000-coordinate dense atom is
the 1.4 µs in the last row. A forced tie with a distinct atom falls back
to the scan and costs what the scan costs (798, 2,091 and 357 ns), and no
real run produced one. Weighted by BPCG's real call rates, the total per
iteration is 0.16-0.21 ns for the certificate against 1.0-3.1 ns for the
index (its insert and delete-repair costs from its own committed sweep)
and 5-40 ns for the scan (`results_certificate_total.csv`). For pairwise
Frank-Wolfe on Birkhoff n=60, the 6.35 s the scan spent across 20,001
calls becomes, at the per-call costs above, about 2 ms: that figure is
arithmetic from measured parts, not an end-to-end run, because reaching it
means changing `pairwise.jl` itself.

## Reproducing the measurements

With Julia 1.10 or later, from the repository root:

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. measurement/run.jl                     # results.csv, results_algorithms.csv
julia --project=. microbenchmark/run_certificate.jl      # the certificate against key, hash and scan
julia --project=. microbenchmark/run.jl                  # the full-hash-vs-scan sweep
julia --project=. microbenchmark/run_prefix.jl           # the value-prefix-hash sweep
julia --project=. microbenchmark/run_lifecycle.jl        # sparse-pattern key + trie, full lifecycle
julia --project=. microbenchmark/run_pattern_key_reps.jl # UInt64/NTuple/Vector{Int} representations
julia --project=. test/test_certificate.jl               # the certificate's own correctness tests
julia --project=. test/test_public_api.jl                # the index's own correctness tests
```

Every script writes a committed `results*.csv` beside itself; `MEASURING.md`
states, once, how every timing here was taken and what it does not claim
(none are asserted by CI, only run end to end; the one assertion is that
the certificate never contradicts the scan). `TESTING.md` covers every
correctness test and what each protects.

## Reading further

How both methods work and why they are correct, including why BPCG never
needed the scan and the signed-zero subtlety of the dense key: `METHOD.md`.
What was tried and refused, with numbers, and how the two implementations
that do solve this elsewhere (`copt`, `linearFW`) do it: `REJECTED.md`,
`references.md`. What each test protects: `TESTING.md`. The machine and
its noise: `MEASURING.md`. Every judgement call and open question:
`DECISIONS.md`. Every file, one line each: `what-is-where.md`.

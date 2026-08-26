# What is where

```
Project.toml               dependencies: FrankWolfe, TimerOutputs, LinearAlgebra
Manifest.toml               the exact resolved versions; gitignored, see .gitignore
src/                        the method itself, usable without reading a benchmark
  ActiveSetLookup.jl          the module: includes the three files below and
                               re-exports their public names in one place
  keys.jl                     computing a key from an atom: the folded sparse
                               pattern key, the dense value key, DEFAULT_K=4
  confirm.jl                  the confirmation step: exact equality, mirroring
                               FrankWolfe.jl's own _unsafe_equal dispatch
  index.jl                    the index structure: build_index, lookup_atom,
                               push_atom!, delete_atom! over a caller-owned
                               atoms Vector, plus bucket_health, whose mean
                               bucket size shows METHOD.md's precondition
                               going unmet
test/                       tests for src/, independent of microbenchmark/'s
                             own test suite (which tests the comparisons that
                             led to this design, not this module's code)
  test_public_api.jl          @test: build_index dispatch, lookup-vs-scan
                               equivalence, push_atom!/delete_atom! lifecycle,
                               the signed-zero fix, a forced fold collision
measurement/                a real BPCG run, instrumented, on problems where the
                             active set grows
  instrumentation.jl          times, counts, and counts hits for find_atom without
                               editing FrankWolfe.jl
  problems.jl                  the three problems: two Birkhoff sizes, one L∞-ball
  run.jl                        runs all three, writes results.csv
  results.csv                    max/mean active-set size, calls, hits, lookup share,
                                  per problem
microbenchmark/              the lookup itself, isolated from any solver
  lookup_methods.jl            linear scan, full-atom Dict, and prefix-hash lookup
                                (dense and sparse-atom variants), copied from find_atom
  timing.jl                     warm up, batch to clear 1 ms, fastest of five, plus a
                                 variant that cycles through a pre-generated query mix
  run.jl                        sweeps size x dimension x atom scenario, writes results.csv
  results.csv                    scan and dict time per point, generic and adversarial atoms
  run_prefix.jl                  sweeps size x k x query type x atom alphabet against
                                  FrankWolfe.jl's own atom shapes, writes the three
                                  results_prefix_*.csv files below
  results_prefix_timing.csv        scan and prefix-hash time per point
  results_prefix_collisions.csv    bucket collision stats per point
  results_prefix_crossover.csv     smallest size where each prefix-hash variant beats
                                    the scan, or none
  sparse_pattern.jl              idea 1: key a Birkhoff atom on SparseMatrixCSC's own
                                  rowval instead of a flattened value prefix
  hash_trie.jl                    idea 2: a recursive hash trie over coordinate blocks,
                                   with four coordinate-selection strategies
  bucket_lifecycle.jl              insert! and delete-repair bookkeeping shared by any
                                    flat key->Vector{Int} bucket index (prefix, pattern)
  run_lifecycle.jl                  sweeps trie (k, max_depth, order) by collision stats,
                                     then times lookup, insert and deletion-repair for
                                     scan/prefix/pattern/trie, writes the four results
                                     files below and the total per-iteration comparison
  results_lifecycle_collisions.csv    the trie sweep's collision/depth stats
  results_lifecycle_timing.csv         lookup/insert/delete-repair ns per structure
  results_lifecycle_total.csv           the same, weighted by each real run's own
                                         call rate, into one total ns per iteration
  results_sparse_pattern_collisions.csv  the sparse-pattern key's own collision sweep,
                                          for the direct 9-buckets-become-N comparison
  test_soundness.jl              @test, not just a script: the signed-zero hazard,
                                  checked against the prefix hash, the pattern key
                                  (immune by construction) and the trie key (closed the
                                  same way); the one file in this repository whose
                                  failure means something is actually wrong, not noisy
  pattern_key_reps.jl             two allocation-free representations of the same
                                   pattern key (UInt64 incremental hash, NTuple{K,Int})
  run_pattern_key_reps.jl          sweeps k x representation, measures allocation
                                    (@allocated/@allocations) and time for lookup,
                                    insert and deletion-repair, writes the four
                                    results files below
  results_pattern_key_reps_timing.csv       lookup/insert/delete-repair ns per representation
  results_pattern_key_reps_allocations.csv   bytes and allocation count per representation
  results_pattern_key_reps_collisions.csv     bucket collision stats per representation
  results_pattern_key_reps_total.csv           weighted by each real run's own call rate,
                                                into one total ns and one total bytes per
                                                iteration
  test_pattern_key_reps.jl        @test: forces a real UInt64 fold collision (two
                                   different real patterns, narrowed to 1 bit) and checks
                                   the structure still answers correctly; confirms the
                                   NTuple key has no equivalent hazard; cross-checks all
                                   three representations agree with each other
  test_atom_generators.jl         real-atom generators for all three alphabets, plus the
                                   sparse/dense dispatch (route_build/lookup/scan) shared
                                   by the four confirm-and-sustain test files below
  test_equivalence.jl             @test: property test, structure vs. scan, seeds x sizes
                                   (0-500) x k (2/4/8/16) x alphabet, plus duplicates
  test_lifecycle.jl               @test: randomised insert/delete sequences, checked
                                   against a fresh scan after every single operation
  test_fold_quality.jl            @test: the UInt64 fold's collision rate at narrowed bit
                                   widths vs. Julia's own hash and the birthday prediction
  test_dispatch.jl                @test: sparse routes to the pattern key, dense to the
                                   value prefix, both agree with the scan
TESTING.md                  what each test protects, how to run the suite, what is not
                             tested yet
METHOD.md                  how the folded sparse-pattern key works and why it is
                             correct, including the signed-zero subtlety
REJECTED.md                 what was tried and refused, with the numbers that
                             killed it; the most useful file for the next person
                             who proposes hashing the whole atom, or a trie
README.md                  the question, the numbers, the answer; leads with the
                             answer, then a usable code sample and how to reproduce
README.fr.md                the same, in French, marked as wanting a native pass
references.md               the papers and the issue, cited never redistributed
MEASURING.md                the machine, the noise, and what a number here does not claim
DECISIONS.md                what is Mohamed's to decide, including the draft issue comment
CITATION.cff                how to cite this repository
LICENSE                     MIT
.github/workflows/ci.yml    runs every script on every push and asserts
                             every test_*.jl file's @tests; times nothing
explain/                    gitignored; HTML explainers for one reader, never shipped
```

Two directories, one house style: `measurement/` asks whether the lookup
costs anything in a real run; `microbenchmark/` asks how much the lookup
itself would cost either way, with no solver around it to hide the answer
in. None of the scripts needs another script's own *code*
(`measurement/run.jl` needs `instrumentation.jl` and `problems.jl` only;
`microbenchmark/run.jl` needs `lookup_methods.jl` and `timing.jl` only;
`microbenchmark/run_prefix.jl` and `run_lifecycle.jl` need the same two
plus whichever of `sparse_pattern.jl`/`hash_trie.jl`/`bucket_lifecycle.jl`
they use, plus the `FrankWolfe` package itself to generate atoms in the
exact shape `measurement/problems.jl` uses, never any of `measurement/`'s
own modules; `run_pattern_key_reps.jl` needs `lookup_methods.jl`,
`sparse_pattern.jl`, `pattern_key_reps.jl`, `bucket_lifecycle.jl` and
`timing.jl`, the same shape one level further), so any script can be
read, or rerun, on its own. One disclosed exception, data rather than
code: `run_lifecycle.jl`'s and `run_pattern_key_reps.jl`'s final steps
each read `measurement/results.csv` as plain CSV, to weight their own
measured lookup/insert/repair costs by the real per-iteration call rates
that file records, rather than reimplementing a second BPCG harness to
get the same numbers.

`test_atom_generators.jl` is the one deliberate exception to "no script
needs another *test* script's code": `test_equivalence.jl`,
`test_lifecycle.jl`, and `test_dispatch.jl` all include it, for the same
atom generators and the same sparse/dense dispatch, so a routing bug
shows up in all three rather than being defined three different ways.
`test_fold_quality.jl` needs no such sharing (it only ever touches
Birkhoff atoms) and stays self-contained like `test_soundness.jl` and
`test_pattern_key_reps.jl` before it.

`src/` is a third, separate concern, not a directory this house style's
dependency rules apply to the same way: it needs nothing from
`measurement/` or `microbenchmark/`, only `SparseArrays` from Julia's own
standard library, and nothing in `measurement/` or `microbenchmark/`
needs anything from it either. That is deliberate, not an oversight: the
research harness's own copies of the key/index/confirm logic
(`sparse_pattern.jl`, `pattern_key_reps.jl`, `bucket_lifecycle.jl`,
`lookup_methods.jl`) exist to compare several structures against each
other and against a scan, which is a different job from `src/`'s, "be the
one structure a caller actually uses," and sharing code between the two
would mean either the benchmark harness importing from the very module
its own numbers are meant to justify, or `src/` depending on files whose
job is to keep changing as new ideas get swept. `test/test_public_api.jl`
depends only on `src/`, for the same reason: it is the test for the
module a stranger would use, not another copy of the research harness's
own equivalence checks (which stay under `microbenchmark/`, and are
listed above).

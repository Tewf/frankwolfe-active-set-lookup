# What is where

```
Project.toml               dependencies: FrankWolfe, TimerOutputs, LinearAlgebra
Manifest.toml               the exact resolved versions; gitignored, see .gitignore
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
README.md                  the question, the numbers, the answer; leads with the answer
README.fr.md                the same, in French, updated for the new answer; the
                             "Prefix hashing" section's detail was not translated,
                             see DECISIONS.md
references.md               the papers and the issue, cited never redistributed
MEASURING.md                the machine, the noise, and what a number here does not claim
DECISIONS.md                what is Mohamed's to decide, including the draft issue comment
CITATION.cff                how to cite this repository
LICENSE                     MIT
.github/workflows/ci.yml    runs every script on every push and asserts
                             test_soundness.jl's @tests; times nothing
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

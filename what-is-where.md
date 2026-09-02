# What is where

```
Project.toml                the package ActiveSetLookup: name, version, its
                             one dependency (SparseArrays), compat bounds, and
                             the test target (Aqua, FrankWolfe for real atoms)
CONTRIBUTING.md             how this repository is written, how to work on it,
                             and what went upstream, for the maintainers
Manifest.toml               the exact resolved versions; gitignored, see .gitignore
src/                        the method itself, usable without reading a benchmark
  ActiveSetLookup.jl          the module: includes the four files below and
                               re-exports their public names in one place
  certificate.jl              the lookup that does not search: certified_lookup
                               proves the LMO vertex absent from two inner
                               products the FW step already holds, compares it
                               with the best atom on a tie, and searches only
                               on a tie between distinct atoms
  keys.jl                     computing a key from an atom: the folded sparse
                               pattern key, the dense value key, DEFAULT_K=4
  confirm.jl                  the confirmation step: exact equality, mirroring
                               FrankWolfe.jl's own _unsafe_equal dispatch
  index.jl                    the index structure: build_index, lookup_atom,
                               push_atom!, delete_atom! over a caller-owned
                               atoms Vector, plus bucket_health, whose mean
                               bucket size shows METHOD.md's precondition
                               going unmet
guide/                      Frank-Wolfe from zero for a reader who has never
                             opened the library: no dependency on FrankWolfe.jl
  README.md                   the walkthrough: the problem, the vocabulary, the
                               two algorithms, the three answers, a worked n=3
                               example, what the real library adds
  birkhoff.jl                 the toy polytope: permutations, vertices, a
                               brute-force LMO, the quadratic objective
  active_set.jl               the active set (atoms, weights, x) and the two
                               moves that never add or remove an atom
  lookups.jl                  the three answers side by side, each owning the
                               append and removal of atoms; CrossChecked runs
                               all three and records any disagreement
  frank_wolfe.jl              plain Frank-Wolfe, and blended pairwise with the
                               library's step rule; every iteration recorded
  run.jl                      both algorithms on one problem, narrated
test/                       tests for src/ and guide/, independent of
                             microbenchmark/'s own test suite (which tests the
                             comparisons that led to this design)
  runtests.jl                 Pkg.test()'s entry: Aqua's package-quality checks,
                               then the three suites below, each in its own module
  test_public_api.jl          @test: build_index dispatch, lookup-vs-scan
                               equivalence, push_atom!/delete_atom! lifecycle,
                               the signed-zero fix, a forced fold collision
  test_certificate.jl         @test: dot is a function of values (equal inputs,
                               either order, any alignment); certified_lookup
                               agrees with a scan; an LMO vertex never reaches
                               the fall-back except on a tie; ties, duplicates,
                               NaN/Inf gradients, signed zero
  test_guide.jl               @test: the guide's oracle minimises, plain FW
                               meets the 1/t bound inside the polytope, the
                               active set is a valid mixture after every step, the three
                               lookups agree, the oracle's vertex is never active
                               in the FW branch, the worked n=3 example
measurement/                real runs, instrumented, on problems where the
                             active set grows
  Project.toml                the harness's own environment: FrankWolfe.jl and
                               TimerOutputs, kept out of the package's dependencies
  instrumentation.jl          times, counts, and counts hits for find_atom without
                               editing FrankWolfe.jl; RecordingLMO keeps the
                               gradient so the certificate can be tallied at
                               every call and checked against the scan
  problems.jl                  the three problems: two Birkhoff sizes, one L∞-ball
  run.jl                        runs BPCG, PFW, lazy PFW and BCG on all three,
                                 writes the two files below
  results.csv                    BPCG: max/mean active-set size, calls, hits,
                                  deletions, lookup share, certificate tally;
                                  the per-iteration rates the sweeps weight by
  results_algorithms.csv         the same columns for PFW, lazy PFW and BCG
  run_end_to_end.jl              PFW and BPCG timed on stock FrankWolfe.jl and on the
                                  certificate branch, two processes, fastest of three
  results_end_to_end.csv         its rows: both variants, times, counts, final primal
microbenchmark/              the lookup itself, isolated from any solver
  Project.toml                 the sweeps' own environment: FrankWolfe.jl for real
                                atoms, kept out of the package's dependencies
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
  certificate.jl                  the research copy of the absence certificate, kept
                                   apart from src/ like every other structure here
  run_certificate.jl              times the certificate on a miss, a hit, a forced
                                   tie and the fingerprint walk, beside the pattern
                                   key, the prefix hash and the scan on the same
                                   atoms, at the three real sizes; writes the two
                                   results files below
  results_certificate_timing.csv    ns per lookup, per structure and case
  results_certificate_total.csv     per-iteration total at the real call rates; the
                                     certificate's insert and repair cost is zero,
                                     the baselines' are taken from their own sweeps
  test_atom_generators.jl         real-atom generators for all three alphabets, plus the
                                   sparse/dense dispatch (route_build/lookup/scan) shared
                                   by the four correctness-test stage test files below
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
README.fr.md                the same, in French; the English version prevails
                             where the two differ
references.md               the papers and the issue, cited never redistributed
MEASURING.md                the machine, the noise, and what a number here does not claim
DECISIONS.md                every judgement call this repository made, and
                            the ones still open
CITATION.cff                how to cite this repository
CHANGELOG.md                what each version added, and what is unreleased
LICENSE                     MIT
.github/workflows/ci.yml    Pkg.test() on Julia 1.10 and the newest release, then
                             every script end to end and every microbenchmark
                             test_*.jl's @tests; times nothing
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
